use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    sequence::{Sequence, SequenceType},
};
use std::io::Write;
use std::path::Path;

pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut lines = text.lines();

    let header = lines.next().unwrap_or("");
    if !header.to_ascii_uppercase().starts_with("CLUSTAL") {
        return Err(FormatError::Parse(
            "Not a ClustalW file: missing CLUSTAL header".into(),
        ));
    }

    // order of first appearance
    let mut order: Vec<String> = Vec::new();
    // accumulated residues per name
    let mut residues: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    for line in lines {
        // blank line — block separator, ignore
        if line.trim().is_empty() {
            continue;
        }
        // conservation/ruler line starts with whitespace
        if line.starts_with(' ') || line.starts_with('\t') {
            continue;
        }

        // sequence line: split on whitespace → (name, residues[, position_count])
        let mut parts = line.split_whitespace();
        let seq_name = match parts.next() {
            Some(n) => n.to_string(),
            None => continue,
        };
        let res_str = match parts.next() {
            Some(r) => r,
            None => continue,
        };
        // third token (position counter) is intentionally ignored

        if !residues.contains_key(&seq_name) {
            order.push(seq_name.clone());
            residues.insert(seq_name.clone(), Vec::new());
        }
        let entry = residues.get_mut(&seq_name).unwrap();
        for b in res_str.bytes() {
            entry.push(b);
        }
    }

    if order.is_empty() {
        return Err(FormatError::Parse(
            "No sequences found in Clustal data".into(),
        ));
    }

    let mut aln = Alignment::new(name);
    for seq_name in &order {
        let res = residues.remove(seq_name).unwrap_or_default();
        let seq_type = guess_type(&res);
        let mut seq = Sequence::new(seq_name.clone(), res);
        seq.seq_type = seq_type;
        aln.push(seq);
    }

    Ok(aln)
}

fn guess_type(residues: &[u8]) -> SequenceType {
    let non_gap: Vec<u8> = residues
        .iter()
        .copied()
        .filter(|&b| !matches!(b, b'-' | b'~' | b'.' | b'*'))
        .collect();

    if non_gap.is_empty() {
        return SequenceType::Unknown;
    }

    let nuc_set: &[u8] = b"ACGTURYN";
    let is_protein = non_gap
        .iter()
        .any(|&b| !nuc_set.contains(&b.to_ascii_uppercase()));

    if is_protein {
        return SequenceType::Protein;
    }

    let has_u = non_gap.iter().any(|&b| b.eq_ignore_ascii_case(&b'U'));
    if has_u {
        SequenceType::Rna
    } else {
        SequenceType::Dna
    }
}

pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    let content = to_string(aln)?;
    std::fs::write(path, content).map_err(FormatError::Io)
}

pub fn to_string(aln: &Alignment) -> Result<String> {
    let mut buf = Vec::new();
    write_writer(aln, &mut buf)?;
    String::from_utf8(buf).map_err(|e| FormatError::Parse(e.to_string()))
}

pub fn write_writer(aln: &Alignment, w: &mut dyn Write) -> Result<()> {
    writeln!(w, "CLUSTAL W (2.1) multiple sequence alignments").map_err(FormatError::Io)?;
    writeln!(w).map_err(FormatError::Io)?;

    if aln.is_empty() {
        return Ok(());
    }

    let max_len = aln
        .sequences
        .iter()
        .map(|s| s.residues.len())
        .max()
        .unwrap_or(0);
    let block_size = 60usize;
    let mut offset = 0;
    while offset < max_len {
        let end = (offset + block_size).min(max_len);
        for seq in &aln.sequences {
            let chunk = if offset < seq.residues.len() {
                &seq.residues[offset..end.min(seq.residues.len())]
            } else {
                &[]
            };
            let name_padded = format!("{:<16}", seq.name);
            write!(w, "{}", name_padded).map_err(FormatError::Io)?;
            w.write_all(chunk).map_err(FormatError::Io)?;
            writeln!(w).map_err(FormatError::Io)?;
        }
        writeln!(w).map_err(FormatError::Io)?;
        offset += block_size;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = "CLUSTAL W (1.83) multiple sequence alignments\n\nseq1  ACGT--ACGT\nseq2  ACGTAAACGT\n      ****  ****\n\nseq1  TTTT\nseq2  TTTT\n      ****\n";

    #[test]
    fn basic_parse() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.seq_count(), 2);
        assert_eq!(aln.sequences[0].residues.len(), 14);
    }

    #[test]
    fn names_correct() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.sequences[0].name, "seq1");
        assert_eq!(aln.sequences[1].name, "seq2");
    }

    #[test]
    fn roundtrip_preserves_residues() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        let out = super::to_string(&aln).unwrap();
        let aln2 = super::parse_str(&out, "test2").unwrap();
        assert_eq!(aln.seq_count(), aln2.seq_count());
        for (a, b) in aln.sequences.iter().zip(aln2.sequences.iter()) {
            assert_eq!(
                a.residues, b.residues,
                "residues differ for '{}' after round-trip",
                a.name
            );
        }
    }

    #[test]
    fn residues_include_gaps() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        // seq1 has "--" at positions 4-5 in the first block
        assert!(aln.sequences[0].residues.contains(&b'-'));
    }
}
