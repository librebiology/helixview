use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    sequence::{Sequence, SequenceType},
};
use std::path::Path;

pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut lines = text.lines();

    // Parse header line: two integers (n_seqs and n_cols)
    let header = loop {
        match lines.next() {
            None => return Err(FormatError::Parse("Empty Phylip file".into())),
            Some(l) if l.trim().is_empty() => continue,
            Some(l) => break l,
        }
    };

    let mut header_parts = header.split_whitespace();
    let n_seqs: usize = header_parts
        .next()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| FormatError::Parse("Invalid Phylip header: expected two integers".into()))?;
    let _n_cols: usize = header_parts
        .next()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| FormatError::Parse("Invalid Phylip header: expected two integers".into()))?;

    if n_seqs == 0 {
        return Err(FormatError::Parse("Phylip header says 0 sequences".into()));
    }

    // Collect all remaining lines (non-blank for interleaved handling later)
    let remaining: Vec<&str> = lines.collect();

    // Parse first n_seqs lines: name (first 10 chars) + residues (rest)
    let mut names: Vec<String> = Vec::with_capacity(n_seqs);
    let mut residues: Vec<Vec<u8>> = Vec::with_capacity(n_seqs);

    // We need to pick exactly n_seqs non-blank lines from remaining for the first block
    let mut line_idx = 0;
    let mut parsed = 0;
    while parsed < n_seqs {
        if line_idx >= remaining.len() {
            return Err(FormatError::Parse(format!(
                "Expected {n_seqs} sequences but file ended early"
            )));
        }
        let line = remaining[line_idx];
        line_idx += 1;
        if line.trim().is_empty() {
            continue;
        }

        // Name: first 10 characters, trimmed
        let name_part = if line.len() >= 10 { &line[..10] } else { line };
        let seq_name = name_part.trim().to_string();

        // Residues: rest of line after the first 10 chars
        let res_part = if line.len() > 10 { &line[10..] } else { "" };
        let res: Vec<u8> = res_part
            .bytes()
            .filter(|b| !b.is_ascii_whitespace() && !b.is_ascii_digit())
            .collect();

        names.push(seq_name);
        residues.push(res);
        parsed += 1;
    }

    // All subsequent non-blank lines append residues cycling through sequences
    let mut seq_idx = 0usize;
    while line_idx < remaining.len() {
        let line = remaining[line_idx];
        line_idx += 1;
        if line.trim().is_empty() {
            continue;
        }
        let res: Vec<u8> = line
            .bytes()
            .filter(|b| !b.is_ascii_whitespace() && !b.is_ascii_digit())
            .collect();
        if !res.is_empty() {
            residues[seq_idx % n_seqs].extend_from_slice(&res);
            seq_idx += 1;
        }
    }

    let mut aln = Alignment::new(name);
    for (seq_name, res) in names.into_iter().zip(residues.into_iter()) {
        let seq_type = infer_type(&res);
        let mut seq = Sequence::new(seq_name, res);
        seq.seq_type = seq_type;
        aln.push(seq);
    }

    Ok(aln)
}

fn infer_type(residues: &[u8]) -> SequenceType {
    let non_gap: Vec<u8> = residues
        .iter()
        .copied()
        .filter(|&b| !matches!(b, b'-' | b'~' | b'.'))
        .collect();

    if non_gap.is_empty() {
        return SequenceType::Unknown;
    }

    let dna_rna_set: &[u8] = b"ACGTURYNacgturyn";
    let is_protein = non_gap.iter().any(|b| !dna_rna_set.contains(b));

    if is_protein {
        return SequenceType::Protein;
    }

    let has_u = non_gap.iter().any(|&b| b.to_ascii_uppercase() == b'U');
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
    if aln.is_empty() {
        return Ok(String::new());
    }

    let n_seqs = aln.seq_count();
    let n_cols = aln
        .sequences
        .iter()
        .map(|s| s.residues.len())
        .max()
        .unwrap_or(0);

    let mut out = String::new();
    out.push_str(&format!(" {} {}\n", n_seqs, n_cols));

    let block_size = 60usize;
    let mut offset = 0usize;
    while offset < n_cols {
        let end = (offset + block_size).min(n_cols);
        for seq in &aln.sequences {
            if offset == 0 {
                // First block: include name padded to 10 chars
                let name_padded = format!("{:<10}", &seq.name[..seq.name.len().min(10)]);
                out.push_str(&name_padded);
            }
            let chunk = if offset < seq.residues.len() {
                let e = end.min(seq.residues.len());
                std::str::from_utf8(&seq.residues[offset..e])
                    .map_err(|e| FormatError::Parse(e.to_string()))?
            } else {
                ""
            };
            out.push_str(chunk);
            out.push('\n');
        }
        if end < n_cols {
            out.push('\n');
        }
        offset += block_size;
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    const SEQUENTIAL: &str =
        " 2 20\nSeq1      ACGTACGTACGTACGTACGT\nSeq2      TTTTCCCCGGGGAAAATTTT\n";
    const INTERLEAVED: &str =
        " 2 20\nSeq1      ACGTACGT\nSeq2      TTTTCCCC\n\nACGTACGT\nGGGGAAAA\n\nACGT\nTTTT\n";

    #[test]
    fn sequential_two_seqs() {
        let a = super::parse_str(SEQUENTIAL, "t").unwrap();
        assert_eq!(a.seq_count(), 2);
    }

    #[test]
    fn sequential_correct_len() {
        let a = super::parse_str(SEQUENTIAL, "t").unwrap();
        assert_eq!(a.sequences[0].residues.len(), 20);
    }

    #[test]
    fn interleaved_two_seqs() {
        let a = super::parse_str(INTERLEAVED, "t").unwrap();
        assert_eq!(a.seq_count(), 2);
    }

    #[test]
    fn interleaved_correct_len() {
        let a = super::parse_str(INTERLEAVED, "t").unwrap();
        assert_eq!(a.sequences[0].residues.len(), 20);
    }
}
