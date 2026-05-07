use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    sequence::{Sequence, SequenceType},
};

pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut order: Vec<String> = Vec::new();
    let mut residues: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    for line in text.lines() {
        let line = line.trim_end();

        // End of block marker
        if line == "//" {
            continue;
        }

        // Skip blank lines and all annotation/markup lines (#=... or leading #)
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Sequence line: "name  residues"
        let mut parts = line.splitn(2, |c: char| c.is_ascii_whitespace());
        let seq_name = match parts.next() {
            Some(n) if !n.is_empty() => n.to_string(),
            _ => continue,
        };
        let res_str = match parts.next() {
            Some(r) => r.trim(),
            None => continue,
        };
        if res_str.is_empty() {
            continue;
        }

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
            "No sequences found in Stockholm data".into(),
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
    if non_gap.iter().any(|&b| b.to_ascii_uppercase() == b'U') {
        SequenceType::Rna
    } else {
        SequenceType::Dna
    }
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = "# STOCKHOLM 1.0\n#=GF ID TEST\nseq1  ACGU--\nseq2  ACGUAA\n//\n";

    #[test]
    fn basic_parse() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.seq_count(), 2);
        assert_eq!(aln.sequences[0].residues.len(), 6);
    }

    #[test]
    fn skips_markup() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert!(aln.sequences.iter().all(|s| !s.name.starts_with('#')));
    }
}
