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

/// Strip NEXUS block comments `[...]` from a line.
fn strip_comments(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut depth = 0usize;
    for ch in s.chars() {
        match ch {
            '[' => depth += 1,
            ']' => {
                depth = depth.saturating_sub(1);
            }
            _ => {
                if depth == 0 {
                    out.push(ch);
                }
            }
        }
    }
    out
}

pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    // Check for #NEXUS header (case-insensitive, skip BOM / leading whitespace)
    let trimmed = text.trim_start_matches('\u{feff}').trim_start();
    if !trimmed[..trimmed.len().min(6)].eq_ignore_ascii_case("#NEXUS") {
        return Err(FormatError::Parse(
            "Not a NEXUS file: missing #NEXUS header".into(),
        ));
    }

    let mut seq_type = SequenceType::NucleicAcid;
    let mut order: Vec<String> = Vec::new();
    let mut residues: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    // State machine
    enum State {
        Scanning,
        InData,
        InMatrix,
    }

    let mut state = State::Scanning;
    let lines: Vec<String> = text
        .lines()
        .map(|l| strip_comments(l).trim().to_string())
        .collect();

    let mut i = 0;
    while i < lines.len() {
        let line = &lines[i];
        let upper = line.to_ascii_uppercase();

        match state {
            State::Scanning => {
                if upper.starts_with("BEGIN DATA") || upper.starts_with("BEGIN CHARACTERS") {
                    state = State::InData;
                }
            }
            State::InData => {
                if upper.starts_with("END") {
                    state = State::Scanning;
                } else if upper.contains("DATATYPE=") {
                    let dt = extract_value(&upper, "DATATYPE");
                    seq_type = match dt.as_deref() {
                        Some("DNA") | Some("NUCLEOTIDE") => SequenceType::Dna,
                        Some("RNA") => SequenceType::Rna,
                        Some("PROTEIN") => SequenceType::Protein,
                        _ => SequenceType::NucleicAcid,
                    };
                } else if upper.trim() == "MATRIX" || upper.starts_with("MATRIX ") {
                    state = State::InMatrix;
                }
            }
            State::InMatrix => {
                // A semicolon alone (or line starting with ';') ends the matrix
                if line.trim() == ";" {
                    state = State::InData;
                } else if line.trim().is_empty() {
                    // blank line inside matrix — interleaved block separator, skip
                } else {
                    // sequence line: split on whitespace → (name, residues)
                    let mut parts = line.split_whitespace();
                    if let (Some(seq_name), Some(res_str)) = (parts.next(), parts.next()) {
                        let key = seq_name.to_string();
                        if !residues.contains_key(&key) {
                            order.push(key.clone());
                            residues.insert(key.clone(), Vec::new());
                        }
                        let entry = residues.get_mut(&key).unwrap();
                        for b in res_str.bytes() {
                            match b {
                                b'?' => entry.push(b'-'),
                                b if !b.is_ascii_whitespace() => entry.push(b),
                                _ => {}
                            }
                        }
                    }
                }
            }
        }

        i += 1;
    }

    if order.is_empty() {
        return Err(FormatError::Parse(
            "No sequences found in NEXUS MATRIX".into(),
        ));
    }

    let mut aln = Alignment::new(name);
    for seq_name in &order {
        let res = residues.remove(seq_name).unwrap_or_default();
        let mut seq = Sequence::new(seq_name.clone(), res);
        seq.seq_type = seq_type.clone();
        aln.push(seq);
    }

    Ok(aln)
}

/// Extract the value for `KEY=value` from an uppercase line.
fn extract_value(upper: &str, key: &str) -> Option<String> {
    let key_eq = format!("{}=", key);
    let pos = upper.find(&key_eq)?;
    let rest = &upper[pos + key_eq.len()..];
    // Value ends at whitespace or ';'
    let end = rest
        .find(|c: char| c.is_ascii_whitespace() || c == ';')
        .unwrap_or(rest.len());
    Some(rest[..end].to_string())
}

pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    let content = to_string(aln)?;
    std::fs::write(path, content).map_err(FormatError::Io)
}

pub fn to_string(aln: &Alignment) -> Result<String> {
    let n_tax = aln.seq_count();
    let n_char = aln
        .sequences
        .iter()
        .map(|s| s.residues.len())
        .max()
        .unwrap_or(0);

    let datatype = if n_tax > 0 {
        match aln.sequences[0].seq_type {
            SequenceType::Dna => "DNA",
            SequenceType::Rna => "RNA",
            SequenceType::Protein => "PROTEIN",
            _ => "DNA",
        }
    } else {
        "DNA"
    };

    let mut out = String::new();
    out.push_str("#NEXUS\n\n");
    out.push_str("BEGIN DATA;\n");
    out.push_str(&format!("  DIMENSIONS NTAX={} NCHAR={};\n", n_tax, n_char));
    out.push_str(&format!("  FORMAT DATATYPE={} GAP=-;\n", datatype));
    out.push_str("  MATRIX\n");

    // Pad names to the longest name length for alignment
    let max_name = aln
        .sequences
        .iter()
        .map(|s| s.name.len())
        .max()
        .unwrap_or(0);

    for seq in &aln.sequences {
        let res_str =
            std::str::from_utf8(&seq.residues).map_err(|e| FormatError::Parse(e.to_string()))?;
        out.push_str(&format!(
            "    {:<width$}  {}\n",
            seq.name,
            res_str,
            width = max_name
        ));
    }

    out.push_str("  ;\n");
    out.push_str("END;\n");

    Ok(out)
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = "#NEXUS\n\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=12;\n  FORMAT DATATYPE=DNA GAP=-;\n  MATRIX\n    seq1  ACGTACGTACGT\n    seq2  TTTTGGGGCCCC\n  ;\nEND;\n";

    #[test]
    fn basic_parse() {
        let a = super::parse_str(SAMPLE, "t").unwrap();
        assert_eq!(a.seq_count(), 2);
    }

    #[test]
    fn correct_residues() {
        let a = super::parse_str(SAMPLE, "t").unwrap();
        assert_eq!(a.sequences[0].residues.len(), 12);
    }

    #[test]
    fn names_parsed() {
        let a = super::parse_str(SAMPLE, "t").unwrap();
        assert_eq!(a.sequences[0].name, "seq1");
    }
}
