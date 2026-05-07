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
    let mut aln = Alignment::new(name);

    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        if lines[i].starts_with('>') {
            let header = lines[i];
            // Parse type prefix and name from ">TYPE;name"
            let after_gt = &header[1..];
            let (type_str, seq_name) = if let Some(semi) = after_gt.find(';') {
                (after_gt[..semi].trim(), after_gt[semi + 1..].trim())
            } else {
                ("", after_gt.trim())
            };

            let seq_type = parse_type_prefix(type_str);

            i += 1;
            // second line = description (skip it)
            if i < lines.len() && !lines[i].starts_with('>') {
                i += 1; // skip description line
            }

            // Accumulate residues until '*' found (may be mid-line)
            let mut residues: Vec<u8> = Vec::new();
            let mut terminated = false;
            while i < lines.len() && !lines[i].starts_with('>') {
                let line = lines[i];
                for b in line.bytes() {
                    if b == b'*' {
                        terminated = true;
                        break;
                    }
                    if !b.is_ascii_whitespace() {
                        residues.push(b);
                    }
                }
                i += 1;
                if terminated {
                    break;
                }
            }

            if !seq_name.is_empty() {
                let mut seq = Sequence::new(seq_name.to_string(), residues);
                seq.seq_type = seq_type;
                aln.push(seq);
            }
        } else {
            i += 1;
        }
    }

    if aln.is_empty() {
        return Err(FormatError::Parse(
            "No sequences found in NBRF/PIR data".into(),
        ));
    }

    Ok(aln)
}

fn parse_type_prefix(prefix: &str) -> SequenceType {
    match prefix.to_ascii_uppercase().as_str() {
        "P1" | "P3" => SequenceType::Protein,
        "DL" | "DC" => SequenceType::Dna,
        "RL" | "RC" => SequenceType::Rna,
        "N3" | "N1" => SequenceType::NucleicAcid,
        _ => SequenceType::Unknown,
    }
}

fn type_to_prefix(seq_type: &SequenceType) -> &'static str {
    match seq_type {
        SequenceType::Protein => "P1",
        SequenceType::Dna => "DL",
        SequenceType::Rna => "RL",
        _ => "N3",
    }
}

pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    let content = to_string(aln)?;
    std::fs::write(path, content).map_err(FormatError::Io)
}

pub fn to_string(aln: &Alignment) -> Result<String> {
    let mut out = String::new();
    for seq in &aln.sequences {
        let prefix = type_to_prefix(&seq.seq_type);
        out.push_str(&format!(">{};{}\n", prefix, seq.name));
        out.push_str(&format!("{}\n", seq.name));
        let res_str =
            std::str::from_utf8(&seq.residues).map_err(|e| FormatError::Parse(e.to_string()))?;
        out.push_str(res_str);
        out.push_str("*\n\n");
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = ">P1;HUMAN_HBB\nHemoglobin subunit beta\nMVHLTPEEKSAVTALWGKVNVDEVGGEAL*\n\n>DL;TESTDNA\nTest DNA sequence\nATCGATCGATCG*\n";

    #[test]
    fn parses_two_records() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.seq_count(), 2);
    }

    #[test]
    fn protein_type_detected() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert!(matches!(
            aln.sequences[0].seq_type,
            helixview_core::SequenceType::Protein
        ));
    }

    #[test]
    fn dna_type_detected() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert!(matches!(
            aln.sequences[1].seq_type,
            helixview_core::SequenceType::Dna
        ));
    }
}
