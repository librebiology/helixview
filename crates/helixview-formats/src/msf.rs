//! MSF (GCG Multiple Sequence Format) parser and writer.
//!
//! Handles:
//!   - `!!AA_MULTIPLE_ALIGNMENT` (protein) and `!!NA_MULTIPLE_ALIGNMENT` (nucleotide) headers
//!   - Name block before `//`, sequence blocks after `//`
//!   - `.` gap characters converted to `-`
//!   - Sequence type determined from file header or info line, with residue heuristic fallback

use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    sequence::{Sequence, SequenceType},
};
use std::path::Path;

/// Parse an MSF byte slice into an `Alignment`.
pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

/// Parse an MSF string into an `Alignment`.
pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut aln = Alignment::new(name);

    // Determine seq type from the `!!` header line if present.
    let mut file_type: Option<SequenceType> = None;
    let mut info_type: Option<SequenceType> = None;

    // Collect ordered sequence names from the Name: block.
    let mut seq_names: Vec<String> = Vec::new();
    // Map from name → accumulated residues.
    let mut seq_map: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    let mut past_separator = false;

    for line in text.lines() {
        let trimmed = line.trim();

        if !past_separator {
            // Check for `!!` type header.
            if trimmed.starts_with("!!AA_MULTIPLE_ALIGNMENT") {
                file_type = Some(SequenceType::Protein);
            } else if trimmed.starts_with("!!NA_MULTIPLE_ALIGNMENT") {
                file_type = Some(SequenceType::Dna);
            }

            // Check info line for `Type: P/N/D`.
            if trimmed.contains("Type:") {
                if let Some(pos) = trimmed.find("Type:") {
                    let after = trimmed[pos + 5..].trim();
                    if after.starts_with('P') {
                        info_type = Some(SequenceType::Protein);
                    } else if after.starts_with('N') || after.starts_with('D') {
                        info_type = Some(SequenceType::Dna);
                    }
                }
            }

            // Collect name lines: line must start with " Name: " (leading space).
            if line.starts_with(" Name:") {
                if let Some(seq_name) = parse_name_line(line) {
                    if !seq_map.contains_key(&seq_name) {
                        seq_map.insert(seq_name.clone(), Vec::new());
                        seq_names.push(seq_name);
                    }
                }
            }

            // `//` separator ends the name block.
            if trimmed == "//" {
                past_separator = true;
            }
        } else {
            // After `//`: parse sequence lines.
            if trimmed.is_empty() {
                continue;
            }
            // First token is the sequence name; remainder is residues.
            let mut tokens = trimmed.splitn(2, |c: char| c.is_ascii_whitespace());
            if let Some(seq_name) = tokens.next() {
                if let Some(entry) = seq_map.get_mut(seq_name) {
                    if let Some(residue_str) = tokens.next() {
                        for b in residue_str.bytes() {
                            if b.is_ascii_alphabetic() {
                                entry.push(b);
                            } else if b == b'.' {
                                // Convert GCG gap character to standard gap.
                                entry.push(b'-');
                            } else if b == b'-' || b == b'~' {
                                entry.push(b);
                            }
                            // Spaces and digits in the residue portion are skipped.
                        }
                    }
                }
            }
        }
    }

    if seq_names.is_empty() {
        return Err(FormatError::Parse("No sequences found in MSF data".into()));
    }

    // Determine the sequence type to assign.
    let resolved_type = file_type.or(info_type).unwrap_or(SequenceType::Unknown);

    for seq_name in seq_names {
        let residues = seq_map.remove(&seq_name).unwrap_or_default();
        let mut seq = Sequence::new(seq_name, residues);
        seq.seq_type = if resolved_type == SequenceType::Unknown {
            crate::fasta::guess_type(&seq.residues)
        } else {
            resolved_type.clone()
        };
        aln.push(seq);
    }

    Ok(aln)
}

/// Parse a ` Name: SEQ1  Len: 48  Check: 5720  Weight: 1.00` line.
fn parse_name_line(line: &str) -> Option<String> {
    // Find "Name:" and extract the first token after it.
    let pos = line.find("Name:")?;
    let after = line[pos + 5..].trim_start();
    let name = after.split_whitespace().next()?.to_string();
    Some(name)
}

/// Write an alignment to a file in MSF format.
pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    let content = to_string(aln)?;
    std::fs::write(path, content).map_err(FormatError::Io)
}

/// Serialize an alignment to an MSF-format string.
///
/// Uses 50 residues per block, 10 residues per group, with padded names.
pub fn to_string(aln: &Alignment) -> Result<String> {
    if aln.is_empty() {
        return Err(FormatError::Parse("Alignment is empty".into()));
    }

    let mut out = String::new();

    // Determine type header.
    let type_char = match aln.sequences[0].seq_type {
        SequenceType::Protein => "AA",
        _ => "NA",
    };
    let mol_type_letter = match aln.sequences[0].seq_type {
        SequenceType::Protein => 'P',
        _ => 'N',
    };

    let seq_len = aln
        .sequences
        .iter()
        .map(|s| s.residues.len())
        .max()
        .unwrap_or(0);

    out.push_str(&format!("!!{}_MULTIPLE_ALIGNMENT 1.0\n\n", type_char));
    out.push_str(&format!(
        " alignment.msf  MSF: {}  Type: {}  Check: 0 ..\n\n",
        seq_len, mol_type_letter
    ));

    // Name block.
    let name_width = aln
        .sequences
        .iter()
        .map(|s| s.name.len())
        .max()
        .unwrap_or(4);
    for seq in &aln.sequences {
        out.push_str(&format!(
            " Name: {:<width$}  Len: {}  Check: 0  Weight: 1.00\n",
            seq.name,
            seq.residues.len(),
            width = name_width
        ));
    }

    out.push_str("\n//\n\n");

    // Sequence blocks: 50 residues per block, 10 per group.
    const BLOCK: usize = 50;
    const GROUP: usize = 10;

    let mut offset = 0usize;
    while offset < seq_len {
        let block_end = (offset + BLOCK).min(seq_len);

        for seq in &aln.sequences {
            let slice = if offset < seq.residues.len() {
                let end = block_end.min(seq.residues.len());
                &seq.residues[offset..end]
            } else {
                b""
            };

            out.push_str(&format!("{:<width$}  ", seq.name, width = name_width));

            let mut col = 0usize;
            while col < slice.len() {
                if col > 0 {
                    out.push(' ');
                }
                let group_end = (col + GROUP).min(slice.len());
                out.push_str(
                    std::str::from_utf8(&slice[col..group_end])
                        .map_err(|e| FormatError::Parse(e.to_string()))?,
                );
                col += GROUP;
            }
            out.push('\n');
        }
        out.push('\n');
        offset += BLOCK;
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = "!!AA_MULTIPLE_ALIGNMENT 1.0\n\n test.msf  MSF: 12  Type: P  Check: 0 ..\n\n Name: SEQ1  Len: 12  Check: 0  Weight: 1.00\n Name: SEQ2  Len: 12  Check: 0  Weight: 1.00\n\n//\n\nSEQ1      ACDEFGHIKL MN\nSEQ2      MNLKIHGFED CA\n\n";

    #[test]
    fn parses_two_seqs() {
        let a = super::parse_str(SAMPLE, "t").unwrap();
        assert_eq!(a.seq_count(), 2);
    }

    #[test]
    fn correct_length() {
        let a = super::parse_str(SAMPLE, "t").unwrap();
        assert_eq!(a.sequences[0].residues.len(), 12);
    }

    #[test]
    fn dot_converted_to_gap() {
        let raw = "!!AA_MULTIPLE_ALIGNMENT 1.0\n\n t.msf  MSF: 3  Type: P  Check: 0 ..\n\n Name: S1  Len: 3  Check: 0  Weight: 1.00\n\n//\n\nS1  A.C\n\n";
        let a = super::parse_str(raw, "t").unwrap();
        assert_eq!(a.sequences[0].residues[1], b'-');
    }
}
