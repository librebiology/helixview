//! FASTA format parser and writer.
//!
//! Handles:
//!   - Standard FASTA (>header\nsequence lines)
//!   - Multi-sequence files (alignments)
//!   - Gaps in residues (-, ~, .) preserved as-is
//!   - Sequence type guessed from residue composition

use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    sequence::{Sequence, SequenceType},
};
use std::io::Write;
use std::path::Path;

/// Parse a FASTA byte slice into an `Alignment`.
pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

/// Parse a FASTA string into an `Alignment`.
pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut aln = Alignment::new(name);
    let mut current_name: Option<String> = None;
    let mut current_residues: Vec<u8> = Vec::new();

    for line in text.lines() {
        let line = line.trim_end();

        if let Some(rest) = line.strip_prefix('>') {
            // Flush previous entry.
            if let Some(seq_name) = current_name.take() {
                aln.push(build_sequence(
                    seq_name,
                    std::mem::take(&mut current_residues),
                ));
            }
            // Start new entry — everything after '>' is the name/description.
            current_name = Some(rest.trim().to_string());
        } else if current_name.is_some() {
            // Accumulate residues, skip whitespace.
            for b in line.bytes() {
                if !b.is_ascii_whitespace() {
                    current_residues.push(b);
                }
            }
        }
        // Lines before the first '>' are silently ignored.
    }

    // Flush the last entry.
    if let Some(seq_name) = current_name {
        aln.push(build_sequence(seq_name, current_residues));
    }

    if aln.is_empty() {
        return Err(FormatError::Parse(
            "No sequences found in FASTA data".into(),
        ));
    }

    Ok(aln)
}

fn build_sequence(name: String, residues: Vec<u8>) -> Sequence {
    let mut seq = Sequence::new(name, residues);
    seq.seq_type = guess_type(&seq.residues);
    seq
}

/// Heuristic: if >85% of non-gap characters are in {A,C,G,T,U,N} → nucleic acid.
pub fn guess_type(residues: &[u8]) -> SequenceType {
    let non_gap: Vec<u8> = residues
        .iter()
        .copied()
        .filter(|&b| !matches!(b, b'-' | b'~' | b'.'))
        .collect();

    if non_gap.is_empty() {
        return SequenceType::Unknown;
    }

    let na_chars: usize = non_gap
        .iter()
        .filter(|&&b| {
            matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C'
                    | b'G'
                    | b'T'
                    | b'U'
                    | b'N'
                    | b'R'
                    | b'Y'
                    | b'S'
                    | b'W'
                    | b'K'
                    | b'M'
                    | b'B'
                    | b'D'
                    | b'H'
                    | b'V'
            )
        })
        .count();

    let fraction = na_chars as f64 / non_gap.len() as f64;
    if fraction >= 0.85 {
        // Distinguish DNA vs RNA by presence of U.
        let has_u = non_gap.iter().any(|&b| b.eq_ignore_ascii_case(&b'U'));
        let has_t = non_gap.iter().any(|&b| b.eq_ignore_ascii_case(&b'T'));
        if has_u && !has_t {
            SequenceType::Rna
        } else {
            SequenceType::Dna
        }
    } else {
        SequenceType::Protein
    }
}

/// Write an alignment to a file in FASTA format.
pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)?;
    write_writer(aln, &mut f)
}

/// Write an alignment to any `Write` implementor.
pub fn write_writer(aln: &Alignment, w: &mut dyn Write) -> Result<()> {
    for seq in &aln.sequences {
        writeln!(w, ">{}", seq.name).map_err(FormatError::Io)?;
        // Write residues 60 per line (standard FASTA convention).
        for chunk in seq.residues.chunks(60) {
            w.write_all(chunk).map_err(FormatError::Io)?;
            writeln!(w).map_err(FormatError::Io)?;
        }
    }
    Ok(())
}

/// Serialize to a `String` (useful for clipboard export).
pub fn to_string(aln: &Alignment) -> Result<String> {
    let mut buf = Vec::new();
    write_writer(aln, &mut buf)?;
    String::from_utf8(buf).map_err(|e| FormatError::Parse(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
>Sequence_1 description here
ATCGATCGATCG
ATCG
>Sequence_2
GGCCTTAA
>Sequence_3_protein
MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNG
";

    #[test]
    fn parse_basic() {
        let aln = parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.seq_count(), 3);
        assert_eq!(aln.sequences[0].name, "Sequence_1 description here");
        assert_eq!(aln.sequences[0].true_len(), 16); // 12 + 4
        assert_eq!(aln.sequences[1].true_len(), 8);
    }

    #[test]
    fn type_guessing() {
        let dna = guess_type(b"ATCGATCG");
        assert_eq!(dna, SequenceType::Dna);

        let rna = guess_type(b"AUCGAUCG");
        assert_eq!(rna, SequenceType::Rna);

        let prot = guess_type(b"MAEGEITTFTAL");
        assert_eq!(prot, SequenceType::Protein);
    }

    #[test]
    fn roundtrip() {
        let aln = parse_str(SAMPLE, "test").unwrap();
        let out = to_string(&aln).unwrap();
        let aln2 = parse_str(&out, "test2").unwrap();
        assert_eq!(aln.seq_count(), aln2.seq_count());
        for (a, b) in aln.sequences.iter().zip(aln2.sequences.iter()) {
            assert_eq!(a.residues, b.residues);
        }
    }

    #[test]
    fn gaps_preserved() {
        let input = ">seq\nAT-CG~TA\n";
        let aln = parse_str(input, "test").unwrap();
        assert_eq!(&aln.sequences[0].residues, b"AT-CG~TA");
    }
}
