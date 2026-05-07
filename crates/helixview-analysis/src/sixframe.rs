//! Six-frame translation and ORF finding for nucleotide sequences.

/// A single ORF found in a sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Orf {
    /// Reading frame: +1, +2, +3, -1, -2, -3.
    pub frame: i8,
    /// 0-based nucleotide position on the original (forward) strand (inclusive).
    pub start: usize,
    /// Exclusive end position on the original (forward) strand.
    pub end: usize,
    /// Translated amino acid bytes (no stop codon included).
    pub protein: Vec<u8>,
}

/// Translate a single nucleotide codon (3 bytes, uppercase) to an amino acid byte.
/// Stop codons return `b'*'`. Gaps or ambiguous bases return `b'X'`.
pub fn translate_codon(codon: &[u8; 3]) -> u8 {
    match codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" => b'L',
        b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" => b'R',
        b"AGA" | b"AGG" => b'R',
        b"AGT" | b"AGC" => b'S',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

/// Reverse-complement a nucleotide byte slice (gaps stripped beforehand is fine,
/// but gaps are preserved here as-is if present).
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Strip gap characters (`b'-'`, `b'~'`, `b'.'`) from a sequence.
fn strip_gaps(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .copied()
        .filter(|&b| !matches!(b, b'-' | b'~' | b'.'))
        .collect()
}

/// Translate a nucleotide sequence in one reading frame.
///
/// `frame` is 0-based:
/// - 0, 1, 2 → forward strand offsets 0, 1, 2
/// - 3, 4, 5 → reverse-complement strand offsets 0, 1, 2
///
/// Gap characters are stripped before translation.
/// Stop codons produce `b'*'`; gaps/ambiguous nucleotides produce `b'X'`.
pub fn translate_frame(seq: &[u8], frame: usize) -> Vec<u8> {
    let stripped = strip_gaps(seq);
    let working: Vec<u8> = if frame < 3 {
        stripped[frame..].to_vec()
    } else {
        let rc = reverse_complement(&stripped);
        let offset = frame - 3;
        rc[offset..].to_vec()
    };

    let mut protein = Vec::with_capacity(working.len() / 3);
    let mut i = 0;
    while i + 3 <= working.len() {
        let codon: &[u8; 3] = working[i..i + 3].try_into().unwrap();
        let codon_upper: [u8; 3] = [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ];
        protein.push(translate_codon(&codon_upper));
        i += 3;
    }
    protein
}

/// Translate a nucleotide sequence in all 6 reading frames.
///
/// Returns `[(frame_label, protein)]` where `frame_label` is +1/+2/+3/-1/-2/-3.
pub fn six_frame_all(seq: &[u8]) -> [(i8, Vec<u8>); 6] {
    let labels: [i8; 6] = [1, 2, 3, -1, -2, -3];
    core::array::from_fn(|i| (labels[i], translate_frame(seq, i)))
}

/// Find all ORFs (ATG … stop) with at least `min_aa_len` amino acids in all 6 frames.
///
/// `start` and `end` are 0-based positions on the **original forward strand** sequence
/// (after gap stripping). For reverse-complement frames the mapping is:
///   `fwd_start = len - rc_end`, `fwd_end = len - rc_start`.
pub fn find_orfs(seq: &[u8], min_aa_len: usize) -> Vec<Orf> {
    let stripped = strip_gaps(seq);
    let len = stripped.len();
    let mut orfs = Vec::new();

    for frame in 0usize..6 {
        // Build the working strand for this frame.
        let (strand_seq, is_rc) = if frame < 3 {
            (stripped.clone(), false)
        } else {
            (reverse_complement(&stripped), true)
        };
        let offset = frame % 3;

        // Translate the whole frame first.
        let protein = translate_frame(seq, frame);

        // Scan for ATG → stop stretches.
        let mut i = 0usize;
        while i < protein.len() {
            if protein[i] == b'M' {
                // Found a start; look for the next stop.
                let start_aa = i;
                let mut j = i + 1;
                while j < protein.len() && protein[j] != b'*' {
                    j += 1;
                }
                let aa_len = j - start_aa;
                if aa_len >= min_aa_len {
                    let prot = protein[start_aa..j].to_vec();

                    // Nucleotide positions in the working (possibly RC) strand.
                    // Include the stop codon in the end position.
                    let nt_start_in_strand = offset + start_aa * 3;
                    let stop_included = if j < protein.len() { 3 } else { 0 };
                    let nt_end_in_strand = offset + j * 3 + stop_included;

                    let (fwd_start, fwd_end) = if is_rc {
                        // Map RC positions back to forward strand.
                        let rc_len = strand_seq.len();
                        let fwd_end = rc_len - nt_start_in_strand;
                        let fwd_start = rc_len - nt_end_in_strand;
                        (fwd_start, fwd_end)
                    } else {
                        (nt_start_in_strand, nt_end_in_strand)
                    };

                    // Clamp to actual sequence length.
                    let fwd_start = fwd_start.min(len);
                    let fwd_end = fwd_end.min(len);

                    let frame_label: i8 = if frame < 3 {
                        (frame as i8) + 1
                    } else {
                        -((frame as i8) - 2)
                    };

                    orfs.push(Orf {
                        frame: frame_label,
                        start: fwd_start,
                        end: fwd_end,
                        protein: prot,
                    });
                }
                // Continue scanning after the stop (or end of sequence).
                i = j + 1;
            } else {
                i += 1;
            }
        }
    }

    orfs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translate_atg_to_met() {
        assert_eq!(translate_codon(b"ATG"), b'M');
    }

    #[test]
    fn translate_stop_codons() {
        assert_eq!(translate_codon(b"TAA"), b'*');
        assert_eq!(translate_codon(b"TAG"), b'*');
        assert_eq!(translate_codon(b"TGA"), b'*');
    }

    #[test]
    fn translate_frame_forward() {
        // ATGAAATAA → M K *
        let seq = b"ATGAAATAA";
        let prot = translate_frame(seq, 0);
        assert_eq!(prot, b"MK*");
    }

    #[test]
    fn translate_frame_strips_gaps() {
        // ATG-AAA-TAA with gaps should give same result as without
        let seq_with_gaps = b"ATG-AAA-TAA";
        let seq_clean = b"ATGAAATAA";
        assert_eq!(
            translate_frame(seq_with_gaps, 0),
            translate_frame(seq_clean, 0)
        );
    }

    #[test]
    fn find_orfs_finds_known_orf() {
        // Frame +1: ATG AAA TAA  → MK (2 aa)
        let seq = b"ATGAAATAA";
        let orfs = find_orfs(seq, 1);
        assert!(!orfs.is_empty(), "should find at least one ORF");
        let orf = orfs
            .iter()
            .find(|o| o.frame == 1)
            .expect("frame +1 ORF missing");
        assert_eq!(orf.protein, b"MK");
        assert_eq!(orf.start, 0);
        assert_eq!(orf.end, 9);
    }

    #[test]
    fn find_orfs_min_len_filter() {
        // Only 2 aa long; requiring 3 should return nothing in frame +1.
        let seq = b"ATGAAATAA";
        let orfs = find_orfs(seq, 3);
        assert!(orfs.iter().all(|o| o.protein.len() >= 3));
    }
}
