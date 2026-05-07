//! Consensus sequence generation from a multiple sequence alignment.

use std::collections::BTreeSet;
use helixview_core::{Alignment, sequence::{is_gap, SequenceType}};

/// Strategy for choosing the consensus residue at each alignment column.
#[derive(Debug, Clone)]
pub enum ConsensusMethod {
    /// Most frequent non-gap residue (plurality vote).
    Plurality,
    /// IUPAC ambiguity codes: include all residues present at ≥ `threshold` frequency.
    Iupac { threshold: f64 },
}

/// Compute a consensus sequence for an alignment.
///
/// Returns a `Vec<u8>` of length `aln.col_count()`.
/// Columns where all residues are gaps return `b'-'`.
pub fn consensus(aln: &Alignment, method: ConsensusMethod) -> Vec<u8> {
    let cols = aln.col_count();
    if cols == 0 {
        return Vec::new();
    }

    // Decide whether this is a DNA or protein alignment.
    let is_dna = !aln.sequences.iter().any(|s| s.seq_type == SequenceType::Protein);

    (0..cols)
        .map(|col| {
            let mut counts = [0u32; 256];
            let mut total_non_gap = 0u32;

            for seq in &aln.sequences {
                if !seq.seq_type.is_sequence() {
                    continue;
                }
                if let Some(&b) = seq.residues.get(col) {
                    let b = b.to_ascii_uppercase();
                    if !is_gap(b) {
                        counts[b as usize] += 1;
                        total_non_gap += 1;
                    }
                }
            }

            if total_non_gap == 0 {
                return b'-';
            }

            match &method {
                ConsensusMethod::Plurality => {
                    counts
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, &v)| v)
                        .map(|(i, _)| i as u8)
                        .unwrap_or(b'-')
                }
                ConsensusMethod::Iupac { threshold } => {
                    let mut present: BTreeSet<u8> = BTreeSet::new();
                    for (byte, &cnt) in counts.iter().enumerate() {
                        if cnt > 0 {
                            let freq = cnt as f64 / total_non_gap as f64;
                            if freq >= *threshold {
                                present.insert(byte as u8);
                            }
                        }
                    }
                    if present.is_empty() {
                        // Shouldn't happen (total_non_gap > 0), but be safe.
                        return b'-';
                    }
                    to_iupac(&present, is_dna)
                }
            }
        })
        .collect()
}

/// Map a set of residues to an IUPAC ambiguity code.
///
/// For nucleotides, uses the standard IUPAC alphabet:
///   R=AG, Y=CT, S=GC, W=AT, K=GT, M=AC, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT.
/// For amino acids, returns `b'X'` for any ambiguous set.
/// A singleton set returns the residue itself.
pub fn to_iupac(residues: &BTreeSet<u8>, is_dna: bool) -> u8 {
    if residues.len() == 1 {
        return *residues.iter().next().unwrap();
    }
    if !is_dna {
        return b'X';
    }

    // Convert set to a sorted Vec for matching.
    let v: Vec<u8> = residues.iter().copied().collect();
    match v.as_slice() {
        [b'A', b'G']             => b'R',
        [b'C', b'T']             => b'Y',
        [b'C', b'G']             => b'S',
        [b'A', b'T']             => b'W',
        [b'G', b'T']             => b'K',
        [b'A', b'C']             => b'M',
        [b'C', b'G', b'T']       => b'B',
        [b'A', b'G', b'T']       => b'D',
        [b'A', b'C', b'T']       => b'H',
        [b'A', b'C', b'G']       => b'V',
        [b'A', b'C', b'G', b'T'] => b'N',
        _                        => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use helixview_core::{Alignment, sequence::Sequence};

    fn make_seq(name: &str, res: &[u8]) -> Sequence {
        Sequence::new(name, res.to_vec())
    }

    #[test]
    fn plurality_consensus_simple() {
        // Two sequences: ACGT and ACGT → consensus is ACGT.
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"ACGT"));
        aln.push(make_seq("b", b"ACGT"));
        let cons = consensus(&aln, ConsensusMethod::Plurality);
        assert_eq!(cons, b"ACGT");
    }

    #[test]
    fn plurality_consensus_majority_wins() {
        // Column 0: A A G → A wins.
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"A"));
        aln.push(make_seq("b", b"A"));
        aln.push(make_seq("c", b"G"));
        let cons = consensus(&aln, ConsensusMethod::Plurality);
        assert_eq!(cons, b"A");
    }

    #[test]
    fn plurality_consensus_all_gaps_returns_dash() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"-"));
        aln.push(make_seq("b", b"~"));
        let cons = consensus(&aln, ConsensusMethod::Plurality);
        assert_eq!(cons, b"-");
    }

    #[test]
    fn iupac_consensus_ag_gives_r() {
        // Two sequences at one column: A and G → R with threshold ≤ 0.5.
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"A"));
        aln.push(make_seq("b", b"G"));
        let cons = consensus(&aln, ConsensusMethod::Iupac { threshold: 0.5 });
        assert_eq!(cons, b"R");
    }

    #[test]
    fn to_iupac_singleton() {
        let mut set = BTreeSet::new();
        set.insert(b'A');
        assert_eq!(to_iupac(&set, true), b'A');
    }

    #[test]
    fn to_iupac_protein_ambiguous() {
        let mut set = BTreeSet::new();
        set.insert(b'A');
        set.insert(b'L');
        assert_eq!(to_iupac(&set, false), b'X');
    }
}
