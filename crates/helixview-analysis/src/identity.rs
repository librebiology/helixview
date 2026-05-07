use rayon::prelude::*;
use helixview_core::Alignment;

/// Fraction of aligned, non-gap columns where both sequences have the same residue.
/// Positions where either sequence has a gap are excluded from denominator.
pub fn pairwise_identity(a: &[u8], b: &[u8]) -> f64 {
    let mut numer = 0u64;
    let mut denom = 0u64;
    for (&ra, &rb) in a.iter().zip(b.iter()) {
        let is_gap = |c: u8| matches!(c, b'-' | b'~' | b'.');
        if is_gap(ra) || is_gap(rb) {
            continue;
        }
        denom += 1;
        if ra.to_ascii_uppercase() == rb.to_ascii_uppercase() {
            numer += 1;
        }
    }
    if denom == 0 { 0.0 } else { numer as f64 / denom as f64 }
}

/// NxN symmetric matrix of pairwise identities (0.0–1.0).
/// Row i, col j = pairwise_identity(seq_i, seq_j).
/// Diagonal is 1.0.
pub fn identity_matrix(aln: &Alignment) -> Vec<Vec<f64>> {
    let n = aln.seq_count();
    let mut matrix = vec![vec![0.0f64; n]; n];
    for row in &mut matrix {
        for (j, val) in row.iter_mut().enumerate() {
            let _ = j;
            let _ = val;
        }
    }

    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect();

    let results: Vec<(usize, usize, f64)> = pairs
        .into_par_iter()
        .map(|(i, j)| {
            let id = pairwise_identity(
                &aln.sequences[i].residues,
                &aln.sequences[j].residues,
            );
            (i, j, id)
        })
        .collect();

    for i in 0..n {
        matrix[i][i] = 1.0;
    }
    for (i, j, id) in results {
        matrix[i][j] = id;
        matrix[j][i] = id;
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use helixview_core::sequence::Sequence;

    fn make_aln(seqs: &[(&str, &[u8])]) -> Alignment {
        let mut aln = Alignment::new("test");
        for (name, res) in seqs {
            aln.push(Sequence::new(*name, res.to_vec()));
        }
        aln
    }

    #[test]
    fn identical_sequences_give_one() {
        assert!((pairwise_identity(b"ACGT", b"ACGT") - 1.0).abs() < 1e-10);
    }

    #[test]
    fn completely_different_gives_zero() {
        assert!((pairwise_identity(b"AAAA", b"CCCC") - 0.0).abs() < 1e-10);
    }

    #[test]
    fn gaps_excluded_from_denominator() {
        // "A-C" vs "A-G": position 1 is gap in both → skip; positions 0 and 2 count
        // pos 0: A==A → match; pos 2: C vs G → no match → 1/2 = 0.5
        let id = pairwise_identity(b"A-C", b"A-G");
        assert!((id - 0.5).abs() < 1e-10);
    }

    #[test]
    fn identity_matrix_two_identical() {
        let aln = make_aln(&[("a", b"ACGT"), ("b", b"ACGT")]);
        let mat = identity_matrix(&aln);
        assert_eq!(mat.len(), 2);
        assert!((mat[0][0] - 1.0).abs() < 1e-10);
        assert!((mat[0][1] - 1.0).abs() < 1e-10);
        assert!((mat[1][0] - 1.0).abs() < 1e-10);
        assert!((mat[1][1] - 1.0).abs() < 1e-10);
    }
}
