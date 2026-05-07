//! Mutual Information (MI) analysis between alignment columns.
//!
//! MI(i,j) = Σ_a Σ_b P(a,b) × ln[ P(a,b) / (P(a) × P(b)) ]
//!
//! Only positions where neither sequence has a gap are counted.
//! Sequences with `SequenceType::Comment / SequenceMask / RnaStructureMask`
//! are excluded from the analysis.

use helixview_core::Alignment as SeqAlignment;
use rayon::prelude::*;
use std::collections::HashMap;

/// A single (column_i, column_j) MI result.
#[derive(Debug, Clone)]
pub struct MiPair {
    pub col_a: usize,
    pub col_b: usize,
    /// Raw mutual information in nats.
    pub mi: f64,
    /// Number of paired observations used.
    pub n_obs: usize,
    /// Number of distinct Watson-Crick pair types observed (0-6).
    /// Canonical: A-U, U-A, G-C, C-G; wobble: G-U, U-G.
    pub wc_types: u8,
    /// Fraction of observations that are Watson-Crick complementary.
    pub wc_frac: f64,
    /// Covariation score = wc_types × MI (BioEdit-style composite).
    pub cov_score: f64,
}

/// Full mutual information result over all column pairs.
#[derive(Debug, Clone)]
pub struct MiResult {
    /// Number of alignment columns.
    pub n_cols: usize,
    /// Upper-triangle MI values. Retrieve with `get(i, j)`.
    /// Stored as flat Vec of length n_cols*(n_cols-1)/2.
    pub values: Vec<f64>,
    /// Top pairs sorted by MI descending.
    pub top_pairs: Vec<MiPair>,
    /// Same pairs sorted by covariation score descending.
    pub top_covarying: Vec<MiPair>,
}

impl MiResult {
    /// Get the MI value for columns i and j (order-independent).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            return 0.0;
        }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        let idx = upper_idx(a, b, self.n_cols);
        self.values.get(idx).copied().unwrap_or(0.0)
    }

    /// Maximum MI value in the matrix (for colour scaling).
    pub fn max_mi(&self) -> f64 {
        self.values.iter().cloned().fold(0.0_f64, f64::max)
    }
}

fn upper_idx(i: usize, j: usize, n: usize) -> usize {
    // i < j guaranteed by caller
    i * n - i * (i + 1) / 2 + (j - i - 1)
}

/// Compute pairwise mutual information for all column pairs in `aln`.
///
/// * `min_obs` – minimum number of non-gap pairs required to compute MI for a
///   column pair (pairs below this threshold get MI = 0).
/// * `top_n`   – how many top pairs to include in `MiResult::top_pairs`.
pub fn mutual_information(aln: &SeqAlignment, min_obs: usize, top_n: usize) -> MiResult {
    let ncols = aln.col_count();
    if ncols < 2 {
        return MiResult {
            n_cols: ncols,
            values: vec![],
            top_pairs: vec![],
            top_covarying: vec![],
        };
    }

    // Only biological sequences.
    let seqs: Vec<&[u8]> = aln
        .sequences
        .iter()
        .filter(|s| s.seq_type.is_sequence())
        .map(|s| s.residues.as_slice())
        .collect();
    let nseqs = seqs.len();
    if nseqs == 0 {
        return MiResult {
            n_cols: ncols,
            values: vec![],
            top_pairs: vec![],
            top_covarying: vec![],
        };
    }

    // Precompute per-column residue vectors (uppercase, gap → None).
    let cols: Vec<Vec<Option<u8>>> = (0..ncols)
        .map(|c| {
            seqs.iter()
                .map(|s| {
                    let b = *s.get(c).unwrap_or(&b'-');
                    if is_gap(b) {
                        None
                    } else {
                        Some(b.to_ascii_uppercase())
                    }
                })
                .collect()
        })
        .collect();

    // Compute MI for every upper-triangle pair in parallel.
    let n_pairs = ncols * (ncols - 1) / 2;
    let mut values = vec![0.0_f64; n_pairs];

    // Build list of (i, j, flat_idx) so we can parallelise.
    let pairs: Vec<(usize, usize, usize)> = (0..ncols)
        .flat_map(|i| (i + 1..ncols).map(move |j| (i, j, upper_idx(i, j, ncols))))
        .collect();

    let mi_vals: Vec<(usize, f64)> = pairs
        .par_iter()
        .map(|&(i, j, flat_idx)| {
            let mi = compute_mi_pair(&cols[i], &cols[j], min_obs);
            (flat_idx, mi)
        })
        .collect();

    for (idx, mi) in mi_vals {
        values[idx] = mi;
    }

    // Collect top pairs with covariation annotation.
    let mut pair_list: Vec<MiPair> = pairs
        .iter()
        .map(|&(i, j, flat_idx)| {
            let mi = values[flat_idx];
            let (n_obs, wc_types, wc_frac) = covariation_stats(&cols[i], &cols[j]);
            let cov_score = (wc_types as f64) * mi;
            MiPair {
                col_a: i,
                col_b: j,
                mi,
                n_obs,
                wc_types,
                wc_frac,
                cov_score,
            }
        })
        .filter(|p| p.mi > 0.0)
        .collect();

    pair_list.sort_by(|a, b| b.mi.partial_cmp(&a.mi).unwrap());
    pair_list.truncate(top_n);

    let mut cov_list = pair_list.clone();
    cov_list.sort_by(|a, b| b.cov_score.partial_cmp(&a.cov_score).unwrap());

    MiResult {
        n_cols: ncols,
        values,
        top_pairs: pair_list,
        top_covarying: cov_list,
    }
}

fn compute_mi_pair(col_a: &[Option<u8>], col_b: &[Option<u8>], min_obs: usize) -> f64 {
    let mut joint: HashMap<(u8, u8), u32> = HashMap::new();
    let mut n = 0u32;
    for (a, b) in col_a.iter().zip(col_b.iter()) {
        if let (Some(ra), Some(rb)) = (a, b) {
            *joint.entry((*ra, *rb)).or_insert(0) += 1;
            n += 1;
        }
    }
    if (n as usize) < min_obs {
        return 0.0;
    }
    let nf = n as f64;

    // Marginals.
    let mut marg_a: HashMap<u8, u32> = HashMap::new();
    let mut marg_b: HashMap<u8, u32> = HashMap::new();
    for (&(ra, rb), &cnt) in &joint {
        *marg_a.entry(ra).or_insert(0) += cnt;
        *marg_b.entry(rb).or_insert(0) += cnt;
    }

    let mut mi = 0.0_f64;
    for (&(ra, rb), &cnt) in &joint {
        if cnt == 0 {
            continue;
        }
        let p_ab = cnt as f64 / nf;
        let p_a = *marg_a.get(&ra).unwrap_or(&0) as f64 / nf;
        let p_b = *marg_b.get(&rb).unwrap_or(&0) as f64 / nf;
        if p_a > 0.0 && p_b > 0.0 {
            mi += p_ab * (p_ab / (p_a * p_b)).ln();
        }
    }
    mi.max(0.0)
}

/// Compute covariation stats for a column pair.
/// Returns (n_obs, distinct_wc_types, wc_fraction).
fn covariation_stats(col_a: &[Option<u8>], col_b: &[Option<u8>]) -> (usize, u8, f64) {
    // WC pair types: AU=0, UA=1, GC=2, CG=3, GU=4, UG=5
    let mut wc_seen = [false; 6];
    let mut n_wc = 0usize;
    let mut n_obs = 0usize;

    for (a, b) in col_a.iter().zip(col_b.iter()) {
        if let (Some(ra), Some(rb)) = (a, b) {
            n_obs += 1;
            if let Some(kind) = wc_kind(*ra, *rb) {
                wc_seen[kind] = true;
                n_wc += 1;
            }
        }
    }

    let wc_types = wc_seen.iter().filter(|&&v| v).count() as u8;
    let wc_frac = if n_obs > 0 {
        n_wc as f64 / n_obs as f64
    } else {
        0.0
    };
    (n_obs, wc_types, wc_frac)
}

/// Return the Watson-Crick pair type index for bases (ra, rb), or None.
/// Canonical: A-U/U-A, G-C/C-G; wobble: G-U/U-G. DNA T treated as U.
#[inline]
fn wc_kind(ra: u8, rb: u8) -> Option<usize> {
    // Normalise T→U
    let a = if ra == b'T' { b'U' } else { ra };
    let b = if rb == b'T' { b'U' } else { rb };
    match (a, b) {
        (b'A', b'U') => Some(0),
        (b'U', b'A') => Some(1),
        (b'G', b'C') => Some(2),
        (b'C', b'G') => Some(3),
        (b'G', b'U') => Some(4),
        (b'U', b'G') => Some(5),
        _ => None,
    }
}

#[inline]
fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'.' | b'~' | b' ')
}

#[cfg(test)]
mod tests {
    use super::*;
    use helixview_core::sequence::Sequence;
    use helixview_core::Alignment;

    fn make_aln(seqs: &[&[u8]]) -> Alignment {
        let sequences = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| Sequence::new(format!("s{i}"), s.to_vec()))
            .collect();
        Alignment {
            sequences,
            ..Default::default()
        }
    }

    #[test]
    fn perfect_covariation() {
        // Col 0 and col 1 always co-vary perfectly: A→C, G→T
        let aln = make_aln(&[b"AC", b"AC", b"GT", b"GT"]);
        let res = mutual_information(&aln, 1, 10);
        let mi = res.get(0, 1);
        assert!(mi > 0.5, "expected high MI, got {mi}");
    }

    #[test]
    fn independent_columns() {
        // Col 0 = AAGG, col 1 = AGAG — independent
        let aln = make_aln(&[b"AA", b"AG", b"GA", b"GG"]);
        let res = mutual_information(&aln, 1, 10);
        let mi = res.get(0, 1);
        assert!(mi < 0.01, "expected near-zero MI, got {mi}");
    }

    #[test]
    fn wc_covariation_canonical_pairs() {
        // Col 0 → Col 1: A→U and G→C — two distinct canonical WC pair types.
        let aln = make_aln(&[b"AU", b"AU", b"GC", b"GC"]);
        let res = mutual_information(&aln, 1, 10);
        let pair = res
            .top_covarying
            .first()
            .expect("expected at least one pair");
        assert_eq!((pair.col_a, pair.col_b), (0, 1));
        assert!(
            pair.wc_types >= 2,
            "expected wc_types≥2, got {}",
            pair.wc_types
        );
        assert!(
            (pair.wc_frac - 1.0).abs() < 1e-6,
            "expected wc_frac=1.0, got {}",
            pair.wc_frac
        );
        assert!(pair.cov_score > 0.0);
    }

    #[test]
    fn non_wc_pairs_have_zero_wc_types() {
        // Perfectly correlated but AA/GG — not Watson-Crick pairs.
        let aln = make_aln(&[b"AA", b"AA", b"GG", b"GG"]);
        let res = mutual_information(&aln, 1, 10);
        if let Some(pair) = res.top_covarying.first() {
            assert_eq!(pair.wc_types, 0, "AA/GG should not score as WC");
        }
    }

    #[test]
    fn cov_score_equals_wc_types_times_mi() {
        let aln = make_aln(&[b"AU", b"AU", b"GC", b"GC"]);
        let res = mutual_information(&aln, 1, 10);
        for p in &res.top_covarying {
            let expected = p.wc_types as f64 * p.mi;
            assert!(
                (p.cov_score - expected).abs() < 1e-10,
                "cov_score {} ≠ wc_types×MI {}",
                p.cov_score,
                expected
            );
        }
    }
}
