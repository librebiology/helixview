/// Default scoring parameters.
pub const DEFAULT_MATCH:    i32 =  2;
pub const DEFAULT_MISMATCH: i32 = -1;
pub const DEFAULT_GAP:      i32 = -2;

/// Result of a pairwise alignment.
#[derive(Debug, Clone)]
pub struct PairwiseResult {
    pub aligned_a: Vec<u8>,
    pub aligned_b: Vec<u8>,
    pub score: i32,
    /// Fraction of aligned (non-double-gap) positions that are identical.
    pub identity: f64,
}

fn compute_identity(aligned_a: &[u8], aligned_b: &[u8]) -> f64 {
    let mut matches = 0usize;
    let mut valid   = 0usize;
    for (&a, &b) in aligned_a.iter().zip(aligned_b.iter()) {
        let a_gap = a == b'-';
        let b_gap = b == b'-';
        if a_gap && b_gap { continue; }
        valid += 1;
        if !a_gap && !b_gap && a.to_ascii_uppercase() == b.to_ascii_uppercase() {
            matches += 1;
        }
    }
    if valid == 0 { 0.0 } else { matches as f64 / valid as f64 }
}

/// Global alignment (Needleman-Wunsch) with linear gap penalty.
/// `match_score` > 0, `mismatch` <= 0, `gap_penalty` <= 0.
pub fn needleman_wunsch(
    a: &[u8],
    b: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_penalty: i32,
) -> PairwiseResult {
    let m = a.len();
    let n = b.len();

    // dp[i][j] = best score aligning a[..i] with b[..j]
    let mut dp = vec![vec![0i32; n + 1]; m + 1];

    // Initialise gap rows/columns
    for i in 0..=m { dp[i][0] = i as i32 * gap_penalty; }
    for j in 0..=n { dp[0][j] = j as i32 * gap_penalty; }

    for i in 1..=m {
        for j in 1..=n {
            let diag = dp[i-1][j-1] + if a[i-1].to_ascii_uppercase() == b[j-1].to_ascii_uppercase() {
                match_score
            } else {
                mismatch
            };
            let from_top  = dp[i-1][j] + gap_penalty;
            let from_left = dp[i][j-1] + gap_penalty;
            dp[i][j] = diag.max(from_top).max(from_left);
        }
    }

    let score = dp[m][n];

    // Traceback
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let mut i = m;
    let mut j = n;
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let diag_score = dp[i-1][j-1] + if a[i-1].to_ascii_uppercase() == b[j-1].to_ascii_uppercase() {
                match_score
            } else {
                mismatch
            };
            if dp[i][j] == diag_score {
                aligned_a.push(a[i-1]);
                aligned_b.push(b[j-1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i-1][j] + gap_penalty {
            aligned_a.push(a[i-1]);
            aligned_b.push(b'-');
            i -= 1;
        } else {
            aligned_a.push(b'-');
            aligned_b.push(b[j-1]);
            j -= 1;
        }
    }

    aligned_a.reverse();
    aligned_b.reverse();

    let identity = compute_identity(&aligned_a, &aligned_b);
    PairwiseResult { aligned_a, aligned_b, score, identity }
}

/// Local alignment (Smith-Waterman) with linear gap penalty.
pub fn smith_waterman(
    a: &[u8],
    b: &[u8],
    match_score: i32,
    mismatch: i32,
    gap_penalty: i32,
) -> PairwiseResult {
    if a.is_empty() || b.is_empty() {
        return PairwiseResult {
            aligned_a: Vec::new(),
            aligned_b: Vec::new(),
            score: 0,
            identity: 0.0,
        };
    }

    let m = a.len();
    let n = b.len();

    let mut dp = vec![vec![0i32; n + 1]; m + 1];
    let mut max_score = 0i32;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    for i in 1..=m {
        for j in 1..=n {
            let diag = dp[i-1][j-1] + if a[i-1].to_ascii_uppercase() == b[j-1].to_ascii_uppercase() {
                match_score
            } else {
                mismatch
            };
            let from_top  = dp[i-1][j] + gap_penalty;
            let from_left = dp[i][j-1] + gap_penalty;
            dp[i][j] = 0.max(diag).max(from_top).max(from_left);
            if dp[i][j] > max_score {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback from max cell until score reaches 0
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let mut i = max_i;
    let mut j = max_j;

    while i > 0 && j > 0 && dp[i][j] > 0 {
        if i > 0 && j > 0 {
            let diag_score = dp[i-1][j-1] + if a[i-1].to_ascii_uppercase() == b[j-1].to_ascii_uppercase() {
                match_score
            } else {
                mismatch
            };
            if dp[i][j] == diag_score {
                aligned_a.push(a[i-1]);
                aligned_b.push(b[j-1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && dp[i][j] == dp[i-1][j] + gap_penalty {
            aligned_a.push(a[i-1]);
            aligned_b.push(b'-');
            i -= 1;
        } else {
            aligned_a.push(b'-');
            aligned_b.push(b[j-1]);
            j -= 1;
        }
    }

    aligned_a.reverse();
    aligned_b.reverse();

    let identity = compute_identity(&aligned_a, &aligned_b);
    PairwiseResult { aligned_a, aligned_b, score: max_score, identity }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nw_identical_sequences() {
        let r = needleman_wunsch(b"ACGT", b"ACGT", 2, -1, -2);
        assert_eq!(r.score, 8);
        assert!((r.identity - 1.0).abs() < 1e-6);
    }

    #[test]
    fn nw_one_gap() {
        // "ACT" vs "ACGT" — optimal: A-CT aligned to ACGT
        let r = needleman_wunsch(b"ACT", b"ACGT", 2, -1, -2);
        assert!(r.score > 0);
        assert_eq!(r.aligned_a.len(), r.aligned_b.len());
    }

    #[test]
    fn sw_local_finds_match() {
        let r = smith_waterman(b"XXXACGTXXX", b"ACGT", 2, -1, -2);
        // Should find the ACGT local match
        assert!(r.score >= 8);
    }

    #[test]
    fn sw_empty_returns_empty() {
        let r = smith_waterman(b"", b"ACGT", 2, -1, -2);
        assert_eq!(r.score, 0);
    }
}
