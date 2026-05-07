//! Conservation region search (§5.8).
//!
//! Scans the alignment for windows of consecutive columns where the average
//! per-column identity meets or exceeds a user-defined threshold.

use helixview_core::Alignment;

/// A conserved region found by `find_conserved_regions`.
#[derive(Debug, Clone)]
pub struct ConservedRegion {
    /// Start column (0-based, inclusive).
    pub start: usize,
    /// End column (0-based, inclusive).
    pub end:   usize,
    /// Average per-column identity across the region.
    pub avg_identity: f64,
    /// Plurality consensus sequence for the region.
    pub consensus: Vec<u8>,
}

impl ConservedRegion {
    pub fn width(&self) -> usize {
        self.end - self.start + 1
    }
}

/// Find all contiguous regions where every column has identity ≥ `threshold`.
///
/// `min_width` is the minimum number of consecutive conserved columns
/// needed to report a region (filter noise / single-column spikes).
///
/// Returns regions sorted by start position.
pub fn find_conserved_regions(
    aln:       &Alignment,
    threshold: f64,
    min_width: usize,
) -> Vec<ConservedRegion> {
    let ncols = aln.col_count();
    if ncols == 0 || aln.seq_count() == 0 {
        return vec![];
    }

    let mut regions: Vec<ConservedRegion> = Vec::new();
    let mut run_start: Option<usize> = None;
    let mut run_sum  = 0.0f64;

    for col in 0..ncols {
        let ident = aln.column_identity(col) as f64;
        if ident >= threshold {
            if run_start.is_none() { run_start = Some(col); run_sum = 0.0; }
            run_sum += ident;
        } else {
            if let Some(start) = run_start.take() {
                let end = col - 1;
                let width = end - start + 1;
                if width >= min_width {
                    regions.push(build_region(aln, start, end, run_sum / width as f64));
                }
                run_sum = 0.0;
            }
        }
    }
    // Close any open run at end
    if let Some(start) = run_start {
        let end = ncols - 1;
        let width = end - start + 1;
        if width >= min_width {
            regions.push(build_region(aln, start, end, run_sum / width as f64));
        }
    }

    regions
}

fn build_region(aln: &Alignment, start: usize, end: usize, avg: f64) -> ConservedRegion {
    let consensus = (start..=end)
        .map(|col| aln.column_consensus_residue(col))
        .collect();
    ConservedRegion { start, end, avg_identity: avg, consensus }
}

#[cfg(test)]
mod tests {
    use super::*;
    use helixview_core::{Alignment, sequence::Sequence};

    fn aln(seqs: &[&[u8]]) -> Alignment {
        let mut a = Alignment::new("t");
        for (i, s) in seqs.iter().enumerate() {
            a.push(Sequence::new(format!("s{i}"), s.to_vec()));
        }
        a
    }

    #[test]
    fn all_identical_one_full_region() {
        let a = aln(&[b"ACGT", b"ACGT", b"ACGT"]);
        let regions = find_conserved_regions(&a, 1.0, 1);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start, 0);
        assert_eq!(regions[0].end, 3);
    }

    #[test]
    fn broken_column_splits_regions() {
        // Cols 0-1 and col 3 match; col 2 breaks.
        let a = aln(&[b"AACG", b"AAGG"]);
        let regions = find_conserved_regions(&a, 1.0, 1);
        // Expect region [0,1] and region [3,3].
        assert!(regions.iter().any(|r| r.start == 0 && r.end == 1));
        assert!(regions.iter().any(|r| r.start == 3 && r.end == 3));
    }

    #[test]
    fn min_width_filters_short_runs() {
        // Cols 0 and 2 are conserved (1-col runs), col 1 is not. min_width=2 → nothing.
        let a = aln(&[b"ACG", b"AGG"]);
        let regions = find_conserved_regions(&a, 1.0, 2);
        assert_eq!(regions.len(), 0);
    }

    #[test]
    fn empty_alignment_returns_empty() {
        let a = Alignment::new("empty");
        assert!(find_conserved_regions(&a, 0.5, 1).is_empty());
    }

    #[test]
    fn consensus_sequence_in_region_is_correct() {
        let a = aln(&[b"AAAA", b"AAAA"]);
        let regions = find_conserved_regions(&a, 1.0, 1);
        assert_eq!(regions[0].consensus, b"AAAA");
    }
}
