use helixview_core::Alignment;

/// Shannon entropy per alignment column.
/// Re-exported from helixview-core for convenience; extended analysis goes here.
pub fn column_entropy(aln: &Alignment) -> Vec<f64> {
    aln.column_entropy()
}

/// Information content = max_entropy − observed_entropy.
/// For nucleotides: max = ln(4) ≈ 1.386; for amino acids: max = ln(20) ≈ 3.0.
pub fn information_content(aln: &Alignment, is_protein: bool) -> Vec<f64> {
    let max_h = if is_protein {
        (20f64).ln()
    } else {
        (4f64).ln()
    };
    aln.column_entropy()
        .into_iter()
        .map(|h| (max_h - h).max(0.0))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use helixview_core::{sequence::Sequence, Alignment};

    fn aln(seqs: &[&[u8]]) -> Alignment {
        let mut a = Alignment::new("t");
        for (i, s) in seqs.iter().enumerate() {
            a.push(Sequence::new(format!("s{i}"), s.to_vec()));
        }
        a
    }

    #[test]
    fn constant_column_zero_entropy() {
        let a = aln(&[b"AAAA", b"AAAA", b"AAAA"]);
        for h in column_entropy(&a) {
            assert!(h.abs() < 1e-10, "expected 0, got {h}");
        }
    }

    #[test]
    fn uniform_four_base_column_max_entropy() {
        // One column with A, C, G, T (equal frequency) → H = ln(4).
        let a = aln(&[b"A", b"C", b"G", b"T"]);
        let h = column_entropy(&a)[0];
        assert!(
            (h - 4f64.ln()).abs() < 1e-10,
            "expected ln(4)≈{:.4}, got {h}",
            4f64.ln()
        );
    }

    #[test]
    fn two_equal_residues_give_ln2() {
        // 50% A, 50% C → H = ln(2).
        let a = aln(&[b"A", b"A", b"C", b"C"]);
        let h = column_entropy(&a)[0];
        assert!(
            (h - 2f64.ln()).abs() < 1e-10,
            "expected ln(2)≈{:.4}, got {h}",
            2f64.ln()
        );
    }

    #[test]
    fn length_matches_column_count() {
        let a = aln(&[b"ACGT", b"ACGT"]);
        assert_eq!(column_entropy(&a).len(), 4);
    }

    #[test]
    fn information_content_conserved_column_is_max() {
        // H = 0 for conserved column → IC = ln(4) for nucleotides.
        let a = aln(&[b"A", b"A", b"A"]);
        let ic = information_content(&a, false);
        assert!((ic[0] - 4f64.ln()).abs() < 1e-10, "got {}", ic[0]);
    }

    #[test]
    fn information_content_uniform_column_is_zero() {
        // H = ln(4) → IC = max_H − ln(4) = 0.
        let a = aln(&[b"A", b"C", b"G", b"T"]);
        let ic = information_content(&a, false);
        assert!(ic[0].abs() < 1e-10, "got {}", ic[0]);
    }
}
