use helixview_core::sequence::Sequence;

/// Kyte-Doolittle hydrophobicity scale for standard amino acids.
/// Returns 0.0 for unknown or gap characters.
pub fn kd_scale(aa: u8) -> f64 {
    match aa.to_ascii_uppercase() {
        b'I' =>  4.5,
        b'V' =>  4.2,
        b'L' =>  3.8,
        b'F' =>  2.8,
        b'C' =>  2.5,
        b'M' =>  1.9,
        b'A' =>  1.8,
        b'G' => -0.4,
        b'T' => -0.7,
        b'W' => -0.9,
        b'S' => -0.8,
        b'Y' => -1.3,
        b'P' => -1.6,
        b'H' => -3.2,
        b'E' => -3.5,
        b'Q' => -3.5,
        b'D' => -3.5,
        b'N' => -3.5,
        b'K' => -3.9,
        b'R' => -4.5,
        _    =>  0.0,
    }
}

/// Kyte-Doolittle hydrophobicity profile for a sequence.
/// Returns mean hydrophobicity for each window of size `window` (sliding by 1).
/// Output length = max(0, seq_len - window + 1).
/// Gaps are skipped (only true residues are used).
pub fn kyte_doolittle_profile(seq: &Sequence, window: usize) -> Vec<f64> {
    let residues: Vec<u8> = seq.residues.iter()
        .copied()
        .filter(|&b| !matches!(b, b'-' | b'~' | b'.'))
        .collect();
    if window == 0 || residues.len() < window { return vec![]; }
    residues.windows(window)
        .map(|w| w.iter().map(|&b| kd_scale(b)).sum::<f64>() / window as f64)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_length_correct() {
        let seq = Sequence::new("test", b"ACDEFGHIKLMNPQRSTVWY".to_vec());
        let p = kyte_doolittle_profile(&seq, 3);
        assert_eq!(p.len(), 18); // 20 - 3 + 1
    }

    #[test]
    fn all_ile_gives_max() {
        let seq = Sequence::new("test", b"IIIIII".to_vec());
        let p = kyte_doolittle_profile(&seq, 3);
        assert!(p.iter().all(|&v| (v - 4.5).abs() < 1e-9));
    }
}
