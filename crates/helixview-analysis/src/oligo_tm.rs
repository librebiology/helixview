//! Oligonucleotide melting temperature (Tm) calculation.
//!
//! Uses the SantaLucia (1998) nearest-neighbor thermodynamic parameters for DNA.
//! Assumptions: 50 mM Na⁺, oligo concentration 250 nM (typical PCR conditions).

/// SantaLucia 1998 nearest-neighbor parameters: (ΔH kcal/mol, ΔS cal/mol·K).
/// Keys are two-letter DNA dinucleotides in 5'→3' direction (uppercase).
const NN_PARAMS: &[(&[u8; 2], f64, f64)] = &[
    (b"AA", -7.9,  -22.2),
    (b"AT", -7.2,  -20.4),
    (b"TA", -7.2,  -21.3),
    (b"CA", -8.5,  -22.7),
    (b"GT", -8.4,  -22.4),
    (b"CT", -7.8,  -21.0),
    (b"GA", -8.2,  -22.2),
    (b"CG", -10.6, -27.2),
    (b"GC", -9.8,  -24.4),
    (b"GG", -8.0,  -19.9),
    // Complements (reverse strand reads 5'→3' so GG complement = CC)
    (b"TT", -7.9,  -22.2), // complement of AA
    (b"AC", -7.8,  -21.0), // complement of GT
    (b"TC", -8.2,  -22.2), // complement of GA
    (b"CC", -8.0,  -19.9), // complement of GG
    (b"AG", -7.8,  -21.0), // complement of CT
    (b"TG", -8.5,  -22.7), // complement of CA
    (b"GT", -8.4,  -22.4),
    (b"TT", -7.9,  -22.2),
];

/// Initiation parameters (SantaLucia 1998, Table 2).
/// Terminal GC: ΔH = 0.1 kcal/mol, ΔS = −2.8 cal/mol·K
/// Terminal AT: ΔH = 2.3 kcal/mol, ΔS = 4.1 cal/mol·K
const INIT_GC: (f64, f64) = (0.1,  -2.8);
const INIT_AT: (f64, f64) = (2.3,   4.1);

/// R in cal/(mol·K)
const R: f64 = 1.987;

/// Result of a Tm calculation.
#[derive(Debug, Clone)]
pub struct TmResult {
    /// Melting temperature in °C (assuming 50 mM Na⁺, 250 nM oligo).
    pub tm_celsius: f64,
    /// Total ΔH (kcal/mol).
    pub delta_h: f64,
    /// Total ΔS (cal/mol·K).
    pub delta_s: f64,
    /// GC content fraction.
    pub gc_fraction: f64,
    /// Length of the sequence used (gap-stripped).
    pub length: usize,
}

/// Calculate Tm for a raw DNA sequence (gaps stripped, case-insensitive).
///
/// Uses:
///   Tm = ΔH / (ΔS + R × ln(CT/4)) − 273.15  (in Kelvin → Celsius)
///
/// where CT = 250 × 10⁻⁹ mol/L (250 nM total strand concentration).
pub fn calculate_tm(seq: &[u8]) -> Option<TmResult> {
    // Strip gaps, uppercase.
    let raw: Vec<u8> = seq.iter()
        .filter(|&&b| !matches!(b, b'-' | b'~' | b'.'))
        .map(|&b| b.to_ascii_uppercase())
        .collect();

    if raw.len() < 2 { return None; }

    let mut dh: f64 = 0.0;
    let mut ds: f64 = 0.0;

    // Sum nearest-neighbor parameters
    for i in 0..raw.len() - 1 {
        let pair = [raw[i], raw[i + 1]];
        if let Some(&(_, h, s)) = NN_PARAMS.iter().find(|(k, _, _)| k[0] == pair[0] && k[1] == pair[1]) {
            dh += h;
            ds += s;
        }
        // Unknown pairs: use average (−8.0, −22.0)
        else {
            dh += -8.0;
            ds += -22.0;
        }
    }

    // Initiation correction
    let add_init = |b: u8, dh: &mut f64, ds: &mut f64| {
        let (h, s) = if matches!(b, b'G' | b'C') { INIT_GC } else { INIT_AT };
        *dh += h; *ds += s;
    };
    add_init(*raw.first().unwrap(), &mut dh, &mut ds);
    add_init(*raw.last().unwrap(),  &mut dh, &mut ds);

    // Convert units: ΔH kcal/mol → cal/mol
    let dh_cal = dh * 1000.0;

    // CT = 250 nM = 250e-9 mol/L
    let ct: f64 = 250e-9;

    // Tm (K) = ΔH / (ΔS + R × ln(CT/4))
    let tm_k = dh_cal / (ds + R * (ct / 4.0_f64).ln());
    let tm_c = tm_k - 273.15;

    let gc = raw.iter().filter(|&&b| matches!(b, b'G' | b'C')).count();
    let gc_fraction = gc as f64 / raw.len() as f64;

    Some(TmResult {
        tm_celsius: tm_c,
        delta_h:    dh,
        delta_s:    ds,
        gc_fraction,
        length:     raw.len(),
    })
}

/// Quick rule-of-thumb Tm for short oligos (< 14 nt): Wallace rule.
/// Tm = 2(A+T) + 4(G+C)
pub fn wallace_tm(seq: &[u8]) -> f64 {
    let raw: Vec<u8> = seq.iter()
        .filter(|&&b| !matches!(b, b'-' | b'~' | b'.'))
        .map(|&b| b.to_ascii_uppercase())
        .collect();
    let at = raw.iter().filter(|&&b| matches!(b, b'A' | b'T' | b'U')).count() as f64;
    let gc = raw.iter().filter(|&&b| matches!(b, b'G' | b'C')).count() as f64;
    2.0 * at + 4.0 * gc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wallace_known_value() {
        // AATTGGCC: AT=4, GC=4 → 2×4 + 4×4 = 24.
        assert!((wallace_tm(b"AATTGGCC") - 24.0).abs() < 1e-10);
    }

    #[test]
    fn wallace_gaps_stripped() {
        // "A-T" strips to "AT": AT=2, GC=0 → 2×2 = 4.
        assert!((wallace_tm(b"A-T") - 4.0).abs() < 1e-10);
    }

    #[test]
    fn wallace_all_gc() {
        // "GCGC": AT=0, GC=4 → 4×4 = 16.
        assert!((wallace_tm(b"GCGC") - 16.0).abs() < 1e-10);
    }

    #[test]
    fn calculate_tm_too_short_returns_none() {
        assert!(calculate_tm(b"").is_none());
        assert!(calculate_tm(b"A").is_none());
    }

    #[test]
    fn calculate_tm_valid_returns_some() {
        let r = calculate_tm(b"ATCGATCG").unwrap();
        // Tm should be a plausible real number (not NaN/inf), somewhere in 10–80 °C range.
        assert!(r.tm_celsius.is_finite());
        assert!(r.tm_celsius > 5.0 && r.tm_celsius < 90.0,
            "Tm out of expected range: {}", r.tm_celsius);
    }

    #[test]
    fn calculate_tm_gc_rich_higher_than_at_rich() {
        let gc = calculate_tm(b"GCGCGCGCGCGCGCGCGCGC").unwrap();
        let at = calculate_tm(b"ATATATATATATATATATATAT").unwrap();
        assert!(gc.tm_celsius > at.tm_celsius,
            "expected GC Tm ({:.1}) > AT Tm ({:.1})", gc.tm_celsius, at.tm_celsius);
    }

    #[test]
    fn calculate_tm_gc_fraction_correct() {
        let r = calculate_tm(b"GCGCGCGC").unwrap();
        assert!((r.gc_fraction - 1.0).abs() < 1e-6, "expected 1.0, got {}", r.gc_fraction);
        let r2 = calculate_tm(b"ATATATATAT").unwrap();
        assert!(r2.gc_fraction.abs() < 1e-6, "expected 0.0, got {}", r2.gc_fraction);
    }

    #[test]
    fn calculate_tm_strips_gaps() {
        let with_gaps    = calculate_tm(b"AT-CG-AT").unwrap();
        let without_gaps = calculate_tm(b"ATCGAT").unwrap();
        assert!((with_gaps.tm_celsius - without_gaps.tm_celsius).abs() < 1e-6);
    }
}
