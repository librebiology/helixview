use helixview_core::sequence::Sequence;

#[derive(Debug, Clone)]
pub struct AminoAcidComposition {
    pub a: u64,
    pub r: u64,
    pub n: u64,
    pub d: u64,
    pub c: u64,
    pub e: u64,
    pub q: u64,
    pub g: u64,
    pub h: u64,
    pub i: u64,
    pub l: u64,
    pub k: u64,
    pub m: u64,
    pub f: u64,
    pub p: u64,
    pub s: u64,
    pub t: u64,
    pub w: u64,
    pub y: u64,
    pub v: u64,
    pub other: u64,
    pub total: u64,
    /// Approximate average molecular weight in Daltons.
    pub mw: f64,
}

pub fn amino_acid_composition(seq: &Sequence) -> AminoAcidComposition {
    let mut a = 0u64;
    let mut r = 0u64;
    let mut n = 0u64;
    let mut d = 0u64;
    let mut c = 0u64;
    let mut e = 0u64;
    let mut q = 0u64;
    let mut g = 0u64;
    let mut h = 0u64;
    let mut i = 0u64;
    let mut l = 0u64;
    let mut k = 0u64;
    let mut m = 0u64;
    let mut f = 0u64;
    let mut p = 0u64;
    let mut s = 0u64;
    let mut t = 0u64;
    let mut w = 0u64;
    let mut y = 0u64;
    let mut v = 0u64;
    let mut other = 0u64;

    for &b in &seq.residues {
        match b.to_ascii_uppercase() {
            b'A' => a += 1,
            b'R' => r += 1,
            b'N' => n += 1,
            b'D' => d += 1,
            b'C' => c += 1,
            b'E' => e += 1,
            b'Q' => q += 1,
            b'G' => g += 1,
            b'H' => h += 1,
            b'I' => i += 1,
            b'L' => l += 1,
            b'K' => k += 1,
            b'M' => m += 1,
            b'F' => f += 1,
            b'P' => p += 1,
            b'S' => s += 1,
            b'T' => t += 1,
            b'W' => w += 1,
            b'Y' => y += 1,
            b'V' => v += 1,
            b'-' | b'~' | b'.' => {}
            _ => other += 1,
        }
    }

    let total =
        a + r + n + d + c + e + q + g + h + i + l + k + m + f + p + s + t + w + y + v + other;

    let residue_mw = a as f64 * 89.09
        + r as f64 * 174.20
        + n as f64 * 132.12
        + d as f64 * 133.10
        + c as f64 * 121.16
        + e as f64 * 147.13
        + q as f64 * 146.15
        + g as f64 * 75.03
        + h as f64 * 155.16
        + i as f64 * 131.17
        + l as f64 * 131.17
        + k as f64 * 146.19
        + m as f64 * 149.21
        + f as f64 * 165.19
        + p as f64 * 115.13
        + s as f64 * 105.09
        + t as f64 * 119.12
        + w as f64 * 204.23
        + y as f64 * 181.19
        + v as f64 * 117.15;

    let aa_count = total - other;
    let mw = if aa_count > 0 {
        residue_mw - (aa_count as f64 - 1.0) * 18.02
    } else {
        0.0
    };

    AminoAcidComposition {
        a,
        r,
        n,
        d,
        c,
        e,
        q,
        g,
        h,
        i,
        l,
        k,
        m,
        f,
        p,
        s,
        t,
        w,
        y,
        v,
        other,
        total,
        mw,
    }
}

#[cfg(test)]
mod aa_tests {
    use super::*;

    #[test]
    fn all_twenty_aas_counted() {
        let seq = Sequence::new("test", b"ACDEFGHIKLMNPQRSTVWY".to_vec());
        let comp = amino_acid_composition(&seq);
        assert_eq!(comp.a, 1);
        assert_eq!(comp.c, 1);
        assert_eq!(comp.d, 1);
        assert_eq!(comp.e, 1);
        assert_eq!(comp.f, 1);
        assert_eq!(comp.g, 1);
        assert_eq!(comp.h, 1);
        assert_eq!(comp.i, 1);
        assert_eq!(comp.k, 1);
        assert_eq!(comp.l, 1);
        assert_eq!(comp.m, 1);
        assert_eq!(comp.n, 1);
        assert_eq!(comp.p, 1);
        assert_eq!(comp.q, 1);
        assert_eq!(comp.r, 1);
        assert_eq!(comp.s, 1);
        assert_eq!(comp.t, 1);
        assert_eq!(comp.v, 1);
        assert_eq!(comp.w, 1);
        assert_eq!(comp.y, 1);
        assert_eq!(comp.other, 0);
    }

    #[test]
    fn mw_reasonable_magnitude() {
        let seq = Sequence::new("test", b"ACDEFGHIKLMNPQRSTVWY".to_vec());
        let comp = amino_acid_composition(&seq);
        assert!(comp.mw > 0.0);
        assert!(comp.mw > 1000.0 && comp.mw < 5000.0);
    }
}

/// Returns (A, C, G, T, other, gc_percent, at_percent, molecular_weight_Da)
pub fn nucleotide_composition(seq: &Sequence) -> NucleotideComposition {
    let mut a = 0u64;
    let mut c = 0u64;
    let mut g = 0u64;
    let mut t = 0u64;
    let mut u = 0u64;
    let mut other = 0u64;

    for &b in &seq.residues {
        match b.to_ascii_uppercase() {
            b'A' => a += 1,
            b'C' => c += 1,
            b'G' => g += 1,
            b'T' => t += 1,
            b'U' => u += 1,
            b'-' | b'~' | b'.' => {}
            _ => other += 1,
        }
    }
    let total = a + c + g + t + u + other;
    let gc = g + c;
    let at = a + t + u;
    let gc_pct = if total > 0 {
        gc as f64 / total as f64 * 100.0
    } else {
        0.0
    };
    let at_pct = if total > 0 {
        at as f64 / total as f64 * 100.0
    } else {
        0.0
    };

    // Average molecular weights (monoisotopic, single-stranded, approximate)
    let mw = a as f64 * 313.21 + c as f64 * 289.18 + g as f64 * 329.21 + (t + u) as f64 * 304.19;

    NucleotideComposition {
        a,
        c,
        g,
        t,
        u,
        other,
        total,
        gc_pct,
        at_pct,
        mw,
    }
}

#[derive(Debug, Clone)]
pub struct NucleotideComposition {
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub u: u64,
    pub other: u64,
    pub total: u64,
    pub gc_pct: f64,
    pub at_pct: f64,
    /// Approximate monoisotopic molecular weight in Daltons.
    pub mw: f64,
}
