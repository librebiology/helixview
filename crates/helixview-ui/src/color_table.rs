//! Per-residue display colors, editable at runtime.
//!
//! Two independent tables: one for nucleotides, one for amino acids.
//! Both are stored as arrays indexed by `byte - b'A'`.

use iced::Color;

/// Mutable display color table.  Stored as `Arc<ColorTable>` in `HelixViewApp`;
/// swapped out on every edit so the canvas version stamp picks up the pointer change.
#[derive(Debug, Clone)]
pub struct ColorTable {
    pub nuc: [Color; 26],
    pub aa:  [Color; 26],
}

impl Default for ColorTable {
    fn default() -> Self {
        Self {
            nuc: default_nuc_table(),
            aa:  default_aa_table(),
        }
    }
}

impl ColorTable {
    /// Look up the nucleotide background color for byte `b` (any case).
    pub fn nuc_color(&self, b: u8) -> Color {
        let idx = (b.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
        if idx < 26 { self.nuc[idx] } else { Color::from_rgb(0.90, 0.90, 0.90) }
    }

    /// Look up the amino-acid background color for byte `b` (any case).
    pub fn aa_color(&self, b: u8) -> Color {
        let idx = (b.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
        if idx < 26 { self.aa[idx] } else { Color::from_rgb(0.90, 0.90, 0.90) }
    }

    /// Set a nucleotide color (by uppercase byte).
    pub fn set_nuc(&mut self, b: u8, color: Color) {
        let idx = (b.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
        if idx < 26 { self.nuc[idx] = color; }
    }

    /// Set an amino-acid color (by uppercase byte).
    pub fn set_aa(&mut self, b: u8, color: Color) {
        let idx = (b.to_ascii_uppercase().wrapping_sub(b'A')) as usize;
        if idx < 26 { self.aa[idx] = color; }
    }
}

// ── Default palettes ──────────────────────────────────────────────────────────

fn default_nuc_table() -> [Color; 26] {
    let grey  = Color::from_rgb(0.90, 0.90, 0.90);
    let mut t = [grey; 26];
    let set = |t: &mut [Color; 26], b: u8, c: Color| { t[(b - b'A') as usize] = c; };

    // Pastel four-color scheme
    set(&mut t, b'A', Color::from_rgb(0.72, 0.94, 0.62)); // soft green
    set(&mut t, b'T', Color::from_rgb(0.98, 0.76, 0.76)); // soft coral
    set(&mut t, b'U', Color::from_rgb(0.98, 0.76, 0.76)); // same as T
    set(&mut t, b'G', Color::from_rgb(0.97, 0.92, 0.62)); // soft gold
    set(&mut t, b'C', Color::from_rgb(0.72, 0.86, 0.98)); // soft blue
    set(&mut t, b'N', Color::from_rgb(0.88, 0.88, 0.88));
    t
}

fn default_aa_table() -> [Color; 26] {
    let grey  = Color::from_rgb(0.90, 0.90, 0.90);
    let mut t = [grey; 26];
    let set = |t: &mut [Color; 26], b: u8, c: Color| { t[(b - b'A') as usize] = c; };

    // ClustalX-inspired pastel palette
    let hydrophobic = Color::from_rgb(0.98, 0.88, 0.68);
    let polar       = Color::from_rgb(0.75, 0.94, 0.80);
    let basic       = Color::from_rgb(0.75, 0.86, 0.98);
    let acidic      = Color::from_rgb(0.98, 0.78, 0.78);
    let amide       = Color::from_rgb(0.80, 0.95, 0.98);

    for b in [b'A', b'V', b'I', b'L', b'M', b'F', b'W', b'P'] {
        set(&mut t, b, hydrophobic);
    }
    for b in [b'S', b'T', b'Y', b'H', b'C'] {
        set(&mut t, b, polar);
    }
    for b in [b'K', b'R'] {
        set(&mut t, b, basic);
    }
    for b in [b'D', b'E'] {
        set(&mut t, b, acidic);
    }
    for b in [b'N', b'Q', b'G'] {
        set(&mut t, b, amide);
    }
    t
}

// ── Display helpers ───────────────────────────────────────────────────────────

/// Residues shown in the nucleotide section of the color editor.
pub const NUC_RESIDUES: &[u8] = b"ACGTU N";

/// Residues shown in the amino-acid section of the color editor.
pub const AA_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
