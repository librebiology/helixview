/// Platform-independent RGBA color with f32 components in [0.0, 1.0].
/// Kept in core so it has no iced dependency; the UI layer converts as needed.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Color {
    pub r: f32,
    pub g: f32,
    pub b: f32,
    pub a: f32,
}

impl Color {
    pub const BLACK: Self = Self {
        r: 0.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    pub const WHITE: Self = Self {
        r: 1.0,
        g: 1.0,
        b: 1.0,
        a: 1.0,
    };
    pub const TRANSPARENT: Self = Self {
        r: 0.0,
        g: 0.0,
        b: 0.0,
        a: 0.0,
    };

    // Standard nucleotide colors (Megan convention)
    pub const DNA_A: Self = Self::from_u8(0, 200, 0);
    pub const DNA_C: Self = Self::from_u8(0, 0, 200);
    pub const DNA_G: Self = Self::from_u8(220, 180, 0);
    pub const DNA_T: Self = Self::from_u8(200, 0, 0);
    pub const DNA_U: Self = Self::from_u8(200, 0, 0);
    pub const DNA_GAP: Self = Self::WHITE;

    pub const fn rgb(r: f32, g: f32, b: f32) -> Self {
        Self { r, g, b, a: 1.0 }
    }

    pub const fn from_u8(r: u8, g: u8, b: u8) -> Self {
        Self::rgb(r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0)
    }

    pub const fn with_alpha(mut self, a: f32) -> Self {
        self.a = a;
        self
    }
}

impl Default for Color {
    fn default() -> Self {
        Self::BLACK
    }
}

/// A mapping from a residue byte value to its display color.
/// Note: large arrays ([T; 256]) don't derive serde — serialization handled separately.
#[derive(Debug, Clone)]
pub struct ColorTable {
    /// Separate tables for nucleotides and amino acids.
    pub nucleotide: [Color; 256],
    pub amino_acid: [Color; 256],
}

impl Default for ColorTable {
    fn default() -> Self {
        let mut nt = [Color::WHITE; 256];
        nt[b'A' as usize] = Color::DNA_A;
        nt[b'a' as usize] = Color::DNA_A;
        nt[b'C' as usize] = Color::DNA_C;
        nt[b'c' as usize] = Color::DNA_C;
        nt[b'G' as usize] = Color::DNA_G;
        nt[b'g' as usize] = Color::DNA_G;
        nt[b'T' as usize] = Color::DNA_T;
        nt[b't' as usize] = Color::DNA_T;
        nt[b'U' as usize] = Color::DNA_U;
        nt[b'u' as usize] = Color::DNA_U;

        // Amino acid defaults: simple CPK-inspired scheme
        let mut aa = [Color::WHITE; 256];
        // Hydrophobic (orange)
        for b in b"AILMFWVailmfwv" {
            aa[*b as usize] = Color::from_u8(255, 140, 0);
        }
        // Positive (blue)
        for b in b"RKHrkh" {
            aa[*b as usize] = Color::from_u8(0, 70, 255);
        }
        // Negative (red)
        for b in b"DEde" {
            aa[*b as usize] = Color::from_u8(220, 0, 0);
        }
        // Polar uncharged (green)
        for b in b"STNQstnq" {
            aa[*b as usize] = Color::from_u8(0, 180, 0);
        }
        // Special
        aa[b'G' as usize] = Color::from_u8(180, 180, 180);
        aa[b'g' as usize] = Color::from_u8(180, 180, 180);
        aa[b'P' as usize] = Color::from_u8(255, 100, 180);
        aa[b'p' as usize] = Color::from_u8(255, 100, 180);
        aa[b'C' as usize] = Color::from_u8(200, 200, 0);
        aa[b'c' as usize] = Color::from_u8(200, 200, 0);
        aa[b'Y' as usize] = Color::from_u8(0, 180, 180);
        aa[b'y' as usize] = Color::from_u8(0, 180, 180);

        Self {
            nucleotide: nt,
            amino_acid: aa,
        }
    }
}
