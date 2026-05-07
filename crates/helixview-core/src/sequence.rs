use crate::{color::Color, feature::Feature};

/// The biological type of a sequence row.
#[derive(Debug, Clone, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum SequenceType {
    Dna,
    Rna,
    /// Nucleic acid, strand unspecified
    NucleicAcid,
    Protein,
    #[default]
    Unknown,
    /// Non-sequence comment row — excluded from all calculations
    Comment,
    /// Mask row defining which positions count as valid residues
    SequenceMask,
    /// RNA secondary structure mask
    RnaStructureMask,
}

impl SequenceType {
    pub fn is_nucleic(&self) -> bool {
        matches!(self, Self::Dna | Self::Rna | Self::NucleicAcid)
    }

    pub fn is_sequence(&self) -> bool {
        !matches!(
            self,
            Self::Comment | Self::SequenceMask | Self::RnaStructureMask
        )
    }
}

impl std::fmt::Display for SequenceType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Dna => write!(f, "DNA"),
            Self::Rna => write!(f, "RNA"),
            Self::NucleicAcid => write!(f, "Nucleic Acid"),
            Self::Protein => write!(f, "Protein"),
            Self::Unknown => write!(f, "Unknown"),
            Self::Comment => write!(f, "Comment"),
            Self::SequenceMask => write!(f, "Sequence Mask"),
            Self::RnaStructureMask => write!(f, "RNA Structure Mask"),
        }
    }
}

/// GenBank flat-file metadata stored per sequence.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct GenBankMeta {
    pub locus: Option<String>,
    pub definition: Option<String>,
    pub accession: Option<String>,
    pub version: Option<String>,
    pub pid: Option<String>,
    pub dbsource: Option<String>,
    pub keywords: Option<String>,
    pub source: Option<String>,
    pub references: Vec<String>,
    pub comment: Option<String>,
    /// Raw FEATURES block text, preserved exactly as parsed.
    pub features_raw: Option<String>,
}

/// Gap character interpretation:
///   b'-'  = locked gap (will not be crunched during sliding)
///   b'~'  = unlocked gap
///   b'.'  = unlocked gap (compatibility alias)
pub const GAP_LOCKED: u8 = b'-';
pub const GAP_UNLOCKED: u8 = b'~';
pub const GAP_PERIOD: u8 = b'.';

pub fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'~' | b'.')
}

pub fn is_locked_gap(b: u8) -> bool {
    b == GAP_LOCKED
}

/// A single sequence row within an alignment document.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Sequence {
    /// Display name shown in the title panel.
    pub name: String,
    /// Residues including gap characters, as raw bytes.
    /// Invariant: len == alignment column count.
    pub residues: Vec<u8>,
    pub seq_type: SequenceType,
    /// If true, hand-alignment edits (slide/grab) are disabled for this row.
    /// Gap insert/delete via right-click is still allowed.
    pub locked: bool,
    /// Per-position gap locks. Index matches residues.
    /// Only meaningful for positions that are gaps.
    pub gap_locks: Vec<bool>,
    pub features: Vec<Feature>,
    pub genbank: GenBankMeta,
    /// Color used to paint this sequence's title label.
    pub label_color: Color,
    /// Family/group index — None = ungrouped.
    pub group: Option<usize>,
}

impl Sequence {
    pub fn new(name: impl Into<String>, residues: Vec<u8>) -> Self {
        let len = residues.len();
        Self {
            name: name.into(),
            residues,
            seq_type: SequenceType::Unknown,
            locked: false,
            gap_locks: vec![false; len],
            features: Vec::new(),
            genbank: GenBankMeta::default(),
            label_color: Color::BLACK,
            group: None,
        }
    }

    /// Length including gaps (= number of alignment columns).
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Length excluding all gap characters (true biological length).
    pub fn true_len(&self) -> usize {
        self.residues.iter().filter(|&&b| !is_gap(b)).count()
    }

    /// Convert an alignment column index to the true (ungapped) position.
    /// Returns None if the column is a gap.
    pub fn true_position(&self, col: usize) -> Option<usize> {
        if col >= self.residues.len() {
            return None;
        }
        if is_gap(self.residues[col]) {
            return None;
        }
        Some(self.residues[..col].iter().filter(|&&b| !is_gap(b)).count())
    }

    /// Convert a true position to the corresponding alignment column.
    pub fn alignment_col(&self, true_pos: usize) -> Option<usize> {
        let mut count = 0usize;
        for (col, &b) in self.residues.iter().enumerate() {
            if !is_gap(b) {
                if count == true_pos {
                    return Some(col);
                }
                count += 1;
            }
        }
        None
    }

    /// Raw sequence without gaps, as a String (lossy UTF-8).
    pub fn raw_sequence(&self) -> String {
        String::from_utf8_lossy(
            &self
                .residues
                .iter()
                .copied()
                .filter(|&b| !is_gap(b))
                .collect::<Vec<_>>(),
        )
        .into_owned()
    }

    /// Uppercase all residues (not gaps).
    pub fn to_uppercase(&mut self) {
        for b in &mut self.residues {
            if !is_gap(*b) {
                b.make_ascii_uppercase();
            }
        }
    }

    /// Lowercase all residues (not gaps).
    pub fn to_lowercase(&mut self) {
        for b in &mut self.residues {
            if !is_gap(*b) {
                b.make_ascii_lowercase();
            }
        }
    }

    /// Reverse the sequence (gaps move too, for raw reversal).
    pub fn reverse(&mut self) {
        self.residues.reverse();
        self.gap_locks.reverse();
    }

    /// DNA complement in place (A↔T, C↔G). No-op on protein/unknown.
    pub fn complement(&mut self) {
        if !self.seq_type.is_nucleic() {
            return;
        }
        for b in &mut self.residues {
            *b = complement_base(*b);
        }
    }

    /// Reverse complement in place.
    pub fn reverse_complement(&mut self) {
        self.complement();
        self.reverse();
    }

    /// Convert T→U (DNA to RNA).
    pub fn dna_to_rna(&mut self) {
        for b in &mut self.residues {
            if *b == b'T' {
                *b = b'U';
            }
            if *b == b't' {
                *b = b'u';
            }
        }
        self.seq_type = SequenceType::Rna;
    }

    /// Convert U→T (RNA to DNA).
    pub fn rna_to_dna(&mut self) {
        for b in &mut self.residues {
            if *b == b'U' {
                *b = b'T';
            }
            if *b == b'u' {
                *b = b't';
            }
        }
        self.seq_type = SequenceType::Dna;
    }
}

/// Public complement helper used by the free `reverse_complement` function.
/// Identical to `complement_base` but pub(crate).
pub(crate) fn complement_base_rc(b: u8) -> u8 {
    complement_base(b)
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'U' => b'A',
        b'u' => b'a',
        b'R' => b'Y',
        b'Y' => b'R',
        b'r' => b'y',
        b'y' => b'r',
        b'S' => b'S',
        b's' => b's',
        b'W' => b'W',
        b'w' => b'w',
        b'K' => b'M',
        b'M' => b'K',
        b'k' => b'm',
        b'm' => b'k',
        b'B' => b'V',
        b'V' => b'B',
        b'b' => b'v',
        b'v' => b'b',
        b'D' => b'H',
        b'H' => b'D',
        b'd' => b'h',
        b'h' => b'd',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn true_len_counts_only_residues() {
        let seq = Sequence::new("test", b"A-T~CG.".to_vec());
        assert_eq!(seq.true_len(), 4); // A, T, C, G
    }

    #[test]
    fn true_position_skips_gaps() {
        let seq = Sequence::new("test", b"A-TG".to_vec());
        assert_eq!(seq.true_position(0), Some(0)); // A → pos 0
        assert_eq!(seq.true_position(1), None); // gap
        assert_eq!(seq.true_position(2), Some(1)); // T → pos 1
        assert_eq!(seq.true_position(3), Some(2)); // G → pos 2
    }

    #[test]
    fn alignment_col_roundtrip() {
        let seq = Sequence::new("test", b"A-TG".to_vec());
        assert_eq!(seq.alignment_col(0), Some(0));
        assert_eq!(seq.alignment_col(1), Some(2));
        assert_eq!(seq.alignment_col(2), Some(3));
    }

    #[test]
    fn complement_dna() {
        let mut seq = Sequence::new("test", b"ATCG".to_vec());
        seq.seq_type = SequenceType::Dna;
        seq.complement();
        assert_eq!(&seq.residues, b"TAGC");
    }

    #[test]
    fn reverse_complement() {
        let mut seq = Sequence::new("test", b"ATCG".to_vec());
        seq.seq_type = SequenceType::Dna;
        seq.reverse_complement();
        assert_eq!(&seq.residues, b"CGAT");
    }
}
