use crate::color::Color;

#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum FeatureShape {
    #[default]
    Rectangle,
    Oval,
    Diamond,
    Arrow,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum FeatureDirection {
    Forward,
    Reverse,
    #[default]
    None,
}

/// All 67 standard GenBank feature key types, plus a catch-all.
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum FeatureType {
    ThreePrimeClip,
    ThreePrimeUtr,
    FivePrimeClip,
    FivePrimeUtr,
    MinusTenSignal,
    MinusThirtyFiveSignal,
    Allele,
    Attenuator,
    CaatSignal,
    CDS,
    Conflict,
    DLoop,
    Enhancer,
    Exon,
    GCSignal,
    Gene,
    ImpropRegion,
    Intron,
    LTR,
    MatPeptide,
    MiscBinding,
    MiscDifference,
    #[default]
    MiscFeature,
    MiscRecomb,
    MiscRNA,
    MiscSignal,
    MiscStructure,
    ModifiedBase,
    MRNA,
    MutationRegion,
    NRegion,
    OldSequence,
    OpCode,
    Operon,
    PolyASignal,
    PolyASite,
    PrecursorRNA,
    PrimBind,
    PrimTranscript,
    Primer,
    PromotionalRegion,
    ProprotPeptide,
    Protein,
    RBS,
    RepeatRegion,
    RepeatUnit,
    RRNA,
    SRegion,
    Satellite,
    ScRNA,
    SecStr,
    SigPeptide,
    SnRNA,
    Source,
    StemLoop,
    STS,
    TATA,
    Terminator,
    TransitPeptide,
    Transposon,
    TRNA,
    Unsure,
    VRegion,
    VarSeq,
    Variation,
    Virion,
    NonStdResidue,
    Custom(String),
}

impl FeatureType {
    /// Default display color for rendering this feature type.
    pub fn default_color(&self) -> Color {
        match self {
            Self::Gene => Color::from_u8(52, 120, 205), // blue
            Self::CDS => Color::from_u8(220, 115, 25),  // orange
            Self::MRNA => Color::from_u8(50, 180, 65),  // green
            Self::MiscRNA | Self::SnRNA | Self::ScRNA => Color::from_u8(80, 180, 80),
            Self::Exon => Color::from_u8(140, 50, 190), // purple
            Self::Intron => Color::from_u8(180, 180, 200),
            Self::TRNA => Color::from_u8(215, 50, 50), // red
            Self::RRNA => Color::from_u8(220, 60, 120), // rose
            Self::RepeatRegion | Self::RepeatUnit => Color::from_u8(210, 200, 50), // yellow
            Self::Source => Color::TRANSPARENT,
            Self::StemLoop => Color::from_u8(200, 100, 100),
            Self::Terminator | Self::PromotionalRegion => Color::from_u8(150, 80, 200),
            _ => Color::from_u8(120, 120, 180), // slate
        }
    }

    pub fn from_genbank_key(key: &str) -> Self {
        match key {
            "3'clip" => Self::ThreePrimeClip,
            "3'UTR" => Self::ThreePrimeUtr,
            "5'clip" => Self::FivePrimeClip,
            "5'UTR" => Self::FivePrimeUtr,
            "-10_signal" => Self::MinusTenSignal,
            "-35_signal" => Self::MinusThirtyFiveSignal,
            "allele" => Self::Allele,
            "attenuator" => Self::Attenuator,
            "CAAT_signal" => Self::CaatSignal,
            "CDS" => Self::CDS,
            "conflict" => Self::Conflict,
            "D-loop" => Self::DLoop,
            "enhancer" => Self::Enhancer,
            "exon" => Self::Exon,
            "GC_signal" => Self::GCSignal,
            "gene" => Self::Gene,
            "intron" => Self::Intron,
            "LTR" => Self::LTR,
            "mat_peptide" => Self::MatPeptide,
            "misc_binding" => Self::MiscBinding,
            "misc_difference" => Self::MiscDifference,
            "misc_feature" => Self::MiscFeature,
            "misc_recomb" => Self::MiscRecomb,
            "misc_RNA" => Self::MiscRNA,
            "misc_signal" => Self::MiscSignal,
            "misc_structure" => Self::MiscStructure,
            "modified_base" => Self::ModifiedBase,
            "mRNA" => Self::MRNA,
            "N_region" => Self::NRegion,
            "old_sequence" => Self::OldSequence,
            "operon" => Self::Operon,
            "polyA_signal" => Self::PolyASignal,
            "polyA_site" => Self::PolyASite,
            "precursor_RNA" => Self::PrecursorRNA,
            "prim_transcript" => Self::PrimTranscript,
            "primer" => Self::Primer,
            "primer_bind" => Self::PrimBind,
            "promoter" => Self::PromotionalRegion,
            "propeptide" => Self::ProprotPeptide,
            "protein_bind" => Self::Protein,
            "RBS" => Self::RBS,
            "repeat_region" => Self::RepeatRegion,
            "repeat_unit" => Self::RepeatUnit,
            "rRNA" => Self::RRNA,
            "S_region" => Self::SRegion,
            "satellite" => Self::Satellite,
            "scRNA" => Self::ScRNA,
            "SecStr" => Self::SecStr,
            "sig_peptide" => Self::SigPeptide,
            "snRNA" => Self::SnRNA,
            "source" => Self::Source,
            "stem_loop" => Self::StemLoop,
            "STS" => Self::STS,
            "TATA_signal" => Self::TATA,
            "terminator" => Self::Terminator,
            "transit_peptide" => Self::TransitPeptide,
            "transposon" => Self::Transposon,
            "tRNA" => Self::TRNA,
            "unsure" => Self::Unsure,
            "V_region" => Self::VRegion,
            "variation" => Self::Variation,
            "virion" => Self::Virion,
            "NonStdResidue" => Self::NonStdResidue,
            other => Self::Custom(other.to_string()),
        }
    }
}

/// A single annotated region on a sequence.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Feature {
    pub name: String,
    /// Raw qualifier block (GenBank format) or user description.
    pub description: String,
    /// True (ungapped) start position, 0-indexed.
    pub start: usize,
    /// True (ungapped) end position, 0-indexed, inclusive.
    pub end: usize,
    pub shape: FeatureShape,
    pub color: Color,
    pub direction: FeatureDirection,
    pub feature_type: FeatureType,
}

impl Feature {
    /// Parse a feature type from a GenBank key string (convenience wrapper).
    pub fn type_from_str(key: &str) -> FeatureType {
        FeatureType::from_genbank_key(key)
    }

    pub fn new(name: impl Into<String>, start: usize, end: usize) -> Self {
        Self {
            name: name.into(),
            description: String::new(),
            start,
            end,
            shape: FeatureShape::default(),
            color: Color::from_u8(0, 120, 255),
            direction: FeatureDirection::default(),
            feature_type: FeatureType::default(),
        }
    }

    /// Length in residues (ungapped).
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start) + 1
    }

    pub fn is_empty(&self) -> bool {
        self.start == self.end
    }
}
