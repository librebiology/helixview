//! Analysis algorithms — Phase 5 onwards.
//! Stub crate: populated as phases are implemented.

pub mod composition;
pub mod consensus;
pub mod conservation;
pub mod entropy;
pub mod hydrophobicity;
pub mod identity;
pub mod mutual_info;
pub mod nj_tree;
pub mod oligo_tm;
pub mod pairwise;
pub mod restriction;
pub mod sixframe;

pub use composition::amino_acid_composition;
pub use composition::nucleotide_composition;
pub use consensus::{consensus, to_iupac, ConsensusMethod};
pub use conservation::{find_conserved_regions, ConservedRegion};
pub use entropy::column_entropy;
pub use hydrophobicity::{kd_scale, kyte_doolittle_profile};
pub use identity::{identity_matrix, pairwise_identity};
pub use mutual_info::{mutual_information, MiPair, MiResult};
pub use nj_tree::nj_tree_from_alignment;
pub use oligo_tm::{calculate_tm, wallace_tm, TmResult};
pub use pairwise::{
    needleman_wunsch, smith_waterman, PairwiseResult, DEFAULT_GAP, DEFAULT_MATCH, DEFAULT_MISMATCH,
};
pub use restriction::{find_sites, RestrictionSite, ENZYMES};
pub use sixframe::{find_orfs, six_frame_all, translate_codon, translate_frame, Orf};
