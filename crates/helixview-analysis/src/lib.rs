//! Analysis algorithms — Phase 5 onwards.
//! Stub crate: populated as phases are implemented.

pub mod composition;
pub mod entropy;
pub mod identity;
pub mod pairwise;
pub mod hydrophobicity;
pub mod mutual_info;
pub mod sixframe;
pub mod consensus;
pub mod restriction;
pub mod oligo_tm;
pub mod conservation;
pub mod nj_tree;

pub use composition::nucleotide_composition;
pub use composition::amino_acid_composition;
pub use entropy::column_entropy;
pub use identity::{pairwise_identity, identity_matrix};
pub use pairwise::{needleman_wunsch, smith_waterman, PairwiseResult, DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP};
pub use hydrophobicity::{kyte_doolittle_profile, kd_scale};
pub use sixframe::{translate_frame, find_orfs, translate_codon, six_frame_all, Orf};
pub use consensus::{consensus, ConsensusMethod, to_iupac};
pub use restriction::{find_sites, RestrictionSite, ENZYMES};
pub use oligo_tm::{calculate_tm, wallace_tm, TmResult};
pub use conservation::{find_conserved_regions, ConservedRegion};
pub use mutual_info::{mutual_information, MiResult, MiPair};
pub use nj_tree::nj_tree_from_alignment;
