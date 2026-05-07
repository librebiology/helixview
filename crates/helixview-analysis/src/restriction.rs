//! Restriction enzyme recognition site finder.
//!
//! Searches a raw DNA sequence (gaps already stripped, uppercase) for the
//! recognition sequences of a built-in set of common restriction enzymes and
//! returns the position (0-based, in the raw sequence) of each cut.

/// A restriction enzyme entry.
#[derive(Debug, Clone)]
pub struct RestrictionEnzyme {
    /// Common name, e.g. "EcoRI".
    pub name:        &'static str,
    /// Recognition sequence (uppercase ASCII, e.g. b"GAATTC").
    pub recognition: &'static [u8],
    /// Colour index (0–19) into the UI palette.
    pub color_idx:   usize,
}

/// Built-in set of common Type-II restriction enzymes.
pub static ENZYMES: &[RestrictionEnzyme] = &[
    RestrictionEnzyme { name: "EcoRI",  recognition: b"GAATTC", color_idx:  0 },
    RestrictionEnzyme { name: "BamHI",  recognition: b"GGATCC", color_idx:  1 },
    RestrictionEnzyme { name: "HindIII",recognition: b"AAGCTT", color_idx:  2 },
    RestrictionEnzyme { name: "XbaI",   recognition: b"TCTAGA", color_idx:  3 },
    RestrictionEnzyme { name: "SalI",   recognition: b"GTCGAC", color_idx:  4 },
    RestrictionEnzyme { name: "KpnI",   recognition: b"GGTACC", color_idx:  5 },
    RestrictionEnzyme { name: "SacI",   recognition: b"GAGCTC", color_idx:  6 },
    RestrictionEnzyme { name: "EcoRV",  recognition: b"GATATC", color_idx:  7 },
    RestrictionEnzyme { name: "NcoI",   recognition: b"CCATGG", color_idx:  8 },
    RestrictionEnzyme { name: "NdeI",   recognition: b"CATATG", color_idx:  9 },
    RestrictionEnzyme { name: "BglII",  recognition: b"AGATCT", color_idx: 10 },
    RestrictionEnzyme { name: "ClaI",   recognition: b"ATCGAT", color_idx: 11 },
    RestrictionEnzyme { name: "NheI",   recognition: b"GCTAGC", color_idx: 12 },
    RestrictionEnzyme { name: "SpeI",   recognition: b"ACTAGT", color_idx: 13 },
    RestrictionEnzyme { name: "XhoI",   recognition: b"CTCGAG", color_idx: 14 },
    RestrictionEnzyme { name: "PstI",   recognition: b"CTGCAG", color_idx: 15 },
    RestrictionEnzyme { name: "SphI",   recognition: b"GCATGC", color_idx: 16 },
    RestrictionEnzyme { name: "SmaI",   recognition: b"CCCGGG", color_idx: 17 },
    RestrictionEnzyme { name: "MluI",   recognition: b"ACGCGT", color_idx: 18 },
    RestrictionEnzyme { name: "AvrII",  recognition: b"CCTAGG", color_idx: 19 },
];

/// A single enzyme hit in a sequence.
#[derive(Debug, Clone)]
pub struct RestrictionSite {
    /// Index into [`ENZYMES`].
    pub enzyme_idx: usize,
    /// 0-based start position in the *ungapped* sequence.
    pub true_pos:   usize,
}

/// Find all restriction sites in a raw (gap-free, uppercase ASCII) DNA sequence.
/// Returns sites sorted by `true_pos`.
pub fn find_sites(raw_seq: &[u8]) -> Vec<RestrictionSite> {
    let mut sites: Vec<RestrictionSite> = Vec::new();

    for (ei, enzyme) in ENZYMES.iter().enumerate() {
        let recog = enzyme.recognition;
        let rlen  = recog.len();
        if raw_seq.len() < rlen { continue; }
        for i in 0..=(raw_seq.len() - rlen) {
            if raw_seq[i..i + rlen] == *recog {
                sites.push(RestrictionSite { enzyme_idx: ei, true_pos: i });
            }
        }
    }

    sites.sort_unstable_by_key(|s| s.true_pos);
    sites
}
