//! EMBL flat-file parser and writer.
//!
//! Handles the EBI/UniProt EMBL flat-file format:
//!   - Records separated by `//`
//!   - `ID` line: entry name and molecule info
//!   - `DE` lines: description
//!   - `FT` lines: feature table (gene, CDS, exon, misc_feature, etc.)
//!   - `SQ` line + indented sequence blocks
//!   - Gap characters (`-`) preserved as-is

use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    feature::{Feature, FeatureDirection, FeatureShape, FeatureType},
    sequence::{Sequence, SequenceType},
};
use std::path::Path;

pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut aln = Alignment::new(name);

    // Split records at lines that are exactly "//" (with optional whitespace).
    // We split on "\n//" which handles both "\n//\n" and a trailing "//".
    let records: Vec<&str> = text.split("\n//").collect();

    for record in records {
        let record = record.trim();
        if record.is_empty() {
            continue;
        }
        // Skip if the record is just "//" leftovers
        if record == "//" {
            continue;
        }
        if let Some(seq) = parse_record(record) {
            aln.push(seq);
        }
    }

    if aln.is_empty() {
        return Err(FormatError::Parse("No sequences found in EMBL data".into()));
    }

    Ok(aln)
}

fn parse_record(record: &str) -> Option<Sequence> {
    let mut seq_name = String::new();
    let mut description_parts: Vec<String> = Vec::new();
    let mut residues: Vec<u8> = Vec::new();
    let mut seq_type = SequenceType::Dna;
    let mut features: Vec<Feature> = Vec::new();
    let mut in_sequence = false;

    // Feature parsing state
    let mut cur_key: Option<String> = None;
    let mut cur_loc = String::new();
    let mut cur_quals: Vec<(String, String)> = Vec::new();
    let mut cur_qname: Option<String> = None;
    let mut cur_qval = String::new();
    let mut in_ft_location = true;

    for line in record.lines() {
        // EMBL lines have a 2-char tag followed by 3 spaces: "XX   content"
        // Special: SQ sequence lines are indented with spaces (no tag).

        if in_sequence {
            // Sequence data lines: "     acgtacgt acgtacgt   60"
            // Just collect all alphabetic characters and '-' gaps.
            if line.starts_with("//") {
                break;
            }
            for b in line.bytes() {
                if b.is_ascii_alphabetic() || b == b'-' {
                    residues.push(b);
                }
            }
            continue;
        }

        if line.len() < 2 {
            continue;
        }

        let tag = &line[..2];

        match tag {
            "ID" => {
                // ID   <entry_name>; SV <ver>; linear; <mol_type>; ...
                // or:  ID   <entry_name>; <status>; <mol_type>; <len> BP.
                let content = line[2..].trim();
                // Entry name is the first token before ';'
                seq_name = content
                    .split(';')
                    .next()
                    .unwrap_or(content)
                    .trim()
                    .to_string();

                // Detect molecule type from any semicolon-separated field
                for field in content.split(';') {
                    let f = field.trim().to_uppercase();
                    if f == "DNA" || f == "GENOMIC DNA" || f == "OTHER DNA" {
                        seq_type = SequenceType::Dna;
                        break;
                    } else if f == "RNA" || f.contains("RNA") {
                        seq_type = SequenceType::Rna;
                        break;
                    } else if f == "AA" || f == "PROTEIN" {
                        seq_type = SequenceType::Protein;
                        break;
                    }
                }
            }

            "DE" => {
                let content = line[2..].trim().to_string();
                if !content.is_empty() {
                    description_parts.push(content);
                }
            }

            "FT" => {
                // Feature table lines.
                // Key lines:   "FT   <key>          <location>"  (key starts at col 5, location at col 21)
                // Qualifier:   "FT                   /key=value"   (starts at col 21)
                // Continuation:"FT                   more_value"   (starts at col 21)

                if line.len() < 5 {
                    continue;
                }

                let ft_content = if line.len() > 5 { &line[5..] } else { "" };

                // Determine indentation within the FT content portion (after "FT   ")
                let leading = ft_content.bytes().take_while(|&b| b == b' ').count();

                if leading == 0 && !ft_content.trim().is_empty() {
                    // New feature key line: "FT   gene            100..200"
                    // Flush previous feature
                    flush_embl_qualifier(&mut cur_qname, &mut cur_qval, &mut cur_quals);
                    if let Some(key) = cur_key.take() {
                        if let Some(f) = build_embl_feature(&key, &cur_loc, &cur_quals) {
                            features.push(f);
                        }
                    }
                    cur_quals.clear();
                    cur_loc.clear();
                    in_ft_location = true;

                    // Split key and location
                    let (key, loc) = split_ft_key_location(ft_content);
                    cur_key = Some(key);
                    cur_loc = loc;
                } else if leading >= 16 {
                    let rest = ft_content[leading..].trim_start();
                    if let Some(body) = rest.strip_prefix('/') {
                        // New qualifier
                        flush_embl_qualifier(&mut cur_qname, &mut cur_qval, &mut cur_quals);
                        in_ft_location = false;
                        if let Some(eq) = body.find('=') {
                            let qkey = body[..eq].to_string();
                            let raw_val = body[eq + 1..].trim_matches('"').to_string();
                            cur_qname = Some(qkey);
                            cur_qval = raw_val;
                        } else {
                            cur_quals.push((body.to_string(), String::new()));
                        }
                    } else if !rest.is_empty() {
                        // Continuation line
                        let c = rest.trim_matches('"');
                        if in_ft_location {
                            cur_loc.push_str(c);
                        } else {
                            if !cur_qval.is_empty() {
                                cur_qval.push(' ');
                            }
                            cur_qval.push_str(c);
                        }
                    }
                }
            }

            "SQ" => {
                // SQ   Sequence 929 BP; 234 A; 237 C; 220 G; 238 T; 0 other;
                // Flush the last FT feature before entering sequence mode
                flush_embl_qualifier(&mut cur_qname, &mut cur_qval, &mut cur_quals);
                if let Some(key) = cur_key.take() {
                    if let Some(f) = build_embl_feature(&key, &cur_loc, &cur_quals) {
                        features.push(f);
                    }
                }
                in_sequence = true;
            }

            _ => {}
        }
    }

    if residues.is_empty() {
        return None;
    }

    if seq_name.is_empty() {
        seq_name = "unnamed".to_string();
    }

    let description = description_parts.join(" ");

    let mut seq = Sequence::new(seq_name, residues);
    seq.seq_type = seq_type;
    seq.features = features;
    // Store description in genbank.definition for reuse in write
    seq.genbank.definition = if description.is_empty() {
        None
    } else {
        Some(description)
    };

    Some(seq)
}

fn split_ft_key_location(s: &str) -> (String, String) {
    let mut it = s.splitn(2, |c: char| c.is_ascii_whitespace());
    let key = it.next().unwrap_or("").to_string();
    let loc = it.next().map(|v| v.trim().to_string()).unwrap_or_default();
    (key, loc)
}

fn flush_embl_qualifier(
    qname: &mut Option<String>,
    qval: &mut String,
    quals: &mut Vec<(String, String)>,
) {
    if let Some(k) = qname.take() {
        let v = std::mem::take(qval).trim_matches('"').trim().to_string();
        quals.push((k, v));
    } else {
        qval.clear();
    }
}

fn build_embl_feature(key: &str, location: &str, quals: &[(String, String)]) -> Option<Feature> {
    let feature_type = FeatureType::from_genbank_key(key);

    if matches!(feature_type, FeatureType::Source) {
        return None;
    }

    let (start, end, direction) = parse_embl_location(location.trim())?;

    let name = quals
        .iter()
        .find(|(k, _)| k == "gene" || k == "locus_tag" || k == "product")
        .map(|(_, v)| v.clone())
        .unwrap_or_else(|| key.to_string());

    let description = quals
        .iter()
        .map(|(k, v)| {
            if v.is_empty() {
                format!("/{k}")
            } else {
                format!("/{k}={v}")
            }
        })
        .collect::<Vec<_>>()
        .join("\n");

    let color = feature_type.default_color();
    let shape = match feature_type {
        FeatureType::Gene | FeatureType::CDS => FeatureShape::Arrow,
        _ => FeatureShape::Rectangle,
    };

    Some(Feature {
        name,
        description,
        start,
        end,
        shape,
        color,
        direction,
        feature_type,
    })
}

/// Parse an EMBL location string. EMBL uses the same syntax as GenBank
/// (1-based, `N..M`, `complement(...)`, `join(...)`).
fn parse_embl_location(loc: &str) -> Option<(usize, usize, FeatureDirection)> {
    if loc.is_empty() {
        return None;
    }

    if loc.starts_with("complement(") {
        let inner = loc.strip_prefix("complement(")?.trim_end_matches(')');
        let (s, e, _) = parse_embl_location(inner.trim())?;
        return Some((s, e, FeatureDirection::Reverse));
    }

    if loc.starts_with("join(") {
        let inner = loc.strip_prefix("join(")?.trim_end_matches(')');
        let mut min_start = usize::MAX;
        let mut max_end = 0usize;
        for part in inner.split(',') {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }
            if let Some((s, e, _)) = parse_embl_location(part) {
                min_start = min_start.min(s);
                max_end = max_end.max(e);
            }
        }
        if min_start == usize::MAX {
            return None;
        }
        return Some((min_start, max_end, FeatureDirection::Forward));
    }

    // Simple range: N..M
    let s = loc.trim_start_matches('<').trim_start_matches('>');
    if let Some(dot_pos) = s.find("..") {
        let s_part = s[..dot_pos]
            .trim_start_matches('<')
            .trim_start_matches('>')
            .trim();
        let e_part = s[dot_pos + 2..]
            .trim_start_matches('<')
            .trim_start_matches('>')
            .trim();
        let start: usize = s_part.parse().ok()?;
        let end: usize = e_part.parse().ok()?;
        // EMBL is 1-based inclusive; convert to 0-based inclusive.
        Some((
            start.saturating_sub(1),
            end.saturating_sub(1),
            FeatureDirection::Forward,
        ))
    } else {
        let pos: usize = s.trim().parse().ok()?;
        Some((
            pos.saturating_sub(1),
            pos.saturating_sub(1),
            FeatureDirection::None,
        ))
    }
}

pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    let content = to_string(aln)?;
    std::fs::write(path, content).map_err(FormatError::Io)
}

pub fn to_string(aln: &Alignment) -> Result<String> {
    let mut out = String::new();
    for seq in &aln.sequences {
        write_record(&mut out, seq)?;
        out.push_str("//\n");
    }
    Ok(out)
}

fn write_record(out: &mut String, seq: &Sequence) -> Result<()> {
    let len = seq.true_len();
    let mol_type = match seq.seq_type {
        SequenceType::Dna => "DNA",
        SequenceType::Rna => "RNA",
        SequenceType::Protein => "AA",
        _ => "DNA",
    };

    // ID line
    out.push_str(&format!(
        "ID   {}; SV 1; linear; {}; STD; UNK; {} BP.\n",
        seq.name, mol_type, len
    ));
    out.push_str("XX\n");

    // DE line (description)
    if let Some(def) = &seq.genbank.definition {
        out.push_str(&format!("DE   {}\n", def));
        out.push_str("XX\n");
    }

    // AC line
    if let Some(acc) = &seq.genbank.accession {
        out.push_str(&format!("AC   {};\n", acc));
        out.push_str("XX\n");
    }

    // FT lines for features
    if !seq.features.is_empty() {
        out.push_str("FH   Key             Location/Qualifiers\n");
        out.push_str("FH\n");
        for feat in &seq.features {
            let key = feature_type_to_embl_key(&feat.feature_type);
            // Location: 1-based
            let start = feat.start + 1;
            let end = feat.end + 1;
            let location = match feat.direction {
                FeatureDirection::Reverse => format!("complement({}..{})", start, end),
                _ => format!("{}..{}", start, end),
            };
            out.push_str(&format!("FT   {:<16}{}\n", key, location));
            if !feat.name.is_empty() && feat.name != key {
                out.push_str(&format!("FT                   /gene=\"{}\"\n", feat.name));
            }
        }
        out.push_str("XX\n");
    }

    // SQ line
    let raw: Vec<u8> = seq
        .residues
        .iter()
        .copied()
        .filter(|&b| !matches!(b, b'~' | b'.'))
        .collect();

    // Count nucleotide composition for SQ header
    let a = raw.iter().filter(|&&b| b == b'A' || b == b'a').count();
    let c = raw.iter().filter(|&&b| b == b'C' || b == b'c').count();
    let g = raw.iter().filter(|&&b| b == b'G' || b == b'g').count();
    let t = raw.iter().filter(|&&b| b == b'T' || b == b't').count();
    let other = raw.len().saturating_sub(a + c + g + t);

    out.push_str(&format!(
        "SQ   Sequence {} BP; {} A; {} C; {} G; {} T; {} other;\n",
        raw.len(),
        a,
        c,
        g,
        t,
        other
    ));

    // Write sequence in blocks of 60 chars (6 groups of 10), with position at end
    for (i, chunk) in raw.chunks(60).enumerate() {
        out.push_str("     ");
        for (j, sub) in chunk.chunks(10).enumerate() {
            if j > 0 {
                out.push(' ');
            }
            let s = std::str::from_utf8(sub).map_err(|e| FormatError::Parse(e.to_string()))?;
            out.push_str(&s.to_lowercase());
        }
        // Pad to 66 chars of sequence (6 groups of 10 + 5 spaces), then position
        let written_chars = chunk.len() + (chunk.len() / 10).saturating_sub(1).min(5);
        let pad = 66usize.saturating_sub(written_chars);
        for _ in 0..pad {
            out.push(' ');
        }
        out.push_str(&format!("{:>9}\n", (i * 60) + chunk.len()));
    }

    Ok(())
}

fn feature_type_to_embl_key(ft: &FeatureType) -> &'static str {
    match ft {
        FeatureType::Gene => "gene",
        FeatureType::CDS => "CDS",
        FeatureType::Exon => "exon",
        FeatureType::Intron => "intron",
        FeatureType::MRNA => "mRNA",
        FeatureType::TRNA => "tRNA",
        FeatureType::RRNA => "rRNA",
        FeatureType::RepeatRegion => "repeat_region",
        FeatureType::Source => "source",
        FeatureType::MiscFeature => "misc_feature",
        FeatureType::StemLoop => "stem_loop",
        FeatureType::Terminator => "terminator",
        FeatureType::PromotionalRegion => "promoter",
        _ => "misc_feature",
    }
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = "\
ID   AJ000001; SV 1; linear; DNA; STD; HUM; 929 BP.
XX
DE   Homo sapiens HBB gene for hemoglobin subunit beta.
XX
FT   source          1..929
FT                   /organism=\"Homo sapiens\"
FT   gene            1..929
FT                   /gene=\"HBB\"
FT   CDS             142..495
FT                   /gene=\"HBB\"
FT                   /product=\"hemoglobin subunit beta\"
XX
SQ   Sequence 20 BP; 5 A; 5 C; 5 G; 5 T; 0 other;
     acgtacgtac gtacgtacgt       20
//
";

    #[test]
    fn parses_one_record() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.seq_count(), 1);
    }

    #[test]
    fn correct_name() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.sequences[0].name, "AJ000001");
    }

    #[test]
    fn correct_description() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(
            aln.sequences[0].genbank.definition.as_deref(),
            Some("Homo sapiens HBB gene for hemoglobin subunit beta.")
        );
    }

    #[test]
    fn sequence_residues_collected() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.sequences[0].residues.len(), 20);
    }

    #[test]
    fn features_parsed() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        // source is skipped; gene and CDS remain
        assert_eq!(aln.sequences[0].features.len(), 2);
    }

    #[test]
    fn cds_location_0based() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        let cds = aln.sequences[0]
            .features
            .iter()
            .find(|f| matches!(f.feature_type, helixview_core::feature::FeatureType::CDS))
            .unwrap();
        assert_eq!(cds.start, 141); // 142 → 141
        assert_eq!(cds.end, 494); // 495 → 494
    }

    #[test]
    fn roundtrip_does_not_panic() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        let _ = super::to_string(&aln).unwrap();
    }
}
