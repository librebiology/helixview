use helixview_core::feature::{Feature, FeatureDirection, FeatureShape, FeatureType};

pub fn parse_features(raw: &str) -> Vec<Feature> {
    let mut features: Vec<Feature> = Vec::new();

    let mut cur_type: Option<String> = None;
    let mut cur_loc = String::new();
    let mut cur_quals: Vec<(String, String)> = Vec::new();
    let mut cur_qname: Option<String> = None;
    let mut cur_qval = String::new();
    let mut in_location = true;

    for line in raw.lines() {
        if line.trim_start().starts_with("Location") {
            continue;
        }

        let leading = line.bytes().take_while(|&b| b == b' ').count();
        if leading >= line.len() {
            continue;
        }

        let content = &line[leading..];

        if leading == 5 && !content.starts_with('/') {
            flush_qualifier(&mut cur_qname, &mut cur_qval, &mut cur_quals);
            if let Some(key) = cur_type.take() {
                if let Some(f) = build_feature(&key, &cur_loc, &cur_quals) {
                    features.push(f);
                }
            }
            cur_quals.clear();
            cur_loc.clear();
            in_location = true;
            let (key, loc) = split_key_location(content);
            cur_type = Some(key);
            cur_loc = loc;
        } else if leading >= 21 && content.starts_with('/') {
            flush_qualifier(&mut cur_qname, &mut cur_qval, &mut cur_quals);
            in_location = false;
            let body = &content[1..];
            if let Some(eq) = body.find('=') {
                let qkey = body[..eq].to_string();
                let val_start = body[eq + 1..].trim_start_matches('"');
                let clean = val_start.trim_end_matches('"').to_string();
                cur_qname = Some(qkey);
                cur_qval = clean;
            } else {
                cur_quals.push((body.to_string(), String::new()));
            }
        } else if leading >= 21 {
            let c = content.trim_end_matches('"').trim_start_matches('"');
            if in_location {
                cur_loc.push_str(c);
            } else {
                if !cur_qval.is_empty() {
                    cur_qval.push(' ');
                }
                cur_qval.push_str(c);
            }
        }
    }

    flush_qualifier(&mut cur_qname, &mut cur_qval, &mut cur_quals);
    if let Some(key) = cur_type.take() {
        if let Some(f) = build_feature(&key, &cur_loc, &cur_quals) {
            features.push(f);
        }
    }

    const CYCLE: [(u8, u8, u8); 10] = [
        (52, 120, 205),
        (220, 115, 25),
        (50, 180, 65),
        (215, 50, 50),
        (140, 50, 190),
        (5, 165, 165),
        (215, 50, 120),
        (189, 166, 12),
        (100, 160, 100),
        (180, 80, 40),
    ];
    let mut ci = 0usize;
    for feat in &mut features {
        if matches!(feat.feature_type, FeatureType::Source) {
            continue;
        }
        let (r, g, b) = CYCLE[ci % CYCLE.len()];
        feat.color = helixview_core::color::Color::from_u8(r, g, b);
        ci += 1;
    }

    features
}

fn flush_qualifier(
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

fn split_key_location(s: &str) -> (String, String) {
    let mut it = s.splitn(2, |c: char| c.is_ascii_whitespace());
    let key = it.next().unwrap_or("").to_string();
    let loc = it.next().map(|v| v.trim().to_string()).unwrap_or_default();
    (key, loc)
}

fn build_feature(key: &str, location: &str, quals: &[(String, String)]) -> Option<Feature> {
    let feature_type = FeatureType::from_genbank_key(key);
    if matches!(feature_type, FeatureType::Source) {
        return None;
    }

    let (start, end, direction) = parse_location(location.trim())?;

    // Prefer /gene > /locus_tag > /product > /note for display name.
    // /note carries restriction-site names (e.g. misc_feature + /note="HindIII").
    let name = quals
        .iter()
        .find(|(k, _)| k == "gene" || k == "locus_tag" || k == "product" || k == "note")
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
        color: feature_type.default_color(),
        direction,
        feature_type,
    })
}

fn parse_location(loc: &str) -> Option<(usize, usize, FeatureDirection)> {
    if loc.is_empty() {
        return None;
    }

    if loc.starts_with("complement(") {
        let inner = loc.strip_prefix("complement(")?.trim_end_matches(')');
        let (s, e, _) = parse_location(inner.trim())?;
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
            if let Some((s, e, _)) = parse_location(part) {
                min_start = min_start.min(s);
                max_end = max_end.max(e);
            }
        }
        if min_start == usize::MAX {
            return None;
        }
        return Some((min_start, max_end, FeatureDirection::Forward));
    }

    parse_simple_range(loc)
}

fn parse_simple_range(s: &str) -> Option<(usize, usize, FeatureDirection)> {
    let s = s.trim_start_matches('<').trim_start_matches('>');
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
        // GenBank is 1-based inclusive; convert to 0-based.
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

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
Location/Qualifiers
     source          1..929
                     /organism=\"Homo sapiens\"
                     /mol_type=\"mRNA\"
     gene            1..929
                     /gene=\"HBB\"
                     /db_xref=\"GeneID:3043\"
     CDS             142..495
                     /gene=\"HBB\"
                     /codon_start=1
                     /product=\"hemoglobin subunit beta\"
                     /translation=\"MVHLTPEEKSAV\"";

    #[test]
    fn basic_parse_count() {
        assert_eq!(parse_features(SAMPLE).len(), 2);
    }

    #[test]
    fn gene_name_from_qualifier() {
        assert!(parse_features(SAMPLE).iter().any(|f| f.name == "HBB"));
    }

    #[test]
    fn location_1based_to_0based() {
        let feats = parse_features(SAMPLE);
        let cds = feats
            .iter()
            .find(|f| matches!(f.feature_type, FeatureType::CDS))
            .unwrap();
        assert_eq!(cds.start, 141);
        assert_eq!(cds.end, 494);
    }

    #[test]
    fn complement_sets_direction() {
        let raw = "Location/Qualifiers\n     CDS             complement(10..90)\n                     /gene=\"TEST\"";
        let feats = parse_features(raw);
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].direction, FeatureDirection::Reverse);
        assert_eq!(feats[0].start, 9);
        assert_eq!(feats[0].end, 89);
    }

    #[test]
    fn join_takes_span() {
        let raw = "Location/Qualifiers\n     CDS             join(10..50,70..90)\n                     /gene=\"X\"";
        let feats = parse_features(raw);
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].start, 9);
        assert_eq!(feats[0].end, 89);
    }

    #[test]
    fn multiline_qualifier_value() {
        let raw = "\
Location/Qualifiers
     CDS             1..30
                     /translation=\"MKLPPTQFGE
                     TAFACR\"";
        let feats = parse_features(raw);
        assert_eq!(feats.len(), 1);
        assert!(feats[0].description.contains("translation"));
    }

    #[test]
    fn feature_type_mapped_correctly() {
        assert!(parse_features(SAMPLE)
            .iter()
            .any(|f| matches!(f.feature_type, FeatureType::CDS)));
    }
}
