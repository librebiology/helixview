use crate::{FormatError, Result};
use helixview_core::{
    alignment::Alignment,
    sequence::{GenBankMeta, Sequence, SequenceType},
};
use std::path::Path;

pub fn parse_bytes(data: &[u8], name: &str) -> Result<Alignment> {
    let text =
        std::str::from_utf8(data).map_err(|e| FormatError::Parse(format!("Invalid UTF-8: {e}")))?;
    parse_str(text, name)
}

pub fn parse_str(text: &str, name: &str) -> Result<Alignment> {
    let mut aln = Alignment::new(name);
    for record in text.split("\n//") {
        let record = record.trim();
        if record.is_empty() {
            continue;
        }
        if let Some(seq) = parse_record(record) {
            aln.push(seq);
        }
    }
    if aln.is_empty() {
        return Err(FormatError::Parse(
            "No sequences found in GenBank data".into(),
        ));
    }
    Ok(aln)
}

fn parse_record(record: &str) -> Option<Sequence> {
    let mut meta = GenBankMeta::default();
    let mut seq_name = String::new();
    let mut residues = Vec::<u8>::new();
    let mut in_origin = false;
    let mut current_field: Option<String> = None;
    let mut current_value = String::new();

    for raw_line in record.lines() {
        // Normalise CRLF from HTTP responses before any index or length operations.
        let line = raw_line.trim_end_matches('\r');

        if line.starts_with("ORIGIN") {
            flush_field(&mut current_field, &mut current_value, &mut meta);
            in_origin = true;
            continue;
        }

        if in_origin {
            for b in line.bytes() {
                if b.is_ascii_alphabetic() {
                    residues.push(b);
                }
            }
            continue;
        }

        if line.len() > 12 && !line.starts_with(' ') {
            flush_field(&mut current_field, &mut current_value, &mut meta);
            let field = line[..12].trim().to_string();
            let value = line[12..].trim().to_string();
            if field == "LOCUS" {
                seq_name = value.split_whitespace().next().unwrap_or("").to_string();
            }
            current_field = Some(field);
            current_value = value;
        } else if line.starts_with(' ') && current_field.is_some() {
            current_value.push('\n');
            // FEATURES needs raw indentation to distinguish keys from qualifiers.
            if current_field.as_deref() == Some("FEATURES") {
                current_value.push_str(line);
            } else {
                current_value.push_str(line.trim());
            }
        }
    }
    flush_field(&mut current_field, &mut current_value, &mut meta);

    if residues.is_empty() {
        return None;
    }
    if seq_name.is_empty() {
        seq_name = "unnamed".to_string();
    }

    let mut seq = Sequence::new(seq_name, residues);
    seq.seq_type = SequenceType::Dna;
    if let Some(raw) = &meta.features_raw {
        seq.features = crate::features::parse_features(raw);
    }
    seq.genbank = meta;
    Some(seq)
}

fn flush_field(field: &mut Option<String>, value: &mut String, meta: &mut GenBankMeta) {
    if let Some(f) = field.take() {
        let v = std::mem::take(value).trim().to_string();
        match f.as_str() {
            "LOCUS" => meta.locus = Some(v),
            "DEFINITION" => meta.definition = Some(v),
            "ACCESSION" => meta.accession = Some(v),
            "VERSION" => meta.version = Some(v),
            "DBSOURCE" => meta.dbsource = Some(v),
            "KEYWORDS" => meta.keywords = Some(v),
            "SOURCE" => meta.source = Some(v),
            "COMMENT" => meta.comment = Some(v),
            "FEATURES" => meta.features_raw = Some(v),
            "REFERENCE" => meta.references.push(v),
            _ => {}
        }
    }
}

pub fn write_file(aln: &Alignment, path: &Path) -> Result<()> {
    std::fs::write(path, to_string(aln)?).map_err(FormatError::Io)
}

pub fn to_string(aln: &Alignment) -> Result<String> {
    let mut out = String::new();
    for seq in &aln.sequences {
        write_record(&mut out, seq);
        out.push_str("//\n");
    }
    Ok(out)
}

fn write_record(out: &mut String, seq: &Sequence) {
    let len = seq.true_len();
    let mol = match seq.seq_type {
        SequenceType::Dna => "DNA",
        SequenceType::Rna => "RNA",
        SequenceType::Protein => "AA ",
        _ => "   ",
    };
    out.push_str(&format!(
        "LOCUS       {:<16} {:>10} bp    {}    linear\n",
        seq.name, len, mol
    ));
    if let Some(def) = &seq.genbank.definition {
        out.push_str(&format!("DEFINITION  {}\n", def));
    }
    if let Some(acc) = &seq.genbank.accession {
        out.push_str(&format!("ACCESSION   {}\n", acc));
    }
    if let Some(ver) = &seq.genbank.version {
        out.push_str(&format!("VERSION     {}\n", ver));
    }
    out.push_str("ORIGIN\n");
    let raw: Vec<u8> = seq
        .residues
        .iter()
        .copied()
        .filter(|&b| !matches!(b, b'-' | b'~' | b'.'))
        .collect();
    for (i, chunk) in raw.chunks(60).enumerate() {
        out.push_str(&format!("{:>9} ", (i * 60) + 1));
        for (j, sub) in chunk.chunks(10).enumerate() {
            if j > 0 {
                out.push(' ');
            }
            out.push_str(std::str::from_utf8(sub).unwrap_or("?"));
        }
        out.push('\n');
    }
}

#[cfg(test)]
mod tests {
    const SAMPLE: &str = "\
LOCUS       seq1              12 bp    DNA    linear
ORIGIN
        1 atcgatcgatcg
//
LOCUS       seq2               8 bp    DNA    linear
ORIGIN
        1 aaaacccc
//";

    #[test]
    fn parse_seq_count_and_names() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.seq_count(), 2);
        assert_eq!(aln.sequences[0].name, "seq1");
        assert_eq!(aln.sequences[1].name, "seq2");
    }

    #[test]
    fn parse_residues_correct() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        assert_eq!(aln.sequences[0].true_len(), 12);
        assert_eq!(aln.sequences[1].true_len(), 8);
        assert!(aln.sequences[0].residues.starts_with(b"atcg"));
    }

    #[test]
    fn roundtrip_preserves_names_and_residues() {
        let aln = super::parse_str(SAMPLE, "test").unwrap();
        let out = super::to_string(&aln).unwrap();
        let aln2 = super::parse_str(&out, "test2").unwrap();
        assert_eq!(aln.seq_count(), aln2.seq_count());
        for (a, b) in aln.sequences.iter().zip(aln2.sequences.iter()) {
            assert_eq!(a.name, b.name);
            assert_eq!(a.residues, b.residues);
        }
    }
}
