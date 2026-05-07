pub mod abi;
pub mod clustal;
pub mod detect;
pub mod embl;
pub mod fasta;
pub mod features;
pub mod genbank;
pub mod msf;
pub mod nbrf;
pub mod newick;
pub mod nexus;
pub mod phylip;
pub mod stockholm;

use helixview_core::Alignment;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FormatError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
    #[error("Unknown or unsupported file format")]
    UnknownFormat,
    #[error("Empty file")]
    EmptyFile,
}

pub type Result<T> = std::result::Result<T, FormatError>;

/// Auto-detect format from file extension and content, then parse.
pub fn open(path: &Path) -> Result<Alignment> {
    let data = std::fs::read(path)?;
    if data.is_empty() {
        return Err(FormatError::EmptyFile);
    }

    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("alignment")
        .to_string();

    match detect::detect_format(&data, path) {
        detect::Format::Fasta => fasta::parse_bytes(&data, &name),
        detect::Format::GenBank => genbank::parse_bytes(&data, &name),
        detect::Format::Clustal => clustal::parse_bytes(&data, &name),
        detect::Format::Nbrf => nbrf::parse_bytes(&data, &name),
        detect::Format::Phylip => phylip::parse_bytes(&data, &name),
        detect::Format::Nexus => nexus::parse_bytes(&data, &name),
        detect::Format::Msf => msf::parse_bytes(&data, &name),
        detect::Format::Embl => embl::parse_bytes(&data, &name),
        detect::Format::Stockholm => stockholm::parse_bytes(&data, &name),
        detect::Format::Unknown => {
            // Try FASTA first (most common), then GenBank.
            fasta::parse_bytes(&data, &name)
                .or_else(|_| genbank::parse_bytes(&data, &name))
                .map_err(|_| FormatError::UnknownFormat)
        }
    }
}

/// Write an alignment to a file, inferring format from extension.
pub fn save(aln: &Alignment, path: &Path) -> Result<()> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();
    match ext.as_str() {
        "fasta" | "fa" | "fna" | "faa" | "ffn" | "frn" => fasta::write_file(aln, path),
        "gb" | "gbk" | "genbank" => genbank::write_file(aln, path),
        "aln" => clustal::write_file(aln, path),
        "pir" | "nbrf" => nbrf::write_file(aln, path),
        "phy" | "phylip" => phylip::write_file(aln, path),
        "nex" | "nexus" | "nxs" => nexus::write_file(aln, path),
        "msf" => msf::write_file(aln, path),
        "embl" | "dat" => embl::write_file(aln, path),
        _ => fasta::write_file(aln, path),
    }
}
