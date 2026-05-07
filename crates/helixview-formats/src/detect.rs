use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Format {
    Fasta,
    GenBank,
    Clustal,
    Nbrf,
    Phylip,
    Nexus,
    Msf,
    Embl,
    Stockholm,
    Unknown,
}

/// Detect the file format from the file extension, falling back to
/// content sniffing.
pub fn detect_format(data: &[u8], path: &Path) -> Format {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "fasta" | "fa" | "fna" | "faa" | "ffn" | "frn" => return Format::Fasta,
        "gb" | "gbk" | "genbank" => return Format::GenBank,
        "bio" => return Format::GenBank, // BioEdit project files use GenBank-like format
        "aln" => return Format::Clustal,
        "pir" | "nbrf" => return Format::Nbrf,
        "phy" | "phylip" => return Format::Phylip,
        "nex" | "nexus" | "nxs" => return Format::Nexus,
        "msf" => return Format::Msf,
        "embl" | "dat" => return Format::Embl,
        "sto" | "stockholm" | "stk" => return Format::Stockholm,
        _ => {}
    }

    // Content sniffing: look at the first non-whitespace byte/line.
    let trimmed = data
        .iter()
        .position(|&b| !b.is_ascii_whitespace())
        .map(|i| &data[i..])
        .unwrap_or(data);

    // Stockholm: starts with "# STOCKHOLM"
    {
        let upper_start: Vec<u8> = trimmed
            .iter()
            .take(11)
            .map(|b| b.to_ascii_uppercase())
            .collect();
        if upper_start == b"# STOCKHOLM" {
            return Format::Stockholm;
        }
    }

    // NEXUS: starts with #NEXUS (case-insensitive)
    {
        let upper_start: Vec<u8> = trimmed
            .iter()
            .take(6)
            .map(|b| b.to_ascii_uppercase())
            .collect();
        if upper_start == b"#NEXUS" {
            return Format::Nexus;
        }
    }

    // Phylip: first line is exactly two integers separated by whitespace
    {
        let first_line_end = trimmed
            .iter()
            .position(|&b| b == b'\n')
            .unwrap_or(trimmed.len());
        let first_line = std::str::from_utf8(&trimmed[..first_line_end]).unwrap_or("");
        let parts: Vec<&str> = first_line.split_whitespace().collect();
        if parts.len() == 2
            && parts[0].parse::<usize>().is_ok()
            && parts[1].parse::<usize>().is_ok()
        {
            return Format::Phylip;
        }
    }

    // ClustalW content sniffing
    let upper_start: Vec<u8> = trimmed
        .iter()
        .take(7)
        .map(|b| b.to_ascii_uppercase())
        .collect();
    if upper_start.starts_with(b"CLUSTAL") {
        return Format::Clustal;
    }

    if trimmed.starts_with(b"LOCUS") || trimmed.starts_with(b"ORIGIN") {
        return Format::GenBank;
    }

    // EMBL: first line starts with "ID   "
    if trimmed.starts_with(b"ID   ") {
        return Format::Embl;
    }

    // MSF: starts with "!!" or an early line contains "MSF:"
    if trimmed.starts_with(b"!!") {
        return Format::Msf;
    }
    {
        let text_sample = std::str::from_utf8(data).unwrap_or("");
        let early: &str = &text_sample[..text_sample.len().min(2048)];
        if early.lines().take(10).any(|l| l.contains("MSF:")) {
            return Format::Msf;
        }
    }

    // NBRF/PIR: starts with '>' and first line contains ';'
    // FASTA also starts with '>' but does not have ';' on the header line
    if trimmed.starts_with(b">") {
        // Check if the first line contains a semicolon (PIR type prefix pattern)
        let first_line_end = trimmed
            .iter()
            .position(|&b| b == b'\n')
            .unwrap_or(trimmed.len());
        let first_line = &trimmed[..first_line_end];
        if first_line.contains(&b';') {
            return Format::Nbrf;
        }
        return Format::Fasta;
    }

    Format::Unknown
}
