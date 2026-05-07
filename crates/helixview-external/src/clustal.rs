//! External multiple-sequence aligner wrappers: ClustalW, ClustalOmega, MUSCLE, MAFFT.
//!
//! Each wrapper:
//!   1. Writes selected sequences to a temp FASTA file.
//!   2. Invokes the external binary with appropriate flags.
//!   3. Parses the resulting ClustalW `.aln` output (or MAFFT stdout).
//!   4. Returns the re-aligned `Alignment` or an error string.

use std::io::Write as _;
use std::process::Command;

use helixview_core::Alignment;
use helixview_formats::clustal;

/// Which external aligner to invoke.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExternalAligner {
    ClustalW,
    ClustalOmega,
    Muscle,
    Mafft,
}

impl ExternalAligner {
    /// Human-readable name for UI messages.
    pub fn display_name(self) -> &'static str {
        match self {
            Self::ClustalW => "ClustalW",
            Self::ClustalOmega => "Clustal Omega",
            Self::Muscle => "MUSCLE",
            Self::Mafft => "MAFFT",
        }
    }

    fn find_binary(self) -> Result<String, String> {
        let candidates: &[&str] = match self {
            Self::ClustalW => &["clustalw2", "clustalw"],
            Self::ClustalOmega => &["clustalo"],
            Self::Muscle => &["muscle"],
            Self::Mafft => &["mafft"],
        };
        for &bin in candidates {
            if which(bin) {
                return Ok(bin.to_owned());
            }
        }
        Err(format!(
            "{} binary not found in PATH. \
             Install it (e.g. `sudo pacman -S mafft` or `sudo apt install mafft`) and try again.",
            self.display_name(),
        ))
    }
}

/// Run an external aligner on the sequences at `indices` (all sequences when empty).
/// Returns the new aligned `Alignment` or an error string.
pub fn run_aligner(
    aln: &Alignment,
    indices: &[usize],
    aligner: ExternalAligner,
) -> Result<Alignment, String> {
    let binary = aligner.find_binary()?;

    let seqs: Vec<(String, Vec<u8>)> = if indices.is_empty() {
        aln.sequences.iter()
    } else {
        // lifetime trick: collect first, then iter
        let _filtered: Vec<_> = indices
            .iter()
            .filter_map(|&i| aln.sequences.get(i))
            .collect();
        // We need owned data; rebuild inline below.
        return run_aligner_seqs(
            indices
                .iter()
                .filter_map(|&i| aln.sequences.get(i))
                .map(|s| {
                    let raw: Vec<u8> = s
                        .residues
                        .iter()
                        .filter(|&&b| !helixview_core::sequence::is_gap(b))
                        .copied()
                        .collect();
                    (s.name.clone(), raw)
                })
                .collect(),
            &binary,
            aligner,
            aln.name.clone(),
        );
    }
    .map(|s| {
        let raw: Vec<u8> = s
            .residues
            .iter()
            .filter(|&&b| !helixview_core::sequence::is_gap(b))
            .copied()
            .collect();
        (s.name.clone(), raw)
    })
    .collect();

    run_aligner_seqs(seqs, &binary, aligner, aln.name.clone())
}

fn run_aligner_seqs(
    seqs: Vec<(String, Vec<u8>)>,
    binary: &str,
    aligner: ExternalAligner,
    aln_name: String,
) -> Result<Alignment, String> {
    if seqs.len() < 2 {
        return Err("Need at least 2 sequences to run a multiple alignment.".to_string());
    }

    let tmp_dir = std::env::temp_dir();
    let in_path = tmp_dir.join("helixview_align_in.fasta");
    let out_path = tmp_dir.join("helixview_align_out.aln");

    {
        let mut f =
            std::fs::File::create(&in_path).map_err(|e| format!("Cannot write temp file: {e}"))?;
        for (name, residues) in &seqs {
            writeln!(f, ">{name}").map_err(|e| e.to_string())?;
            // Write sequence in 70-char lines
            for chunk in residues.chunks(70) {
                writeln!(f, "{}", std::str::from_utf8(chunk).unwrap_or(""))
                    .map_err(|e| e.to_string())?;
            }
        }
    }

    let mut cmd = Command::new(binary);
    match aligner {
        ExternalAligner::ClustalW => {
            cmd.args([
                &format!("-INFILE={}", in_path.display()),
                &format!("-OUTFILE={}", out_path.display()),
                "-OUTPUT=CLUSTAL",
                "-OUTORDER=INPUT",
                "-QUIET",
            ]);
        }
        ExternalAligner::ClustalOmega => {
            cmd.args([
                "-i",
                in_path.to_str().unwrap_or(""),
                "-o",
                out_path.to_str().unwrap_or(""),
                "--outfmt=clu",
                "--force",
                "--quiet",
            ]);
        }
        ExternalAligner::Muscle => {
            // MUSCLE v5 syntax; fall back to v3 syntax on failure
            cmd.args([
                "-align",
                in_path.to_str().unwrap_or(""),
                "-output",
                out_path.to_str().unwrap_or(""),
            ]);
        }
        ExternalAligner::Mafft => {
            // MAFFT writes the aligned FASTA to stdout; stderr carries progress.
            // --auto: choose strategy based on data size
            // --clustalout: emit ClustalW format so we reuse the same parser
            // --quiet: suppress progress to stderr
            // --thread -1: use all available cores
            cmd.args([
                "--auto",
                "--clustalout",
                "--quiet",
                "--thread",
                "-1",
                in_path.to_str().unwrap_or(""),
            ]);
        }
    }

    // MAFFT writes to stdout; all others write to out_path.
    if aligner == ExternalAligner::Mafft {
        let output = cmd
            .output()
            .map_err(|e| format!("Failed to run MAFFT: {e}"))?;
        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(format!("MAFFT failed:\n{}", stderr.trim()));
        }
        let aln_str = std::str::from_utf8(&output.stdout)
            .map_err(|_| "MAFFT output is not valid UTF-8".to_string())?;
        return clustal::parse_str(aln_str, &aln_name)
            .map_err(|e| format!("Cannot parse MAFFT output: {e}"));
    }

    let output = cmd
        .output()
        .map_err(|e| format!("Failed to run {}: {e}", aligner.display_name()))?;

    if !output.status.success() {
        // MUSCLE v3 retry with old-style flags
        if aligner == ExternalAligner::Muscle {
            let output2 = Command::new(binary)
                .args([
                    "-in",
                    in_path.to_str().unwrap_or(""),
                    "-out",
                    out_path.to_str().unwrap_or(""),
                    "-clw",
                    "-quiet",
                ])
                .output()
                .map_err(|e| format!("Failed to run MUSCLE: {e}"))?;
            if !output2.status.success() {
                let stderr = String::from_utf8_lossy(&output2.stderr);
                return Err(format!("MUSCLE failed:\n{}", stderr.trim()));
            }
        } else {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(format!(
                "{} failed:\n{}",
                aligner.display_name(),
                stderr.trim(),
            ));
        }
    }

    let aln_bytes =
        std::fs::read(&out_path).map_err(|e| format!("Cannot read aligner output: {e}"))?;
    let aln_str = std::str::from_utf8(&aln_bytes)
        .map_err(|_| "Aligner output is not valid UTF-8".to_string())?;

    let result = clustal::parse_str(aln_str, &aln_name)
        .map_err(|e| format!("Cannot parse aligner output: {e}"))?;
    Ok(result)
}

fn which(bin: &str) -> bool {
    Command::new("which")
        .arg(bin)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}
