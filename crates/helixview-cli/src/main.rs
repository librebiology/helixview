use std::path::PathBuf;
use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "helixview-cli", version, about = "Headless bioinformatics analysis")]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Convert an alignment between formats (detected by file extension).
    Convert {
        input:  PathBuf,
        output: PathBuf,
    },
    /// Print alignment statistics: sequence count, column count, mean pairwise identity.
    Stats {
        input: PathBuf,
    },
    /// Output per-column Shannon entropy as two-column TSV (col, entropy).
    Entropy {
        input: PathBuf,
        /// Skip columns with entropy below this threshold.
        #[arg(long, default_value_t = 0.0)]
        min_entropy: f64,
    },
    /// Print consensus sequence in FASTA format.
    Consensus {
        input: PathBuf,
        /// Consensus method: plurality or iupac.
        #[arg(long, default_value = "plurality")]
        method: String,
    },
    /// Run mutual information over all column pairs and print top N as TSV.
    Mi {
        input: PathBuf,
        /// Number of top pairs to print (0 = all).
        #[arg(long, default_value_t = 50)]
        top: usize,
        /// Minimum covariation score (wc_types × MI) to include.
        #[arg(long, default_value_t = 0.0)]
        min_cov: f64,
    },
    /// Print pairwise identity matrix as CSV.
    Identity {
        input: PathBuf,
        /// Output as percent (0–100) instead of fraction (0–1).
        #[arg(long, default_value_t = false)]
        percent: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Cmd::Convert { input, output } => cmd_convert(&input, &output),
        Cmd::Stats   { input }        => cmd_stats(&input),
        Cmd::Entropy { input, min_entropy } => cmd_entropy(&input, min_entropy),
        Cmd::Consensus { input, method }    => cmd_consensus(&input, &method),
        Cmd::Mi { input, top, min_cov }     => cmd_mi(&input, top, min_cov),
        Cmd::Identity { input, percent }    => cmd_identity(&input, percent),
    }
}

// ── convert ───────────────────────────────────────────────────────────────────

fn cmd_convert(input: &PathBuf, output: &PathBuf) -> Result<()> {
    let aln = helixview_formats::open(input)
        .with_context(|| format!("reading {}", input.display()))?;
    helixview_formats::save(&aln, output)
        .with_context(|| format!("writing {}", output.display()))?;
    eprintln!(
        "Converted {} sequences ({} cols) → {}",
        aln.seq_count(),
        aln.col_count(),
        output.display()
    );
    Ok(())
}

// ── stats ─────────────────────────────────────────────────────────────────────

fn cmd_stats(input: &PathBuf) -> Result<()> {
    let aln = helixview_formats::open(input)
        .with_context(|| format!("reading {}", input.display()))?;

    let n_seqs = aln.seq_count();
    let n_cols = aln.col_count();

    // Mean pairwise identity over all upper-triangle pairs.
    let mat = helixview_analysis::identity_matrix(&aln);
    let (sum, count) = mat.iter().enumerate()
        .flat_map(|(i, row)| row.iter().enumerate().map(move |(j, &v)| (i, j, v)))
        .filter(|(i, j, _)| i < j)
        .fold((0.0f64, 0usize), |(s, c), (_, _, v)| (s + v, c + 1));
    let mean_id = if count > 0 { sum / count as f64 } else { 1.0 };

    // Count gap-only columns.
    let gap_cols = (0..n_cols).filter(|&c| {
        aln.sequences.iter().all(|s| {
            s.residues.get(c).map_or(true, |&b| matches!(b, b'-' | b'.' | b'~'))
        })
    }).count();

    println!("File:          {}", input.display());
    println!("Sequences:     {}", n_seqs);
    println!("Columns:       {}", n_cols);
    println!("Gap-only cols: {}", gap_cols);
    println!("Mean identity: {:.1}%", mean_id * 100.0);
    Ok(())
}

// ── entropy ───────────────────────────────────────────────────────────────────

fn cmd_entropy(input: &PathBuf, min_entropy: f64) -> Result<()> {
    let aln = helixview_formats::open(input)
        .with_context(|| format!("reading {}", input.display()))?;

    let entropies = helixview_analysis::column_entropy(&aln);
    println!("col\tentropy");
    for (col, h) in entropies.iter().enumerate() {
        if *h >= min_entropy {
            println!("{}\t{:.6}", col + 1, h);
        }
    }
    Ok(())
}

// ── consensus ─────────────────────────────────────────────────────────────────

fn cmd_consensus(input: &PathBuf, method: &str) -> Result<()> {
    use helixview_analysis::ConsensusMethod;

    let aln = helixview_formats::open(input)
        .with_context(|| format!("reading {}", input.display()))?;

    let m = match method {
        "iupac"    => ConsensusMethod::Iupac { threshold: 0.25 },
        "plurality" => ConsensusMethod::Plurality,
        other      => anyhow::bail!("Unknown consensus method '{}' (use plurality or iupac)", other),
    };

    let cons = helixview_analysis::consensus(&aln, m);
    let name = input.file_stem().and_then(|s| s.to_str()).unwrap_or("consensus");
    println!(">{}_{}", name, method);
    // Wrap at 80 characters.
    for chunk in cons.chunks(80) {
        println!("{}", std::str::from_utf8(chunk).unwrap_or("?"));
    }
    Ok(())
}

// ── mi ────────────────────────────────────────────────────────────────────────

fn cmd_mi(input: &PathBuf, top: usize, min_cov: f64) -> Result<()> {
    let aln = helixview_formats::open(input)
        .with_context(|| format!("reading {}", input.display()))?;

    eprintln!("Running mutual information on {} sequences × {} columns…", aln.seq_count(), aln.col_count());
    // If min_cov filtering is active we can't know how many pairs will pass,
    // so request everything; otherwise requesting top*4 gives enough headroom.
    let fetch_n = if min_cov > 0.0 || top == 0 { usize::MAX } else { top * 4 };
    let result = helixview_analysis::mutual_information(&aln, 5, fetch_n);

    println!("col_a\tcol_b\tmi\twc_types\twc_frac\tcov_score");
    let pairs = result.top_covarying.iter()
        .filter(|p| p.cov_score >= min_cov);
    let pairs: Box<dyn Iterator<Item = _>> = if top == 0 {
        Box::new(pairs)
    } else {
        Box::new(pairs.take(top))
    };
    for p in pairs {
        println!(
            "{}\t{}\t{:.6}\t{}\t{:.4}\t{:.4}",
            p.col_a + 1, p.col_b + 1, p.mi, p.wc_types, p.wc_frac, p.cov_score
        );
    }
    Ok(())
}

// ── identity ─────────────────────────────────────────────────────────────────

fn cmd_identity(input: &PathBuf, percent: bool) -> Result<()> {
    let aln = helixview_formats::open(input)
        .with_context(|| format!("reading {}", input.display()))?;

    let mat = helixview_analysis::identity_matrix(&aln);

    // Header row: sequence names.
    let names: Vec<&str> = aln.sequences.iter().map(|s| s.name.as_str()).collect();
    print!("seq");
    for n in &names { print!("\t{}", n); }
    println!();

    for (i, row) in mat.iter().enumerate() {
        print!("{}", names[i]);
        for &v in row.iter() {
            if percent {
                print!("\t{:.1}", v * 100.0);
            } else {
                print!("\t{:.4}", v);
            }
        }
        println!();
    }
    Ok(())
}
