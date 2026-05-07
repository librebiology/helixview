//! Formatted alignment text export.

use helixview_core::Alignment as SeqAlignment;

#[derive(Debug, Clone)]
pub struct TextExportOptions {
    pub residues_per_row: usize,
    pub title_chars: usize,
    pub show_ruler: bool,
    pub number_lines: bool,
}

impl Default for TextExportOptions {
    fn default() -> Self {
        Self {
            residues_per_row: 60,
            title_chars: 16,
            show_ruler: true,
            number_lines: true,
        }
    }
}

/// Format the full alignment as a plain-text block string.
pub fn format_alignment(aln: &SeqAlignment, opts: &TextExportOptions) -> String {
    if aln.seq_count() == 0 {
        return String::new();
    }

    let rpr = opts.residues_per_row.max(10);
    let tc = opts.title_chars.max(4);
    let ncols = aln.col_count();

    // Pre-truncate/pad titles.
    let titles: Vec<String> = aln
        .sequences
        .iter()
        .map(|s| {
            if s.name.len() > tc {
                format!("{:.prec$}", s.name, prec = tc)
            } else {
                format!("{:<width$}", s.name, width = tc)
            }
        })
        .collect();

    let num_width = if opts.number_lines {
        ncols.to_string().len() + 1
    } else {
        0
    };

    let mut out = String::with_capacity(ncols * aln.seq_count() * 2);
    let mut col = 0usize;

    while col < ncols {
        let end = (col + rpr).min(ncols);

        // Ruler row.
        if opts.show_ruler {
            // Left padding (title width + 1 space).
            for _ in 0..(tc + 1) {
                out.push(' ');
            }
            let mut ruler_pos = col;
            let mut ruler = String::new();
            while ruler_pos < end {
                let label = (ruler_pos + 1).to_string();
                // Pad to next 10-unit boundary.
                let next_tick = ((ruler_pos / 10) + 1) * 10;
                let gap = next_tick.min(end) - ruler_pos;
                if ruler_pos % 10 == 0 {
                    // Position label, right-padded to fill until next tick.
                    ruler.push_str(&label);
                    let pad = gap.saturating_sub(label.len());
                    for _ in 0..pad {
                        ruler.push(' ');
                    }
                } else {
                    let remaining = 10 - (ruler_pos % 10);
                    let to_fill = remaining.min(end - ruler_pos);
                    for _ in 0..to_fill {
                        ruler.push(' ');
                    }
                }
                ruler_pos = next_tick.min(end);
            }
            out.push_str(&ruler);
            out.push('\n');
        }

        // Sequence rows.
        for (seq_idx, seq) in aln.sequences.iter().enumerate() {
            out.push_str(&titles[seq_idx]);
            out.push(' ');
            // Residues for this block.
            for res in &seq.residues[col..end] {
                out.push(*res as char);
            }
            // End position.
            if opts.number_lines {
                let true_pos = seq.residues[..end]
                    .iter()
                    .filter(|&&b| !matches!(b, b'-' | b'.' | b'~'))
                    .count();
                out.push_str(&format!("  {:>width$}", true_pos, width = num_width));
            }
            out.push('\n');
        }

        out.push('\n');
        col = end;
    }

    out
}
