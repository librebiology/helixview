//! Publication-quality shaded alignment graphic export (SVG).

use std::collections::BTreeSet;
use std::fmt::Write as FmtWrite;
use std::sync::Arc;

use helixview_core::Alignment as SeqAlignment;
use crate::app::ColorScheme;
use crate::color_table::ColorTable;

const CELL_W:     f32   = 10.0;
const CELL_H:     f32   = 14.0;
const NAME_W:     f32   = 140.0;
const FONT_SZ:    f32   = 10.0;
const RULER_H:    f32   = 14.0;
const PAD_L:      f32   = 4.0;
const PAD_TOP:    f32   = 4.0;
const RULER_TICK: usize = 10;

pub fn shaded_alignment_svg(
    aln:            &Arc<SeqAlignment>,
    rows:           &BTreeSet<usize>,
    col_start:      usize,
    col_end:        usize,
    color_table:    &ColorTable,
    scheme:         &ColorScheme,
    threshold:      f32,
    show_ruler:     bool,
    show_consensus: bool,
) -> String {
    let ncols   = aln.col_count();
    let c_start = col_start.min(ncols);
    let c_end   = col_end.min(ncols);
    if c_start >= c_end {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".into();
    }
    let seq_indices: Vec<usize> = if rows.is_empty() {
        (0..aln.seq_count()).collect()
    } else {
        rows.iter().copied().filter(|&i| i < aln.seq_count()).collect()
    };
    if seq_indices.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".into();
    }

    let num_cols = c_end - c_start;
    let num_rows = seq_indices.len();
    let is_nuc   = aln.sequences.first().map(|s| s.seq_type.is_nucleic()).unwrap_or(true);

    let entropies: Vec<f32> = aln.column_entropy().iter().map(|&e| e as f32).collect();

    let ruler_offset    = if show_ruler    { RULER_H } else { 0.0 };
    let consensus_extra = if show_consensus { 1     } else { 0   };
    let total_rows_f    = (num_rows + consensus_extra) as f32;

    let width  = PAD_L + NAME_W + num_cols as f32 * CELL_W + PAD_L;
    let height = PAD_TOP + ruler_offset + total_rows_f * CELL_H + PAD_TOP;

    let mut svg = String::new();
    let _ = writeln!(svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w:.0}\" height=\"{h:.0}\" \
         font-family=\"monospace\" font-size=\"{fs}\">",
        w = width, h = height, fs = FONT_SZ,
    );
    let _ = writeln!(svg,
        "  <rect width=\"{w:.0}\" height=\"{h:.0}\" fill=\"white\"/>",
        w = width, h = height,
    );

    // ── Ruler ──────────────────────────────────────────────────────────────────
    if show_ruler {
        for ci in 0..num_cols {
            let col = c_start + ci;
            if col % RULER_TICK == 0 {
                let x  = PAD_L + NAME_W + ci as f32 * CELL_W;
                let y0 = PAD_TOP;
                let y1 = PAD_TOP + RULER_H * 0.5;
                let yt = PAD_TOP + RULER_H - 2.0;
                let _ = writeln!(svg,
                    "  <line x1=\"{x:.1}\" y1=\"{y0:.1}\" x2=\"{x:.1}\" y2=\"{y1:.1}\" \
                     stroke=\"#999999\" stroke-width=\"0.5\"/>",
                );
                let _ = writeln!(svg,
                    "  <text x=\"{xt:.1}\" y=\"{yt:.1}\" fill=\"#555555\" font-size=\"8\">{pos}</text>",
                    xt = x + 1.0, yt = yt, pos = col + 1,
                );
            }
        }
    }

    // ── Sequence rows ──────────────────────────────────────────────────────────
    for (row_i, &seq_idx) in seq_indices.iter().enumerate() {
        let seq    = match aln.sequences.get(seq_idx) { Some(s) => s, None => continue };
        let y_top  = PAD_TOP + ruler_offset + row_i as f32 * CELL_H;
        let text_y = y_top + CELL_H * 0.78;

        let _ = writeln!(svg,
            "  <text x=\"{x:.1}\" y=\"{y:.1}\" fill=\"#222222\">{n}</text>",
            x = PAD_L, y = text_y, n = escape_xml(&seq.name),
        );

        for ci in 0..num_cols {
            let col     = c_start + ci;
            let residue = seq.residues.get(col).copied().unwrap_or(b'-');
            if matches!(residue, b'-' | b'.' | b'~') { continue; }

            let x = PAD_L + NAME_W + ci as f32 * CELL_W;
            let base = if is_nuc { color_table.nuc_color(residue) }
                       else      { color_table.aa_color(residue)  };

            let bg = cell_color(base, col, scheme, threshold, &entropies, aln);
            if let Some((r8, g8, b8)) = bg {
                let _ = writeln!(svg,
                    "  <rect x=\"{x:.1}\" y=\"{yt:.1}\" width=\"{cw:.1}\" height=\"{ch:.1}\" \
                     fill=\"#{r:02X}{g:02X}{b:02X}\"/>",
                    x = x, yt = y_top, cw = CELL_W, ch = CELL_H,
                    r = r8, g = g8, b = b8,
                );
                let ch  = char::from(residue.to_ascii_uppercase());
                let lum = 0.299 * r8 as f32 / 255.0
                        + 0.587 * g8 as f32 / 255.0
                        + 0.114 * b8 as f32 / 255.0;
                let fg  = if lum < 0.5 { "white" } else { "#111111" };
                let _ = writeln!(svg,
                    "  <text x=\"{tx:.1}\" y=\"{ty:.1}\" text-anchor=\"middle\" fill=\"{fg}\">{ch}</text>",
                    tx = x + CELL_W * 0.5, ty = text_y, fg = fg, ch = ch,
                );
            }
        }
    }

    // ── Consensus row ──────────────────────────────────────────────────────────
    if show_consensus {
        let y_top  = PAD_TOP + ruler_offset + num_rows as f32 * CELL_H;
        let text_y = y_top + CELL_H * 0.78;
        let _ = writeln!(svg,
            "  <text x=\"{x:.1}\" y=\"{y:.1}\" fill=\"#222222\" font-weight=\"bold\">Consensus</text>",
            x = PAD_L, y = text_y,
        );
        for ci in 0..num_cols {
            let col     = c_start + ci;
            let residue = aln.column_consensus_residue(col);
            if residue == b'-' { continue; }
            let x    = PAD_L + NAME_W + ci as f32 * CELL_W;
            let base = if is_nuc { color_table.nuc_color(residue) }
                       else      { color_table.aa_color(residue)  };
            let alpha = aln.column_identity(col).max(0.3);
            let (r8, g8, b8) = blend_white(base, alpha);
            let _ = writeln!(svg,
                "  <rect x=\"{x:.1}\" y=\"{yt:.1}\" width=\"{cw:.1}\" height=\"{ch:.1}\" \
                 fill=\"#{r:02X}{g:02X}{b:02X}\"/>",
                x = x, yt = y_top, cw = CELL_W, ch = CELL_H,
                r = r8, g = g8, b = b8,
            );
            let ch  = char::from(residue.to_ascii_uppercase());
            let lum = 0.299 * r8 as f32 / 255.0
                    + 0.587 * g8 as f32 / 255.0
                    + 0.114 * b8 as f32 / 255.0;
            let fg  = if lum < 0.5 { "white" } else { "#111111" };
            let _ = writeln!(svg,
                "  <text x=\"{tx:.1}\" y=\"{ty:.1}\" text-anchor=\"middle\" fill=\"{fg}\">{ch}</text>",
                tx = x + CELL_W * 0.5, ty = text_y, fg = fg, ch = ch,
            );
        }
    }

    svg.push_str("</svg>\n");
    svg
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn cell_color(
    base:      iced::Color,
    col:       usize,
    scheme:    &ColorScheme,
    threshold: f32,
    entropies: &[f32],
    aln:       &SeqAlignment,
) -> Option<(u8, u8, u8)> {
    match scheme {
        ColorScheme::Plain => None,
        ColorScheme::ResidueType => Some(to_u8(base)),
        ColorScheme::Identity => {
            if aln.column_identity(col) >= threshold { Some(to_u8(base)) } else { None }
        }
        ColorScheme::Strength => {
            let max_e = if aln.sequences.first()
                .map(|s| s.seq_type.is_nucleic()).unwrap_or(true)
            { 5_f32.ln() } else { 21_f32.ln() };
            let e     = entropies.get(col).copied().unwrap_or(0.0).max(0.0);
            let alpha = (1.0 - e / max_e).max(0.0).powf(0.7);
            if alpha < 0.05 { None } else { Some(blend_white(base, alpha)) }
        }
    }
}

fn to_u8(c: iced::Color) -> (u8, u8, u8) {
    (
        (c.r.clamp(0.0, 1.0) * 255.0) as u8,
        (c.g.clamp(0.0, 1.0) * 255.0) as u8,
        (c.b.clamp(0.0, 1.0) * 255.0) as u8,
    )
}

fn blend_white(c: iced::Color, alpha: f32) -> (u8, u8, u8) {
    let r = c.r * alpha + (1.0 - alpha);
    let g = c.g * alpha + (1.0 - alpha);
    let b = c.b * alpha + (1.0 - alpha);
    (
        (r.clamp(0.0, 1.0) * 255.0) as u8,
        (g.clamp(0.0, 1.0) * 255.0) as u8,
        (b.clamp(0.0, 1.0) * 255.0) as u8,
    )
}

fn escape_xml(s: &str) -> String {
    s.chars().flat_map(|c| match c {
        '&'  => "&amp;".chars().collect::<Vec<_>>(),
        '<'  => "&lt;".chars().collect(),
        '>'  => "&gt;".chars().collect(),
        '"'  => "&quot;".chars().collect(),
        '\'' => "&#39;".chars().collect(),
        c    => vec![c],
    }).collect()
}
