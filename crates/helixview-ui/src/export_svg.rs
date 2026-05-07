//! SVG export of an alignment view.

use crate::color_table::ColorTable;
use helixview_core::Alignment;
use std::collections::BTreeSet;
use std::fmt::Write as FmtWrite;

const CELL_W: f32 = 12.0;
const CELL_H: f32 = 16.0;
const NAME_W: f32 = 150.0;
const FONT_SZ: f32 = 11.0;
const PAD_TOP: f32 = 4.0;

/// Generate an SVG string for a subset of an alignment.
///
/// `rows` and `cols` define which sequences/columns to include.
/// If empty, all rows/cols are used.
pub fn alignment_to_svg(
    aln: &Alignment,
    rows: &BTreeSet<usize>,
    col_start: usize,
    col_end: usize,
    color_table: &ColorTable,
) -> String {
    let seq_indices: Vec<usize> = if rows.is_empty() {
        (0..aln.seq_count()).collect()
    } else {
        rows.iter().copied().collect()
    };

    let ncols = aln.col_count();
    let c_start = col_start.min(ncols);
    let c_end = col_end.min(ncols);
    if c_start >= c_end || seq_indices.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string();
    }
    let num_cols = c_end - c_start;
    let num_rows = seq_indices.len();

    let width = NAME_W + num_cols as f32 * CELL_W;
    let height = PAD_TOP + num_rows as f32 * CELL_H;

    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{:.0}" height="{:.0}" font-family="monospace" font-size="{}">"#,
        width, height, FONT_SZ
    );

    // White background
    let _ = writeln!(
        svg,
        r#"  <rect width="{:.0}" height="{:.0}" fill="white"/>"#,
        width, height
    );

    let is_nuc = aln
        .sequences
        .first()
        .map(|s| s.seq_type.is_nucleic())
        .unwrap_or(true);

    for (row_i, &seq_idx) in seq_indices.iter().enumerate() {
        let seq = match aln.sequences.get(seq_idx) {
            Some(s) => s,
            None => continue,
        };
        let y_top = PAD_TOP + row_i as f32 * CELL_H;
        let text_y = y_top + CELL_H * 0.75; // baseline inside cell

        // Sequence name
        let name = escape_xml(&seq.name);
        let _ = writeln!(
            svg,
            "  <text x=\"2\" y=\"{:.1}\" fill=\"#333333\" font-size=\"{}\">{}</text>",
            text_y, FONT_SZ, name
        );

        // Residue cells
        for ci in 0..num_cols {
            let col = c_start + ci;
            let residue = seq.residues.get(col).copied().unwrap_or(b'-');
            let x = NAME_W + ci as f32 * CELL_W;

            let color = if is_nuc {
                color_table.nuc_color(residue)
            } else {
                color_table.aa_color(residue)
            };
            let (r8, g8, b8) = (
                (color.r * 255.0) as u8,
                (color.g * 255.0) as u8,
                (color.b * 255.0) as u8,
            );

            let is_gap = matches!(residue, b'-' | b'.' | b'~');
            if !is_gap {
                let _ = writeln!(
                    svg,
                    "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#{:02X}{:02X}{:02X}\"/>",
                    x, y_top, CELL_W, CELL_H, r8, g8, b8
                );
            }

            // Character label (skip gaps and very narrow cells)
            let ch = char::from(residue.to_ascii_uppercase());
            if !is_gap {
                let text_x = x + CELL_W * 0.5;
                // Use white text on dark backgrounds, dark on light
                let brightness = 0.299 * color.r + 0.587 * color.g + 0.114 * color.b;
                let text_color = if brightness < 0.5 { "white" } else { "#111111" };
                let _ = writeln!(
                    svg,
                    r#"  <text x="{:.1}" y="{:.1}" text-anchor="middle" fill="{}">{}</text>"#,
                    text_x, text_y, text_color, ch
                );
            }
        }
    }

    svg.push_str("</svg>\n");
    svg
}

fn escape_xml(s: &str) -> String {
    s.chars()
        .flat_map(|c| match c {
            '&' => "&amp;".chars().collect::<Vec<_>>(),
            '<' => "&lt;".chars().collect(),
            '>' => "&gt;".chars().collect(),
            '"' => "&quot;".chars().collect(),
            '\'' => "&#39;".chars().collect(),
            c => vec![c],
        })
        .collect()
}
