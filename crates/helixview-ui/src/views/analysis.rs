//! Analysis side panel — sequence stats table and column identity bar chart.

use std::sync::Arc;

use iced::{
    alignment, mouse,
    widget::{button, canvas, column, container, row, scrollable, text},
    Background, Border, Color, Element, Font, Length, Pixels, Point, Rectangle, Size,
};
use iced::widget::canvas::{Frame, Text as CanvasText};
use helixview_core::{self, Alignment};
use helixview_core::sequence::{is_gap, SequenceType};
use helixview_analysis::{kyte_doolittle_profile, calculate_tm, wallace_tm};

use crate::app::Message;
use crate::theme::palette;

// ── Public entry point ────────────────────────────────────────────────────────

/// Returns a fixed-width (280 px) analysis panel showing per-sequence stats
/// and a column-identity bar chart.
pub fn analysis_panel<'a>(
    aln: &'a Arc<Alignment>,
) -> Element<'a, Message> {

    // Detect molecule type from the first sequence.
    // If seq_type is explicitly Protein, or residues contain amino-acid-only
    // letters (F, L, I, P, Q, E, W, Y, H, K, R, D, N, M), treat as protein.
    let is_protein = aln.sequences.first().map(|s| {
        s.seq_type == SequenceType::Protein
            || s.residues.iter().any(|&b| matches!(
                b.to_ascii_uppercase(),
                b'F' | b'L' | b'I' | b'P' | b'Q' | b'E' | b'W' | b'Y'
                    | b'H' | b'K' | b'R' | b'D' | b'N' | b'M'
            ))
    }).unwrap_or(false);

    // Compute overall GC% for the alignment
    let (mut gc, mut atotal) = (0u64, 0u64);
    for seq in &aln.sequences {
        for &b in &seq.residues {
            let b = b.to_ascii_uppercase();
            if !is_gap(b) {
                atotal += 1;
                if matches!(b, b'G' | b'C') { gc += 1; }
            }
        }
    }
    let gc_pct = if atotal > 0 { gc as f32 / atotal as f32 * 100.0 } else { 0.0 };

    // Hydrophobicity: show profile for first sequence.
    let hydro_seq: Option<helixview_core::Sequence> = aln.sequences.first().cloned();

    let content = column![
        section_header("Sequence Stats"),
        stats_section(aln),
        section_header("Column Conservation"),
        identity_canvas(aln),
        section_header("Alignment Summary"),
        summary_section(aln, gc_pct, is_protein),
        section_header("Hydrophobicity (seq 1)"),
        hydro_canvas(hydro_seq.clone()),
        section_header("Oligo Tm (seq 1)"),
        tm_section(&hydro_seq),
        identity_matrix_btn(),
    ]
    .width(Length::Fill);

    container(content)
        .width(280)
        .height(Length::Fill)
        .style(|_| container::Style {
            background: Some(Background::Color(palette::BG_PANEL)),
            border: Border {
                color: palette::BORDER,
                width: 1.0,
                radius: 0.0.into(),
            },
            ..Default::default()
        })
        .into()
}

// ── Section header ────────────────────────────────────────────────────────────

fn section_header<'a>(label: &'a str) -> Element<'a, Message> {
    container(
        text(label)
            .size(12)
            .color(palette::TEXT)
            .font(Font::MONOSPACE),
    )
    .padding([5, 8])
    .width(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        border: Border {
            color: palette::BORDER,
            width: 0.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .into()
}

// ── Sequence stats table ──────────────────────────────────────────────────────

fn stats_section<'a>(aln: &'a Arc<Alignment>) -> Element<'a, Message> {
    const MAX_VISIBLE: usize = 20;

    let total_cols = aln.col_count();
    let seqs = &aln.sequences;
    let shown = seqs.len().min(MAX_VISIBLE);

    let mut rows: Vec<Element<'a, Message>> = Vec::with_capacity(shown + 1);

    for seq in seqs.iter().take(MAX_VISIBLE) {
        let gap_count = seq.residues.iter().filter(|&&b| is_gap(b)).count();
        let residue_len = seq.residues.len() - gap_count;
        let gap_pct = if total_cols > 0 {
            gap_count as f32 / total_cols as f32 * 100.0
        } else {
            0.0
        };

        // Truncate name to 18 chars
        let name_display: String = if seq.name.chars().count() > 18 {
            let truncated: String = seq.name.chars().take(17).collect();
            format!("{truncated}\u{2026}") // …
        } else {
            seq.name.clone()
        };

        let name_cell = text(name_display)
            .size(11)
            .color(palette::TEXT)
            .font(Font::MONOSPACE)
            .width(Length::Fill);

        let len_cell = text(format!("{}aa", residue_len))
            .size(11)
            .color(palette::TEXT_DIM)
            .font(Font::MONOSPACE);

        let gap_cell = text(format!("{:.0}%", gap_pct))
            .size(11)
            .color(palette::TEXT_DIM)
            .font(Font::MONOSPACE);

        let seq_row = container(
            row![name_cell, len_cell, gap_cell]
                .spacing(6)
                .align_y(alignment::Vertical::Center),
        )
        .padding([3, 8])
        .width(Length::Fill)
        .into();

        rows.push(seq_row);
    }

    // "… and N more" note
    if seqs.len() > MAX_VISIBLE {
        let extra = seqs.len() - MAX_VISIBLE;
        let note = container(
            text(format!("\u{2026} and {} more", extra))
                .size(11)
                .color(palette::TEXT_DIM)
                .font(Font::MONOSPACE),
        )
        .padding([3, 8])
        .width(Length::Fill)
        .into();
        rows.push(note);
    }

    let inner = column(rows).width(Length::Fill);

    scrollable(inner)
        .width(Length::Fill)
        .height(Length::Shrink)
        .into()
}

// ── Identity bar chart ────────────────────────────────────────────────────────

fn identity_canvas<'a>(aln: &'a Arc<Alignment>) -> Element<'a, Message> {
    canvas(IdentityBar { aln: Arc::clone(aln) })
        .width(280)
        .height(80)
        .into()
}

// ── Alignment summary ─────────────────────────────────────────────────────────

fn summary_section<'a>(aln: &'a Arc<Alignment>, gc_pct: f32, is_protein: bool) -> Element<'a, Message> {
    let seqs = aln.seq_count();
    let cols = aln.col_count();
    let gc_label = if is_protein {
        format!("Gly+Cys: {gc_pct:.1}%")
    } else {
        format!("GC: {gc_pct:.1}%")
    };
    let info = text(format!(
        "{seqs} seqs × {cols} cols\n{gc_label}"
    ))
    .size(11)
    .color(palette::TEXT)
    .font(Font::MONOSPACE);

    container(info)
        .padding([6, 8])
        .width(Length::Fill)
        .into()
}

// ── Hydrophobicity profile canvas ─────────────────────────────────────────────

fn hydro_canvas<'a>(seq: Option<helixview_core::Sequence>) -> Element<'a, Message> {
    canvas(HydroPlot { seq })
        .width(280)
        .height(60)
        .into()
}

struct HydroPlot {
    seq: Option<helixview_core::Sequence>,
}

impl canvas::Program<Message> for HydroPlot {
    type State = ();

    fn draw(
        &self,
        _state:   &(),
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let mut frame = Frame::new(renderer, bounds.size());
        frame.fill_rectangle(Point::ORIGIN, bounds.size(), palette::BG_PANEL);

        let seq = match &self.seq {
            Some(s) => s,
            None    => return vec![frame.into_geometry()],
        };

        const WINDOW: usize = 7;
        let profile = kyte_doolittle_profile(seq, WINDOW);
        if profile.is_empty() {
            return vec![frame.into_geometry()];
        }

        let w = bounds.width as f64;
        let h = (bounds.height - 14.0) as f64;
        let min_v = profile.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_v = profile.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = (max_v - min_v).max(0.001);

        let n = profile.len();
        let bar_w = (w / n as f64).max(1.0);

        for (i, &v) in profile.iter().enumerate() {
            let norm = (v - min_v) / range;
            let bar_h = (norm * h).max(1.0);
            let color = if v > 0.0 {
                Color::from_rgb(0.10, 0.55, 0.85)
            } else {
                Color::from_rgb(0.85, 0.28, 0.18)
            };
            frame.fill_rectangle(
                Point::new((i as f64 * bar_w) as f32, (h - bar_h) as f32),
                Size::new(((bar_w - 0.5).max(0.5)) as f32, bar_h as f32),
                color,
            );
        }

        // Zero line
        let zero_y = ((0.0_f64 - min_v) / range * h).clamp(0.0, h);
        let mut x = 0.0f32;
        while (x as f64) < w {
            frame.fill_rectangle(
                Point::new(x, (h - zero_y) as f32),
                Size::new(3.0_f32.min(bounds.width - x), 1.0),
                Color { r: 0.4, g: 0.4, b: 0.4, a: 0.6 },
            );
            x += 6.0;
        }

        frame.fill_text(CanvasText {
            content:  format!("KD window={WINDOW}"),
            position: Point::new(4.0, 2.0),
            color:    palette::TEXT_DIM,
            size:     Pixels(9.5),
            font:     Font::MONOSPACE,
            horizontal_alignment: alignment::Horizontal::Left,
            vertical_alignment:   alignment::Vertical::Top,
            ..CanvasText::default()
        });

        vec![frame.into_geometry()]
    }
}

// ── Identity matrix button ────────────────────────────────────────────────────

fn identity_matrix_btn<'a>() -> Element<'a, Message> {
    container(
        button(text("Show Identity Matrix…").size(12))
            .padding([5, 10])
            .on_press(Message::ShowIdentityMatrix),
    )
    .padding([8, 8])
    .width(Length::Fill)
    .into()
}

// ── Oligo Tm section ──────────────────────────────────────────────────────────

fn tm_section<'a>(seq: &Option<helixview_core::Sequence>) -> Element<'a, Message> {
    let info = if let Some(s) = seq {
        let raw = &s.residues;
        if raw.len() < 2 {
            "Sequence too short for Tm".to_string()
        } else if raw.len() <= 13 {
            // Wallace rule for very short oligos
            let tm = wallace_tm(raw);
            format!("Wallace rule: {tm:.1} °C\n(length {}, rule of 2+4)", raw.len())
        } else {
            match calculate_tm(raw) {
                Some(result) => format!(
                    "SantaLucia Tm: {:.1} °C\nGC: {:.1}%  |  len: {} nt\nΔH: {:.1} kcal/mol",
                    result.tm_celsius,
                    result.gc_fraction * 100.0,
                    result.length,
                    result.delta_h,
                ),
                None => "Tm: N/A".to_string(),
            }
        }
    } else {
        "No sequence loaded".to_string()
    };

    container(
        text(info)
            .size(11)
            .color(palette::TEXT)
            .font(Font::MONOSPACE),
    )
    .padding([6, 8])
    .width(Length::Fill)
    .into()
}

/// Canvas program that draws a per-column identity bar histogram.
struct IdentityBar {
    aln: Arc<Alignment>,
}

impl canvas::Program<Message> for IdentityBar {
    type State = ();

    fn draw(
        &self,
        _state:   &(),
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let mut frame = Frame::new(renderer, bounds.size());
        let w = bounds.width;
        let h = bounds.height;

        // Background
        frame.fill_rectangle(Point::ORIGIN, bounds.size(), palette::BG_PANEL);

        let n_cols = self.aln.col_count();
        if n_cols == 0 {
            return vec![frame.into_geometry()];
        }

        // Sample columns to fit within the canvas width (max 280 bars).
        let max_bars = w.floor() as usize;
        let step = if n_cols > max_bars {
            (n_cols + max_bars - 1) / max_bars   // ceiling division
        } else {
            1
        };

        let num_bars = (n_cols + step - 1) / step;
        let bar_w = (w / num_bars as f32).max(1.0);
        let chart_h = h - 10.0; // leave room for label

        for i in 0..num_bars {
            let col = i * step;
            if col >= n_cols { break; }
            let identity = self.aln.column_identity(col);

            let color = if identity >= 0.8 {
                Color::from_rgb(0.15, 0.72, 0.33)
            } else if identity >= 0.5 {
                Color::from_rgb(0.85, 0.72, 0.10)
            } else {
                Color::from_rgb(0.82, 0.22, 0.18)
            };

            let bar_h = (identity * chart_h).max(1.0);
            let x = i as f32 * bar_w;
            let y = chart_h - bar_h;

            frame.fill_rectangle(
                Point::new(x, y),
                Size::new((bar_w - 0.5).max(0.5), bar_h),
                color,
            );
        }

        // 50% reference line (dashed approximation via short segments)
        let ref_y = chart_h * 0.5;
        let dash_len = 4.0f32;
        let gap_len  = 3.0f32;
        let mut x = 0.0f32;
        while x < w {
            let seg_w = dash_len.min(w - x);
            frame.fill_rectangle(
                Point::new(x, ref_y),
                Size::new(seg_w, 1.0),
                Color { r: 0.55, g: 0.55, b: 0.60, a: 0.7 },
            );
            x += dash_len + gap_len;
        }

        // "Identity" label in the top-left
        frame.fill_text(CanvasText {
            content:  "Identity".to_string(),
            position: Point::new(4.0, 3.0),
            color:    palette::TEXT_DIM,
            size:     Pixels(10.0),
            font:     Font::MONOSPACE,
            horizontal_alignment: alignment::Horizontal::Left,
            vertical_alignment:   alignment::Vertical::Top,
            ..CanvasText::default()
        });

        vec![frame.into_geometry()]
    }
}
