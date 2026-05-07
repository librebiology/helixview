//! Dot-plot view — sequence vs sequence comparison matrix.

use std::sync::Arc;
use std::cell::Cell;

use iced::{
    alignment, mouse,
    widget::{button, canvas, column, container, row, slider, text},
    Alignment as IAlignment, Background, Color, Element, Font,
    Length, Pixels, Point, Rectangle, Size,
};
use iced::widget::canvas::{Frame, Text as CanvasText};
use helixview_core::Alignment as SeqAlignment;
use helixview_core::sequence::is_gap;

use crate::app::Message;
use crate::theme::palette;

// ── Public entry point ────────────────────────────────────────────────────────

pub fn dot_plot_view<'a>(
    aln:       &'a Arc<SeqAlignment>,
    seq_a_idx: usize,
    seq_b_idx: usize,
    window:    usize,
    min_match: usize,
) -> Element<'a, Message> {
    let n = aln.seq_count();

    let name_a = aln.sequences.get(seq_a_idx)
        .map(|s| truncate(&s.name, 28))
        .unwrap_or_else(|| "(none)".to_string());
    let name_b = aln.sequences.get(seq_b_idx)
        .map(|s| truncate(&s.name, 28))
        .unwrap_or_else(|| "(none)".to_string());

    // Length display
    let len_a = aln.sequences.get(seq_a_idx)
        .map(|s| s.residues.iter().filter(|&&b| !is_gap(b)).count())
        .unwrap_or(0);
    let len_b = aln.sequences.get(seq_b_idx)
        .map(|s| s.residues.iter().filter(|&&b| !is_gap(b)).count())
        .unwrap_or(0);

    // ── Controls ──────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 10])
        .on_press(Message::CloseDotPlot);

    let title_lbl = text("Dot Plot").size(14).color(palette::TEXT);

    // Sequence A selector
    let prev_a = if seq_a_idx > 0 {
        button(text("<").size(11)).padding([3, 6])
            .on_press(Message::DotPlotSeqA(seq_a_idx - 1))
    } else {
        button(text("<").size(11)).padding([3, 6])
    };
    let next_a = if seq_a_idx + 1 < n {
        button(text(">").size(11)).padding([3, 6])
            .on_press(Message::DotPlotSeqA(seq_a_idx + 1))
    } else {
        button(text(">").size(11)).padding([3, 6])
    };
    let seq_a_ctrl = row![
        text("X:").size(12).color(palette::TEXT_DIM),
        prev_a,
        text(format!("{} ({}bp)", name_a, len_a)).size(12).color(palette::TEXT),
        next_a,
    ].spacing(4).align_y(IAlignment::Center);

    // Sequence B selector
    let prev_b = if seq_b_idx > 0 {
        button(text("<").size(11)).padding([3, 6])
            .on_press(Message::DotPlotSeqB(seq_b_idx - 1))
    } else {
        button(text("<").size(11)).padding([3, 6])
    };
    let next_b = if seq_b_idx + 1 < n {
        button(text(">").size(11)).padding([3, 6])
            .on_press(Message::DotPlotSeqB(seq_b_idx + 1))
    } else {
        button(text(">").size(11)).padding([3, 6])
    };
    let seq_b_ctrl = row![
        text("Y:").size(12).color(palette::TEXT_DIM),
        prev_b,
        text(format!("{} ({}bp)", name_b, len_b)).size(12).color(palette::TEXT),
        next_b,
    ].spacing(4).align_y(IAlignment::Center);

    // Window size slider (1–25)
    let win_slider = row![
        text("Window:").size(12).color(palette::TEXT_DIM),
        slider(1u32..=25, window as u32, |v| Message::DotPlotWindow(v as usize))
            .width(80),
        text(format!("{}", window)).size(12).color(palette::TEXT),
    ].spacing(6).align_y(IAlignment::Center);

    // Min-match slider (1..=window)
    let min_max = (window as u32).max(1);
    let cur_min = (min_match as u32).min(min_max);
    let min_slider = row![
        text("Min match:").size(12).color(palette::TEXT_DIM),
        slider(1u32..=min_max, cur_min, |v| Message::DotPlotMinMatch(v as usize))
            .width(80),
        text(format!("{}", min_match)).size(12).color(palette::TEXT),
    ].spacing(6).align_y(IAlignment::Center);

    let toolbar = container(
        row![
            back_btn,
            title_lbl,
            iced::widget::horizontal_space(),
            seq_a_ctrl,
            seq_b_ctrl,
            win_slider,
            min_slider,
        ]
        .spacing(12)
        .align_y(IAlignment::Center)
        .padding([5, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        ..Default::default()
    })
    .width(Length::Fill);

    // ── Dot-plot canvas ───────────────────────────────────────────────────────
    let dp = DotPlot {
        aln:       Arc::clone(aln),
        seq_a_idx,
        seq_b_idx,
        window,
        min_match,
    };

    let dp_canvas = canvas(dp)
        .width(Length::Fill)
        .height(Length::Fill);

    column![toolbar, dp_canvas]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Canvas program ────────────────────────────────────────────────────────────

struct DotPlot {
    aln:       Arc<SeqAlignment>,
    seq_a_idx: usize,
    seq_b_idx: usize,
    window:    usize,
    min_match: usize,
}

pub struct DotPlotState {
    cache:    canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for DotPlotState {
    fn default() -> Self {
        Self { cache: canvas::Cache::default(), last_key: Cell::new(u64::MAX) }
    }
}

impl canvas::Program<Message> for DotPlot {
    type State = DotPlotState;

    fn draw(
        &self,
        state:    &DotPlotState,
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let key = (Arc::as_ptr(&self.aln) as u64)
            .wrapping_add((self.seq_a_idx as u64).wrapping_mul(0x9e3779b97f4a7c15))
            .wrapping_add((self.seq_b_idx as u64).wrapping_mul(0x6c62272e07bb0142))
            .wrapping_add((self.window    as u64).wrapping_mul(0x517cc1b727220a95))
            .wrapping_add((self.min_match as u64).wrapping_mul(0x6b8b4567))
            .wrapping_add(bounds.width.to_bits()  as u64)
            .wrapping_add(bounds.height.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let geom = state.cache.draw(renderer, bounds.size(), |frame| {
            draw_dot_plot(
                frame, bounds.size(), &self.aln,
                self.seq_a_idx, self.seq_b_idx,
                self.window, self.min_match,
            );
        });

        vec![geom]
    }
}

// ── Drawing logic ─────────────────────────────────────────────────────────────

fn draw_dot_plot(
    frame:     &mut Frame<iced::Renderer>,
    size:      Size,
    aln:       &SeqAlignment,
    a_idx:     usize,
    b_idx:     usize,
    window:    usize,
    min_match: usize,
) {
    frame.fill_rectangle(Point::ORIGIN, size, Color::WHITE);

    let seq_a = match aln.sequences.get(a_idx) {
        Some(s) => s.residues.iter().filter(|&&b| !is_gap(b))
            .map(|&b| b.to_ascii_uppercase())
            .collect::<Vec<u8>>(),
        None => { draw_placeholder(frame, size, "No sequence A"); return; }
    };
    let seq_b = match aln.sequences.get(b_idx) {
        Some(s) => s.residues.iter().filter(|&&b| !is_gap(b))
            .map(|&b| b.to_ascii_uppercase())
            .collect::<Vec<u8>>(),
        None => { draw_placeholder(frame, size, "No sequence B"); return; }
    };

    if seq_a.is_empty() || seq_b.is_empty() {
        draw_placeholder(frame, size, "Empty sequence");
        return;
    }

    // Matrix plotting area (leave room for axis labels)
    const AXIS_PAD: f32 = 24.0;
    let plot_w = (size.width  - AXIS_PAD).max(1.0);
    let plot_h = (size.height - AXIS_PAD).max(1.0);

    let len_a = seq_a.len();
    let len_b = seq_b.len();

    // Resolution: cap at canvas pixels so we never draw more dots than pixels.
    let res_x = (plot_w as usize).min(len_a);
    let res_y = (plot_h as usize).min(len_b);
    let step_x = (len_a + res_x - 1) / res_x;
    let step_y = (len_b + res_y - 1) / res_y;
    let cell_w = plot_w / res_x as f32;
    let cell_h = plot_h / res_y as f32;

    let dot_color = Color::from_rgb(0.05, 0.30, 0.70);

    for pi in 0..res_x {
        for pj in 0..res_y {
            let i = pi * step_x;
            let j = pj * step_y;
            if window_matches(&seq_a, &seq_b, i, j, window * step_x, min_match) {
                let x = AXIS_PAD + pi as f32 * cell_w;
                let y = pj as f32 * cell_h;
                frame.fill_rectangle(
                    Point::new(x, y),
                    Size::new(cell_w.max(1.0), cell_h.max(1.0)),
                    dot_color,
                );
            }
        }
    }

    // Axis border
    frame.fill_rectangle(
        Point::new(AXIS_PAD, 0.0),
        Size::new(1.0, size.height),
        palette::BORDER,
    );
    frame.fill_rectangle(
        Point::new(AXIS_PAD, plot_h),
        Size::new(size.width, 1.0),
        palette::BORDER,
    );

    // Axis labels (bp counts)
    let lbl_y = format!("{} bp", len_b);
    frame.fill_text(CanvasText {
        content:  lbl_y,
        position: Point::new(2.0, plot_h * 0.5),
        color:    palette::TEXT_DIM,
        size:     Pixels(10.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Left,
        vertical_alignment:   alignment::Vertical::Center,
        ..CanvasText::default()
    });

    let lbl_x = format!("{} bp", len_a);
    frame.fill_text(CanvasText {
        content:  lbl_x,
        position: Point::new(AXIS_PAD + plot_w * 0.5, plot_h + 4.0),
        color:    palette::TEXT_DIM,
        size:     Pixels(10.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment:   alignment::Vertical::Top,
        ..CanvasText::default()
    });


}

/// Check whether seq_a[i..] and seq_b[j..] have at least `min_match` matches
/// in a window of length `window`.
#[inline]
fn window_matches(a: &[u8], b: &[u8], i: usize, j: usize, window: usize, min_match: usize) -> bool {
    let avail = window
        .min(a.len().saturating_sub(i))
        .min(b.len().saturating_sub(j));
    if avail == 0 { return false; }
    let needed = min_match.min(avail);
    let mut matches = 0usize;
    for k in 0..avail {
        if a[i + k] == b[j + k] {
            matches += 1;
            if matches >= needed { return true; }
        }
    }
    false
}

fn draw_placeholder(frame: &mut Frame<iced::Renderer>, size: Size, msg: &str) {
    frame.fill_text(CanvasText {
        content:  msg.to_string(),
        position: Point::new(size.width * 0.5, size.height * 0.5),
        color:    palette::TEXT_DIM,
        size:     Pixels(13.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment:   alignment::Vertical::Center,
        ..CanvasText::default()
    });
}

fn truncate(s: &str, max: usize) -> String {
    if s.chars().count() > max {
        let t: String = s.chars().take(max - 1).collect();
        format!("{t}\u{2026}")
    } else {
        s.to_owned()
    }
}
