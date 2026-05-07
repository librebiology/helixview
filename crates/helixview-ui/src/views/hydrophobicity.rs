//! Kyte-Doolittle hydrophobicity profile view — one line plot per sequence.

use std::cell::Cell;
use std::sync::Arc;

use iced::{
    mouse,
    widget::{button, canvas, column, container, row, slider, text},
    Alignment, Background, Color, Element, Length, Point, Rectangle, Size,
};
use iced::widget::canvas::Frame;
use helixview_core::Alignment as SeqAlignment;
use helixview_analysis::kyte_doolittle_profile;

use crate::app::Message;
use crate::theme::palette;

// ── Public entry point ────────────────────────────────────────────────────────

pub fn hydrophobicity_view<'a>(
    aln:    &'a Arc<SeqAlignment>,
    window: usize,
) -> Element<'a, Message> {
    let back_btn = button(text("< Back").size(13))
        .padding([4, 10])
        .on_press(Message::CloseHydrophobicity);

    let win_ctrl = row![
        text("Window:").size(12).color(palette::TEXT_DIM),
        slider(3u32..=25, window as u32, |v| Message::HydroWindow(v as usize)).width(110),
        text(format!("{window}")).size(12).color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let n_seq = aln.seq_count();
    let toolbar = container(
        row![
            back_btn,
            text("Hydrophobicity (Kyte-Doolittle)").size(14).color(palette::TEXT),
            iced::widget::horizontal_space(),
            win_ctrl,
            text(format!("{n_seq} sequences")).size(12).color(palette::TEXT_DIM),
        ]
        .spacing(12)
        .align_y(Alignment::Center)
        .padding([5, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        ..Default::default()
    })
    .width(Length::Fill);

    let plot: Element<Message> = canvas(HydroPlot {
        aln: Arc::clone(aln),
        window,
    })
    .width(Length::Fill)
    .height(Length::Fill)
    .into();

    column![toolbar, plot]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Canvas program ────────────────────────────────────────────────────────────

struct HydroPlot {
    aln:    Arc<SeqAlignment>,
    window: usize,
}

pub struct HydroPlotState {
    cache:    canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for HydroPlotState {
    fn default() -> Self {
        Self { cache: canvas::Cache::default(), last_key: Cell::new(u64::MAX) }
    }
}

impl canvas::Program<Message> for HydroPlot {
    type State = HydroPlotState;

    fn draw(
        &self,
        state:    &HydroPlotState,
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let key = (Arc::as_ptr(&self.aln) as u64)
            .wrapping_add(self.window as u64)
            .wrapping_add(bounds.width.to_bits() as u64)
            .wrapping_add(bounds.height.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let geom = state.cache.draw(renderer, bounds.size(), |frame| {
            draw_hydro(frame, bounds.size(), &self.aln, self.window);
        });
        vec![geom]
    }
}

// ── Drawing ───────────────────────────────────────────────────────────────────

const PAD_L: f32 = 52.0;
const PAD_B: f32 = 28.0;
const PAD_T: f32 = 14.0;
const PAD_R: f32 = 140.0; // legend

// Y range for KD scale (absolute min/max ≈ -4.5 to +4.5, give headroom).
const Y_MIN: f32 = -5.0;
const Y_MAX: f32 =  5.0;

const SEQ_COLORS: &[Color] = &[
    Color { r: 0.18, g: 0.45, b: 0.84, a: 1.0 }, // blue
    Color { r: 0.84, g: 0.18, b: 0.18, a: 1.0 }, // red
    Color { r: 0.12, g: 0.62, b: 0.28, a: 1.0 }, // green
    Color { r: 0.62, g: 0.18, b: 0.78, a: 1.0 }, // purple
    Color { r: 0.90, g: 0.55, b: 0.05, a: 1.0 }, // orange
    Color { r: 0.10, g: 0.62, b: 0.72, a: 1.0 }, // teal
    Color { r: 0.80, g: 0.35, b: 0.60, a: 1.0 }, // pink
    Color { r: 0.50, g: 0.40, b: 0.20, a: 1.0 }, // brown
];

fn seq_color(i: usize) -> Color {
    SEQ_COLORS[i % SEQ_COLORS.len()]
}

fn draw_hydro(frame: &mut Frame<iced::Renderer>, size: Size, aln: &SeqAlignment, window: usize) {
    frame.fill_rectangle(Point::ORIGIN, size, Color::WHITE);

    let pw = (size.width  - PAD_L - PAD_R).max(1.0);
    let ph = (size.height - PAD_T - PAD_B).max(1.0);

    // Compute profiles for all sequences.
    let profiles: Vec<Vec<f64>> = aln.sequences.iter()
        .map(|s| kyte_doolittle_profile(s, window))
        .collect();

    let max_len = profiles.iter().map(|p| p.len()).max().unwrap_or(1).max(1);

    // Grid lines and Y-axis labels.
    draw_axes(frame, pw, ph, max_len);

    // Sequence lines.
    for (i, profile) in profiles.iter().enumerate() {
        if profile.is_empty() { continue; }
        draw_line(frame, pw, ph, profile, seq_color(i));
    }

    // Zero line (bold).
    let y0 = value_to_y(0.0, ph);
    let path = canvas::Path::line(
        Point::new(PAD_L, PAD_T + y0),
        Point::new(PAD_L + pw, PAD_T + y0),
    );
    frame.stroke(&path, canvas::Stroke::default()
        .with_color(Color::from_rgb(0.4, 0.4, 0.4))
        .with_width(1.0));

    // Legend (right side).
    draw_legend(frame, size, aln);
}

fn draw_axes(frame: &mut Frame<iced::Renderer>, pw: f32, ph: f32, max_len: usize) {
    let grid_color = Color::from_rgba(0.80, 0.80, 0.84, 1.0);
    let text_color = Color::from_rgb(0.35, 0.35, 0.40);

    // Horizontal grid lines at -4, -2, 0, +2, +4.
    for &y_val in &[-4.0_f32, -2.0, 0.0, 2.0, 4.0] {
        let y = value_to_y(y_val, ph);
        let path = canvas::Path::line(
            Point::new(PAD_L, PAD_T + y),
            Point::new(PAD_L + pw, PAD_T + y),
        );
        let w = if y_val == 0.0 { 1.0 } else { 0.5 };
        frame.stroke(&path, canvas::Stroke::default().with_color(grid_color).with_width(w));

        // Y-axis label.
        frame.fill_text(canvas::Text {
            content: format!("{:+.0}", y_val),
            position: Point::new(PAD_L - 6.0, PAD_T + y),
            size: iced::Pixels(10.0),
            color: text_color,
            horizontal_alignment: iced::alignment::Horizontal::Right,
            vertical_alignment:   iced::alignment::Vertical::Center,
            ..Default::default()
        });
    }

    // Axis border.
    let border = canvas::Path::rectangle(
        Point::new(PAD_L, PAD_T),
        Size::new(pw, ph),
    );
    frame.stroke(&border, canvas::Stroke::default()
        .with_color(Color::from_rgb(0.60, 0.60, 0.65))
        .with_width(1.0));

    // X-axis tick labels (roughly 5 ticks).
    let n_ticks = 5.min(max_len);
    for i in 0..=n_ticks {
        let pos = i * max_len / n_ticks.max(1);
        let x = pos as f32 / max_len as f32 * pw + PAD_L;
        frame.fill_text(canvas::Text {
            content: format!("{}", pos + 1),
            position: Point::new(x, PAD_T + ph + 4.0),
            size: iced::Pixels(10.0),
            color: text_color,
            horizontal_alignment: iced::alignment::Horizontal::Center,
            vertical_alignment:   iced::alignment::Vertical::Top,
            ..Default::default()
        });
    }

    // Y-axis title.
    frame.fill_text(canvas::Text {
        content: "KD score".to_string(),
        position: Point::new(8.0, PAD_T + ph / 2.0),
        size: iced::Pixels(10.0),
        color: text_color,
        horizontal_alignment: iced::alignment::Horizontal::Left,
        vertical_alignment:   iced::alignment::Vertical::Center,
        ..Default::default()
    });
}

fn draw_line(frame: &mut Frame<iced::Renderer>, pw: f32, ph: f32, profile: &[f64], color: Color) {
    if profile.len() < 2 { return; }
    let n = profile.len();
    let mut builder = canvas::path::Builder::new();
    for (i, &val) in profile.iter().enumerate() {
        let x = PAD_L + i as f32 / (n - 1) as f32 * pw;
        let y = PAD_T + value_to_y(val as f32, ph);
        if i == 0 { builder.move_to(Point::new(x, y)); }
        else       { builder.line_to(Point::new(x, y)); }
    }
    let path = builder.build();
    frame.stroke(&path, canvas::Stroke::default().with_color(color).with_width(1.5));
}

fn draw_legend(frame: &mut Frame<iced::Renderer>, size: Size, aln: &SeqAlignment) {
    let x0 = size.width - PAD_R + 10.0;
    let mut y = PAD_T + 10.0;
    for (i, seq) in aln.sequences.iter().enumerate().take(12) {
        let color = seq_color(i);
        // Color swatch.
        frame.fill_rectangle(
            Point::new(x0, y),
            Size::new(14.0, 3.0),
            color,
        );
        // Name.
        let name = if seq.name.len() > 18 {
            format!("{}…", &seq.name[..17])
        } else {
            seq.name.clone()
        };
        frame.fill_text(canvas::Text {
            content: name,
            position: Point::new(x0 + 18.0, y + 1.5),
            size: iced::Pixels(10.0),
            color: Color::from_rgb(0.2, 0.2, 0.25),
            horizontal_alignment: iced::alignment::Horizontal::Left,
            vertical_alignment:   iced::alignment::Vertical::Center,
            ..Default::default()
        });
        y += 16.0;
        if y > size.height - PAD_B { break; }
    }
    if aln.seq_count() > 12 {
        frame.fill_text(canvas::Text {
            content: format!("… +{} more", aln.seq_count() - 12),
            position: Point::new(x0, y),
            size: iced::Pixels(10.0),
            color: Color::from_rgb(0.5, 0.5, 0.5),
            horizontal_alignment: iced::alignment::Horizontal::Left,
            vertical_alignment:   iced::alignment::Vertical::Top,
            ..Default::default()
        });
    }
}

#[inline]
fn value_to_y(v: f32, ph: f32) -> f32 {
    (1.0 - (v - Y_MIN) / (Y_MAX - Y_MIN)) * ph
}
