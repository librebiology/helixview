//! Pairing-arcs strip — draws bracket arcs above alignment columns for top covarying pairs.
//!
//! Each arc connects two columns (col_a, col_b) with a quadratic bezier curve.
//! Color encodes the covariation score: blue (weak) → red (strong).

use std::cell::Cell;

use iced::{
    alignment, mouse,
    widget::canvas::{self, Frame, Path, Stroke, Text},
    Color, Font, Pixels, Point, Rectangle, Size,
};

use crate::app::Message;
use crate::widgets::grid::TITLE_W;

pub const ARCS_H: f32 = 60.0;

/// (col_a, col_b, cov_score) — both columns are 0-based alignment indices.
pub type PairingPair = (usize, usize, f64);

pub struct PairingArcsStrip {
    /// Pairs to draw, pre-filtered and sorted by score descending.
    pub pairs:      Vec<PairingPair>,
    pub scroll_col: usize,
    pub cell_w:     f32,
    /// Maximum cov_score in the full set (for normalised coloring).
    pub max_score:  f64,
}

pub struct PairingArcsState {
    cache:    canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for PairingArcsState {
    fn default() -> Self {
        Self { cache: canvas::Cache::default(), last_key: Cell::new(u64::MAX) }
    }
}

impl canvas::Program<Message> for PairingArcsStrip {
    type State = PairingArcsState;

    fn draw(
        &self,
        state:    &PairingArcsState,
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let key = (self.scroll_col as u64)
            .wrapping_add(self.cell_w.to_bits() as u64)
            .wrapping_add(self.pairs.len() as u64 * 0x9e3779b9)
            .wrapping_add(self.max_score.to_bits())
            .wrapping_add(bounds.width.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let pairs      = &self.pairs;
        let scroll_col = self.scroll_col;
        let cell_w     = self.cell_w;
        let max_score  = self.max_score;

        let geom = state.cache.draw(renderer, bounds.size(), |frame| {
            draw_arcs(frame, bounds.size(), pairs, scroll_col, cell_w, max_score);
        });

        vec![geom]
    }
}

fn draw_arcs(
    frame:      &mut Frame<iced::Renderer>,
    size:       Size,
    pairs:      &[PairingPair],
    scroll_col: usize,
    cell_w:     f32,
    max_score:  f64,
) {
    // White background with a faint bottom border.
    frame.fill_rectangle(Point::ORIGIN, size, Color::from_rgb(0.97, 0.97, 1.0));

    let base_y   = size.height - 4.0;   // arcs bottom out here
    let max_arc_h = size.height - 12.0; // tallest arc height
    let visible_w = size.width - TITLE_W;
    let n_visible = (visible_w / cell_w).ceil() as usize + 2;
    let max_s = max_score.max(1e-9);

    // Label: strip title on the left in the title margin.
    frame.fill_text(Text {
        content:  "Covariation".to_string(),
        position: Point::new(2.0, size.height * 0.5),
        color:    Color::from_rgb(0.40, 0.30, 0.55),
        size:     Pixels(9.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Left,
        vertical_alignment:   alignment::Vertical::Center,
        ..Text::default()
    });

    for &(col_a, col_b, score) in pairs {
        // Both endpoints must be visible (or close enough to clip nicely).
        let lo = scroll_col.saturating_sub(1);
        let hi = scroll_col + n_visible;
        if col_b < lo || col_a > hi { continue; }

        let x_a = TITLE_W + (col_a as f32 - scroll_col as f32 + 0.5) * cell_w;
        let x_b = TITLE_W + (col_b as f32 - scroll_col as f32 + 0.5) * cell_w;

        // Clamp to visible area so arcs partially off-screen look correct.
        let x_a = x_a.clamp(TITLE_W, size.width);
        let x_b = x_b.clamp(TITLE_W, size.width);

        let span  = (x_b - x_a).abs();
        // Arc height proportional to column distance, capped.
        let arc_h = (span * 0.4).min(max_arc_h);
        if arc_h < 3.0 { continue; }

        let t = (score / max_s) as f32;
        let color = arc_color(t);

        // Quadratic bezier: p0 → control → p1
        let p0  = Point::new(x_a, base_y);
        let p1  = Point::new(x_b, base_y);
        let cp  = Point::new((x_a + x_b) * 0.5, base_y - arc_h);

        let arc = Path::new(|b| {
            b.move_to(p0);
            b.quadratic_curve_to(cp, p1);
        });
        let stroke_w = (1.0 + t * 1.5).min(3.0);
        frame.stroke(&arc, Stroke::default().with_color(color).with_width(stroke_w));

        // Short vertical tick at each endpoint.
        let tick_h = 4.0;
        for &x in &[x_a, x_b] {
            let tick = Path::new(|b| {
                b.move_to(Point::new(x, base_y));
                b.line_to(Point::new(x, base_y - tick_h));
            });
            frame.stroke(&tick, Stroke::default().with_color(color).with_width(stroke_w));
        }

        // Score label near the arc peak (only if span is wide enough).
        if span >= cell_w * 4.0 {
            frame.fill_text(Text {
                content:  format!("{:.1}", score),
                position: Point::new((x_a + x_b) * 0.5, base_y - arc_h - 2.0),
                color:    darken(color, 0.20),
                size:     Pixels(8.0),
                font:     Font::MONOSPACE,
                horizontal_alignment: alignment::Horizontal::Center,
                vertical_alignment:   alignment::Vertical::Bottom,
                ..Text::default()
            });
        }
    }

    // Bottom baseline.
    let baseline = Path::new(|b| {
        b.move_to(Point::new(TITLE_W, base_y));
        b.line_to(Point::new(size.width, base_y));
    });
    frame.stroke(&baseline, Stroke::default()
        .with_color(Color::from_rgb(0.75, 0.72, 0.85))
        .with_width(1.0));
}

/// Blue → purple → red color scale for covariation score.
fn arc_color(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    if t < 0.5 {
        let s = t * 2.0;
        Color::from_rgb(0.20 + s * 0.30, 0.20 - s * 0.15, 0.85 - s * 0.15)
    } else {
        let s = (t - 0.5) * 2.0;
        Color::from_rgb(0.50 + s * 0.45, 0.05, 0.70 - s * 0.60)
    }
}

fn darken(c: Color, f: f32) -> Color {
    Color { r: (c.r * (1.0 - f)).max(0.0), g: (c.g * (1.0 - f)).max(0.0), b: (c.b * (1.0 - f)).max(0.0), a: 1.0 }
}
