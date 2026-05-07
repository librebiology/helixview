//! Entropy ruler — a thin canvas strip showing Shannon entropy per alignment column.

use std::cell::Cell;
use std::sync::Arc;

use helixview_core::Alignment as SeqAlignment;
use iced::{
    alignment, mouse,
    widget::canvas::{self, Frame, Text},
    Color, Font, Pixels, Point, Rectangle, Size,
};

use crate::app::Message;
use crate::widgets::grid::TITLE_W;

const MAX_ENTROPY: f64 = 3.0;
pub const STRIP_H: f32 = 48.0;

pub struct EntropyStrip {
    pub aln: Arc<SeqAlignment>,
    pub scroll_col: usize,
    pub cell_w: f32,
}

/// Persistent state for the entropy strip — holds the draw cache.
pub struct EntropyState {
    cache: canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for EntropyState {
    fn default() -> Self {
        Self {
            cache: canvas::Cache::default(),
            last_key: Cell::new(u64::MAX),
        }
    }
}

impl canvas::Program<Message> for EntropyStrip {
    type State = EntropyState;

    fn draw(
        &self,
        state: &EntropyState,
        renderer: &iced::Renderer,
        _theme: &iced::Theme,
        bounds: Rectangle,
        _cursor: mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        // Invalidate if scroll position, alignment, or cell width changed.
        let key = (Arc::as_ptr(&self.aln) as u64)
            .wrapping_add(self.scroll_col as u64 * 0x9e3779b9)
            .wrapping_add(self.cell_w.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let aln = Arc::clone(&self.aln);
        let scroll_col = self.scroll_col;
        let cell_w = self.cell_w;
        const SNAP: f32 = 32.0;
        let snapped = Size {
            width: (bounds.width / SNAP).ceil() * SNAP,
            height: (bounds.height / SNAP).ceil() * SNAP,
        };
        let visible = bounds.size();
        let geom = state.cache.draw(renderer, snapped, |frame| {
            draw_entropy(frame, visible, &aln, scroll_col, cell_w);
        });

        vec![geom]
    }
}

fn draw_entropy(
    frame: &mut Frame<iced::Renderer>,
    size: Size,
    aln: &SeqAlignment,
    scroll_col: usize,
    cell_w: f32,
) {
    use crate::theme::palette;

    frame.fill_rectangle(Point::ORIGIN, size, palette::BG_PANEL);

    // Top separator — visually divides the entropy strip from the sequence grid above.
    frame.fill_rectangle(
        Point::ORIGIN,
        Size::new(size.width, 1.0),
        palette::BORDER_STRONG,
    );

    frame.fill_text(Text {
        content: "entropy".to_string(),
        position: Point::new(TITLE_W * 0.5, size.height * 0.5),
        color: palette::TEXT_DIM,
        size: Pixels(10.0),
        font: Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment: alignment::Vertical::Center,
        ..Text::default()
    });

    let ncols = aln.col_count();
    let bar_area = size.height - 4.0;
    let vis_cols = ((size.width - TITLE_W) / cell_w).ceil() as usize + 1;

    // Compute entropy only for the visible column range — O(n_seqs × vis_cols)
    // instead of O(n_seqs × total_cols) which freezes on large alignments.
    let end_col = (scroll_col + vis_cols).min(ncols);
    let n_seqs = aln
        .sequences
        .iter()
        .filter(|s| s.seq_type.is_sequence())
        .count() as f64;
    let vis_entropy: Vec<f64> = if n_seqs == 0.0 {
        vec![0.0; vis_cols]
    } else {
        (scroll_col..end_col)
            .map(|col| {
                let mut counts = [0u32; 256];
                for seq in &aln.sequences {
                    if !seq.seq_type.is_sequence() {
                        continue;
                    }
                    if let Some(&b) = seq.residues.get(col) {
                        counts[b as usize] += 1;
                    }
                }
                let mut h = 0.0f64;
                for &c in counts.iter() {
                    if c > 0 {
                        let p = c as f64 / n_seqs;
                        h -= p * p.ln();
                    }
                }
                h
            })
            .collect()
    };

    for ci in 0..vis_cols {
        let col = scroll_col + ci;
        if col >= ncols {
            break;
        }

        let h_val = vis_entropy.get(ci).copied().unwrap_or(0.0);
        let frac = (h_val / MAX_ENTROPY).clamp(0.0, 1.0) as f32;

        let color = lerp_color(
            Color::from_rgb(0.18, 0.72, 0.58),
            Color::from_rgb(0.93, 0.45, 0.12),
            frac,
        );

        let bar_h = (frac * bar_area).max(1.0);
        let x = TITLE_W + ci as f32 * cell_w;
        let y = size.height - 2.0 - bar_h;

        frame.fill_rectangle(
            Point::new(x + 0.5, y),
            Size::new((cell_w - 1.0).max(1.0), bar_h),
            color,
        );
    }

    frame.fill_rectangle(
        Point::new(0.0, size.height - 1.0),
        Size::new(size.width, 1.0),
        palette::BORDER,
    );
}

fn lerp_color(a: Color, b: Color, t: f32) -> Color {
    Color {
        r: a.r + (b.r - a.r) * t,
        g: a.g + (b.g - a.g) * t,
        b: a.b + (b.b - a.b) * t,
        a: 1.0,
    }
}
