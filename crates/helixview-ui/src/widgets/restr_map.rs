//! Restriction map strip — a thin canvas strip showing enzyme cut sites.
//!
//! Only meaningful for DNA sequences; call sites once per alignment load and
//! cache them via `RestrMapState`.

use std::cell::Cell;
use std::sync::Arc;

use iced::{
    alignment,
    mouse,
    widget::canvas::{self, Frame, Text},
    Color, Font, Pixels, Point, Rectangle, Size,
};
use helixview_core::Alignment as SeqAlignment;
use helixview_core::sequence::is_gap;
use helixview_analysis::{find_sites, ENZYMES};

use crate::app::Message;
use crate::widgets::grid::TITLE_W;

pub const RESTR_H: f32 = 42.0;

/// Twenty visually distinct colours for enzyme labels / ticks.
const PALETTE: [Color; 20] = [
    Color { r: 0.85, g: 0.15, b: 0.15, a: 1.0 }, // 0 EcoRI   red
    Color { r: 0.10, g: 0.60, b: 0.10, a: 1.0 }, // 1 BamHI   green
    Color { r: 0.10, g: 0.30, b: 0.85, a: 1.0 }, // 2 HindIII blue
    Color { r: 0.90, g: 0.50, b: 0.05, a: 1.0 }, // 3 XbaI    orange
    Color { r: 0.55, g: 0.10, b: 0.70, a: 1.0 }, // 4 SalI    purple
    Color { r: 0.05, g: 0.65, b: 0.65, a: 1.0 }, // 5 KpnI    teal
    Color { r: 0.65, g: 0.35, b: 0.05, a: 1.0 }, // 6 SacI    brown
    Color { r: 0.65, g: 0.05, b: 0.40, a: 1.0 }, // 7 EcoRV   magenta
    Color { r: 0.20, g: 0.55, b: 0.25, a: 1.0 }, // 8 NcoI    forest
    Color { r: 0.75, g: 0.65, b: 0.05, a: 1.0 }, // 9 NdeI    gold
    Color { r: 0.15, g: 0.50, b: 0.75, a: 1.0 }, // 10 BglII  steel
    Color { r: 0.85, g: 0.30, b: 0.60, a: 1.0 }, // 11 ClaI   pink
    Color { r: 0.35, g: 0.65, b: 0.35, a: 1.0 }, // 12 NheI   sage
    Color { r: 0.60, g: 0.25, b: 0.75, a: 1.0 }, // 13 SpeI   violet
    Color { r: 0.25, g: 0.70, b: 0.55, a: 1.0 }, // 14 XhoI   cyan-green
    Color { r: 0.75, g: 0.40, b: 0.20, a: 1.0 }, // 15 PstI   terra
    Color { r: 0.45, g: 0.20, b: 0.60, a: 1.0 }, // 16 SphI   plum
    Color { r: 0.20, g: 0.55, b: 0.50, a: 1.0 }, // 17 SmaI   jade
    Color { r: 0.70, g: 0.15, b: 0.15, a: 1.0 }, // 18 MluI   crimson
    Color { r: 0.30, g: 0.30, b: 0.70, a: 1.0 }, // 19 AvrII  slate
];

pub struct RestrMapStrip {
    pub aln:        Arc<SeqAlignment>,
    pub scroll_col: usize,
    pub cell_w:     f32,
    /// Sequence row to scan (default 0; first selected).
    pub seq_idx:    usize,
}

pub struct RestrMapState {
    cache:    canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for RestrMapState {
    fn default() -> Self {
        Self { cache: canvas::Cache::default(), last_key: Cell::new(u64::MAX) }
    }
}

impl canvas::Program<Message> for RestrMapStrip {
    type State = RestrMapState;

    fn draw(
        &self,
        state:    &RestrMapState,
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let key = (Arc::as_ptr(&self.aln) as u64)
            .wrapping_add((self.seq_idx   as u64).wrapping_mul(0x9e3779b97f4a7c15))
            .wrapping_add((self.scroll_col as u64).wrapping_mul(0x6c62272e07bb0142))
            .wrapping_add(self.cell_w.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let aln        = Arc::clone(&self.aln);
        let seq_idx    = self.seq_idx;
        let scroll_col = self.scroll_col;
        let cell_w     = self.cell_w;

        const SNAP: f32 = 32.0;
        let snapped = Size {
            width:  (bounds.width  / SNAP).ceil() * SNAP,
            height: (bounds.height / SNAP).ceil() * SNAP,
        };
        let geom = state.cache.draw(renderer, snapped, |frame| {
            draw_restr_map(frame, snapped, &aln, seq_idx, scroll_col, cell_w);
        });

        vec![geom]
    }
}

fn draw_restr_map(
    frame:      &mut Frame<iced::Renderer>,
    size:       Size,
    aln:        &SeqAlignment,
    seq_idx:    usize,
    scroll_col: usize,
    cell_w:     f32,
) {
    use crate::theme::palette;

    frame.fill_rectangle(Point::ORIGIN, size, palette::BG_PANEL);

    // "RE map" label on the left gutter
    frame.fill_text(Text {
        content:  "RE map".to_string(),
        position: Point::new(TITLE_W * 0.5, size.height * 0.5),
        color:    palette::TEXT_DIM,
        size:     Pixels(10.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment:   alignment::Vertical::Center,
        ..Text::default()
    });

    let seq = match aln.sequences.get(seq_idx) {
        Some(s) => s,
        None    => return,
    };

    // Build a map: alignment_col → true_pos, and its inverse true_pos → alignment_col.
    // We need true_pos → alignment_col to go from site position to pixel x.
    let mut true_to_col: Vec<usize> = Vec::with_capacity(seq.residues.len());
    for (col, &b) in seq.residues.iter().enumerate() {
        if !is_gap(b) {
            true_to_col.push(col);
        }
    }

    // Build raw (ungapped) sequence for scanning
    let raw: Vec<u8> = seq.residues.iter()
        .filter(|&&b| !is_gap(b))
        .map(|&b| b.to_ascii_uppercase())
        .collect();

    let sites = find_sites(&raw);
    if sites.is_empty() {
        // No sites found — show a note
        frame.fill_text(Text {
            content:  "no sites".to_string(),
            position: Point::new(TITLE_W + 10.0, size.height * 0.5),
            color:    palette::TEXT_DIM,
            size:     Pixels(9.0),
            font:     Font::MONOSPACE,
            horizontal_alignment: alignment::Horizontal::Left,
            vertical_alignment:   alignment::Vertical::Center,
            ..Text::default()
        });
        frame.fill_rectangle(
            Point::new(0.0, size.height - 1.0),
            Size::new(size.width, 1.0),
            palette::BORDER,
        );
        return;
    }

    let vis_cols  = ((size.width - TITLE_W) / cell_w).ceil() as usize + 1;
    let col_start = scroll_col;
    let col_end   = scroll_col + vis_cols;

    // Tick geometry
    let tick_top    = 4.0f32;
    let tick_bottom = size.height - 14.0;
    let tick_h      = tick_bottom - tick_top;

    for site in &sites {
        let aln_col = match true_to_col.get(site.true_pos) {
            Some(&c) => c,
            None     => continue,
        };
        if aln_col < col_start || aln_col >= col_end { continue; }

        let vis_idx = aln_col - col_start;
        let x = TITLE_W + vis_idx as f32 * cell_w + cell_w * 0.5;

        let color = PALETTE[ENZYMES[site.enzyme_idx].color_idx % PALETTE.len()];

        // Vertical tick line
        frame.fill_rectangle(
            Point::new(x - 0.5, tick_top),
            Size::new(1.5, tick_h),
            color,
        );

        // Enzyme name label
        let name = ENZYMES[site.enzyme_idx].name;
        frame.fill_text(Text {
            content:  name.to_string(),
            position: Point::new(x, tick_bottom + 2.0),
            color,
            size:     Pixels(8.5),
            font:     Font::MONOSPACE,
            horizontal_alignment: alignment::Horizontal::Center,
            vertical_alignment:   alignment::Vertical::Top,
            ..Text::default()
        });
    }

    // Bottom border line
    frame.fill_rectangle(
        Point::new(0.0, size.height - 1.0),
        Size::new(size.width, 1.0),
        palette::BORDER,
    );
}
