use std::cell::Cell;
use std::collections::{BTreeSet, HashSet};
use std::sync::Arc;

use helixview_core::color::Color as CoreColor;
use helixview_core::feature::FeatureType;
use helixview_core::{Alignment as SeqAlignment, SequenceType};
use iced::{
    alignment, mouse,
    widget::canvas::{self, Frame, Path, Text},
    Color, Font, Pixels, Point, Rectangle, Size,
};

use crate::app::{ColorScheme, Message};
use crate::color_table::ColorTable;

// ── Layout constants ──────────────────────────────────────────────────────────

pub const TITLE_W: f32 = 200.0;
pub const RULER_H: f32 = 22.0;
pub const CELL_W: f32 = 10.0;
pub const CELL_H: f32 = 22.0;
pub const FONT_SIZE: f32 = 12.5;
/// Each feature tier is 4 px tall; 3 tiers = 12 px total stripe area.
const FEATURE_LAYER_H: f32 = 4.0;
const FEATURE_TIERS: usize = 3;

/// Identity threshold below which a column is rendered un-colored.
const IDENTITY_THRESHOLD: f32 = 0.50;

// ── AlignmentGrid ─────────────────────────────────────────────────────────────

/// Per-interaction state for the alignment grid canvas.
pub struct GridState {
    drag_from: Option<usize>,
    drag_to: Option<usize>,
    is_dragging: bool,
    col_drag_start: Option<usize>,
    /// Cached geometry for the main grid (residues, ruler, titles, features).
    /// Cleared when `cached_version` differs from `AlignmentGrid::version`.
    grid_cache: canvas::Cache,
    /// Version stamp seen when the cache was last filled.
    cached_version: Cell<u64>,
    /// For double-click detection: (row_idx, click_time).
    last_title_click: Option<(usize, std::time::Instant)>,
}

impl Default for GridState {
    fn default() -> Self {
        Self {
            drag_from: None,
            drag_to: None,
            is_dragging: false,
            col_drag_start: None,
            grid_cache: canvas::Cache::default(),
            cached_version: Cell::new(u64::MAX),
            last_title_click: None,
        }
    }
}

impl std::fmt::Debug for GridState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GridState")
            .field("drag_from", &self.drag_from)
            .field("drag_to", &self.drag_to)
            .field("is_dragging", &self.is_dragging)
            .field("col_drag_start", &self.col_drag_start)
            .finish_non_exhaustive()
    }
}

pub struct AlignmentGrid {
    pub aln: Arc<SeqAlignment>,
    pub scroll_col: usize,
    pub scroll_row: usize,
    pub selected: BTreeSet<usize>,
    pub scheme: ColorScheme,
    pub ctrl_held: bool,
    pub shift_held: bool,
    /// Currently selected cell for residue editing (row, col).
    pub edit_cell: Option<(usize, usize)>,
    /// Dynamic cell width (CELL_W * zoom).
    pub cell_w: f32,
    /// Dynamic cell height (CELL_H * zoom).
    pub cell_h: f32,
    /// Dynamic font size (FONT_SIZE * zoom).
    pub font_size: f32,
    /// Currently selected columns (for highlight and bulk delete).
    pub selected_cols: HashSet<usize>,
    /// Whether to draw feature annotation stripes.
    pub show_features: bool,
    /// Pattern search hits (seq_idx, col) — highlighted with a yellow border.
    pub pattern_matches: Vec<(usize, usize)>,
    /// Index of the currently focused pattern match.
    pub pattern_match_idx: usize,
    /// Opaque stamp — when this differs from `GridState::cached_version`, the
    /// grid cache is cleared and everything is redrawn.  Computed in the view
    /// function from scroll position, alignment pointer, zoom, etc.
    pub version: u64,
    /// Active residue color table.
    pub color_table: Arc<ColorTable>,
}

impl canvas::Program<Message> for AlignmentGrid {
    type State = GridState;

    fn draw(
        &self,
        state: &GridState,
        renderer: &iced::Renderer,
        _theme: &iced::Theme,
        bounds: Rectangle,
        _cursor: mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        // ── Layer 1: main grid (cached) ───────────────────────────────────
        // Invalidate whenever the visual data changes (scroll, edit, zoom…).
        // Hover-only updates do NOT change `self.version`, so they reuse the cache.
        if self.version != state.cached_version.get() {
            state.grid_cache.clear();
            state.cached_version.set(self.version);
        }

        // Snap the cache allocation to 32-pixel steps so that live window-resize only
        // invalidates the cache every 32 px instead of every pixel (~30× fewer
        // re-renders during drag). Drawing uses actual bounds to prevent overflow
        // into adjacent widgets (e.g. the entropy strip below).
        const SNAP: f32 = 32.0;
        let snapped = Size {
            width: (bounds.width / SNAP).ceil() * SNAP,
            height: (bounds.height / SNAP).ceil() * SNAP,
        };

        let grid = state.grid_cache.draw(renderer, snapped, |frame| {
            self.draw_grid(frame, bounds.size());
        });

        // ── Layer 2: drag indicator (live, cheap) ─────────────────────────
        // Drawn fresh every frame so the drop-line follows the cursor during
        // a row-drag without forcing a full grid redraw.
        let mut overlay = Frame::new(renderer, bounds.size());
        if state.is_dragging {
            self.draw_drag_overlay(&mut overlay, bounds.size(), state);
        }

        vec![grid, overlay.into_geometry()]
    }

    fn update(
        &self,
        state: &mut GridState,
        event: canvas::Event,
        bounds: Rectangle,
        cursor: mouse::Cursor,
    ) -> (canvas::event::Status, Option<Message>) {
        match event {
            canvas::Event::Mouse(mouse::Event::WheelScrolled { delta }) => {
                let (dx, dy) = match delta {
                    mouse::ScrollDelta::Lines { x, y } => {
                        if self.ctrl_held {
                            // Ctrl+scroll: force horizontal
                            (-(y as isize) * 5, 0)
                        } else {
                            (x as isize, -(y as isize) * 3)
                        }
                    }
                    mouse::ScrollDelta::Pixels { x, y } => {
                        if self.ctrl_held {
                            (-(y / self.cell_w) as isize, 0)
                        } else {
                            ((x / self.cell_w) as isize, -(y / self.cell_h) as isize)
                        }
                    }
                };
                if dx != 0 || dy != 0 {
                    return (
                        canvas::event::Status::Captured,
                        Some(Message::GridScrolled { dx, dy }),
                    );
                }
            }

            canvas::Event::Mouse(mouse::Event::CursorMoved { position }) => {
                if let Some(pos) = cursor.position_in(bounds) {
                    // Update drag target if we are currently dragging in the title area.
                    if let Some(from) = state.drag_from {
                        if pos.x < TITLE_W && pos.y >= RULER_H {
                            let ri = ((pos.y - RULER_H) / self.cell_h) as usize;
                            let target =
                                (self.scroll_row + ri).min(self.aln.seq_count().saturating_sub(1));
                            state.drag_to = Some(target);
                            if target != from {
                                state.is_dragging = true;
                            }
                            return (canvas::event::Status::Captured, None);
                        }
                    }

                    // Normal hover update over the residue grid.
                    if pos.x > TITLE_W && pos.y > RULER_H {
                        let ci = ((pos.x - TITLE_W) / self.cell_w) as usize;
                        let ri = ((pos.y - RULER_H) / self.cell_h) as usize;
                        let col = self.scroll_col + ci;
                        let row = self.scroll_row + ri;
                        if col < self.aln.col_count() && row < self.aln.seq_count() {
                            return (
                                canvas::event::Status::Ignored,
                                Some(Message::HoverPosition {
                                    col: Some(col),
                                    row: Some(row),
                                }),
                            );
                        }
                    }
                    return (
                        canvas::event::Status::Ignored,
                        Some(Message::HoverPosition {
                            col: None,
                            row: None,
                        }),
                    );
                }
                let _ = position;
            }

            canvas::Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Left)) => {
                if let Some(pos) = cursor.position_in(bounds) {
                    if pos.x >= TITLE_W && pos.y < RULER_H {
                        // Ruler area → start column selection drag.
                        let ci = ((pos.x - TITLE_W) / self.cell_w) as usize;
                        let col = self.scroll_col + ci;
                        if col < self.aln.col_count() {
                            state.col_drag_start = Some(col);
                            return (canvas::event::Status::Captured, None);
                        }
                    } else if pos.x < TITLE_W && pos.y >= RULER_H {
                        // Title area → start row drag/select; double-click opens editor.
                        let ri = ((pos.y - RULER_H) / self.cell_h) as usize;
                        let actual_row = self.scroll_row + ri;
                        if actual_row < self.aln.seq_count() {
                            // Double-click detection (within 400 ms, same row).
                            let now = std::time::Instant::now();
                            let is_double = state
                                .last_title_click
                                .as_ref()
                                .map(|(r, t)| *r == actual_row && t.elapsed().as_millis() < 400)
                                .unwrap_or(false);
                            if is_double {
                                state.last_title_click = None;
                                return (
                                    canvas::event::Status::Captured,
                                    Some(Message::EditSequenceRaw { idx: actual_row }),
                                );
                            }
                            state.last_title_click = Some((actual_row, now));
                            state.drag_from = Some(actual_row);
                            state.drag_to = Some(actual_row);
                            state.is_dragging = false;
                            return (canvas::event::Status::Captured, None);
                        }
                    } else if pos.x >= TITLE_W && pos.y >= RULER_H {
                        // Residue grid → select cell for editing.
                        let ci = ((pos.x - TITLE_W) / self.cell_w) as usize;
                        let ri = ((pos.y - RULER_H) / self.cell_h) as usize;
                        let col = self.scroll_col + ci;
                        let row = self.scroll_row + ri;
                        if col < self.aln.col_count() && row < self.aln.seq_count() {
                            return (
                                canvas::event::Status::Captured,
                                Some(Message::CellClicked { row, col }),
                            );
                        }
                    }
                }
            }

            canvas::Event::Mouse(mouse::Event::ButtonReleased(mouse::Button::Left)) => {
                // Finalise column selection drag.
                if let Some(start) = state.col_drag_start.take() {
                    if let Some(pos) = cursor.position_in(bounds) {
                        let ci = ((pos.x - TITLE_W).max(0.0) / self.cell_w) as usize;
                        let end =
                            (self.scroll_col + ci).min(self.aln.col_count().saturating_sub(1));
                        return (
                            canvas::event::Status::Captured,
                            Some(Message::ColumnsSelected { start, end }),
                        );
                    }
                }

                // Finalise row drag/select.
                let from = state.drag_from.take();
                let to = state.drag_to.take();
                let was_dragging = state.is_dragging;
                state.is_dragging = false;

                if let Some(from_row) = from {
                    if was_dragging {
                        let to_row = to.unwrap_or(from_row);
                        if from_row != to_row {
                            return (
                                canvas::event::Status::Captured,
                                Some(Message::SequenceMoved {
                                    from: from_row,
                                    to: to_row,
                                }),
                            );
                        }
                    } else {
                        return (
                            canvas::event::Status::Captured,
                            Some(Message::SequenceClicked {
                                idx: from_row,
                                ctrl: self.ctrl_held,
                                shift: self.shift_held,
                            }),
                        );
                    }
                }
            }

            canvas::Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Right)) => {
                if let Some(pos) = cursor.position_in(bounds) {
                    if pos.x < TITLE_W && pos.y >= RULER_H {
                        // Right-click on sequence title → context menu.
                        let ri = ((pos.y - RULER_H) / self.cell_h) as usize;
                        let actual_row = self.scroll_row + ri;
                        if actual_row < self.aln.seq_count() {
                            return (
                                canvas::event::Status::Captured,
                                Some(Message::TitleRightClick { idx: actual_row }),
                            );
                        }
                    } else if pos.x > TITLE_W && pos.y > RULER_H {
                        let ci = ((pos.x - TITLE_W) / self.cell_w) as usize;
                        let col = self.scroll_col + ci;
                        if col < self.aln.col_count() {
                            return (
                                canvas::event::Status::Captured,
                                Some(Message::GridRightClick { col }),
                            );
                        }
                    }
                }
            }

            _ => {}
        }
        (canvas::event::Status::Ignored, None)
    }

    fn mouse_interaction(
        &self,
        state: &GridState,
        bounds: Rectangle,
        cursor: mouse::Cursor,
    ) -> mouse::Interaction {
        if let Some(pos) = cursor.position_in(bounds) {
            if pos.x < TITLE_W && pos.y >= RULER_H {
                return if state.is_dragging {
                    mouse::Interaction::Grabbing
                } else {
                    mouse::Interaction::Grab
                };
            }
            return mouse::Interaction::Crosshair;
        }
        mouse::Interaction::default()
    }
}

// ── Drawing ───────────────────────────────────────────────────────────────────

impl AlignmentGrid {
    fn draw_grid(&self, frame: &mut Frame<iced::Renderer>, size: Size) {
        let nrows = self.aln.seq_count();
        let ncols = self.aln.col_count();

        let vis_cols = ((size.width - TITLE_W) / self.cell_w).ceil() as usize + 2;
        let vis_rows = ((size.height - RULER_H) / self.cell_h).ceil() as usize + 1;

        // ── Background ───────────────────────────────────────────────────
        frame.fill_rectangle(Point::ORIGIN, size, Color::from_rgb(0.98, 0.98, 0.98));

        // ── Precompute per-column entropy for Strength scheme ────────────
        let entropy_scores: Vec<f32> = if self.scheme == ColorScheme::Strength {
            (0..vis_cols)
                .map(|ci| {
                    let c = self.scroll_col + ci;
                    if c < ncols {
                        col_entropy_at(&self.aln, c) as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        } else {
            vec![]
        };

        // ── Precompute per-column identity for Identity scheme ────────────
        let (identity_scores, consensus_residues) = if self.scheme == ColorScheme::Identity {
            let scores: Vec<f32> = (0..vis_cols)
                .map(|ci| {
                    let c = self.scroll_col + ci;
                    if c < ncols {
                        self.aln.column_identity(c)
                    } else {
                        0.0
                    }
                })
                .collect();
            let cons: Vec<u8> = (0..vis_cols)
                .map(|ci| {
                    let c = self.scroll_col + ci;
                    if c < ncols {
                        self.aln.column_consensus_residue(c)
                    } else {
                        b'-'
                    }
                })
                .collect();
            (scores, cons)
        } else {
            (vec![], vec![])
        };

        // ── Residue rows ─────────────────────────────────────────────────
        for ri in 0..vis_rows {
            let actual_row = self.scroll_row + ri;
            if actual_row >= nrows {
                break;
            }

            let seq = &self.aln.sequences[actual_row];
            let y = RULER_H + ri as f32 * self.cell_h;

            // Whole-row background (grid area)
            let row_bg = if self.selected.contains(&actual_row) {
                Color::from_rgb(0.78, 0.88, 1.00)
            } else if ri % 2 == 0 {
                Color::WHITE
            } else {
                Color::from_rgb(0.96, 0.97, 1.00)
            };
            frame.fill_rectangle(
                Point::new(TITLE_W, y),
                Size::new(size.width - TITLE_W, self.cell_h),
                row_bg,
            );

            // ── Selected column highlight ─────────────────────────────────
            for ci in 0..vis_cols {
                let actual_col = self.scroll_col + ci;
                if self.selected_cols.contains(&actual_col) {
                    let x = TITLE_W + ci as f32 * self.cell_w;
                    frame.fill_rectangle(
                        Point::new(x, y),
                        Size::new(self.cell_w, self.cell_h),
                        Color {
                            r: 0.30,
                            g: 0.55,
                            b: 1.00,
                            a: 0.25,
                        },
                    );
                }
            }

            // Per-cell coloring
            for ci in 0..vis_cols {
                let actual_col = self.scroll_col + ci;
                if actual_col >= ncols {
                    break;
                }

                let x = TITLE_W + ci as f32 * self.cell_w;
                let residue = seq.residues[actual_col];

                let colored = match self.scheme {
                    ColorScheme::Plain => false,
                    ColorScheme::ResidueType => !is_gap(residue),
                    ColorScheme::Identity => {
                        if is_gap(residue) {
                            false
                        } else {
                            let ident = *identity_scores.get(ci).unwrap_or(&0.0);
                            let cons = *consensus_residues.get(ci).unwrap_or(&b'-');
                            ident >= IDENTITY_THRESHOLD && residue.to_ascii_uppercase() == cons
                        }
                    }
                    ColorScheme::Strength => !is_gap(residue),
                };

                if colored {
                    const MAX_ENTROPY: f32 = 3.0;
                    let bg = if self.scheme == ColorScheme::Strength {
                        let ent = entropy_scores.get(ci).copied().unwrap_or(0.0);
                        let conservation = (1.0 - ent / MAX_ENTROPY).clamp(0.0, 1.0);
                        let base = residue_bg(residue, seq.seq_type.clone(), &self.color_table);
                        Color {
                            r: 1.0 - conservation * (1.0 - base.r),
                            g: 1.0 - conservation * (1.0 - base.g),
                            b: 1.0 - conservation * (1.0 - base.b),
                            a: 1.0,
                        }
                    } else {
                        residue_bg(residue, seq.seq_type.clone(), &self.color_table)
                    };
                    frame.fill_rectangle(
                        Point::new(x, y),
                        Size::new(self.cell_w - 0.5, self.cell_h - 0.5),
                        bg,
                    );
                }

                // Skip text rendering for gaps (they show as plain background) and
                // for cells too small to read — both eliminate expensive glyph calls.
                if !is_gap(residue) && self.cell_w >= 5.5 {
                    frame.fill_text(Text {
                        content: (residue as char).to_string(),
                        position: Point::new(x + self.cell_w * 0.5, y + self.cell_h * 0.5),
                        color: Color::from_rgb(0.08, 0.08, 0.10),
                        size: Pixels(self.font_size),
                        font: Font::MONOSPACE,
                        horizontal_alignment: alignment::Horizontal::Center,
                        vertical_alignment: alignment::Vertical::Center,
                        ..Text::default()
                    });
                }
            }

            // ── Edit-cell cursor highlight ────────────────────────────────
            if let Some((edit_row, edit_col)) = self.edit_cell {
                if edit_row == actual_row {
                    if let Some(ci) = edit_col.checked_sub(self.scroll_col) {
                        if ci < vis_cols {
                            let x = TITLE_W + ci as f32 * self.cell_w;
                            frame.stroke(
                                &Path::rectangle(
                                    Point::new(x + 0.5, y + 0.5),
                                    Size::new(self.cell_w - 1.0, self.cell_h - 1.0),
                                ),
                                canvas::Stroke::default()
                                    .with_color(Color::from_rgb(0.05, 0.45, 0.85))
                                    .with_width(2.0),
                            );
                        }
                    }
                }
            }

            // ── Pattern match highlights ──────────────────────────────────
            for ci in 0..vis_cols {
                let actual_col = self.scroll_col + ci;
                if self.pattern_matches.contains(&(actual_row, actual_col)) {
                    let x = TITLE_W + ci as f32 * self.cell_w;
                    let is_focused = self
                        .pattern_matches
                        .get(self.pattern_match_idx)
                        .map(|&hit| hit == (actual_row, actual_col))
                        .unwrap_or(false);
                    let outline_color = if is_focused {
                        Color::from_rgb(1.0, 0.55, 0.0) // orange = focused hit
                    } else {
                        Color::from_rgb(0.95, 0.85, 0.0) // yellow = other hit
                    };
                    frame.fill_rectangle(
                        Point::new(x, y),
                        Size::new(self.cell_w, self.cell_h),
                        Color {
                            a: 0.35,
                            ..outline_color
                        },
                    );
                }
            }

            // ── Feature stripes for this row ──────────────────────────────
            // Three fixed tiers drawn as 2 px bands stacked at the cell bottom:
            //   tier 0 (bottom): sub-features — exon, intron, peptides, misc
            //   tier 1 (middle): functional   — CDS, mRNA, tRNA, rRNA
            //   tier 2 (top):    gene-level   — gene, operon
            // Each tier independently picks the first feature that covers the column.
            if self.show_features && !seq.features.is_empty() {
                // Build three per-tier band lists. alignment_col() is O(seq_len);
                // call it once per feature here rather than per cell.
                let mut tier_bands: [Vec<(usize, usize, Color)>; FEATURE_TIERS] =
                    [Vec::new(), Vec::new(), Vec::new()];

                for f in &seq.features {
                    let fc = f.feature_type.default_color();
                    if fc.a < 0.01 {
                        continue;
                    } // transparent → skip (Source)
                    let Some(aln_start) = seq.alignment_col(f.start) else {
                        continue;
                    };
                    let aln_end = seq
                        .alignment_col(f.end)
                        .unwrap_or_else(|| seq.len().saturating_sub(1));
                    let tier = feature_tier(&f.feature_type);
                    tier_bands[tier].push((aln_start, aln_end, core_to_iced(fc)));
                }

                let any = tier_bands.iter().any(|v| !v.is_empty());
                if any {
                    for ci in 0..vis_cols {
                        let actual_col = self.scroll_col + ci;
                        if actual_col >= ncols {
                            break;
                        }
                        let x = TITLE_W + ci as f32 * self.cell_w;

                        for tier in 0..FEATURE_TIERS {
                            if let Some(&(_, _, fc)) = tier_bands[tier]
                                .iter()
                                .find(|&&(fs, fe, _)| actual_col >= fs && actual_col <= fe)
                            {
                                // Tier 0 = bottom-most stripe; tier 2 = top stripe.
                                let stripe_y =
                                    y + self.cell_h - (tier as f32 + 1.0) * FEATURE_LAYER_H;
                                frame.fill_rectangle(
                                    Point::new(x, stripe_y),
                                    Size::new(self.cell_w, FEATURE_LAYER_H),
                                    fc,
                                );
                            }
                        }
                    }
                }
            }
        }

        // ── Title panel background ───────────────────────────────────────
        frame.fill_rectangle(
            Point::new(0.0, RULER_H),
            Size::new(TITLE_W, size.height - RULER_H),
            Color::from_rgb(0.94, 0.94, 0.96),
        );

        // ── Sequence names ───────────────────────────────────────────────
        for ri in 0..vis_rows {
            let actual_row = self.scroll_row + ri;
            if actual_row >= nrows {
                break;
            }

            let seq = &self.aln.sequences[actual_row];
            let y = RULER_H + ri as f32 * self.cell_h;

            if self.selected.contains(&actual_row) {
                frame.fill_rectangle(
                    Point::new(0.0, y),
                    Size::new(TITLE_W - 1.0, self.cell_h),
                    Color::from_rgb(0.78, 0.88, 1.00),
                );
            } else if ri % 2 == 1 {
                frame.fill_rectangle(
                    Point::new(0.0, y),
                    Size::new(TITLE_W - 1.0, self.cell_h),
                    Color::from_rgb(0.91, 0.92, 0.95),
                );
            }

            frame.fill_text(Text {
                content: truncate_name(&seq.name, 26),
                position: Point::new(6.0, y + self.cell_h * 0.5),
                color: Color::from_rgb(0.10, 0.10, 0.12),
                size: Pixels(11.5),
                font: Font::MONOSPACE,
                horizontal_alignment: alignment::Horizontal::Left,
                vertical_alignment: alignment::Vertical::Center,
                ..Text::default()
            });
        }

        // ── Ruler background ─────────────────────────────────────────────
        frame.fill_rectangle(
            Point::new(TITLE_W, 0.0),
            Size::new(size.width - TITLE_W, RULER_H),
            Color::from_rgb(0.88, 0.91, 0.97),
        );

        // ── Selected column highlight in ruler ────────────────────────────
        for ci in 0..vis_cols {
            let actual_col = self.scroll_col + ci;
            if self.selected_cols.contains(&actual_col) {
                let x = TITLE_W + ci as f32 * self.cell_w;
                frame.fill_rectangle(
                    Point::new(x, 0.0),
                    Size::new(self.cell_w, RULER_H),
                    Color {
                        r: 0.05,
                        g: 0.45,
                        b: 0.85,
                        a: 0.35,
                    },
                );
            }
        }

        // ── Ruler ticks every 10 columns ─────────────────────────────────
        for ci in 0..vis_cols {
            let actual_col = self.scroll_col + ci;
            if actual_col >= ncols {
                break;
            }
            if (actual_col + 1) % 10 == 0 || actual_col == 0 {
                let x = TITLE_W + ci as f32 * self.cell_w;
                frame.fill_text(Text {
                    content: (actual_col + 1).to_string(),
                    position: Point::new(x + 1.0, RULER_H * 0.5),
                    color: Color::from_rgb(0.10, 0.10, 0.12),
                    size: Pixels(10.0),
                    font: Font::MONOSPACE,
                    horizontal_alignment: alignment::Horizontal::Left,
                    vertical_alignment: alignment::Vertical::Center,
                    ..Text::default()
                });
            }
        }

        // ── Corner ───────────────────────────────────────────────────────
        frame.fill_rectangle(
            Point::ORIGIN,
            Size::new(TITLE_W, RULER_H),
            Color::from_rgb(0.88, 0.91, 0.97),
        );

        // ── Dividers ─────────────────────────────────────────────────────
        let border = Color::from_rgb(0.78, 0.78, 0.82);
        let stroke = canvas::Stroke::default().with_color(border).with_width(1.0);

        frame.stroke(
            &Path::line(Point::new(TITLE_W, 0.0), Point::new(TITLE_W, size.height)),
            stroke.clone(),
        );
        frame.stroke(
            &Path::line(Point::new(0.0, RULER_H), Point::new(size.width, RULER_H)),
            stroke,
        );
    }

    /// Draws the drag-reorder indicator into a separate (non-cached) overlay frame.
    /// Called every render frame while a row drag is active, which is why it lives
    /// outside the main cache — it only draws 1–2 cheap rectangles.
    fn draw_drag_overlay(&self, frame: &mut Frame<iced::Renderer>, size: Size, state: &GridState) {
        if let (Some(from), Some(to)) = (state.drag_from, state.drag_to) {
            let vis_rows = ((size.height - RULER_H) / self.cell_h).ceil() as usize + 1;

            // Dim the source row label.
            let from_vis = from.wrapping_sub(self.scroll_row);
            if from_vis < vis_rows {
                let y_src = RULER_H + from_vis as f32 * self.cell_h;
                frame.fill_rectangle(
                    Point::new(0.0, y_src),
                    Size::new(TITLE_W - 1.0, self.cell_h),
                    Color {
                        r: 0.05,
                        g: 0.45,
                        b: 0.85,
                        a: 0.25,
                    },
                );
            }

            // Drop-target line.
            let to_vis = to.wrapping_sub(self.scroll_row);
            let line_y = if to >= from {
                RULER_H + (to_vis + 1) as f32 * self.cell_h
            } else {
                RULER_H + to_vis as f32 * self.cell_h
            };
            if line_y >= RULER_H && line_y <= size.height {
                frame.fill_rectangle(
                    Point::new(0.0, line_y - 1.5),
                    Size::new(TITLE_W - 1.0, 3.0),
                    Color::from_rgb(0.05, 0.45, 0.85),
                );
            }
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn col_entropy_at(aln: &SeqAlignment, col: usize) -> f64 {
    let mut counts = [0u32; 256];
    let mut total = 0u32;
    for seq in &aln.sequences {
        if col < seq.residues.len() {
            let b = seq.residues[col].to_ascii_uppercase();
            if !matches!(b, b'-' | b'~' | b'.') {
                counts[b as usize] += 1;
                total += 1;
            }
        }
    }
    if total == 0 {
        return 0.0;
    }
    let mut h = 0.0f64;
    for &c in counts.iter() {
        if c > 0 {
            let p = c as f64 / total as f64;
            h -= p * p.ln();
        }
    }
    h
}

fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'~' | b'.')
}

fn truncate_name(name: &str, max_chars: usize) -> String {
    let count = name.chars().count();
    if count <= max_chars {
        name.to_owned()
    } else {
        let s: String = name.chars().take(max_chars - 1).collect();
        format!("{s}…")
    }
}

fn residue_bg(b: u8, seq_type: SequenceType, table: &ColorTable) -> Color {
    match seq_type {
        SequenceType::Dna | SequenceType::Rna | SequenceType::NucleicAcid => table.nuc_color(b),
        SequenceType::Protein => table.aa_color(b),
        // Unknown: try nucleotide table first; if grey (default), try AA table.
        _ => {
            let c = table.nuc_color(b);
            if c.r > 0.89 && c.g > 0.89 && c.b > 0.89 {
                // Likely a default grey — try amino acid table instead
                let ca = table.aa_color(b);
                if ca.r < 0.89 || ca.g < 0.89 || ca.b < 0.89 {
                    ca
                } else {
                    c
                }
            } else {
                c
            }
        }
    }
}

/// Convert a helixview-core Color (f32 RGBA) to an iced Color.
#[inline]
fn core_to_iced(c: CoreColor) -> Color {
    Color {
        r: c.r,
        g: c.g,
        b: c.b,
        a: c.a,
    }
}

/// Map a feature type to its display tier (0 = bottom/specific, 2 = top/broad).
///
/// Tier 0 — sub-features: exon, intron, peptide signals, misc annotations
/// Tier 1 — functional:   CDS, mRNA, tRNA, rRNA, misc_RNA
/// Tier 2 — gene-level:   gene, operon
fn feature_tier(ft: &FeatureType) -> usize {
    match ft {
        FeatureType::Gene | FeatureType::Operon => 2,
        FeatureType::CDS
        | FeatureType::MRNA
        | FeatureType::TRNA
        | FeatureType::RRNA
        | FeatureType::PrecursorRNA
        | FeatureType::MiscRNA
        | FeatureType::SnRNA
        | FeatureType::ScRNA
        | FeatureType::RBS => 1,
        _ => 0, // Exon, Intron, SigPeptide, MiscFeature, RepeatRegion, etc.
    }
}
