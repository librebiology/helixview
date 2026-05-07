//! Phylogenetic tree viewer — cladogram canvas with leaf selection.

use std::cell::Cell;
use std::collections::BTreeSet;
use std::sync::Arc;

use helixview_core::{layout_cladogram, LayoutNode, PhyloTree};
use iced::{
    alignment, mouse,
    widget::{
        button,
        canvas::{self, Frame, Path, Stroke, Text},
        column, horizontal_space, row, text,
    },
    Color, Element, Font, Length, Pixels, Point, Rectangle, Size,
};

use crate::app::Message;

// ── Public entry ──────────────────────────────────────────────────────────────

pub fn tree_view<'a>(tree: &Arc<PhyloTree>, selected: &BTreeSet<String>) -> Element<'a, Message> {
    let leaf_count = tree.leaf_count();

    let toolbar = {
        let btn = |label: &str, msg: Message| -> iced::widget::Button<'static, Message> {
            button(text(label.to_string()).size(12))
                .padding([3, 10])
                .style(button::secondary)
                .on_press(msg)
        };

        let sel_label = if selected.is_empty() {
            "No selection".to_string()
        } else {
            format!("{} selected", selected.len())
        };

        row![
            btn("< Back", Message::CloseTree),
            horizontal_space(),
            text(format!("{leaf_count} leaves")).size(12),
            horizontal_space(),
            text(sel_label).size(12),
            btn("Reorder Alignment", Message::TreeReorderAlignment),
            btn("Clear Selection", Message::TreeClearSelection),
        ]
        .spacing(6)
        .padding([4, 8])
    };

    let canvas_widget = iced::widget::canvas(TreeCanvas {
        tree: Arc::clone(tree),
        selected: selected.clone(),
    })
    .width(Length::Fill)
    .height(Length::Fill);

    column![toolbar, canvas_widget]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Canvas ────────────────────────────────────────────────────────────────────

struct TreeCanvas {
    tree: Arc<PhyloTree>,
    selected: BTreeSet<String>,
}

pub struct TreeCanvasState {
    cache: canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for TreeCanvasState {
    fn default() -> Self {
        Self {
            cache: canvas::Cache::default(),
            last_key: Cell::new(u64::MAX),
        }
    }
}

impl canvas::Program<Message> for TreeCanvas {
    type State = TreeCanvasState;

    fn update(
        &self,
        _state: &mut TreeCanvasState,
        event: canvas::Event,
        bounds: Rectangle,
        cursor: mouse::Cursor,
    ) -> (canvas::event::Status, Option<Message>) {
        if let canvas::Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Left)) = event {
            if let mouse::Cursor::Available(pt) = cursor {
                if bounds.contains(pt) {
                    let rel = Point::new(pt.x - bounds.x, pt.y - bounds.y);
                    let w = bounds.width;
                    let h = bounds.height;
                    const MARGIN_X: f32 = 20.0;
                    const LABEL_W: f32 = 160.0;
                    const MARGIN_Y: f32 = 10.0;
                    let plot_w = (w - MARGIN_X - LABEL_W).max(1.0);
                    let plot_h = (h - MARGIN_Y * 2.0).max(1.0);

                    // Compute layout fresh — state.nodes is never populated from draw()
                    // because draw() takes &State (immutable). Layout is fast for
                    // typical tree sizes so recomputing on click is fine.
                    let nodes = layout_cladogram(&self.tree.root);

                    for node in &nodes {
                        if !node.is_leaf {
                            continue;
                        }
                        let nx = MARGIN_X + node.x * plot_w;
                        let ny = MARGIN_Y + node.y * plot_h;
                        let dist = ((rel.x - nx).powi(2) + (rel.y - ny).powi(2)).sqrt();
                        if dist < 10.0 {
                            return (
                                canvas::event::Status::Captured,
                                Some(Message::TreeNodeClicked(node.name.clone())),
                            );
                        }
                    }
                }
            }
        }
        (canvas::event::Status::Ignored, None)
    }

    fn draw(
        &self,
        state: &TreeCanvasState,
        renderer: &iced::Renderer,
        _theme: &iced::Theme,
        bounds: Rectangle,
        _cursor: mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let sel_len = self.selected.len() as u64;
        let key = (Arc::as_ptr(&self.tree) as u64)
            .wrapping_add(sel_len.wrapping_mul(0x9e3779b9))
            .wrapping_add(bounds.width.to_bits() as u64)
            .wrapping_add(bounds.height.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
            // Safety: canvas::Program::draw takes &State not &mut State, but we
            // need to rebuild the node list for hit-testing. Use interior mutability
            // via the cache key path — we rebuild nodes inside the closure below.
        }

        let tree = Arc::clone(&self.tree);
        let selected = self.selected.clone();

        const SNAP: f32 = 32.0;
        let snapped = Size {
            width: (bounds.width / SNAP).ceil() * SNAP,
            height: (bounds.height / SNAP).ceil() * SNAP,
        };

        let geom = state.cache.draw(renderer, snapped, |frame| {
            let nodes = layout_cladogram(&tree.root);
            draw_tree(frame, snapped, &nodes, &selected);
        });

        vec![geom]
    }
}

// ── Drawing ───────────────────────────────────────────────────────────────────

fn draw_tree(
    frame: &mut Frame<iced::Renderer>,
    size: Size,
    nodes: &[LayoutNode],
    selected: &BTreeSet<String>,
) {
    use crate::theme::palette;

    frame.fill_rectangle(Point::ORIGIN, size, palette::BG_PANEL);

    let margin_x = 20.0f32;
    let label_area = 160.0f32;
    let margin_y = 10.0f32;
    let plot_w = (size.width - margin_x - label_area).max(1.0);
    let plot_h = (size.height - margin_y * 2.0).max(1.0);

    let to_pt =
        |n: &LayoutNode| -> Point { Point::new(margin_x + n.x * plot_w, margin_y + n.y * plot_h) };

    // Draw branches
    for node in nodes {
        let pt = to_pt(node);

        if let Some(pidx) = node.parent_idx {
            let parent = &nodes[pidx];
            let ppt = to_pt(parent);
            // Cladogram: horizontal line from parent_x to node_x, then vertical to node_y
            let corner = Point::new(ppt.x, pt.y);
            let branch = Path::new(|b| {
                b.move_to(ppt);
                b.line_to(corner);
                b.line_to(pt);
            });
            let color = if node.is_leaf && selected.contains(&node.name) {
                Color::from_rgb(0.2, 0.5, 1.0)
            } else {
                palette::TEXT_DIM
            };
            frame.stroke(&branch, Stroke::default().with_color(color).with_width(1.5));
        }

        // Support value on internal nodes
        if !node.is_leaf {
            if let Some(sup) = node.support {
                if sup >= 50.0 {
                    frame.fill_text(Text {
                        content: format!("{:.0}", sup),
                        position: Point::new(pt.x + 2.0, pt.y - 6.0),
                        color: Color::from_rgb(0.45, 0.45, 0.45),
                        size: Pixels(7.5),
                        font: Font::MONOSPACE,
                        horizontal_alignment: alignment::Horizontal::Left,
                        vertical_alignment: alignment::Vertical::Bottom,
                        ..Text::default()
                    });
                }
            }
        }
    }

    // Draw leaf tips and labels (on top of branches)
    for node in nodes {
        if !node.is_leaf {
            continue;
        }
        let pt = to_pt(node);
        let is_sel = selected.contains(&node.name);

        // Dot
        let dot = Path::circle(pt, if is_sel { 4.0 } else { 2.5 });
        let dot_color = if is_sel {
            Color::from_rgb(0.2, 0.5, 1.0)
        } else {
            Color::from_rgb(0.35, 0.35, 0.35)
        };
        frame.fill(&dot, dot_color);

        // Label
        let label_x = pt.x + 6.0;
        let label = if node.name.len() > 22 {
            format!("{}…", &node.name[..21])
        } else {
            node.name.clone()
        };
        frame.fill_text(Text {
            content: label,
            position: Point::new(label_x, pt.y),
            color: if is_sel {
                Color::from_rgb(0.1, 0.4, 0.9)
            } else {
                palette::TEXT_DIM
            },
            size: Pixels(10.0),
            font: Font::MONOSPACE,
            horizontal_alignment: alignment::Horizontal::Left,
            vertical_alignment: alignment::Vertical::Center,
            ..Text::default()
        });
    }
}
