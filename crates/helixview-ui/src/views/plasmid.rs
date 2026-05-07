//! Circular plasmid map viewer with feature list panel and inline editing.

use std::cell::Cell;
use std::f32::consts::{PI, TAU};
use std::sync::Arc;

use helixview_analysis::{find_sites, RestrictionSite, ENZYMES};
use helixview_core::color::Color as CoreColor;
use helixview_core::feature::{FeatureDirection, FeatureType};
use helixview_core::Alignment as SeqAlignment;
use iced::{
    alignment, mouse,
    widget::{
        button,
        canvas::{self, Frame, Path, Stroke, Text},
        column, container, horizontal_space, row, scrollable, text, text_input,
    },
    Background, Border, Color, Element, Font, Length, Pixels, Point, Rectangle, Size,
};

use crate::app::Message;
use crate::theme::palette;

// ── Palette for RE sites ──────────────────────────────────────────────────────

const RE_PALETTE: [Color; 20] = [
    Color {
        r: 0.85,
        g: 0.15,
        b: 0.15,
        a: 1.0,
    },
    Color {
        r: 0.10,
        g: 0.60,
        b: 0.10,
        a: 1.0,
    },
    Color {
        r: 0.10,
        g: 0.30,
        b: 0.85,
        a: 1.0,
    },
    Color {
        r: 0.90,
        g: 0.50,
        b: 0.05,
        a: 1.0,
    },
    Color {
        r: 0.55,
        g: 0.10,
        b: 0.70,
        a: 1.0,
    },
    Color {
        r: 0.05,
        g: 0.65,
        b: 0.65,
        a: 1.0,
    },
    Color {
        r: 0.65,
        g: 0.35,
        b: 0.05,
        a: 1.0,
    },
    Color {
        r: 0.65,
        g: 0.05,
        b: 0.40,
        a: 1.0,
    },
    Color {
        r: 0.20,
        g: 0.55,
        b: 0.25,
        a: 1.0,
    },
    Color {
        r: 0.75,
        g: 0.65,
        b: 0.05,
        a: 1.0,
    },
    Color {
        r: 0.15,
        g: 0.50,
        b: 0.75,
        a: 1.0,
    },
    Color {
        r: 0.85,
        g: 0.30,
        b: 0.60,
        a: 1.0,
    },
    Color {
        r: 0.35,
        g: 0.65,
        b: 0.35,
        a: 1.0,
    },
    Color {
        r: 0.60,
        g: 0.25,
        b: 0.75,
        a: 1.0,
    },
    Color {
        r: 0.25,
        g: 0.70,
        b: 0.55,
        a: 1.0,
    },
    Color {
        r: 0.75,
        g: 0.40,
        b: 0.20,
        a: 1.0,
    },
    Color {
        r: 0.45,
        g: 0.20,
        b: 0.60,
        a: 1.0,
    },
    Color {
        r: 0.20,
        g: 0.55,
        b: 0.50,
        a: 1.0,
    },
    Color {
        r: 0.70,
        g: 0.15,
        b: 0.15,
        a: 1.0,
    },
    Color {
        r: 0.30,
        g: 0.30,
        b: 0.70,
        a: 1.0,
    },
];

/// 16-color palette for user-assignable feature colors.
const FEAT_PALETTE: [CoreColor; 16] = [
    CoreColor::rgb(0.10, 0.60, 0.10), // green
    CoreColor::rgb(0.10, 0.30, 0.85), // blue
    CoreColor::rgb(0.85, 0.15, 0.15), // red
    CoreColor::rgb(0.90, 0.50, 0.05), // orange
    CoreColor::rgb(0.55, 0.10, 0.70), // purple
    CoreColor::rgb(0.05, 0.65, 0.65), // teal
    CoreColor::rgb(0.85, 0.30, 0.60), // pink
    CoreColor::rgb(0.75, 0.65, 0.05), // gold
    CoreColor::rgb(0.15, 0.50, 0.75), // sky
    CoreColor::rgb(0.65, 0.35, 0.05), // brown
    CoreColor::rgb(0.35, 0.65, 0.35), // sage
    CoreColor::rgb(0.60, 0.25, 0.75), // violet
    CoreColor::rgb(0.25, 0.70, 0.55), // seafoam
    CoreColor::rgb(0.75, 0.40, 0.20), // sienna
    CoreColor::rgb(0.45, 0.45, 0.45), // grey
    CoreColor::rgb(0.20, 0.20, 0.65), // indigo
];

// ── Public entry ──────────────────────────────────────────────────────────────

pub fn plasmid_view<'a>(
    aln: &'a Arc<SeqAlignment>,
    seq_idx: usize,
    show_features: bool,
    show_re: bool,
    selected_feat: Option<usize>,
    edit_name: &'a str,
    edit_start: &'a str,
    edit_end: &'a str,
    edit_color: CoreColor,
    edit_hex: &'a str,
    zoom: f32,
    pan_x: f32,
    pan_y: f32,
    re_sites_cache: Option<Arc<Vec<RestrictionSite>>>, // None while computing or RE off
) -> Element<'a, Message> {
    let n_seqs = aln.seq_count();
    let seq = aln.sequences.get(seq_idx);
    let seq_name = seq.map(|s| s.name.as_str()).unwrap_or("?");
    let label = if seq_name.len() > 22 {
        format!("{}…", &seq_name[..21])
    } else {
        seq_name.to_string()
    };

    let btn = |s: &str, m: Message| {
        button(text(s.to_string()).size(12))
            .padding([3, 8])
            .style(button::secondary)
            .on_press(m)
    };
    let feat_label = if show_features {
        "Feat ●"
    } else {
        "Feat ○"
    };
    let re_label = if show_re { "RE ●" } else { "RE ○" };

    let toolbar = container(
        row![
            btn("< Back", Message::ClosePlasmid),
            horizontal_space(),
            btn("<", Message::PlasmidPrevSeq),
            text(format!("{label}  ({}/{})", seq_idx + 1, n_seqs)).size(12),
            btn(">", Message::PlasmidNextSeq),
            horizontal_space(),
            btn(feat_label, Message::PlasmidToggleFeatures),
            btn(re_label, Message::PlasmidToggleRe),
            iced::widget::vertical_rule(1),
            btn("−", Message::PlasmidZoomOut),
            text(format!("{:.0}%", zoom * 100.0))
                .size(11)
                .color(crate::theme::palette::TEXT_DIM),
            btn("+", Message::PlasmidZoomIn),
            btn("⟳", Message::PlasmidZoomReset),
            btn("SVG", Message::ExportPlasmidSvg),
        ]
        .spacing(6)
        .padding([4, 8])
        .align_y(iced::Alignment::Center),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .width(Length::Fill)
    .height(Length::Fixed(36.0));

    let canvas_widget = iced::widget::canvas(PlasmidCanvas {
        aln: Arc::clone(aln),
        seq_idx,
        show_features,
        show_re,
        selected_feat,
        zoom,
        pan_x,
        pan_y,
        re_sites: re_sites_cache.clone(),
    })
    .width(Length::Fill)
    .height(Length::Fill);

    // Build RE sites with colors from the pre-computed cache.
    let re_computing = show_re && re_sites_cache.is_none();
    let re_sites: Vec<(String, usize, Color)> = if show_re {
        match &re_sites_cache {
            Some(cached) => cached
                .iter()
                .map(|s| {
                    let color = RE_PALETTE[ENZYMES[s.enzyme_idx].color_idx % RE_PALETTE.len()];
                    (
                        ENZYMES[s.enzyme_idx].name.to_string(),
                        s.true_pos + 1,
                        color,
                    )
                })
                .collect(),
            None => vec![], // still computing — sidebar shows a loading hint
        }
    } else {
        vec![]
    };
    // Feature list panel (right side)
    let feat_panel = feature_panel(
        aln,
        seq_idx,
        show_re,
        selected_feat,
        edit_name,
        edit_start,
        edit_end,
        edit_color,
        edit_hex,
        re_sites,
        re_computing,
    );

    let body = row![canvas_widget, feat_panel,]
        .width(Length::Fill)
        .height(Length::Fill);

    column![toolbar, body]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Feature list panel ────────────────────────────────────────────────────────

fn feature_panel<'a>(
    aln: &'a Arc<SeqAlignment>,
    seq_idx: usize,
    show_re: bool,
    selected_feat: Option<usize>,
    edit_name: &'a str,
    edit_start: &'a str,
    edit_end: &'a str,
    edit_color: CoreColor,
    edit_hex: &'a str,
    re_sites: Vec<(String, usize, Color)>,
    re_computing: bool,
) -> Element<'a, Message> {
    let features = aln
        .sequences
        .get(seq_idx)
        .map(|s| s.features.as_slice())
        .unwrap_or(&[]);

    // Panel header
    let header = container(
        row![
            text("Features").size(13).color(palette::TEXT),
            horizontal_space(),
            button(text("+ New").size(11))
                .padding([2, 6])
                .on_press(Message::PlasmidAddFeature),
        ]
        .align_y(iced::Alignment::Center)
        .padding([4, 8]),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .width(Length::Fill);

    // Feature rows
    let mut rows: Vec<Element<Message>> = Vec::new();
    for (i, feat) in features.iter().enumerate() {
        if matches!(feat.feature_type, FeatureType::Source) {
            continue;
        }
        let is_sel = selected_feat == Some(i);
        let bg = if is_sel {
            Color::from_rgb(0.20, 0.45, 0.80)
        } else if i % 2 == 0 {
            Color::WHITE
        } else {
            Color::from_rgb(0.97, 0.97, 0.99)
        };
        let txt_col = if is_sel { Color::WHITE } else { palette::TEXT };

        let fc = feat.color;
        let swatch: Element<Message> = container(iced::widget::Space::with_width(10))
            .width(10)
            .height(10)
            .style(move |_| container::Style {
                background: Some(Background::Color(Color {
                    r: fc.r,
                    g: fc.g,
                    b: fc.b,
                    a: 1.0,
                })),
                border: Border {
                    color: Color::from_rgb(0.5, 0.5, 0.5),
                    width: 0.5,
                    radius: 2.0.into(),
                },
                ..Default::default()
            })
            .into();

        let name_label = if feat.name.len() > 18 {
            format!("{}…", &feat.name[..17])
        } else {
            feat.name.clone()
        };
        let dir_sym = match feat.direction {
            FeatureDirection::Forward => "→",
            FeatureDirection::Reverse => "<",
            FeatureDirection::None => "—",
        };
        let pos_label = format!("{}..{}", feat.start + 1, feat.end + 1);

        let del_btn: Element<Message> = if is_sel {
            button(text("✕").size(10).color(Color::WHITE))
                .padding([1, 4])
                .style(|_, _| button::Style {
                    background: Some(Background::Color(Color::from_rgb(0.75, 0.20, 0.20))),
                    border: Border {
                        radius: 3.0.into(),
                        ..Default::default()
                    },
                    ..Default::default()
                })
                .on_press(Message::PlasmidDeleteFeature(i))
                .into()
        } else {
            iced::widget::Space::with_width(20).into()
        };

        rows.push(
            button(
                container(
                    row![
                        swatch,
                        text(dir_sym).size(10).color(txt_col),
                        text(name_label).size(11).color(txt_col),
                        horizontal_space(),
                        text(pos_label).size(10).color(if is_sel {
                            Color::from_rgba(1.0, 1.0, 1.0, 0.75)
                        } else {
                            palette::TEXT_DIM
                        }),
                        del_btn,
                    ]
                    .spacing(4)
                    .align_y(iced::Alignment::Center)
                    .padding([4, 6]),
                )
                .width(Length::Fill)
                .style(move |_| container::Style {
                    background: Some(Background::Color(bg)),
                    ..Default::default()
                }),
            )
            .width(Length::Fill)
            .padding(0)
            .style(|_, _| button::Style {
                background: None,
                border: Border::default(),
                ..Default::default()
            })
            .on_press(Message::PlasmidSelectFeature(Some(i)))
            .into(),
        );
    }

    if rows.is_empty() {
        rows.push(
            container(
                text("No features. Click \"+ New\" to add one.")
                    .size(11)
                    .color(palette::TEXT_DIM),
            )
            .padding([10, 8])
            .into(),
        );
    }

    let list = scrollable(column(rows).width(Length::Fill).spacing(0))
        .height(Length::Fill)
        .width(Length::Fill);

    // Inline edit panel for selected feature
    let edit_panel: Element<Message> = if selected_feat.is_some() {
        container(
            column![
                container(text("Edit feature").size(12).color(palette::TEXT))
                    .padding([4, 8])
                    .style(|_| container::Style {
                        background: Some(Background::Color(palette::HEADER_BG)),
                        border: Border {
                            color: palette::BORDER,
                            width: 1.0,
                            radius: 0.0.into()
                        },
                        ..Default::default()
                    })
                    .width(Length::Fill),
                container(
                    column![
                        row![
                            text("Name")
                                .size(11)
                                .color(palette::TEXT_DIM)
                                .width(Length::Fixed(45.0)),
                            text_input("Feature name", edit_name)
                                .on_input(Message::PlasmidEditName)
                                .size(11)
                                .padding([3, 6])
                                .width(Length::Fill),
                        ]
                        .spacing(4)
                        .align_y(iced::Alignment::Center),
                        row![
                            text("Start")
                                .size(11)
                                .color(palette::TEXT_DIM)
                                .width(Length::Fixed(45.0)),
                            text_input("1", edit_start)
                                .on_input(Message::PlasmidEditStart)
                                .size(11)
                                .padding([3, 6])
                                .width(Length::Fill),
                        ]
                        .spacing(4)
                        .align_y(iced::Alignment::Center),
                        row![
                            text("End")
                                .size(11)
                                .color(palette::TEXT_DIM)
                                .width(Length::Fixed(45.0)),
                            text_input("100", edit_end)
                                .on_input(Message::PlasmidEditEnd)
                                .size(11)
                                .padding([3, 6])
                                .width(Length::Fill),
                        ]
                        .spacing(4)
                        .align_y(iced::Alignment::Center),
                        {
                            // Color picker: two rows of 8 swatches.
                            let make_swatch = |c: CoreColor| -> Element<'static, Message> {
                                let ic = Color {
                                    r: c.r,
                                    g: c.g,
                                    b: c.b,
                                    a: 1.0,
                                };
                                let selected = (c.r - edit_color.r).abs() < 0.01
                                    && (c.g - edit_color.g).abs() < 0.01
                                    && (c.b - edit_color.b).abs() < 0.01;
                                button(iced::widget::Space::with_width(Length::Fill))
                                    .width(16)
                                    .height(16)
                                    .padding(0)
                                    .style(move |_, _| button::Style {
                                        background: Some(Background::Color(ic)),
                                        border: Border {
                                            color: if selected {
                                                Color::WHITE
                                            } else {
                                                Color::from_rgba(0.0, 0.0, 0.0, 0.35)
                                            },
                                            width: if selected { 2.5 } else { 1.0 },
                                            radius: 3.0.into(),
                                        },
                                        shadow: if selected {
                                            iced::Shadow {
                                                color: Color::from_rgba(0.0, 0.0, 0.0, 0.4),
                                                offset: iced::Vector::new(0.0, 0.0),
                                                blur_radius: 3.0,
                                            }
                                        } else {
                                            iced::Shadow::default()
                                        },
                                        ..Default::default()
                                    })
                                    .on_press(Message::PlasmidEditColor(c))
                                    .into()
                            };
                            let row1: Vec<Element<Message>> = FEAT_PALETTE[..8]
                                .iter()
                                .copied()
                                .map(&make_swatch)
                                .collect();
                            let row2: Vec<Element<Message>> = FEAT_PALETTE[8..]
                                .iter()
                                .copied()
                                .map(&make_swatch)
                                .collect();
                            row![
                                text("Color")
                                    .size(11)
                                    .color(palette::TEXT_DIM)
                                    .width(Length::Fixed(45.0)),
                                column![
                                    iced::widget::Row::from_vec(row1).spacing(2),
                                    iced::widget::Row::from_vec(row2).spacing(2),
                                ]
                                .spacing(2),
                            ]
                            .spacing(4)
                            .align_y(iced::Alignment::Center)
                        },
                        {
                            // Current color preview swatch + hex input field.
                            let preview_ic = Color {
                                r: edit_color.r,
                                g: edit_color.g,
                                b: edit_color.b,
                                a: 1.0,
                            };
                            row![
                                text("Hex")
                                    .size(11)
                                    .color(palette::TEXT_DIM)
                                    .width(Length::Fixed(45.0)),
                                container(iced::widget::Space::with_width(Length::Fill))
                                    .width(18)
                                    .height(18)
                                    .style(move |_| container::Style {
                                        background: Some(Background::Color(preview_ic)),
                                        border: Border {
                                            color: Color::from_rgba(0.0, 0.0, 0.0, 0.3),
                                            width: 1.0,
                                            radius: 3.0.into()
                                        },
                                        ..Default::default()
                                    }),
                                text_input("#RRGGBB", edit_hex)
                                    .on_input(Message::PlasmidEditHex)
                                    .size(11)
                                    .padding([3, 6])
                                    .width(Length::Fill),
                            ]
                            .spacing(4)
                            .align_y(iced::Alignment::Center)
                        },
                        row![
                            button(text("Apply").size(11))
                                .padding([3, 10])
                                .on_press(Message::PlasmidApplyEdit),
                            button(text("Deselect").size(11))
                                .padding([3, 8])
                                .style(button::secondary)
                                .on_press(Message::PlasmidSelectFeature(None)),
                        ]
                        .spacing(6),
                    ]
                    .spacing(6)
                    .padding([8, 8]),
                )
                .style(|_| container::Style {
                    background: Some(Background::Color(Color::from_rgb(0.97, 0.97, 1.0))),
                    border: Border {
                        color: palette::BORDER,
                        width: 1.0,
                        radius: 0.0.into()
                    },
                    ..Default::default()
                })
                .width(Length::Fill),
            ]
            .spacing(0),
        )
        .width(Length::Fill)
        .into()
    } else {
        iced::widget::Space::with_height(0).into()
    };

    // ── RE sites section ──────────────────────────────────────────────────────
    let n_re = re_sites.len();
    let max_h = (n_re as f32 * 24.0 + 4.0).min(180.0);
    let re_section: Element<Message> = if re_computing {
        container(
            text("⏳ Computing RE sites…")
                .size(11)
                .color(palette::TEXT_DIM),
        )
        .padding([6, 8])
        .width(Length::Fill)
        .into()
    } else if n_re > 0 {
        let re_header = container(
            text(format!("RE Sites  ({n_re})"))
                .size(12)
                .color(palette::TEXT),
        )
        .padding([4, 8])
        .width(Length::Fill)
        .style(|_| container::Style {
            background: Some(Background::Color(palette::HEADER_BG)),
            border: Border {
                color: palette::BORDER,
                width: 1.0,
                radius: 0.0.into(),
            },
            ..Default::default()
        });

        // Consume re_sites so widget closures own the strings.
        let mut re_rows: Vec<Element<Message>> = Vec::new();
        for (i, (name, pos, color)) in re_sites.into_iter().enumerate() {
            let bg = if i % 2 == 0 {
                Color::WHITE
            } else {
                Color::from_rgb(0.97, 0.97, 0.99)
            };
            let dot: Element<Message> = container(iced::widget::Space::with_width(8))
                .width(8)
                .height(8)
                .style(move |_| container::Style {
                    background: Some(Background::Color(color)),
                    border: Border {
                        color: Color::from_rgba(0.0, 0.0, 0.0, 0.2),
                        width: 0.5,
                        radius: 2.0.into(),
                    },
                    ..Default::default()
                })
                .into();
            let pos_str = pos.to_string();
            re_rows.push(
                container(
                    row![
                        dot,
                        text(name).size(11).color(palette::TEXT),
                        horizontal_space(),
                        text(pos_str).size(10).color(palette::TEXT_DIM),
                    ]
                    .spacing(5)
                    .align_y(iced::Alignment::Center)
                    .padding([3, 8]),
                )
                .width(Length::Fill)
                .style(move |_| container::Style {
                    background: Some(Background::Color(bg)),
                    ..Default::default()
                })
                .into(),
            );
        }

        let re_list = scrollable(column(re_rows).width(Length::Fill).spacing(0))
            .height(Length::Fixed(max_h))
            .width(Length::Fill);

        column![re_header, re_list]
            .width(Length::Fill)
            .spacing(0)
            .into()
    } else {
        iced::widget::Space::with_height(0).into()
    };

    container(
        column![header, list, edit_panel, re_section]
            .width(Length::Fill)
            .height(Length::Fill)
            .spacing(0),
    )
    .width(Length::Fixed(260.0))
    .height(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(Color::WHITE)),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .into()
}

// ── Canvas ────────────────────────────────────────────────────────────────────

struct PlasmidCanvas {
    aln: Arc<SeqAlignment>,
    seq_idx: usize,
    show_features: bool,
    show_re: bool,
    selected_feat: Option<usize>,
    zoom: f32,
    pan_x: f32,
    pan_y: f32,
    /// Pre-computed RE sites (background thread). `None` = still computing or RE off.
    re_sites: Option<Arc<Vec<RestrictionSite>>>,
}

pub struct PlasmidCanvasState {
    cache: canvas::Cache,
    last_key: Cell<u64>,
    /// Window-space position where left button was pressed (None = not pressed).
    drag_origin: Cell<Option<Point>>,
    /// Last cursor position during an active drag.
    drag_last: Cell<Option<Point>>,
}

impl Default for PlasmidCanvasState {
    fn default() -> Self {
        Self {
            cache: canvas::Cache::default(),
            last_key: Cell::new(u64::MAX),
            drag_origin: Cell::new(None),
            drag_last: Cell::new(None),
        }
    }
}

impl canvas::Program<Message> for PlasmidCanvas {
    type State = PlasmidCanvasState;

    fn update(
        &self,
        state: &mut PlasmidCanvasState,
        event: canvas::Event,
        bounds: Rectangle,
        cursor: mouse::Cursor,
    ) -> (canvas::event::Status, Option<Message>) {
        let cursor_pos = match cursor {
            mouse::Cursor::Available(p) => Some(p),
            _ => None,
        };
        let in_bounds = cursor_pos.map_or(false, |p| bounds.contains(p));

        match event {
            canvas::Event::Mouse(mouse::Event::ButtonPressed(mouse::Button::Left)) => {
                if let Some(pos) = cursor_pos {
                    if in_bounds {
                        state.drag_origin.set(Some(pos));
                        state.drag_last.set(Some(pos));
                        return (canvas::event::Status::Captured, None);
                    }
                }
            }
            canvas::Event::Mouse(mouse::Event::ButtonReleased(mouse::Button::Left)) => {
                let origin = state.drag_origin.get();
                state.drag_origin.set(None);
                state.drag_last.set(None);
                // Treat as a click only if the total movement was < 6 px
                if let (Some(orig), Some(pos)) = (origin, cursor_pos) {
                    let d = ((pos.x - orig.x).powi(2) + (pos.y - orig.y).powi(2)).sqrt();
                    if d < 6.0 && in_bounds {
                        let rel = Point::new(pos.x - bounds.x, pos.y - bounds.y);
                        let msg = self.hit_test(bounds.size(), rel);
                        return (canvas::event::Status::Captured, msg);
                    }
                }
            }
            canvas::Event::Mouse(mouse::Event::CursorMoved { position }) => {
                if let Some(last) = state.drag_last.get() {
                    let dx = position.x - last.x;
                    let dy = position.y - last.y;
                    state.drag_last.set(Some(position));
                    return (
                        canvas::event::Status::Captured,
                        Some(Message::PlasmidPanDelta(dx, dy)),
                    );
                }
            }
            canvas::Event::Mouse(mouse::Event::WheelScrolled { delta }) => {
                if in_bounds {
                    let step = match delta {
                        mouse::ScrollDelta::Lines { y, .. } => y * 0.12,
                        mouse::ScrollDelta::Pixels { y, .. } => y * 0.003,
                    };
                    return (
                        canvas::event::Status::Captured,
                        Some(Message::PlasmidZoomDelta(step)),
                    );
                }
            }
            _ => {}
        }
        (canvas::event::Status::Ignored, None)
    }

    fn draw(
        &self,
        state: &PlasmidCanvasState,
        renderer: &iced::Renderer,
        _theme: &iced::Theme,
        bounds: Rectangle,
        _cursor: mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let re_ptr = self.re_sites.as_ref().map_or(0, |a| Arc::as_ptr(a) as u64);
        let key = (Arc::as_ptr(&self.aln) as u64)
            .wrapping_add(self.seq_idx as u64 * 0x9e3779b9)
            .wrapping_add(self.show_features as u64 * 4)
            .wrapping_add(self.show_re as u64 * 2)
            .wrapping_add(self.selected_feat.map_or(u64::MAX, |i| i as u64))
            .wrapping_add(bounds.width.to_bits() as u64)
            .wrapping_add(bounds.height.to_bits() as u64)
            .wrapping_add(self.zoom.to_bits() as u64)
            .wrapping_add(self.pan_x.to_bits() as u64)
            .wrapping_add(self.pan_y.to_bits() as u64)
            .wrapping_add(re_ptr);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let aln = Arc::clone(&self.aln);
        let seq_idx = self.seq_idx;
        let show_features = self.show_features;
        let show_re = self.show_re;
        let selected_feat = self.selected_feat;
        let re_sites = self.re_sites.clone();

        // Use actual bounds size for cx/cy so the circle is visually centered.
        // The frame is snapped to 32px multiples for cache efficiency; any extra
        // pixels are clipped by iced so they're invisible.
        let visible = bounds.size();
        const SNAP: f32 = 32.0;
        let snapped = Size {
            width: (visible.width / SNAP).ceil() * SNAP,
            height: (visible.height / SNAP).ceil() * SNAP,
        };

        let zoom = self.zoom;
        let pan_x = self.pan_x;
        let pan_y = self.pan_y;
        let geom = state.cache.draw(renderer, snapped, |frame| {
            draw_plasmid(
                frame,
                snapped,
                visible,
                &aln,
                seq_idx,
                show_features,
                show_re,
                selected_feat,
                zoom,
                pan_x,
                pan_y,
                re_sites.as_ref().map(|v| v.as_slice()),
            );
        });

        vec![geom]
    }
}

impl PlasmidCanvas {
    /// Hit-test a click position. Returns `PlasmidSelectFeature` if a feature arc was clicked.
    fn hit_test(&self, size: Size, pos: Point) -> Option<Message> {
        if !self.show_features {
            return None;
        }
        let seq = self.aln.sequences.get(self.seq_idx)?;

        let raw_len: usize = seq
            .residues
            .iter()
            .filter(|&&b| !matches!(b, b'-' | b'.' | b'~'))
            .count()
            .max(1);

        let cx = size.width * 0.5 + self.pan_x;
        let cy = size.height * 0.5 + self.pan_y;
        let r = (size.width.min(size.height) * 0.5 * 0.52 * self.zoom).max(60.0);
        let feat_r = r + 18.0;
        let feat_half = 14.0 / 2.0 + 4.0;

        let dx = pos.x - cx;
        let dy = pos.y - cy;
        let dist = (dx * dx + dy * dy).sqrt();

        if (dist - feat_r).abs() > feat_half + 8.0 {
            return None;
        }

        let click_a = dy.atan2(dx);

        for (i, feat) in seq.features.iter().enumerate() {
            if matches!(feat.feature_type, FeatureType::Source) {
                continue;
            }
            let a_start = angle_of(feat.start.min(raw_len - 1), raw_len);
            let a_end = angle_of(feat.end.min(raw_len - 1), raw_len);
            if arc_contains(a_start, a_end, click_a) {
                let msg = if self.selected_feat == Some(i) {
                    Message::PlasmidSelectFeature(None)
                } else {
                    Message::PlasmidSelectFeature(Some(i))
                };
                return Some(msg);
            }
        }
        None
    }
}

/// Check if `a` falls within the arc from `start` to `end` (both in radians, CCW).
fn arc_contains(start: f32, end: f32, a: f32) -> bool {
    // Normalise all to [0, TAU)
    let norm = |x: f32| ((x % TAU) + TAU) % TAU;
    let s = norm(start);
    let e = norm(end);
    let p = norm(a);
    if s <= e {
        p >= s && p <= e
    } else {
        p >= s || p <= e
    }
}

// ── Drawing ───────────────────────────────────────────────────────────────────

fn angle_of(pos: usize, len: usize) -> f32 {
    (pos as f32 / len as f32) * TAU - PI / 2.0
}

fn pt(cx: f32, cy: f32, r: f32, a: f32) -> Point {
    Point::new(cx + r * a.cos(), cy + r * a.sin())
}

fn core_to_iced(c: helixview_core::color::Color) -> Color {
    Color {
        r: c.r,
        g: c.g,
        b: c.b,
        a: c.a,
    }
}

fn draw_plasmid(
    frame: &mut Frame<iced::Renderer>,
    _frame_size: Size, // snapped frame (for background fill)
    visible: Size,    // actual canvas bounds (for centering)
    aln: &SeqAlignment,
    seq_idx: usize,
    show_features: bool,
    show_re: bool,
    selected_feat: Option<usize>,
    zoom: f32,
    pan_x: f32,
    pan_y: f32,
    re_sites: Option<&[RestrictionSite]>,
) {
    frame.fill_rectangle(Point::ORIGIN, visible, Color::WHITE);

    let seq = match aln.sequences.get(seq_idx) {
        Some(s) => s,
        None => return,
    };

    let raw: Vec<u8> = seq
        .residues
        .iter()
        .filter(|&&b| !matches!(b, b'-' | b'.' | b'~'))
        .map(|&b| b.to_ascii_uppercase())
        .collect();
    let seq_len = raw.len().max(1);

    // Center in the VISIBLE area, then offset by pan.
    let cx = visible.width * 0.5 + pan_x;
    let cy = visible.height * 0.5 + pan_y;
    let r = (visible.width.min(visible.height) * 0.5 * 0.52 * zoom).max(60.0);

    // ── Ruler ticks ───────────────────────────────────────────────────────────
    let (major_step, minor_step) = tick_steps(seq_len);
    {
        let mut pos = 0usize;
        let mut last_major_a = -999.0f32;
        while pos <= seq_len {
            let a = angle_of(pos, seq_len);
            let is_major = pos % major_step == 0;
            let (r_inner, r_outer) = if is_major {
                (r - 10.0, r + 10.0)
            } else {
                (r - 5.0, r + 5.0)
            };
            let p0 = pt(cx, cy, r_inner, a);
            let p1 = pt(cx, cy, r_outer, a);
            let tick = Path::new(|b| {
                b.move_to(p0);
                b.line_to(p1);
            });
            let color = if is_major {
                Color::from_rgb(0.35, 0.35, 0.40)
            } else {
                Color::from_rgb(0.70, 0.70, 0.74)
            };
            frame.stroke(
                &tick,
                Stroke::default()
                    .with_color(color)
                    .with_width(if is_major { 1.5 } else { 0.8 }),
            );

            if is_major && pos > 0 && (a - last_major_a).abs() > 0.10 {
                last_major_a = a;
                let label_r = r + 18.0;
                let lp = pt(cx, cy, label_r, a);
                let h_align = if a.cos() >= 0.0 {
                    alignment::Horizontal::Left
                } else {
                    alignment::Horizontal::Right
                };
                let label = if pos >= 1_000_000 {
                    format!("{:.0}M", pos as f64 / 1_000_000.0)
                } else if major_step >= 1_000 {
                    format!("{}k", pos / 1_000)
                } else {
                    pos.to_string()
                };
                frame.fill_text(Text {
                    content: label,
                    position: lp,
                    color: Color::from_rgb(0.30, 0.30, 0.35),
                    size: Pixels((9.0 * zoom).max(6.0)),
                    font: Font::MONOSPACE,
                    horizontal_alignment: h_align,
                    vertical_alignment: alignment::Vertical::Center,
                    ..Text::default()
                });
            }
            pos += minor_step;
            if pos > seq_len && pos - minor_step < seq_len {
                pos = seq_len;
            } else if pos > seq_len {
                break;
            }
        }
    }

    // ── Backbone circle ───────────────────────────────────────────────────────
    let backbone = Path::circle(Point::new(cx, cy), r);
    frame.stroke(
        &backbone,
        Stroke::default()
            .with_color(Color::from_rgb(0.25, 0.25, 0.30))
            .with_width(2.5),
    );

    // ── Feature arcs ─────────────────────────────────────────────────────────
    if show_features {
        let feat_r = r + 18.0;
        let feat_w = 14.0;
        let arrow_l = 9.0;
        let arrow_w = 5.5;

        for (feat_idx, feat) in seq.features.iter().enumerate() {
            if matches!(feat.feature_type, FeatureType::Source) {
                continue;
            }
            if feat.start > feat.end && feat.end == 0 {
                continue;
            }

            let is_sel = selected_feat == Some(feat_idx);
            let fc = core_to_iced(feat.color);

            // Thicker stroke + glow outline for selected feature
            let stroke_w = if is_sel { feat_w + 4.0 } else { feat_w };
            if is_sel {
                let a_start = angle_of(feat.start.min(seq_len - 1), seq_len);
                let a_end = angle_of(feat.end.min(seq_len - 1), seq_len);
                let halo = Path::new(|b| {
                    b.arc(canvas::path::Arc {
                        center: Point::new(cx, cy),
                        radius: feat_r,
                        start_angle: iced::Radians(a_start),
                        end_angle: iced::Radians(a_end),
                    });
                });
                frame.stroke(
                    &halo,
                    Stroke::default()
                        .with_color(Color::from_rgba(1.0, 0.85, 0.0, 0.5))
                        .with_width(stroke_w + 6.0),
                );
            }

            let a_start = angle_of(feat.start.min(seq_len - 1), seq_len);
            let a_end = angle_of(feat.end.min(seq_len - 1), seq_len);

            let arc = Path::new(|b| {
                b.arc(canvas::path::Arc {
                    center: Point::new(cx, cy),
                    radius: feat_r,
                    start_angle: iced::Radians(a_start),
                    end_angle: iced::Radians(a_end),
                });
            });
            frame.stroke(&arc, Stroke::default().with_color(fc).with_width(stroke_w));

            let fwd = !matches!(feat.direction, FeatureDirection::Reverse);
            if !matches!(feat.direction, FeatureDirection::None) {
                let tip_a = if fwd { a_end } else { a_start };
                let tip = pt(cx, cy, feat_r, tip_a);
                let (tx, ty) = if fwd {
                    (-tip_a.sin(), tip_a.cos())
                } else {
                    (tip_a.sin(), -tip_a.cos())
                };
                let (nx, ny) = (tip_a.cos(), tip_a.sin());
                let tip_pt = Point::new(tip.x + tx * arrow_l, tip.y + ty * arrow_l);
                let b1 = Point::new(tip.x + nx * arrow_w, tip.y + ny * arrow_w);
                let b2 = Point::new(tip.x - nx * arrow_w, tip.y - ny * arrow_w);
                let arrow = Path::new(|b| {
                    b.move_to(tip_pt);
                    b.line_to(b1);
                    b.line_to(b2);
                    b.close();
                });
                frame.fill(&arrow, fc);
            }

            if !feat.name.is_empty() {
                let mid_a = (a_start + a_end) * 0.5;
                let label_r = feat_r + feat_w * 0.5 + 12.0;
                let lp = pt(cx, cy, label_r, mid_a);
                let on_right = mid_a.cos() >= 0.0;
                let h_align = if on_right {
                    alignment::Horizontal::Left
                } else {
                    alignment::Horizontal::Right
                };
                let label = feat.name.clone();
                let label_color = if is_sel {
                    Color::from_rgb(0.60, 0.42, 0.0)
                } else {
                    darken(fc, 0.35)
                };
                let label_size = ((if is_sel { 12.0_f32 } else { 11.0_f32 }) * zoom).max(8.0);
                // White background chip so text is readable on the white canvas.
                let char_w = label_size * 0.62;
                let chip_w = label.len() as f32 * char_w + 6.0;
                let chip_h = label_size + 4.0;
                let chip_x = if on_right {
                    lp.x - 3.0
                } else {
                    lp.x - chip_w + 3.0
                };
                let chip_y = lp.y - chip_h * 0.5;
                let chip_path = Path::new(|b| {
                    let r2 = (chip_h * 0.3).min(4.0);
                    b.rounded_rectangle(
                        Point::new(chip_x, chip_y),
                        Size::new(chip_w, chip_h),
                        iced::border::Radius::new(r2),
                    );
                });
                frame.fill(&chip_path, Color::from_rgba(1.0, 1.0, 1.0, 0.88));
                frame.stroke(
                    &chip_path,
                    Stroke::default()
                        .with_color(Color::from_rgba(fc.r, fc.g, fc.b, 0.55))
                        .with_width(0.8),
                );
                frame.fill_text(Text {
                    content: label,
                    position: lp,
                    color: label_color,
                    size: Pixels(label_size),
                    font: Font::MONOSPACE,
                    horizontal_alignment: h_align,
                    vertical_alignment: alignment::Vertical::Center,
                    ..Text::default()
                });
            }
        }
    }

    // ── RE sites ─────────────────────────────────────────────────────────────
    if show_re {
        let sites = re_sites.unwrap_or(&[]);
        let re_r_inner = r - 6.0;
        let re_r_outer = r + 6.0;
        // Pixel-constant thresholds: dividing by zoom keeps screen density uniform.
        let min_tick_gap = 0.003 / zoom.max(0.5);
        let label_gap = 0.15 / zoom.max(0.5);
        let mut last_tick_a = -999.0f32;
        let mut last_label_a: Vec<f32> = vec![-999.0; ENZYMES.len()];

        for site in sites {
            let a = angle_of(site.true_pos, seq_len);
            let color = RE_PALETTE[ENZYMES[site.enzyme_idx].color_idx % RE_PALETTE.len()];

            // Skip ticks that would be sub-pixel-close to the previous drawn tick.
            if (a - last_tick_a).abs() >= min_tick_gap {
                last_tick_a = a;
                let p0 = pt(cx, cy, re_r_inner, a);
                let p1 = pt(cx, cy, re_r_outer, a);
                let tick = Path::new(|b| {
                    b.move_to(p0);
                    b.line_to(p1);
                });
                frame.stroke(&tick, Stroke::default().with_color(color).with_width(1.5));
            }

            let prev_a = last_label_a[site.enzyme_idx];
            if (a - prev_a).abs() > label_gap {
                last_label_a[site.enzyme_idx] = a;
                let label_r = re_r_inner - 8.0;
                let lp = pt(cx, cy, label_r, a);
                let h_align = if a.cos() >= 0.0 {
                    alignment::Horizontal::Right
                } else {
                    alignment::Horizontal::Left
                };
                frame.fill_text(Text {
                    content: ENZYMES[site.enzyme_idx].name.to_string(),
                    position: lp,
                    color,
                    size: Pixels((8.0 * zoom).max(6.0)),
                    font: Font::MONOSPACE,
                    horizontal_alignment: h_align,
                    vertical_alignment: alignment::Vertical::Center,
                    ..Text::default()
                });
            }
        }
    }

    // ── Center label ─────────────────────────────────────────────────────────
    let name_label = if seq.name.len() > 20 {
        format!("{}…", &seq.name[..19])
    } else {
        seq.name.clone()
    };
    let size_label = format_bp(seq_len);

    frame.fill_text(Text {
        content: name_label,
        position: Point::new(cx, cy - 10.0 * zoom),
        color: Color::from_rgb(0.10, 0.10, 0.15),
        size: Pixels((13.0 * zoom).max(6.0)),
        font: Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment: alignment::Vertical::Center,
        ..Text::default()
    });
    frame.fill_text(Text {
        content: size_label,
        position: Point::new(cx, cy + 8.0 * zoom),
        color: Color::from_rgb(0.40, 0.40, 0.45),
        size: Pixels((11.0 * zoom).max(6.0)),
        font: Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment: alignment::Vertical::Center,
        ..Text::default()
    });
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn tick_steps(len: usize) -> (usize, usize) {
    if len < 1_000 {
        (100, 20)
    } else if len < 5_000 {
        (500, 100)
    } else if len < 20_000 {
        (1_000, 200)
    } else if len < 100_000 {
        (5_000, 1_000)
    } else if len < 500_000 {
        (50_000, 10_000)
    } else if len < 2_000_000 {
        (200_000, 50_000)
    } else if len < 10_000_000 {
        (500_000, 100_000)
    } else {
        (1_000_000, 200_000)
    }
}

fn format_bp(len: usize) -> String {
    if len < 1_000 {
        format!("{} bp", len)
    } else if len < 1_000_000 {
        format!("{:.2} kbp", len as f64 / 1_000.0)
    } else {
        format!("{:.2} Mbp", len as f64 / 1_000_000.0)
    }
}

fn darken(c: Color, amount: f32) -> Color {
    Color {
        r: (c.r * (1.0 - amount)).max(0.0),
        g: (c.g * (1.0 - amount)).max(0.0),
        b: (c.b * (1.0 - amount)).max(0.0),
        a: 1.0,
    }
}

// ── SVG export ────────────────────────────────────────────────────────────────

pub fn plasmid_to_svg(
    aln: &SeqAlignment,
    seq_idx: usize,
    show_features: bool,
    show_re: bool,
) -> String {
    let seq = match aln.sequences.get(seq_idx) {
        Some(s) => s,
        None => return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string(),
    };

    let raw: Vec<u8> = seq
        .residues
        .iter()
        .filter(|&&b| !matches!(b, b'-' | b'.' | b'~'))
        .map(|&b| b.to_ascii_uppercase())
        .collect();
    let seq_len = raw.len().max(1);

    let size = 800.0f32;
    let cx = size * 0.5;
    let cy = size * 0.5;
    let r = size * 0.28;
    let feat_r = r + 18.0;
    let feat_w = 14.0;

    let mut s = String::new();
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{0}\" height=\"{0}\" font-family=\"monospace\">\n",
        size as u32
    ));
    s.push_str(&format!(
        "  <rect width=\"{0}\" height=\"{0}\" fill=\"white\"/>\n",
        size as u32
    ));

    // Backbone
    s.push_str(&format!(
        "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#404048\" stroke-width=\"2.5\"/>\n",
        cx, cy, r
    ));

    // Ruler ticks
    let (major_step, minor_step) = tick_steps(seq_len);
    let mut pos = 0usize;
    loop {
        let a = angle_of(pos, seq_len);
        let is_major = pos % major_step == 0;
        let (ri, ro) = if is_major {
            (r - 10.0, r + 10.0)
        } else {
            (r - 5.0, r + 5.0)
        };
        let p0 = pt(cx, cy, ri, a);
        let p1 = pt(cx, cy, ro, a);
        let color = if is_major { "#585860" } else { "#b0b0b8" };
        let width = if is_major { 1.5 } else { 0.8 };
        s.push_str(&format!(
            "  <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"{}\"/>\n",
            p0.x, p0.y, p1.x, p1.y, color, width
        ));
        if is_major && pos > 0 {
            let lp = pt(cx, cy, r + 22.0, a);
            let anchor = if a.cos() >= 0.0 { "start" } else { "end" };
            let label = if pos >= 1_000_000 {
                format!("{:.0}M", pos as f64 / 1_000_000.0)
            } else if major_step >= 1_000 {
                format!("{}k", pos / 1_000)
            } else {
                pos.to_string()
            };
            s.push_str(&format!(
                "  <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"9\" text-anchor=\"{}\" dominant-baseline=\"middle\" fill=\"#505058\">{}</text>\n",
                lp.x, lp.y, anchor, label
            ));
        }
        pos += minor_step;
        if pos > seq_len && pos - minor_step < seq_len {
            pos = seq_len;
        } else if pos > seq_len {
            break;
        }
    }

    // Features
    if show_features {
        for feat in &seq.features {
            if matches!(feat.feature_type, FeatureType::Source) {
                continue;
            }
            let fc = feat.color;
            let hex = format!(
                "#{:02X}{:02X}{:02X}",
                (fc.r * 255.0) as u8,
                (fc.g * 255.0) as u8,
                (fc.b * 255.0) as u8
            );
            let a0 = angle_of(feat.start.min(seq_len - 1), seq_len);
            let a1 = angle_of(feat.end.min(seq_len - 1), seq_len);
            let p0 = pt(cx, cy, feat_r, a0);
            let p1 = pt(cx, cy, feat_r, a1);
            let large = if (a1 - a0).abs() > PI { 1 } else { 0 };
            s.push_str(&format!(
                "  <path d=\"M {:.1} {:.1} A {:.1} {:.1} 0 {} 1 {:.1} {:.1}\" fill=\"none\" stroke=\"{}\" stroke-width=\"{}\"/>\n",
                p0.x, p0.y, feat_r, feat_r, large, p1.x, p1.y, hex, feat_w
            ));
            if !feat.name.is_empty() {
                let mid_a = (a0 + a1) * 0.5;
                let lp = pt(cx, cy, feat_r + feat_w * 0.5 + 10.0, mid_a);
                let anchor = if mid_a.cos() >= 0.0 { "start" } else { "end" };
                let label = if feat.name.len() > 16 {
                    format!("{}…", &feat.name[..15])
                } else {
                    feat.name.clone()
                };
                let dark = format!(
                    "#{:02X}{:02X}{:02X}",
                    ((fc.r * 0.75) * 255.0) as u8,
                    ((fc.g * 0.75) * 255.0) as u8,
                    ((fc.b * 0.75) * 255.0) as u8
                );
                s.push_str(&format!(
                    "  <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"9.5\" text-anchor=\"{}\" dominant-baseline=\"middle\" fill=\"{}\">{}</text>\n",
                    lp.x, lp.y, anchor, dark, label
                ));
            }
        }
    }

    // RE sites
    if show_re {
        let sites = find_sites(&raw);
        let mut last_label_a: Vec<f32> = vec![-999.0; ENZYMES.len()];
        for site in &sites {
            let a = angle_of(site.true_pos, seq_len);
            let color = RE_PALETTE[ENZYMES[site.enzyme_idx].color_idx % RE_PALETTE.len()];
            let hex = format!(
                "#{:02X}{:02X}{:02X}",
                (color.r * 255.0) as u8,
                (color.g * 255.0) as u8,
                (color.b * 255.0) as u8
            );
            let p0 = pt(cx, cy, r - 6.0, a);
            let p1 = pt(cx, cy, r + 6.0, a);
            s.push_str(&format!(
                "  <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                p0.x, p0.y, p1.x, p1.y, hex
            ));
            let prev_a = last_label_a[site.enzyme_idx];
            if (a - prev_a).abs() > 0.18 {
                last_label_a[site.enzyme_idx] = a;
                let lp = pt(cx, cy, r - 14.0, a);
                let anchor = if a.cos() >= 0.0 { "end" } else { "start" };
                s.push_str(&format!(
                    "  <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"8\" text-anchor=\"{}\" dominant-baseline=\"middle\" fill=\"{}\">{}</text>\n",
                    lp.x, lp.y, anchor, hex, ENZYMES[site.enzyme_idx].name
                ));
            }
        }
    }

    // Center label
    s.push_str(&format!(
        "  <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"13\" text-anchor=\"middle\" dominant-baseline=\"middle\" fill=\"#1a1a26\">{}</text>\n",
        cx, cy - 10.0,
        if seq.name.len() > 20 { format!("{}…", &seq.name[..19]) } else { seq.name.clone() }
    ));
    s.push_str(&format!(
        "  <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"11\" text-anchor=\"middle\" dominant-baseline=\"middle\" fill=\"#666672\">{}</text>\n",
        cx, cy + 8.0, format_bp(seq_len)
    ));

    s.push_str("</svg>\n");
    s
}
