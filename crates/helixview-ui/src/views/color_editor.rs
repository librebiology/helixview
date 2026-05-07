//! Color table editor — lets the user remap residue background colors.

use std::sync::Arc;

use iced::widget::horizontal_space;
use iced::{
    widget::{button, column, container, row, slider, text},
    Alignment, Background, Border, Color, Element, Length,
};

use crate::app::Message;
use crate::color_table::{ColorTable, AA_RESIDUES, NUC_RESIDUES};
use crate::theme::palette;

// ── Public entry point ────────────────────────────────────────────────────────

pub fn color_editor_view<'a>(
    table: &'a Arc<ColorTable>,
    selected: Option<(u8, bool)>, // (residue_byte, is_nuc)
) -> Element<'a, Message> {
    // ── Toolbar ───────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 10])
        .on_press(Message::CloseColorEditor);

    let reset_btn = button(text("Reset to Defaults").size(12))
        .padding([4, 10])
        .on_press(Message::ResetColorTable);

    let toolbar = container(
        row![
            back_btn,
            text("Color Table Editor").size(14).color(palette::TEXT),
            horizontal_space(),
            reset_btn,
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

    // ── Swatch sections ───────────────────────────────────────────────────────
    let nuc_section = swatch_section("Nucleotides", NUC_RESIDUES, true, table, selected);
    let aa_section = swatch_section("Amino Acids", AA_RESIDUES, false, table, selected);

    // ── RGB sliders for selected residue ──────────────────────────────────────
    let editor_panel: Element<Message> = match selected {
        Some((b, is_nuc)) => rgb_editor(table, b, is_nuc),
        None => container(
            text("Click a residue swatch to edit its color.")
                .size(12)
                .color(palette::TEXT_DIM),
        )
        .padding([12, 16])
        .into(),
    };

    let hint = container(
        text("Changes apply immediately.  'Reset to Defaults' restores the original palette.")
            .size(11)
            .color(palette::TEXT_DIM),
    )
    .padding([4, 10])
    .width(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(palette::BG_PANEL)),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    });

    column![toolbar, nuc_section, aa_section, editor_panel, hint]
        .width(Length::Fill)
        .height(Length::Fill)
        .spacing(0)
        .into()
}

// ── Swatch row ────────────────────────────────────────────────────────────────

fn swatch_section<'a>(
    label: &'a str,
    residues: &'static [u8],
    is_nuc: bool,
    table: &'a Arc<ColorTable>,
    selected: Option<(u8, bool)>,
) -> Element<'a, Message> {
    let header = container(text(label).size(12).color(palette::TEXT))
        .padding([5, 10])
        .width(Length::Fill)
        .style(|_| container::Style {
            background: Some(Background::Color(palette::HEADER_BG)),
            ..Default::default()
        });

    let mut swatches: Vec<Element<Message>> = Vec::new();
    for &b in residues {
        if b == b' ' {
            swatches.push(iced::widget::Space::with_width(8).into());
            continue;
        }
        let color = if is_nuc {
            table.nuc_color(b)
        } else {
            table.aa_color(b)
        };
        let is_selected = selected == Some((b, is_nuc));

        let border_color = if is_selected {
            Color::from_rgb(0.05, 0.35, 0.85)
        } else {
            Color::from_rgb(0.65, 0.65, 0.65)
        };
        let border_w = if is_selected { 2.5 } else { 1.0 };

        let swatch = button(
            container(
                text(b as char)
                    .size(14)
                    .color(Color::from_rgb(0.1, 0.1, 0.1)),
            )
            .width(34)
            .height(34)
            .align_x(Alignment::Center)
            .align_y(Alignment::Center)
            .style(move |_| container::Style {
                background: Some(Background::Color(color)),
                border: Border {
                    color: border_color,
                    width: border_w,
                    radius: 4.0.into(),
                },
                ..Default::default()
            }),
        )
        .padding(0)
        .on_press(Message::SelectColorResidue(b, is_nuc))
        .style(|_theme, _status| iced::widget::button::Style {
            background: None,
            ..Default::default()
        });

        swatches.push(swatch.into());
    }

    let swatch_row = container(
        row(swatches)
            .spacing(4)
            .align_y(Alignment::Center)
            .padding([8, 10]),
    )
    .width(Length::Fill);

    column![header, swatch_row].width(Length::Fill).into()
}

// ── RGB editor ────────────────────────────────────────────────────────────────

fn rgb_editor<'a>(table: &'a Arc<ColorTable>, b: u8, is_nuc: bool) -> Element<'a, Message> {
    let color = if is_nuc {
        table.nuc_color(b)
    } else {
        table.aa_color(b)
    };

    let r = color.r;
    let g = color.g;
    let bl = color.b;

    let hex = format!(
        "#{:02X}{:02X}{:02X}",
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (bl * 255.0) as u8,
    );

    let label = format!(
        "Residue: {}  ({})",
        b as char,
        if is_nuc { "nucleotide" } else { "amino acid" }
    );

    let preview = container(iced::widget::Space::new(60, 36)).style(move |_| container::Style {
        background: Some(Background::Color(color)),
        border: Border {
            color: Color::from_rgb(0.4, 0.4, 0.4),
            width: 1.0,
            radius: 4.0.into(),
        },
        ..Default::default()
    });

    let mk_slider = move |channel: u8, val: f32| {
        let lbl = match channel {
            0 => "R",
            1 => "G",
            _ => "B",
        };
        let msg_fn = move |v: f32| Message::SetResidueColor {
            residue: b,
            is_nuc,
            channel,
            value: v,
        };
        row![
            text(lbl).size(12).color(palette::TEXT_DIM).width(16),
            slider(0.0f32..=1.0, val, msg_fn).step(0.004).width(200),
            text(format!("{:.2}", val))
                .size(11)
                .color(palette::TEXT_DIM)
                .width(40),
        ]
        .spacing(8)
        .align_y(Alignment::Center)
    };

    container(
        row![
            preview,
            column![
                text(label).size(12).color(palette::TEXT),
                mk_slider(0, r),
                mk_slider(1, g),
                mk_slider(2, bl),
                text(hex).size(11).color(palette::TEXT_DIM),
            ]
            .spacing(6),
        ]
        .spacing(16)
        .align_y(Alignment::Center)
        .padding([12, 16]),
    )
    .width(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgb(0.97, 0.97, 1.0))),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .into()
}
