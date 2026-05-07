//! Text alignment export options dialog.

use crate::app::Message;
use crate::text_export::TextExportOptions;
use crate::theme::palette;
use iced::{
    widget::{button, checkbox, column, container, row, slider, text},
    Alignment, Background, Color, Element, Length,
};

pub fn text_export_dialog<'a>(opts: &TextExportOptions, preview: &str) -> Element<'a, Message> {
    // ── Controls toolbar ──────────────────────────────────────────────────────
    let rpr_ctrl = row![
        text("Residues/row:").size(12).color(palette::TEXT_DIM),
        slider(10u32..=120, opts.residues_per_row as u32, |v| {
            Message::TextExportResPerRow(v as usize)
        })
        .width(100),
        text(format!("{}", opts.residues_per_row))
            .size(12)
            .color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let tc_ctrl = row![
        text("Title width:").size(12).color(palette::TEXT_DIM),
        slider(4u32..=40, opts.title_chars as u32, |v| {
            Message::TextExportTitleChars(v as usize)
        })
        .width(100),
        text(format!("{}", opts.title_chars))
            .size(12)
            .color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let ruler_cb = checkbox("Ruler", opts.show_ruler)
        .on_toggle(Message::TextExportShowRuler)
        .size(14)
        .text_size(12);

    let num_cb = checkbox("Position numbers", opts.number_lines)
        .on_toggle(Message::TextExportNumberLines)
        .size(14)
        .text_size(12);

    let save_btn = button(text("Save as .txt…").size(12))
        .padding([4, 12])
        .on_press(Message::TextExportSave);

    let close_btn = button(text("< Back").size(12))
        .padding([4, 10])
        .on_press(Message::CloseTextExport);

    let toolbar = container(
        row![
            close_btn,
            text("Export Alignment as Text")
                .size(14)
                .color(palette::TEXT),
            iced::widget::horizontal_space(),
            rpr_ctrl,
            tc_ctrl,
            ruler_cb,
            num_cb,
            save_btn,
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

    // ── Preview ───────────────────────────────────────────────────────────────
    // Show first ~80 lines of the formatted output so user can tune options.
    let preview_text: String = preview.lines().take(80).collect::<Vec<_>>().join("\n");

    let preview_widget: Element<Message> = iced::widget::scrollable(
        container(
            text(preview_text)
                .size(11)
                .font(iced::Font::MONOSPACE)
                .color(Color::from_rgb(0.15, 0.15, 0.20)),
        )
        .padding([8, 12])
        .width(Length::Fill),
    )
    .width(Length::Fill)
    .height(Length::Fill)
    .into();

    column![toolbar, preview_widget]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
