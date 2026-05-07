//! Options dialog shown before exporting a shaded alignment graphic.

use iced::{
    widget::{button, checkbox, column, container, row, slider, text, horizontal_space},
    Alignment, Element, Length,
};
use crate::app::{ColorScheme, Message, ShadedExportOptions};
use crate::theme::palette;

pub fn shaded_export_dialog<'a>(opts: &ShadedExportOptions) -> Element<'a, Message> {
    let title = text("Export Shaded Graphic").size(15);

    // ── Shading mode buttons ──────────────────────────────────────────────────
    let scheme_btn = |label: &'static str, scheme: ColorScheme| {
        let active = opts.scheme == scheme;
        button(text(label).size(12))
            .padding([4, 10])
            .style(if active { button::primary } else { button::secondary })
            .on_press(Message::ShadedExportSetScheme(scheme))
    };
    let scheme_row = row![
        scheme_btn("Residue Type", ColorScheme::ResidueType),
        scheme_btn("Identity",     ColorScheme::Identity),
        scheme_btn("Strength",     ColorScheme::Strength),
        scheme_btn("Plain",        ColorScheme::Plain),
    ].spacing(6).align_y(Alignment::Center);

    // ── Identity threshold ────────────────────────────────────────────────────
    let threshold_row = row![
        text("Identity threshold:").size(12),
        slider(0.0..=1.0, opts.threshold, Message::ShadedExportThreshold)
            .step(0.05)
            .width(Length::Fixed(160.0)),
        text(format!("{:.0}%", opts.threshold * 100.0)).size(12).color(palette::TEXT_DIM),
    ].spacing(8).align_y(Alignment::Center);

    // ── Options ───────────────────────────────────────────────────────────────
    let ruler_cb = checkbox("Ruler (position numbers)", opts.show_ruler)
        .on_toggle(Message::ShadedExportToggleRuler)
        .size(14)
        .text_size(12);
    let consensus_cb = checkbox("Append consensus row", opts.show_consensus)
        .on_toggle(Message::ShadedExportToggleConsensus)
        .size(14)
        .text_size(12);
    let selected_note = text(
        "Exports selected sequences and column range (or all if nothing selected)."
    ).size(11).color(palette::TEXT_DIM);

    // ── Action buttons ────────────────────────────────────────────────────────
    let export_btn = button(text("Export SVG…").size(12))
        .padding([5, 14])
        .style(button::primary)
        .on_press(Message::ShadedExportRun);
    let cancel_btn = button(text("Cancel").size(12))
        .padding([5, 14])
        .style(button::secondary)
        .on_press(Message::CloseShadedExportDialog);

    let btn_row = row![cancel_btn, horizontal_space(), export_btn]
        .align_y(Alignment::Center);

    let body = column![
        title,
        text("Shading mode:").size(11).color(palette::TEXT_DIM),
        scheme_row,
        threshold_row,
        ruler_cb,
        consensus_cb,
        selected_note,
        btn_row,
    ]
    .spacing(12)
    .padding(20)
    .width(Length::Fixed(480.0));

    container(body)
        .width(Length::Fill)
        .height(Length::Fill)
        .center_x(Length::Fill)
        .center_y(Length::Fill)
        .into()
}
