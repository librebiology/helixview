//! Conservation region search results view.

use std::sync::Arc;

use iced::{
    widget::{button, column, container, row, scrollable, slider, text},
    Alignment, Background, Border, Color, Element, Font, Length,
};
use helixview_core::Alignment as SeqAlignment;
use helixview_analysis::find_conserved_regions;

use crate::app::Message;
use crate::theme::palette;

pub fn conservation_view<'a>(
    aln:       &'a Arc<SeqAlignment>,
    threshold: f64,
    min_width: usize,
) -> Element<'a, Message> {
    let regions = find_conserved_regions(aln, threshold, min_width);

    // ── Toolbar ───────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 10])
        .on_press(Message::CloseConservation);

    let thresh_pct = (threshold * 100.0).round() as u32;
    let thresh_ctrl = row![
        text("Identity ≥").size(12).color(palette::TEXT_DIM),
        slider(50u32..=100, thresh_pct,
            |v| Message::ConservationThreshold(v as f64 / 100.0))
            .width(120),
        text(format!("{thresh_pct}%")).size(12).color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let width_ctrl = row![
        text("Min width").size(12).color(palette::TEXT_DIM),
        slider(1u32..=50, min_width as u32,
            |v| Message::ConservationMinWidth(v as usize))
            .width(80),
        text(format!("{min_width} cols")).size(12).color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let toolbar = container(
        row![
            back_btn,
            text("Conserved Regions").size(14).color(palette::TEXT),
            iced::widget::horizontal_space(),
            thresh_ctrl,
            width_ctrl,
            text(format!("{} region(s) found", regions.len()))
                .size(12).color(palette::TEXT_DIM),
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

    // ── Results table ─────────────────────────────────────────────────────────
    let content: Element<Message> = if regions.is_empty() {
        container(
            text(format!(
                "No conserved regions found with identity ≥ {thresh_pct}% and width ≥ {min_width} columns.\n\
                 Try lowering the threshold or minimum width."
            ))
            .size(13)
            .color(palette::TEXT_DIM),
        )
        .padding([20, 20])
        .into()
    } else {
        let mut rows: Vec<Element<Message>> = Vec::new();

        // Header
        rows.push(
            container(
                row![
                    header_cell("Start",    70),
                    header_cell("End",      70),
                    header_cell("Width",    60),
                    header_cell("Identity", 80),
                    header_cell("Consensus (first 60 chars)", 500),
                ]
                .spacing(0),
            )
            .padding([4, 8])
            .width(Length::Fill)
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(0.88, 0.90, 0.95))),
                ..Default::default()
            })
            .into(),
        );

        for (i, region) in regions.iter().enumerate() {
            let cons_str: String = region.consensus.iter()
                .take(60)
                .map(|&b| b as char)
                .collect();
            let suffix = if region.consensus.len() > 60 { "…" } else { "" };

            let identity_color = if region.avg_identity >= 0.95 {
                Color::from_rgb(0.10, 0.55, 0.20)
            } else if region.avg_identity >= 0.80 {
                Color::from_rgb(0.75, 0.55, 0.05)
            } else {
                palette::TEXT
            };

            let jump_btn: Element<Message> = button(text("→").size(11))
                .padding([2, 5])
                .on_press(Message::JumpToColumn(region.start))
                .into();

            let bg = if i % 2 == 0 {
                Color::WHITE
            } else {
                Color::from_rgb(0.96, 0.97, 1.0)
            };

            rows.push(
                container(
                    row![
                        data_cell(format!("{}", region.start + 1), 70),
                        data_cell(format!("{}", region.end + 1),   70),
                        data_cell(format!("{}", region.width()),    60),
                        text(format!("{:.1}%", region.avg_identity * 100.0))
                            .size(12)
                            .color(identity_color)
                            .width(80),
                        text(format!("{cons_str}{suffix}"))
                            .size(12)
                            .color(palette::TEXT)
                            .font(Font::MONOSPACE)
                            .width(Length::Fill),
                        jump_btn,
                    ]
                    .spacing(0)
                    .align_y(Alignment::Center),
                )
                .padding([4, 8])
                .width(Length::Fill)
                .style(move |_| container::Style {
                    background: Some(Background::Color(bg)),
                    border: Border {
                        color: Color::from_rgb(0.88, 0.88, 0.90),
                        width: 0.0,
                        radius: 0.0.into(),
                    },
                    ..Default::default()
                })
                .into(),
            );
        }

        scrollable(column(rows).width(Length::Fill))
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    };

    column![toolbar, content]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

fn header_cell(label: &str, w: u16) -> Element<'_, Message> {
    text(label)
        .size(12)
        .color(palette::TEXT)
        .font(Font::MONOSPACE)
        .width(w)
        .into()
}

fn data_cell(label: String, w: u16) -> Element<'static, Message> {
    text(label)
        .size(12)
        .color(palette::TEXT)
        .font(Font::MONOSPACE)
        .width(w)
        .into()
}
