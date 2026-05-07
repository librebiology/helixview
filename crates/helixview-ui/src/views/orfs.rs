//! ORF (Open Reading Frame) finder results view.

use helixview_analysis::Orf;
use iced::{
    widget::{button, column, container, row, scrollable, text},
    Background, Border, Element, Font, Length,
};

use crate::app::Message;
use crate::theme::palette;

/// Display ORF finder results for a given sequence.
pub fn orfs_view<'a>(seq_name: &'a str, orfs: &'a [Orf]) -> Element<'a, Message> {
    // ── Header bar ────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 12])
        .on_press(Message::CloseOrfs);

    let header = container(
        row![
            back_btn,
            iced::widget::horizontal_space(),
            text(format!("ORF Finder — {seq_name}"))
                .size(14)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
            iced::widget::horizontal_space(),
        ]
        .align_y(iced::Alignment::Center)
        .padding([4, 10])
        .width(Length::Fill),
    )
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

    // ── Summary line ──────────────────────────────────────────────────────────
    let n = orfs.len();
    let summary = container(
        text(format!("{n} ORF{} found", if n == 1 { "" } else { "s" }))
            .size(12)
            .color(palette::TEXT_DIM)
            .font(Font::MONOSPACE),
    )
    .padding([6, 10])
    .width(Length::Fill);

    // ── Body: empty state or table ────────────────────────────────────────────
    let body: Element<'a, Message> = if orfs.is_empty() {
        container(
            text("No ORFs found with the current minimum length (25 aa)")
                .size(12)
                .color(palette::TEXT_DIM)
                .font(Font::MONOSPACE),
        )
        .width(Length::Fill)
        .padding([20, 10])
        .align_x(iced::alignment::Horizontal::Center)
        .into()
    } else {
        // Table header row
        let col_frame = text("Frame ")
            .size(11)
            .color(palette::ACCENT)
            .font(Font::MONOSPACE);
        let col_start = text(format!("{:<10}", "Start"))
            .size(11)
            .color(palette::ACCENT)
            .font(Font::MONOSPACE);
        let col_end = text(format!("{:<10}", "End"))
            .size(11)
            .color(palette::ACCENT)
            .font(Font::MONOSPACE);
        let col_len = text(format!("{:<10}", "Length"))
            .size(11)
            .color(palette::ACCENT)
            .font(Font::MONOSPACE);
        let col_prot = text("Protein (first 40 aa)")
            .size(11)
            .color(palette::ACCENT)
            .font(Font::MONOSPACE);

        let table_header = container(
            row![col_frame, col_start, col_end, col_len, col_prot]
                .spacing(8)
                .padding([4, 10])
                .width(Length::Fill),
        )
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

        let mut rows: Vec<Element<'a, Message>> = vec![table_header.into()];

        for (i, orf) in orfs.iter().enumerate() {
            let frame_str = format!("{:+}  ", orf.frame,);
            let start_str = format!("{:<10}", orf.start + 1);
            let end_str = format!("{:<10}", orf.end);
            let aa_len = orf.protein.len();
            let len_str = format!("{:<10}", format!("{aa_len} aa"));
            let preview_len = 40.min(orf.protein.len());
            let prot_str = String::from_utf8_lossy(&orf.protein[..preview_len]).into_owned();

            let bg_color = if i % 2 == 0 {
                palette::BG_SEQ
            } else {
                palette::BG_SEQ_ALT
            };

            let row_el = container(
                row![
                    text(frame_str)
                        .size(12)
                        .color(palette::TEXT)
                        .font(Font::MONOSPACE),
                    text(start_str)
                        .size(12)
                        .color(palette::TEXT)
                        .font(Font::MONOSPACE),
                    text(end_str)
                        .size(12)
                        .color(palette::TEXT)
                        .font(Font::MONOSPACE),
                    text(len_str)
                        .size(12)
                        .color(palette::TEXT)
                        .font(Font::MONOSPACE),
                    text(prot_str)
                        .size(12)
                        .color(palette::TEXT)
                        .font(Font::MONOSPACE),
                ]
                .spacing(8)
                .padding([3, 10])
                .width(Length::Fill),
            )
            .width(Length::Fill)
            .style(move |_| container::Style {
                background: Some(Background::Color(bg_color)),
                ..Default::default()
            });

            rows.push(row_el.into());
        }

        let table = column(rows).width(Length::Fill).spacing(0);
        scrollable(table)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    };

    // ── Assemble ──────────────────────────────────────────────────────────────
    column![header, summary, body]
        .width(Length::Fill)
        .height(Length::Fill)
        .spacing(0)
        .into()
}
