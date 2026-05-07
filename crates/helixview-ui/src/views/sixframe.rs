//! Six-frame translation view.
//!
//! Displays all 6 reading-frame translations for a single nucleotide sequence.
//! Forward frames (+1, +2, +3) are shown first, then reverse-complement frames (-1, -2, -3).
//! Protein sequences are chunked into lines of 60 amino acids with position numbers.

use iced::{
    Background, Border, Element, Font, Length,
    widget::{button, column, container, row, scrollable, text},
};

use crate::app::Message;
use crate::theme::palette;

const LINE_WIDTH: usize = 60;

/// Render the six-frame translation panel.
pub fn sixframe_view<'a>(
    seq_name: &'a str,
    nt_len: usize,
    frames: &'a [(i8, Vec<u8>); 6],
) -> Element<'a, Message> {
    // ── Header ────────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 12])
        .on_press(Message::CloseSixFrame);

    let header = container(
        row![
            back_btn,
            iced::widget::horizontal_space(),
            text(format!("6-Frame Translation — {seq_name}  ({nt_len} nt)"))
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

    // ── Body ──────────────────────────────────────────────────────────────────
    let mut body_rows: Vec<Element<'a, Message>> = Vec::new();

    for (frame_label, protein) in frames.iter() {
        let is_forward = *frame_label > 0;
        let label_color = if is_forward { palette::ACCENT } else { palette::TEXT_DIM };

        let aa_count = protein.iter().filter(|&&b| b != b'*').count();
        let stop_count = protein.iter().filter(|&&b| b == b'*').count();

        // Frame header row
        let frame_header = container(
            row![
                text(format!(
                    "Frame {:+2}   {} aa   {} stop codon{}",
                    frame_label, aa_count, stop_count,
                    if stop_count == 1 { "" } else { "s" },
                ))
                .size(12)
                .color(label_color)
                .font(Font::MONOSPACE),
            ]
            .padding([3, 10]),
        )
        .width(Length::Fill)
        .style(|_| container::Style {
            background: Some(Background::Color(palette::BG_PANEL)),
            border: Border { color: palette::BORDER, width: 1.0, radius: 0.0.into() },
            ..Default::default()
        });

        body_rows.push(frame_header.into());

        if protein.is_empty() {
            body_rows.push(
                container(
                    text("  (sequence too short for this frame)")
                        .size(11)
                        .color(palette::TEXT_DIM)
                        .font(Font::MONOSPACE),
                )
                .padding([2, 10])
                .width(Length::Fill)
                .into(),
            );
        } else {
            for (chunk_idx, chunk) in protein.chunks(LINE_WIDTH).enumerate() {
                let pos = chunk_idx * LINE_WIDTH + 1;
                let seq_str = String::from_utf8_lossy(chunk).into_owned();
                let line = format!("  {:>6}  {}", pos, seq_str);
                body_rows.push(
                    text(line)
                        .size(12)
                        .color(palette::TEXT)
                        .font(Font::MONOSPACE)
                        .into(),
                );
            }
        }

        // Spacer between frames
        body_rows.push(
            container(iced::widget::horizontal_space())
                .height(4)
                .width(Length::Fill)
                .into(),
        );
    }

    let body_col = column(body_rows).width(Length::Fill).spacing(0).padding([6, 0]);
    let body = scrollable(body_col)
        .width(Length::Fill)
        .height(Length::Fill);

    column![header, body]
        .width(Length::Fill)
        .height(Length::Fill)
        .spacing(0)
        .into()
}
