//! Pairwise alignment result view.

use helixview_analysis::PairwiseResult;
use iced::{
    widget::{button, column, container, row, scrollable, text},
    Background, Border, Element, Font, Length,
};

use crate::app::Message;
use crate::theme::palette;

const LINE_WIDTH: usize = 60;

/// Display a pairwise alignment result in a clean, scrollable layout.
pub fn pairwise_view<'a>(result: &'a PairwiseResult) -> Element<'a, Message> {
    let total = result.aligned_a.len();

    // Header
    let back_btn = button(text("< Back").size(13))
        .padding([4, 12])
        .on_press(Message::ClosePairwise);

    let header = container(
        row![
            back_btn,
            iced::widget::horizontal_space(),
            text("Pairwise Alignment")
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

    // Score / identity summary row
    let summary = container(
        row![
            text(format!("Score: {}", result.score))
                .size(12)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
            text("   "),
            text(format!("Identity: {:.1}%", result.identity * 100.0))
                .size(12)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
        ]
        .spacing(0),
    )
    .padding([6, 10])
    .width(Length::Fill);

    // Build alignment block lines
    let mut blocks: Vec<Element<'a, Message>> = Vec::new();

    let n_chunks = if total == 0 {
        0
    } else {
        (total + LINE_WIDTH - 1) / LINE_WIDTH
    };

    for chunk in 0..n_chunks {
        let start = chunk * LINE_WIDTH;
        let end = (start + LINE_WIDTH).min(total);

        let chunk_a = &result.aligned_a[start..end];
        let chunk_b = &result.aligned_b[start..end];

        // Build conservation row
        let conservation: String = chunk_a
            .iter()
            .zip(chunk_b.iter())
            .map(|(&a, &b)| {
                let a_gap = a == b'-';
                let b_gap = b == b'-';
                if a_gap && b_gap {
                    ' '
                } else if a_gap || b_gap {
                    ' '
                } else if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
                    '|'
                } else {
                    '.'
                }
            })
            .collect();

        let seq_a_str = String::from_utf8_lossy(chunk_a).into_owned();
        let seq_b_str = String::from_utf8_lossy(chunk_b).into_owned();

        // Position labels: 1-based
        let pos_start = start + 1;
        let pos_end = end;

        let label_start = text(format!("{:<6}", pos_start))
            .size(11)
            .color(palette::TEXT_DIM)
            .font(Font::MONOSPACE);
        let label_end = text(format!("{:>6}", pos_end))
            .size(11)
            .color(palette::TEXT_DIM)
            .font(Font::MONOSPACE);

        // Seq A line
        let line_a = row![
            label_start,
            text("  "),
            text(seq_a_str)
                .size(12)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
            text("  "),
            label_end,
        ]
        .spacing(0);

        // Conservation line (indented to align with sequence)
        let conservation_line = row![
            text("      "), // 6 chars for start label
            text("  "),
            text(conservation)
                .size(12)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
        ]
        .spacing(0);

        // Seq B line (no end label needed — same range)
        let line_b = row![
            text(format!("{:<6}", pos_start))
                .size(11)
                .color(palette::TEXT_DIM)
                .font(Font::MONOSPACE),
            text("  "),
            text(seq_b_str)
                .size(12)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
        ]
        .spacing(0);

        let block = container(column![line_a, conservation_line, line_b].spacing(1))
            .padding([4, 10])
            .width(Length::Fill);

        blocks.push(block.into());
    }

    if total == 0 {
        blocks.push(
            container(
                text("(no alignment)")
                    .size(12)
                    .color(palette::TEXT_DIM)
                    .font(Font::MONOSPACE),
            )
            .padding([8, 10])
            .into(),
        );
    }

    let body = column(blocks).width(Length::Fill).spacing(2);

    let scroll_content = column![header, summary, body]
        .width(Length::Fill)
        .spacing(0);

    scrollable(scroll_content)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
