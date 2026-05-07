//! BLAST results view — shows hits table with import buttons.

use crate::app::{BlastState, Message};
use crate::blast::BlastHit;
use crate::theme::palette;
use iced::{
    widget::{button, column, container, row, scrollable, text},
    Background, Border, Element, Font, Length,
};

/// Render the BLAST panel as a full-screen results view.
pub fn blast_view<'a>(state: &'a BlastState, _rid: Option<&'a str>) -> Element<'a, Message> {
    let back_btn = button(text("< Back").size(13))
        .padding([4, 12])
        .on_press(Message::CloseBlast);

    let header_text = match state {
        BlastState::Idle => "BLAST — Not submitted".to_string(),
        BlastState::Submitted(r) => format!("BLAST — Running  (RID: {r})"),
        BlastState::Complete(hits) => format!("BLAST — {} hit(s) found", hits.len()),
        BlastState::Failed(e) => format!("BLAST — Error: {e}"),
    };

    let check_btn: Option<Element<'a, Message>> = if matches!(state, BlastState::Submitted(_)) {
        Some(
            button(text("Check Results").size(12))
                .padding([3, 10])
                .on_press(Message::BlastCheck)
                .into(),
        )
    } else {
        None
    };

    let header = container(
        row(std::iter::once(back_btn.into())
            .chain(std::iter::once(iced::widget::horizontal_space().into()))
            .chain(std::iter::once(
                text(header_text)
                    .size(14)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE)
                    .into(),
            ))
            .chain(std::iter::once(iced::widget::horizontal_space().into()))
            .chain(check_btn.into_iter())
            .collect::<Vec<Element<'a, Message>>>())
        .spacing(8)
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

    let body: Element<'a, Message> = match state {
        BlastState::Idle => container(
            text("No BLAST search submitted yet.")
                .size(12)
                .color(palette::TEXT_DIM)
                .font(Font::MONOSPACE),
        )
        .padding([20, 10])
        .into(),
        BlastState::Submitted(r) => container(
            column![
                text(format!("Waiting for NCBI results…  RID: {r}"))
                    .size(12)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
                text("Click 'Check Results' to poll for completion (usually 20–60 seconds).")
                    .size(11)
                    .color(palette::TEXT_DIM)
                    .font(Font::MONOSPACE),
            ]
            .spacing(6),
        )
        .padding([20, 10])
        .into(),
        BlastState::Failed(e) => container(
            text(format!("Error: {e}"))
                .size(12)
                .color(iced::Color::from_rgb(0.8, 0.1, 0.1))
                .font(Font::MONOSPACE),
        )
        .padding([20, 10])
        .into(),
        BlastState::Complete(hits) if hits.is_empty() => container(
            text("No significant hits found.")
                .size(12)
                .color(palette::TEXT_DIM)
                .font(Font::MONOSPACE),
        )
        .padding([20, 10])
        .into(),
        BlastState::Complete(hits) => hits_table(hits),
    };

    column![header, body]
        .width(Length::Fill)
        .height(Length::Fill)
        .spacing(0)
        .into()
}

fn hits_table<'a>(hits: &'a [BlastHit]) -> Element<'a, Message> {
    // Column header
    let table_header = container(
        row![
            text(format!("{:<4}", "#"))
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
            text(format!("{:<22}", "Accession"))
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
            text(format!("{:<10}", "Score"))
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
            text(format!("{:<12}", "E-value"))
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
            text(format!("{:<8}", "%ID"))
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
            text(format!("{:<10}", "Len"))
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
            text("Import")
                .size(10)
                .color(palette::ACCENT)
                .font(Font::MONOSPACE),
        ]
        .spacing(4)
        .padding([3, 10]),
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

    for (i, hit) in hits.iter().enumerate() {
        let pct_id = if hit.align_len > 0 {
            format!("{:.0}%", hit.identity as f64 / hit.align_len as f64 * 100.0)
        } else {
            "—".to_string()
        };
        let evalue_str = if hit.evalue < 1e-99 {
            format!("{:.0e}", hit.evalue)
        } else {
            format!("{:.2e}", hit.evalue)
        };

        let acc = hit.accession.clone();
        let bg = if i % 2 == 0 {
            palette::BG_SEQ
        } else {
            palette::BG_SEQ_ALT
        };

        let row_el = container(
            row![
                text(format!("{:<4}", hit.num))
                    .size(11)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
                text(format!(
                    "{:<22}",
                    &hit.accession[..hit.accession.len().min(20)]
                ))
                .size(11)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
                text(format!("{:<10.0}", hit.bit_score))
                    .size(11)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
                text(format!("{:<12}", evalue_str))
                    .size(11)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
                text(format!("{:<8}", pct_id))
                    .size(11)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
                text(format!("{:<10}", hit.align_len))
                    .size(11)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
                button(text("Import").size(10))
                    .padding([1, 6])
                    .on_press(Message::BlastImportHit(acc)),
            ]
            .spacing(4)
            .padding([3, 10])
            .align_y(iced::Alignment::Center)
            .width(Length::Fill),
        )
        .width(Length::Fill)
        .style(move |_| container::Style {
            background: Some(Background::Color(bg)),
            ..Default::default()
        });

        rows.push(row_el.into());
    }

    scrollable(column(rows).spacing(0).width(Length::Fill))
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
