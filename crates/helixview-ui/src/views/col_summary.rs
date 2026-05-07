//! Positional nucleotide/amino-acid summary — per-column frequency table (§5.16).

use std::sync::Arc;

use iced::{
    Background, Border, Color, Element, Font, Length,
    widget::{button, column, container, row, scrollable, text},
};
use helixview_core::Alignment;

use crate::app::Message;
use crate::theme::palette;

const MAX_COLS: usize = 200;

/// Scrollable table of per-column residue frequencies.
pub fn col_summary_view<'a>(aln: &'a Arc<Alignment>) -> Element<'a, Message> {
    let back_btn = button(text("< Back").size(13))
        .padding([4, 12])
        .on_press(Message::CloseColSummary);

    let header = container(
        row![
            back_btn,
            iced::widget::horizontal_space(),
            text("Positional Residue Summary")
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

    let ncols = aln.col_count();
    let nrows = aln.seq_count();

    if ncols == 0 {
        return column![
            header,
            container(text("No alignment loaded.").size(12).color(palette::TEXT_DIM))
                .padding([12, 10]),
        ]
        .into();
    }

    let show_cols = ncols.min(MAX_COLS);

    // Column header
    let mut head_cells: Vec<Element<'a, Message>> = Vec::new();
    head_cells.push(
        container(text("Residue").size(10).color(palette::TEXT_DIM).font(Font::MONOSPACE))
            .width(60)
            .padding([2, 4])
            .style(header_style())
            .into()
    );
    for c in 0..show_cols {
        head_cells.push(
            container(text(format!("{}", c + 1)).size(9).color(palette::TEXT_DIM).font(Font::MONOSPACE))
                .width(32)
                .padding([2, 2])
                .style(header_style())
                .into()
        );
    }
    let col_header: Element<'a, Message> = row(head_cells).into();

    // Collect residue set across visible columns
    let mut all_residues: Vec<u8> = Vec::new();
    for c in 0..show_cols {
        for r in 0..nrows {
            let b = aln.sequences[r].residues.get(c).copied().unwrap_or(b'-');
            let b = b.to_ascii_uppercase();
            if !matches!(b, b'-' | b'~' | b'.') && !all_residues.contains(&b) {
                all_residues.push(b);
            }
        }
    }
    all_residues.sort();

    // Build one row per residue
    let mut table_rows: Vec<Element<'a, Message>> = vec![col_header];

    for &res in &all_residues {
        let mut cells: Vec<Element<'a, Message>> = Vec::new();
        cells.push(
            container(
                text(format!("{}", res as char))
                    .size(10)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
            )
            .width(60)
            .padding([2, 4])
            .style(header_style())
            .into()
        );

        for c in 0..show_cols {
            let count = (0..nrows).filter(|&r| {
                aln.sequences[r].residues.get(c)
                    .map(|&b| b.to_ascii_uppercase() == res)
                    .unwrap_or(false)
            }).count();

            let pct = if nrows > 0 { count as f32 / nrows as f32 } else { 0.0 };
            let intensity = pct.clamp(0.0, 1.0);
            let bg = Color {
                r: 1.0 - intensity * 0.7,
                g: 1.0 - intensity * 0.3,
                b: 1.0 - intensity * 0.8,
                a: 1.0,
            };

            let label = if count > 0 { format!("{}", count) } else { String::new() };
            cells.push(
                container(
                    text(label)
                        .size(9)
                        .color(if pct > 0.5 { Color::WHITE } else { Color::BLACK })
                        .font(Font::MONOSPACE),
                )
                .width(32)
                .padding([2, 2])
                .style(move |_| container::Style {
                    background: Some(Background::Color(bg)),
                    ..Default::default()
                })
                .into()
            );
        }

        table_rows.push(row(cells).into());
    }

    // Gap row
    let mut gap_cells: Vec<Element<'a, Message>> = Vec::new();
    gap_cells.push(
        container(text("gap").size(10).color(palette::TEXT_DIM).font(Font::MONOSPACE))
            .width(60).padding([2, 4]).style(header_style()).into()
    );
    for c in 0..show_cols {
        let count = (0..nrows).filter(|&r| {
            aln.sequences[r].residues.get(c)
                .map(|&b| matches!(b, b'-'|b'~'|b'.'))
                .unwrap_or(true)
        }).count();
        let label = if count > 0 { format!("{}", count) } else { String::new() };
        gap_cells.push(
            container(text(label).size(9).color(palette::TEXT_DIM).font(Font::MONOSPACE))
                .width(32).padding([2, 2]).into()
        );
    }
    table_rows.push(row(gap_cells).into());

    let note: Element<'a, Message> = if ncols > MAX_COLS {
        container(
            text(format!("(showing first {MAX_COLS} of {ncols} columns)"))
                .size(11).color(palette::TEXT_DIM),
        )
        .padding([4, 10])
        .into()
    } else {
        iced::widget::Space::with_height(0).into()
    };

    let body = column(table_rows).spacing(0).width(Length::Shrink);

    column![
        header,
        scrollable(
            column![note, body].spacing(4).padding([4, 8]).width(Length::Fill)
        )
        .width(Length::Fill)
        .height(Length::Fill),
    ]
    .width(Length::Fill)
    .height(Length::Fill)
    .into()
}

fn header_style() -> impl Fn(&iced::Theme) -> container::Style {
    move |_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        ..Default::default()
    }
}
