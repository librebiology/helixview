//! Identity matrix view — NxN pairwise identity table.

use std::sync::Arc;

use crate::app::Message;
use crate::theme::palette;
use helixview_core::Alignment;
use iced::{
    widget::{button, column, container, row, scrollable, text},
    Background, Border, Color, Element, Font, Length,
};

/// Display a scrollable NxN pairwise identity matrix.
/// `precomputed` – if None, shows a "Computing…" message while the async task runs.
pub fn identity_matrix_view<'a>(
    aln: &'a Arc<Alignment>,
    precomputed: Option<&'a Vec<Vec<f64>>>,
) -> Element<'a, Message> {
    let n = aln.seq_count();

    // Header bar
    let back_btn = button(text("< Back").size(13))
        .padding([4, 12])
        .on_press(Message::CloseIdentityMatrix);

    let csv_btn = button(text("Export CSV…").size(13))
        .padding([4, 10])
        .on_press(Message::ExportIdentityCsv);

    let header = container(
        row![
            back_btn,
            iced::widget::horizontal_space(),
            text("Pairwise Identity Matrix")
                .size(14)
                .color(palette::TEXT)
                .font(Font::MONOSPACE),
            iced::widget::horizontal_space(),
            csv_btn,
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

    if n == 0 {
        return column![
            header,
            container(
                text("No sequences loaded.")
                    .size(12)
                    .color(palette::TEXT_DIM)
            )
            .padding([16, 10])
        ]
        .into();
    }

    // Show spinner while async computation runs.
    let matrix = match precomputed {
        None => {
            return column![
                header,
                container(
                    text("Computing identity matrix…")
                        .size(12)
                        .color(palette::TEXT_DIM)
                )
                .padding([20, 16])
            ]
            .into();
        }
        Some(m) => m,
    };

    const MAX_N: usize = 200;
    let show_n = n.min(MAX_N).min(matrix.len());

    // Build the table.
    // Column header row: short seq names
    const COL_W: f32 = 72.0;
    const NAME_W: f32 = 160.0;

    let mut header_cells: Vec<Element<'a, Message>> = Vec::with_capacity(show_n + 1);
    header_cells.push(
        container(text("").size(10))
            .width(NAME_W)
            .padding([2, 4])
            .into(),
    );
    for j in 0..show_n {
        let label = truncate(&aln.sequences[j].name, 8);
        header_cells.push(
            container(
                text(label)
                    .size(12)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
            )
            .width(COL_W)
            .padding([3, 4])
            .style(|_| container::Style {
                background: Some(Background::Color(palette::HEADER_BG)),
                border: Border {
                    color: palette::BORDER,
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .into(),
        );
    }
    let col_header_row: Element<'a, Message> = row(header_cells).into();

    // Data rows
    let mut table_rows: Vec<Element<'a, Message>> = Vec::with_capacity(show_n + 1);
    table_rows.push(col_header_row);

    for i in 0..show_n {
        let mut cells: Vec<Element<'a, Message>> = Vec::with_capacity(show_n + 1);

        // Row name
        cells.push(
            container(
                text(truncate(&aln.sequences[i].name, 22))
                    .size(12)
                    .color(palette::TEXT)
                    .font(Font::MONOSPACE),
            )
            .width(NAME_W)
            .padding([4, 6])
            .style(|_| container::Style {
                background: Some(Background::Color(palette::HEADER_BG)),
                border: Border {
                    color: palette::BORDER,
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .into(),
        );

        for j in 0..show_n {
            let val = matrix[i][j];
            let pct = val * 100.0;
            let bg = identity_color(val as f32);
            let label = if i == j {
                "—".to_string()
            } else {
                format!("{:.0}%", pct)
            };
            // White text on dark colors, dark text on light/medium colors.
            let txt = if val >= 0.80 || val < 0.25 {
                Color::WHITE
            } else {
                Color::from_rgb(0.08, 0.08, 0.12)
            };
            cells.push(
                container(text(label).size(12).color(txt).font(Font::MONOSPACE))
                    .width(COL_W)
                    .padding([4, 4])
                    .style(move |_| container::Style {
                        background: Some(Background::Color(bg)),
                        ..Default::default()
                    })
                    .into(),
            );
        }

        table_rows.push(row(cells).into());
    }

    let note: Element<'a, Message> = if n > MAX_N {
        container(
            text(format!("(showing first {MAX_N} of {n} sequences)"))
                .size(11)
                .color(palette::TEXT_DIM),
        )
        .padding([4, 10])
        .into()
    } else {
        iced::widget::Space::with_height(0).into()
    };

    let body = column(table_rows).spacing(0).width(Length::Fill);

    let content = column![
        header,
        scrollable(
            column![note, body]
                .spacing(4)
                .padding([4, 8])
                .width(Length::Fill)
        )
        .width(Length::Fill)
        .height(Length::Fill),
    ]
    .width(Length::Fill)
    .height(Length::Fill);

    content.into()
}

/// Heat-map color: green (high identity) → yellow → red (low).
fn identity_color(pct: f32) -> Color {
    if pct >= 0.95 {
        return Color::from_rgb(0.08, 0.55, 0.20);
    }
    if pct >= 0.80 {
        return Color::from_rgb(0.22, 0.70, 0.35);
    }
    if pct >= 0.60 {
        return Color::from_rgb(0.55, 0.80, 0.30);
    }
    if pct >= 0.40 {
        return Color::from_rgb(0.90, 0.75, 0.20);
    }
    if pct >= 0.25 {
        return Color::from_rgb(0.90, 0.48, 0.15);
    }
    Color::from_rgb(0.78, 0.20, 0.15)
}

fn truncate(s: &str, max: usize) -> String {
    if s.chars().count() <= max {
        return s.to_owned();
    }
    let t: String = s.chars().take(max - 1).collect();
    format!("{t}…")
}
