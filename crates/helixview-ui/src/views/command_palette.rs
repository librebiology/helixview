//! Command palette overlay — Ctrl+P fuzzy search over all app commands.

use iced::{
    widget::{button, column, container, row, text, text_input},
    Background, Border, Color, Element, Length,
};
use crate::app::Message;

// All palette entries: (display label, message to fire on selection).
// Only simple, parameter-free messages are listed here.
const COMMANDS: &[(&str, fn() -> Message)] = &[
    // File
    ("Open file",                   || Message::OpenFileDialog),
    ("Save file",                   || Message::SaveFileDialog),
    ("New alignment",               || Message::NewAlignment),
    ("New tab",                     || Message::NewTab),
    ("Save project",                || Message::SaveProjectDialog),
    ("Open project",                || Message::OpenProjectDialog),
    // Edit
    ("Undo",                        || Message::Undo),
    ("Redo",                        || Message::Redo),
    ("Select all",                  || Message::SelectAll),
    ("Deselect all",                || Message::SelectNone),
    // Analysis
    ("Pairwise alignment",          || Message::RunPairwise),
    ("ORF finder",                  || Message::RunOrfFinder),
    ("Six-frame translation",       || Message::RunSixFrame),
    ("Composition analysis",        || Message::ShowComposition),
    ("Hydrophobicity profile",      || Message::ShowHydrophobicity),
    ("Oligo Tm calculator",         || Message::ShowOligoTm),
    ("Identity matrix",             || Message::ShowIdentityMatrix),
    ("Dot plot",                    || Message::ShowDotPlot),
    ("Conserved regions",           || Message::ShowConservation),
    ("Column statistics",           || Message::ShowColSummary),
    ("Taxonomy table",              || Message::ShowTaxonomy),
    ("Plasmid map",                 || Message::ShowPlasmid),
    ("BLAST",                       || Message::ToggleBlastPanel),
    // View
    ("Toggle analysis panel",       || Message::ToggleAnalysis),
    ("Toggle features",             || Message::ToggleFeatures),
    ("Color editor",                || Message::OpenColorEditor),
    ("Zoom in",                     || Message::ZoomIn),
    ("Zoom out",                    || Message::ZoomOut),
    ("Reset zoom",                  || Message::ZoomReset),
    ("Toggle search bar",           || Message::ToggleSearch),
    // Export
    ("Export SVG",                  || Message::ExportSvgDialog),
    ("Export PDF",                  || Message::ExportPdfDialog),
    ("Export shaded figure (SVG)",  || Message::ShowShadedExportDialog),
    ("Export text alignment",       || Message::ShowTextExport),
    // External tools
    ("Accessory apps",              || Message::ShowAccessories),
    ("Create BLAST database",       || Message::CreateBlastDbDialog),
    // Settings
    ("Preferences",                 || Message::OpenPrefs),
    ("Keyboard shortcuts",          || Message::ShowHelp),
];

pub fn command_palette<'a>(query: &str, selected: usize) -> Element<'a, Message> {
    // Filter commands by query (case-insensitive substring match).
    let q = query.to_lowercase();
    let matches: Vec<(usize, &(&str, fn() -> Message))> = COMMANDS.iter()
        .enumerate()
        .filter(|(_, (label, _))| q.is_empty() || label.to_lowercase().contains(&q))
        .collect();

    // Text input.
    let input: Element<Message> = text_input("Search commands…", query)
        .on_input(Message::PaletteQuery)
        .on_submit(Message::PaletteConfirm)
        .size(14)
        .padding([8, 12])
        .width(Length::Fill)
        .id(text_input::Id::new("palette_input"))
        .into();

    // Result list.
    let mut list = column![].spacing(0);
    for (i, &(orig_idx, (label, _msg_fn))) in matches.iter().enumerate() {
        let is_sel = i == selected;
        let bg = if is_sel {
            Color::from_rgb(0.20, 0.45, 0.80)
        } else if i % 2 == 0 {
            Color::from_rgb(0.97, 0.97, 0.99)
        } else {
            Color::WHITE
        };
        let txt_color = if is_sel { Color::WHITE } else { Color::from_rgb(0.15, 0.15, 0.20) };
        let orig = orig_idx;
        list = list.push(
            button(
                container(
                    text(label.to_string()).size(13).color(txt_color),
                )
                .padding([6, 14])
                .width(Length::Fill)
                .style(move |_| container::Style {
                    background: Some(Background::Color(bg)),
                    ..Default::default()
                }),
            )
            .width(Length::Fill)
            .style(|_, _| button::Style {
                background: None,
                border: Border::default(),
                ..Default::default()
            })
            .on_press(Message::PaletteSelect(orig))
            .padding(0),
        );
    }

    if matches.is_empty() {
        let empty: Element<Message> = container(
            text("No commands match").size(12).color(Color::from_rgb(0.5, 0.5, 0.5))
        )
        .padding([10, 14])
        .into();
        list = list.push(empty);
    }

    let hint = text("↑↓ navigate · Enter confirm · Esc close")
        .size(10)
        .color(Color::from_rgb(0.55, 0.55, 0.60));

    let inner = container(
        column![
            input,
            container(iced::widget::horizontal_space())
                .width(Length::Fill)
                .height(Length::Fixed(1.0))
                .style(|_| container::Style {
                    background: Some(Background::Color(Color::from_rgb(0.85, 0.86, 0.90))),
                    ..Default::default()
                }),
            iced::widget::scrollable(list)
                .height(Length::Fixed(320.0))
                .width(Length::Fill),
            container(hint)
                .padding([4, 14])
                .width(Length::Fill),
        ]
        .spacing(0),
    )
    .width(Length::Fixed(500.0))
    .style(|_| container::Style {
        background: Some(Background::Color(Color::WHITE)),
        border: Border {
            color: Color::from_rgb(0.55, 0.58, 0.70),
            width: 1.5,
            radius: 6.0.into(),
        },
        shadow: iced::Shadow {
            color: Color::from_rgba(0.0, 0.0, 0.0, 0.25),
            offset: iced::Vector::new(0.0, 4.0),
            blur_radius: 12.0,
        },
        ..Default::default()
    });

    // Dim overlay + centered inner.
    container(
        container(inner)
            .align_x(iced::alignment::Horizontal::Center)
            .width(Length::Fill)
            .padding(iced::Padding { top: 60.0, right: 0.0, bottom: 0.0, left: 0.0 }),
    )
    .width(Length::Fill)
    .height(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgba(0.0, 0.0, 0.0, 0.40))),
        ..Default::default()
    })
    .into()
}

/// Fire the command at index `orig_idx` in the COMMANDS list.
pub fn fire(orig_idx: usize) -> Option<Message> {
    COMMANDS.get(orig_idx).map(|(_, f)| f())
}

/// Return the original COMMANDS index of the `visible_idx`-th item matching `query`.
pub fn resolve_selected(query: &str, visible_idx: usize) -> Option<usize> {
    let q = query.to_lowercase();
    COMMANDS.iter()
        .enumerate()
        .filter(|(_, (label, _))| q.is_empty() || label.to_lowercase().contains(&q))
        .nth(visible_idx)
        .map(|(i, _)| i)
}

/// Count how many commands match `query`.
pub fn match_count(query: &str) -> usize {
    let q = query.to_lowercase();
    COMMANDS.iter()
        .filter(|(label, _)| q.is_empty() || label.to_lowercase().contains(&q))
        .count()
}
