//! Keyboard shortcut reference panel.

use crate::app::Message;
use crate::theme::palette;
use iced::{
    widget::{button, column, container, row, scrollable, text},
    Alignment, Color, Element, Length,
};

const SHORTCUTS: &[(&str, &str, &str)] = &[
    // Category, Key, Action
    ("File", "Ctrl+O", "Open file"),
    ("File", "Ctrl+S", "Save file"),
    ("File", "Ctrl+Shift+S", "Save project"),
    ("File", "Ctrl+Shift+O", "Open project"),
    ("File", "Ctrl+T", "New tab"),
    ("File", "Ctrl+W", "Close tab"),
    ("Edit", "Ctrl+Z", "Undo"),
    ("Edit", "Ctrl+Y", "Redo"),
    ("Edit", "Space", "Insert gap at cursor"),
    ("Edit", "Shift+Space", "Insert gap in all other sequences"),
    ("Edit", "Backspace", "Delete gap before cursor"),
    (
        "Edit",
        "Shift+Backspace",
        "Delete gap before cursor in all other sequences",
    ),
    (
        "Navigate",
        "Arrow keys",
        "Scroll alignment / move edit cursor",
    ),
    ("Navigate", "Ctrl+End", "Scroll to last column"),
    ("Navigate", "Ctrl+F", "Toggle search bar"),
    ("Navigate", "F3 / Enter", "Jump to next search match"),
    ("Navigate", "Shift+F3", "Jump to previous search match"),
    ("Navigate", "Escape", "Close search / cancel selection"),
    ("Selection", "Click", "Select sequence row"),
    ("Selection", "Shift+Click", "Range select rows"),
    ("Selection", "Ctrl+Click", "Toggle row in selection"),
    ("Selection", "Ctrl+A", "Select all sequences"),
    ("Selection", "Ctrl+D", "Deselect all"),
    ("Selection", "Click ruler", "Select column range"),
    ("Zoom", "Ctrl++", "Zoom in"),
    ("Zoom", "Ctrl+−", "Zoom out"),
    ("Zoom", "Ctrl+0", "Reset zoom"),
    (
        "Analysis",
        "Toolbar → Analysis",
        "Composition, Tm, ORFs, BLAST, Tree…",
    ),
    (
        "View",
        "Toolbar → View",
        "Color scheme, Export SVG, Shaded Fig…",
    ),
];

pub fn help_view<'a>() -> Element<'a, Message> {
    let close_btn = button(text("< Back").size(12))
        .padding([4, 10])
        .on_press(Message::CloseHelp);
    let title = text("Keyboard Shortcuts").size(15);
    let toolbar = row![close_btn, title]
        .spacing(12)
        .align_y(Alignment::Center)
        .padding([4, 8]);

    // Header
    let hdr = row![
        cell_hdr("Category", 110.0),
        cell_hdr("Shortcut", 190.0),
        cell_hdr("Action", 340.0),
    ]
    .spacing(1);

    let mut rows_col = column![hdr].spacing(0);
    let mut last_cat = "";
    for (i, (cat, key, action)) in SHORTCUTS.iter().enumerate() {
        let shade = i % 2 == 1;
        let cat_cell = if *cat != last_cat {
            last_cat = cat;
            cell_data(cat, 110.0, shade, true)
        } else {
            cell_data("", 110.0, shade, false)
        };
        rows_col = rows_col.push(
            row![
                cat_cell,
                cell_data(key, 190.0, shade, false),
                cell_data(action, 340.0, shade, false),
            ]
            .spacing(1),
        );
    }

    let content: Element<'a, Message> = scrollable(column![rows_col].padding(8))
        .height(Length::Fill)
        .into();

    column![toolbar, content]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

fn cell_hdr<'a>(label: &str, w: f32) -> Element<'a, Message> {
    container(text(label.to_string()).size(11))
        .width(Length::Fixed(w))
        .padding([3, 6])
        .style(|_| container::Style {
            background: Some(iced::Background::Color(Color::from_rgb(0.88, 0.89, 0.92))),
            ..Default::default()
        })
        .into()
}

fn cell_data<'a>(label: &str, w: f32, shade: bool, bold: bool) -> Element<'a, Message> {
    let bg = if shade {
        Color::from_rgb(0.97, 0.97, 0.99)
    } else {
        Color::WHITE
    };
    let t = if bold {
        text(label.to_string()).size(11).color(palette::TEXT)
    } else {
        text(label.to_string()).size(11)
    };
    container(t)
        .width(Length::Fixed(w))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}
