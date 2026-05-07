//! Raw sequence text editor — opened by double-clicking a sequence title.

use iced::{
    widget::{button, column, container, row, text, text_editor},
    Alignment, Background, Border, Element, Font, Length,
};
use helixview_core::sequence::is_gap;

use crate::app::Message;
use crate::theme::palette;

/// Full-window raw sequence editor panel.
///
/// `content` is the live `text_editor::Content` stored in `HelixViewApp`.
pub fn seq_editor_view<'a>(
    seq_name:  &'a str,
    content:   &'a text_editor::Content,
    seq_idx:   usize,
) -> Element<'a, Message> {

    // Stats derived from the current editor text
    let raw_text = content.text();
    let residue_count = raw_text.bytes()
        .filter(|&b| !b.is_ascii_whitespace() && !is_gap(b))
        .count();
    let gap_count = raw_text.bytes()
        .filter(|&b| is_gap(b))
        .count();
    let total_chars = raw_text.bytes()
        .filter(|&b| !b.is_ascii_whitespace())
        .count();

    // ── Toolbar ───────────────────────────────────────────────────────────────
    let close_btn = button(text("✕ Cancel").size(13))
        .padding([4, 10])
        .on_press(Message::SeqEditorClose);

    let commit_btn = button(text("✓ Commit").size(13))
        .padding([4, 10])
        .on_press(Message::SeqEditorCommit);

    let title_lbl = text(format!("Edit sequence: {seq_name}"))
        .size(14)
        .color(palette::TEXT);

    let stats_lbl = text(format!(
        "{residue_count} residues  ·  {gap_count} gaps  ·  {total_chars} total chars"
    ))
    .size(12)
    .color(palette::TEXT_DIM);

    let hint_lbl = text("Paste or type raw sequence (gaps kept). Whitespace and newlines are stripped on Commit.")
        .size(11)
        .color(palette::TEXT_DIM);

    let toolbar = container(
        row![
            close_btn,
            title_lbl,
            iced::widget::horizontal_space(),
            stats_lbl,
            commit_btn,
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

    // ── Text editor ───────────────────────────────────────────────────────────
    let editor = text_editor(content)
        .on_action(Message::SeqEditorAction)
        .font(Font::MONOSPACE)
        .size(13)
        .height(Length::Fill)
        .padding([8, 10]);

    // ── Status bar ────────────────────────────────────────────────────────────
    let line_col = {
        let cursor = content.cursor_position();
        format!("Line {}, Col {}", cursor.0 + 1, cursor.1 + 1)
    };
    let status_bar = container(
        row![
            text(line_col).size(11).color(palette::TEXT_DIM),
            iced::widget::horizontal_space(),
            hint_lbl,
        ]
        .spacing(10)
        .align_y(Alignment::Center)
        .padding([3, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::BG_PANEL)),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .width(Length::Fill);

    let _ = seq_idx; // may be used for future validation

    column![toolbar, editor, status_bar]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
