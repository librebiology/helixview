//! Global preferences dialog.

use crate::app::Message;
use crate::prefs::Preferences;
use crate::theme::palette;
use iced::{
    widget::{
        button, checkbox, column, container, horizontal_space, row, scrollable, slider, text,
    },
    Alignment, Element, Length,
};

pub fn prefs_dialog<'a>(prefs: &Preferences) -> Element<'a, Message> {
    let title = text("Preferences").size(16);

    // ── Section helper ────────────────────────────────────────────────────────
    let section = |label: &'static str| -> Element<'a, Message> {
        text(label).size(12).color(palette::TEXT_DIM).into()
    };

    // ── Editing ───────────────────────────────────────────────────────────────
    let undo_row = row![
        text("Undo depth:").size(12).width(Length::Fixed(160.0)),
        slider(10..=500, prefs.undo_depth as u32, |v| {
            Message::PrefsUndoDepth(v as usize)
        })
        .width(Length::Fixed(180.0)),
        text(format!("{}", prefs.undo_depth))
            .size(12)
            .color(palette::TEXT_DIM),
    ]
    .spacing(8)
    .align_y(Alignment::Center);

    // ── View defaults ─────────────────────────────────────────────────────────
    let zoom_row = row![
        text("Default zoom:").size(12).width(Length::Fixed(160.0)),
        slider(40..=300, (prefs.default_zoom * 100.0) as u32, |v| {
            Message::PrefsDefaultZoom(v as f32 / 100.0)
        })
        .width(Length::Fixed(180.0)),
        text(format!("{:.0}%", prefs.default_zoom * 100.0))
            .size(12)
            .color(palette::TEXT_DIM),
    ]
    .spacing(8)
    .align_y(Alignment::Center);

    let font_row = row![
        text("Grid font size:").size(12).width(Length::Fixed(160.0)),
        slider(8..=20, prefs.font_size as u32, |v| Message::PrefsFontSize(
            v as f32
        ))
        .width(Length::Fixed(180.0)),
        text(format!("{}pt", prefs.font_size as u32))
            .size(12)
            .color(palette::TEXT_DIM),
    ]
    .spacing(8)
    .align_y(Alignment::Center);

    let features_cb = checkbox("Show feature stripes by default", prefs.show_features)
        .on_toggle(Message::PrefsShowFeatures)
        .size(14)
        .text_size(12);

    let restr_cb = checkbox("Show restriction map by default", prefs.show_restr_map)
        .on_toggle(Message::PrefsShowRestrMap)
        .size(14)
        .text_size(12);

    // ── Identity threshold ────────────────────────────────────────────────────
    let thresh_row = row![
        text("Identity threshold:")
            .size(12)
            .width(Length::Fixed(160.0)),
        slider(
            0.0..=1.0,
            prefs.identity_threshold,
            Message::PrefsIdentityThreshold
        )
        .step(0.05)
        .width(Length::Fixed(180.0)),
        text(format!("{:.0}%", prefs.identity_threshold * 100.0))
            .size(12)
            .color(palette::TEXT_DIM),
    ]
    .spacing(8)
    .align_y(Alignment::Center);

    // ── Recent files ──────────────────────────────────────────────────────────
    let recent_row = row![
        text("Recent files to remember:")
            .size(12)
            .width(Length::Fixed(160.0)),
        slider(0..=20, prefs.recent_files_max as u32, |v| {
            Message::PrefsRecentFilesMax(v as usize)
        })
        .width(Length::Fixed(180.0)),
        text(format!("{}", prefs.recent_files_max))
            .size(12)
            .color(palette::TEXT_DIM),
    ]
    .spacing(8)
    .align_y(Alignment::Center);

    // ── Recent files list ─────────────────────────────────────────────────────
    let mut recent_col = column![].spacing(2);
    if prefs.recent_files.is_empty() {
        recent_col = recent_col.push(text("No recent files.").size(11).color(palette::TEXT_DIM));
    } else {
        for (i, path) in prefs.recent_files.iter().enumerate() {
            let label = path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("?")
                .to_string();
            let full = path.display().to_string();
            let open_btn = button(text(format!("{}. {}", i + 1, label)).size(11))
                .padding([2, 8])
                .style(button::secondary)
                .on_press(Message::OpenRecentFile(path.clone()));
            recent_col = recent_col.push(
                row![open_btn, text(full).size(10).color(palette::TEXT_DIM)]
                    .spacing(6)
                    .align_y(Alignment::Center),
            );
        }
    }
    let clear_recent = button(text("Clear Recent Files").size(11))
        .padding([3, 8])
        .style(button::secondary)
        .on_press(Message::PrefsClearRecent);

    // ── Buttons ───────────────────────────────────────────────────────────────
    let close_btn = button(text("Close").size(12))
        .padding([5, 16])
        .style(button::secondary)
        .on_press(Message::ClosePrefs);
    let save_btn = button(text("Save & Close").size(12))
        .padding([5, 16])
        .style(button::primary)
        .on_press(Message::SavePrefs);
    let btn_row = row![close_btn, horizontal_space(), save_btn].align_y(Alignment::Center);

    let body = column![
        title,
        section("— Editing —"),
        undo_row,
        section("— View Defaults —"),
        zoom_row,
        font_row,
        features_cb,
        restr_cb,
        thresh_row,
        section("— Recent Files —"),
        recent_row,
        scrollable(recent_col).height(Length::Fixed(120.0)),
        clear_recent,
        btn_row,
    ]
    .spacing(12)
    .padding(24)
    .width(Length::Fixed(520.0));

    container(body)
        .width(Length::Fill)
        .height(Length::Fill)
        .center_x(Length::Fill)
        .center_y(Length::Fill)
        .into()
}
