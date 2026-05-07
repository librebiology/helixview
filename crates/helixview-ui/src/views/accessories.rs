//! Accessory app manager — configure and run external CLI tools.

use iced::{
    widget::{button, column, container, row, scrollable, text, text_input},
    Alignment, Background, Border, Color, Element, Length,
};
use crate::app::{AccessoryApp, Message};
use crate::theme::palette;

pub fn accessories_view<'a>(
    apps:             &'a [AccessoryApp],
    editing:          Option<usize>,
    edit_name:        &'a str,
    edit_cmd:         &'a str,
    has_document:     bool,
) -> Element<'a, Message> {
    let back_btn = button(text("< Back").size(12))
        .padding([4, 10])
        .on_press(Message::CloseAccessories);

    let toolbar = container(
        row![
            back_btn,
            text("Accessory Apps").size(14).color(palette::TEXT),
            iced::widget::horizontal_space(),
            button(text("+ New").size(12)).padding([4, 10])
                .on_press(Message::NewAccessory),
        ]
        .spacing(10)
        .align_y(Alignment::Center)
        .padding([5, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        ..Default::default()
    })
    .width(Length::Fill);

    // ── Installed apps list ───────────────────────────────────────────────────
    let mut list = column![].spacing(4).padding(iced::Padding { top: 8.0, right: 8.0, bottom: 0.0, left: 8.0 });

    if apps.is_empty() {
        list = list.push(
            text("No accessory apps configured. Click \"+ New\" to add one.")
                .size(12)
                .color(palette::TEXT_DIM),
        );
    }

    for (i, app_cfg) in apps.iter().enumerate() {
        let name_label = text(app_cfg.name.clone()).size(13).color(palette::TEXT);
        let cmd_label  = text(
            if app_cfg.command.len() > 60 {
                format!("{}…", &app_cfg.command[..59])
            } else {
                app_cfg.command.clone()
            }
        )
        .size(11)
        .color(palette::TEXT_DIM)
        .font(iced::Font::MONOSPACE);

        let run_btn = if has_document {
            button(text("Run").size(11)).padding([3, 8])
                .on_press(Message::RunAccessory(i))
        } else {
            button(text("Run").size(11)).padding([3, 8])
        };

        let edit_btn = button(text("Edit").size(11)).padding([3, 8])
            .style(button::secondary)
            .on_press(Message::EditAccessory(i));

        let del_btn = button(text("Delete").size(11)).padding([3, 8])
            .style(button::secondary)
            .on_press(Message::DeleteAccessory(i));

        list = list.push(
            container(
                row![
                    column![name_label, cmd_label].spacing(2).width(Length::Fill),
                    run_btn,
                    edit_btn,
                    del_btn,
                ]
                .spacing(6)
                .align_y(Alignment::Center),
            )
            .padding([6, 10])
            .width(Length::Fill)
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(0.97, 0.97, 0.99))),
                border: Border {
                    color: Color::from_rgb(0.82, 0.84, 0.90),
                    width: 1.0,
                    radius: 4.0.into(),
                },
                ..Default::default()
            }),
        );
    }

    // ── Edit panel ────────────────────────────────────────────────────────────
    let edit_panel: Element<Message> = if editing.is_some() || edit_name != "" || edit_cmd != "" {
        let title_txt = if editing.is_some() { "Edit Accessory" } else { "New Accessory" };
        let save_btn = button(text("Save").size(12)).padding([4, 12])
            .on_press(Message::SaveAccessory);
        let cancel_btn = button(text("Cancel").size(12)).padding([4, 10])
            .style(button::secondary)
            .on_press(Message::CloseAccessories);

        let name_input = text_input("Display name…", edit_name)
            .on_input(Message::AccessoryFieldName)
            .size(13)
            .padding([6, 10])
            .width(Length::Fill);

        let cmd_input = text_input("Command (use {input} for FASTA path)…", edit_cmd)
            .on_input(Message::AccessoryFieldCmd)
            .size(12)
            .padding([6, 10])
            .width(Length::Fill)
            .font(iced::Font::MONOSPACE);

        container(
            column![
                text(title_txt).size(13).color(palette::TEXT),
                text("Name:").size(11).color(palette::TEXT_DIM),
                name_input,
                text("Command template:").size(11).color(palette::TEXT_DIM),
                cmd_input,
                text("{input} = temp FASTA of selected (or all) sequences")
                    .size(10).color(palette::TEXT_DIM),
                row![save_btn, cancel_btn].spacing(8),
            ]
            .spacing(6),
        )
        .padding([10, 12])
        .width(Length::Fill)
        .style(|_| container::Style {
            background: Some(Background::Color(Color::from_rgb(0.96, 0.97, 1.0))),
            border: Border {
                color: Color::from_rgb(0.70, 0.75, 0.90),
                width: 1.0,
                radius: 4.0.into(),
            },
            ..Default::default()
        })
        .into()
    } else {
        iced::widget::horizontal_space().into()
    };

    let content = scrollable(
        column![list, container(edit_panel).padding([8, 8])].spacing(0),
    )
    .width(Length::Fill)
    .height(Length::Fill);

    column![toolbar, content]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
