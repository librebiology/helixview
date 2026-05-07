use iced::{
    widget::{button, column, container, row, text, vertical_space, horizontal_space},
    Alignment, Background, Border, Color, Element, Length,
};
use crate::app::Message;
use crate::theme::{palette, buttons as btn_style};

pub fn welcome_view() -> Element<'static, Message> {
    // ── Hero ─────────────────────────────────────────────────────────────────
    let title = text("HelixView")
        .size(54)
        .color(palette::ACCENT);

    let subtitle = text("Biological Sequence Analysis")
        .size(18)
        .color(palette::TEXT);

    let dedication = text("Named in honor of Rosalind Franklin (1920–1958)")
        .size(12)
        .color(palette::TEXT_DIM);

    let open_btn = button(text("  Open File…  ").size(15))
        .padding([10, 28])
        .style(btn_style::primary)
        .on_press(Message::OpenFileDialog);

    let new_btn = button(text("  New Alignment  ").size(15))
        .padding([10, 28])
        .style(btn_style::secondary)
        .on_press(Message::NewAlignment);

    let fetch_btn = button(text("  Fetch from NCBI…  ").size(15))
        .padding([10, 28])
        .style(btn_style::secondary)
        .on_press(Message::ToggleFetchBar);

    let cta = row![open_btn, new_btn, fetch_btn]
        .spacing(12)
        .align_y(Alignment::Center);

    // ── Feature cards ────────────────────────────────────────────────────────
    let card = |title: &'static str, body: &'static str| -> Element<'static, Message> {
        container(
            column![
                text(title).size(13).color(palette::ACCENT),
                text(body).size(11).color(palette::TEXT_DIM),
            ]
            .spacing(4)
            .padding([10, 12]),
        )
        .style(|_| container::Style {
            background: Some(Background::Color(Color::WHITE)),
            border: Border {
                color: palette::BORDER,
                width: 1.0,
                radius: 6.0.into(),
            },
            ..Default::default()
        })
        .width(Length::Fixed(200.0))
        .into()
    };

    let cards = row![
        card("Alignment Editor",   "FASTA · GenBank · ClustalW\nPhylip · NEXUS · EMBL · PIR"),
        card("Sequence Analysis",  "Mutual info · Entropy\nORFs · 6-Frame · Identity"),
        card("Visualization",      "Plasmid map · Phylo tree\nABI trace · Dot plot"),
        card("External Tools",     "BLAST · ClustalW · MUSCLE\nMAFFT · ClustalΩ"),
    ]
    .spacing(12)
    .align_y(Alignment::Start);

    // ── Shortcuts reference ──────────────────────────────────────────────────
    let shortcut_row = |keys: &'static str, desc: &'static str| -> Element<'static, Message> {
        row![
            container(text(keys).size(11).color(palette::ACCENT))
                .width(Length::Fixed(130.0)),
            text(desc).size(11).color(palette::TEXT_DIM),
        ]
        .spacing(8)
        .align_y(Alignment::Center)
        .into()
    };

    let shortcuts = container(
        column![
            text("Keyboard shortcuts").size(13).color(palette::TEXT),
            vertical_space().height(6),
            shortcut_row("Ctrl+O / Ctrl+S",  "Open / Save"),
            shortcut_row("Ctrl+Z / Ctrl+Y",  "Undo / Redo"),
            shortcut_row("Ctrl+F",            "Find sequence pattern"),
            shortcut_row("Ctrl+P",            "Command palette"),
            shortcut_row("Click + type",      "Edit residue in cell"),
            shortcut_row("Right-click column","Insert / delete gap column"),
        ]
        .spacing(4)
        .padding([12, 14]),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(Color::WHITE)),
        border: Border { color: palette::BORDER, width: 1.0, radius: 6.0.into() },
        ..Default::default()
    })
    .width(Length::Fixed(320.0));

    // ── Layout ───────────────────────────────────────────────────────────────
    container(
        column![
            vertical_space(),
            title,
            vertical_space().height(2),
            subtitle,
            dedication,
            vertical_space().height(28),
            cta,
            vertical_space().height(32),
            cards,
            vertical_space().height(24),
            row![horizontal_space(), shortcuts, horizontal_space()]
                .width(Length::Fill),
            vertical_space(),
        ]
        .spacing(6)
        .align_x(Alignment::Center)
        .width(Length::Fill),
    )
    .width(Length::Fill)
    .height(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgb(0.96, 0.97, 0.99))),
        ..Default::default()
    })
    .center(Length::Fill)
    .into()
}
