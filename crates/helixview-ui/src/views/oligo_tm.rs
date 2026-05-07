//! Oligonucleotide Tm calculator view.

use helixview_analysis::{calculate_tm, wallace_tm};
use helixview_core::Alignment as SeqAlignment;
use helixview_core::SequenceType;
use iced::{
    widget::{button, column, container, row, scrollable, text},
    Alignment, Color, Element, Length,
};

use crate::app::Message;

pub fn oligo_tm_view<'a>(aln: &SeqAlignment) -> Element<'a, Message> {
    let toolbar = row![
        button(text("< Back").size(12))
            .padding([3, 10])
            .style(button::secondary)
            .on_press(Message::CloseOligoTm),
        text("Oligonucleotide Tm Calculator").size(13),
        iced::widget::horizontal_space(),
        text("Conditions: 250 nM oligo, 50 mM Na⁺  ·  NN = SantaLucia 1998")
            .size(10)
            .color(Color::from_rgb(0.5, 0.5, 0.5)),
    ]
    .spacing(10)
    .padding([4, 8])
    .align_y(Alignment::Center);

    let header = row![
        cell_hdr("Sequence", 200),
        cell_hdr("Length (nt)", 90),
        cell_hdr("GC %", 65),
        cell_hdr("Tm NN (°C)", 90),
        cell_hdr("Tm Wallace (°C)", 110),
        cell_hdr("ΔH (kcal/mol)", 110),
        cell_hdr("ΔS (cal/mol·K)", 110),
    ]
    .spacing(0);

    let mut rows: Vec<Element<Message>> = vec![header.into(), hsep()];

    let dna_seqs: Vec<_> = aln
        .sequences
        .iter()
        .filter(|s| s.seq_type != SequenceType::Protein)
        .collect();

    if dna_seqs.is_empty() {
        rows.push(
            container(text("No nucleotide sequences in this alignment.").size(12))
                .padding(16)
                .into(),
        );
    }

    for (i, seq) in dna_seqs.iter().enumerate() {
        let bg = if i % 2 == 0 {
            Color::from_rgb(0.98, 0.98, 0.99)
        } else {
            Color::WHITE
        };

        let tm_res = calculate_tm(&seq.residues);
        let wallace = wallace_tm(&seq.residues);

        let (len_s, gc_s, tm_nn_s, tm_w_s, dh_s, ds_s) = match &tm_res {
            Some(r) => (
                r.length.to_string(),
                format!("{:.1}", r.gc_fraction * 100.0),
                format!("{:.1}", r.tm_celsius),
                format!("{:.1}", wallace),
                format!("{:.1}", r.delta_h),
                format!("{:.1}", r.delta_s),
            ),
            None => (
                "< 2".to_string(),
                "—".to_string(),
                "—".to_string(),
                "—".to_string(),
                "—".to_string(),
                "—".to_string(),
            ),
        };

        // Color-code Tm: blue < 55, green 55–65, orange 65–75, red > 75
        let tm_color = tm_res
            .as_ref()
            .map(|r| {
                let t = r.tm_celsius;
                if t < 55.0 {
                    Color::from_rgb(0.2, 0.4, 0.9)
                } else if t < 65.0 {
                    Color::from_rgb(0.1, 0.6, 0.3)
                } else if t < 75.0 {
                    Color::from_rgb(0.85, 0.5, 0.1)
                } else {
                    Color::from_rgb(0.85, 0.15, 0.15)
                }
            })
            .unwrap_or(Color::from_rgb(0.5, 0.5, 0.5));

        let r = row![
            cell_name(&seq.name, 200, bg),
            cell_str(len_s, 90, bg, None),
            cell_str(gc_s, 65, bg, None),
            cell_str(tm_nn_s, 90, bg, Some(tm_color)),
            cell_str(tm_w_s, 110, bg, None),
            cell_str(dh_s, 110, bg, None),
            cell_str(ds_s, 110, bg, None),
        ]
        .spacing(0);
        rows.push(r.into());
    }

    rows.push(
        container(
            text("Tm (NN): nearest-neighbor model  ·  Tm (Wallace): 2(A+T) + 4(G+C), useful for oligos < 14 nt")
                .size(10)
                .color(Color::from_rgb(0.5, 0.5, 0.5))
        )
        .padding([6, 8])
        .into()
    );

    column![
        toolbar,
        sep(),
        scrollable(column(rows).width(Length::Fill).spacing(0))
            .width(Length::Fill)
            .height(Length::Fill),
    ]
    .width(Length::Fill)
    .height(Length::Fill)
    .into()
}

// ── Cell helpers ──────────────────────────────────────────────────────────────

fn cell_hdr<'a>(label: &str, w: u16) -> Element<'a, Message> {
    container(text(label.to_string()).size(11))
        .width(Length::Fixed(w as f32))
        .padding([4, 6])
        .style(|_| container::Style {
            background: Some(iced::Background::Color(Color::from_rgb(0.92, 0.93, 0.95))),
            ..Default::default()
        })
        .into()
}

fn cell_name<'a>(name: &str, w: u16, bg: Color) -> Element<'a, Message> {
    let label = if name.len() > 24 {
        format!("{}…", &name[..23])
    } else {
        name.to_string()
    };
    container(text(label).size(11))
        .width(Length::Fixed(w as f32))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}

fn cell_str<'a>(s: String, w: u16, bg: Color, color: Option<Color>) -> Element<'a, Message> {
    let t = match color {
        Some(c) => text(s).size(11).color(c),
        None => text(s).size(11),
    };
    container(t)
        .width(Length::Fixed(w as f32))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}

fn sep<'a>() -> Element<'a, Message> {
    container(iced::widget::horizontal_space())
        .width(Length::Fill)
        .height(Length::Fixed(1.0))
        .style(|_| container::Style {
            background: Some(iced::Background::Color(Color::from_rgb(0.75, 0.76, 0.80))),
            ..Default::default()
        })
        .into()
}

fn hsep<'a>() -> Element<'a, Message> {
    container(iced::widget::horizontal_space())
        .width(Length::Fill)
        .height(Length::Fixed(1.0))
        .style(|_| container::Style {
            background: Some(iced::Background::Color(Color::from_rgb(0.85, 0.86, 0.88))),
            ..Default::default()
        })
        .into()
}
