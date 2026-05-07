//! Sequence composition table view.

use iced::{
    widget::{button, column, container, row, scrollable, text},
    Alignment, Color, Element, Length,
};
use helixview_core::Alignment as SeqAlignment;
use helixview_analysis::{amino_acid_composition, nucleotide_composition};
use helixview_core::SequenceType;

use crate::app::Message;

pub fn composition_view<'a>(aln: &SeqAlignment) -> Element<'a, Message> {
    let toolbar = row![
        button(text("< Back").size(12)).padding([3, 10])
            .style(button::secondary)
            .on_press(Message::CloseComposition),
        text(format!("Composition — {} sequences", aln.seq_count())).size(13),
    ]
    .spacing(10)
    .padding([4, 8])
    .align_y(Alignment::Center);

    // Detect majority sequence type.
    let is_protein = aln.sequences.iter()
        .any(|s| s.seq_type == SequenceType::Protein);

    let table: Element<Message> = if is_protein {
        protein_table(aln)
    } else {
        nucleotide_table(aln)
    };

    column![toolbar, sep(), table]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Nucleotide table ──────────────────────────────────────────────────────────

fn nucleotide_table<'a>(aln: &SeqAlignment) -> Element<'a, Message> {
    let header = row![
        cell_hdr("Sequence",  200),
        cell_hdr("A",          55),
        cell_hdr("T/U",        55),
        cell_hdr("G",          55),
        cell_hdr("C",          55),
        cell_hdr("Other",      55),
        cell_hdr("Total",      65),
        cell_hdr("GC %",       65),
        cell_hdr("AT %",       65),
        cell_hdr("MW (Da)",   100),
    ]
    .spacing(0);

    let mut rows: Vec<Element<Message>> = vec![header.into(), hsep()];

    for (i, seq) in aln.sequences.iter().enumerate() {
        let comp = nucleotide_composition(seq);
        let bg = if i % 2 == 0 {
            Color::from_rgb(0.98, 0.98, 0.99)
        } else {
            Color::WHITE
        };
        let r = row![
            cell_name(&seq.name, 200, bg),
            cell_num(comp.a,         55, bg),
            cell_num(comp.t + comp.u, 55, bg),
            cell_num(comp.g,         55, bg),
            cell_num(comp.c,         55, bg),
            cell_num(comp.other,     55, bg),
            cell_num(comp.total,     65, bg),
            cell_pct(comp.gc_pct,    65, bg),
            cell_pct(comp.at_pct,    65, bg),
            cell_f(comp.mw,         100, bg),
        ]
        .spacing(0);
        rows.push(r.into());
    }

    scrollable(
        column(rows).width(Length::Fill).spacing(0)
    )
    .width(Length::Fill)
    .height(Length::Fill)
    .into()
}

// ── Amino acid table ──────────────────────────────────────────────────────────

const AA_ORDER: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

fn protein_table<'a>(aln: &SeqAlignment) -> Element<'a, Message> {
    let mut hdr = row![cell_hdr("Sequence", 180)].spacing(0);
    for &aa in AA_ORDER {
        hdr = hdr.push(cell_hdr(std::str::from_utf8(&[aa]).unwrap(), 38));
    }
    hdr = hdr.push(cell_hdr("Other", 50))
             .push(cell_hdr("Total", 65))
             .push(cell_hdr("MW (kDa)", 90));

    let mut rows: Vec<Element<Message>> = vec![hdr.into(), hsep()];

    for (i, seq) in aln.sequences.iter().enumerate() {
        let comp = amino_acid_composition(seq);
        let bg = if i % 2 == 0 {
            Color::from_rgb(0.98, 0.98, 0.99)
        } else {
            Color::WHITE
        };
        let counts: [u64; 20] = [
            comp.a, comp.c, comp.d, comp.e, comp.f, comp.g, comp.h,
            comp.i, comp.k, comp.l, comp.m, comp.n, comp.p, comp.q,
            comp.r, comp.s, comp.t, comp.v, comp.w, comp.y,
        ];
        let mut r = row![cell_name(&seq.name, 180, bg)].spacing(0);
        for &cnt in &counts {
            r = r.push(cell_num(cnt, 38, bg));
        }
        r = r.push(cell_num(comp.other, 50, bg))
             .push(cell_num(comp.total, 65, bg))
             .push(cell_f(comp.mw / 1000.0, 90, bg));
        rows.push(r.into());
    }

    scrollable(
        column(rows).width(Length::Fill).spacing(0)
    )
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
    let label = if name.len() > 24 { format!("{}…", &name[..23]) } else { name.to_string() };
    container(text(label).size(11))
        .width(Length::Fixed(w as f32))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}

fn cell_num<'a>(n: u64, w: u16, bg: Color) -> Element<'a, Message> {
    container(text(n.to_string()).size(11))
        .width(Length::Fixed(w as f32))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}

fn cell_pct<'a>(v: f64, w: u16, bg: Color) -> Element<'a, Message> {
    container(text(format!("{:.1}", v)).size(11))
        .width(Length::Fixed(w as f32))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}

fn cell_f<'a>(v: f64, w: u16, bg: Color) -> Element<'a, Message> {
    container(text(format!("{:.1}", v)).size(11))
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
