//! Taxonomy table: shows NCBI lineage parsed from GenBank SOURCE/ORGANISM blocks.

use crate::app::Message;
use crate::theme::palette;
use helixview_core::Alignment as SeqAlignment;
use iced::{
    widget::{button, column, container, horizontal_space, row, scrollable, text},
    Alignment, Color, Element, Length,
};
use std::sync::Arc;

// ── Lineage parsing ───────────────────────────────────────────────────────────

/// Parse organism name and lineage from a GenBank SOURCE block string.
/// The block looks like:
///   "Mus musculus (house mouse)\nORGANISM  Mus musculus\nEukaryota; Metazoa; ..."
fn parse_source(source: &str) -> (String, Vec<String>) {
    let mut organism = String::new();
    let mut lineage_parts: Vec<String> = Vec::new();
    let mut in_lineage = false;
    let mut lineage_buf = String::new();

    for line in source.lines() {
        if let Some(rest) = line.strip_prefix("ORGANISM") {
            organism = rest.trim().to_string();
            in_lineage = true;
        } else if in_lineage {
            if !lineage_buf.is_empty() {
                lineage_buf.push(' ');
            }
            lineage_buf.push_str(line.trim());
        }
    }

    if lineage_buf.is_empty() {
        // Fallback: treat the first line of source as organism if ORGANISM absent
        if let Some(first) = source.lines().next() {
            organism = first.trim().to_string();
        }
    } else {
        lineage_parts = lineage_buf
            .split(';')
            .map(|s| s.trim().trim_end_matches('.').trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
    }

    (organism, lineage_parts)
}

/// Per-sequence parsed taxonomy data.
struct TaxRow {
    seq_idx: usize,
    seq_name: String,
    organism: String,
    ranks: Vec<String>,
}

fn build_rows(aln: &SeqAlignment) -> (Vec<TaxRow>, usize) {
    let rows: Vec<TaxRow> = aln
        .sequences
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            let (organism, ranks) = seq
                .genbank
                .source
                .as_deref()
                .map(parse_source)
                .unwrap_or_default();
            TaxRow {
                seq_idx: i,
                seq_name: seq.name.clone(),
                organism,
                ranks,
            }
        })
        .collect();

    let max_ranks = rows.iter().map(|r| r.ranks.len()).max().unwrap_or(0);
    (rows, max_ranks)
}

// ── Column header labels ──────────────────────────────────────────────────────

/// Best-effort rank labels for the first N levels of a typical NCBI bacterial lineage.
/// The real lineage has variable depth, so beyond the presets we use level numbers.
const RANK_LABELS: &[&str] = &[
    "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species",
];

fn rank_label(i: usize) -> String {
    RANK_LABELS
        .get(i)
        .map(|s| s.to_string())
        .unwrap_or_else(|| format!("Level {}", i + 1))
}

// ── View ──────────────────────────────────────────────────────────────────────

const COL_IDX: f32 = 36.0;
const COL_NAME: f32 = 180.0;
const COL_ORG: f32 = 200.0;
const COL_RANK: f32 = 130.0;

fn hdr_cell<'a>(label: impl Into<String>, width: f32) -> Element<'a, Message> {
    container(text(label.into()).size(11))
        .width(Length::Fixed(width))
        .padding([3, 6])
        .style(|_| container::Style {
            background: Some(iced::Background::Color(Color::from_rgb(0.88, 0.89, 0.92))),
            ..Default::default()
        })
        .into()
}

fn sort_hdr_cell<'a>(
    label: impl Into<String>,
    width: f32,
    rank_idx: usize,
) -> Element<'a, Message> {
    button(text(label.into()).size(11))
        .width(Length::Fixed(width))
        .padding([3, 6])
        .style(button::secondary)
        .on_press(Message::SortByTaxonomy(rank_idx))
        .into()
}

fn data_cell<'a>(label: impl Into<String>, width: f32, shade: bool) -> Element<'a, Message> {
    let bg = if shade {
        Color::from_rgb(0.97, 0.97, 0.99)
    } else {
        Color::WHITE
    };
    container(text(label.into()).size(11))
        .width(Length::Fixed(width))
        .padding([3, 6])
        .style(move |_| container::Style {
            background: Some(iced::Background::Color(bg)),
            ..Default::default()
        })
        .into()
}

pub fn taxonomy_view<'a>(aln: &Arc<SeqAlignment>) -> Element<'a, Message> {
    let (rows, max_ranks) = build_rows(aln);
    let num_rank_cols = max_ranks.min(8);

    let close_btn = button(text("< Back").size(12))
        .padding([4, 10])
        .on_press(Message::CloseTaxonomy);
    let sort_name_btn = button(text("Sort by Organism").size(12))
        .padding([4, 10])
        .style(button::secondary)
        .on_press(Message::SortByTaxonomy(usize::MAX));
    let reorder_note = text("Click a rank column header to reorder alignment sequences")
        .size(11)
        .color(palette::TEXT_DIM);
    let toolbar = row![close_btn, sort_name_btn, horizontal_space(), reorder_note]
        .spacing(8)
        .align_y(Alignment::Center)
        .padding([4, 8]);

    // ── Header row ────────────────────────────────────────────────────────────
    let mut hdr = row![
        hdr_cell("#", COL_IDX),
        hdr_cell("Sequence", COL_NAME),
        hdr_cell("Organism", COL_ORG),
    ]
    .spacing(1);
    for i in 0..num_rank_cols {
        hdr = hdr.push(sort_hdr_cell(rank_label(i), COL_RANK, i));
    }
    let hdr: Element<'a, Message> = hdr.spacing(1).into();

    // ── Data rows ─────────────────────────────────────────────────────────────
    let mut body = column![].spacing(0);
    for (row_i, tr) in rows.iter().enumerate() {
        let shade = row_i % 2 == 1;
        let mut r = row![
            data_cell(format!("{}", tr.seq_idx + 1), COL_IDX, shade),
            data_cell(tr.seq_name.clone(), COL_NAME, shade),
            data_cell(tr.organism.clone(), COL_ORG, shade),
        ]
        .spacing(1);
        for i in 0..num_rank_cols {
            let val = tr.ranks.get(i).cloned().unwrap_or_default();
            r = r.push(data_cell(val, COL_RANK, shade));
        }
        body = body.push(r);
    }

    let content: Element<'a, Message> = scrollable(
        scrollable(column![hdr, body].spacing(1).padding(8))
            .direction(scrollable::Direction::Vertical(scrollable::Scrollbar::new())),
    )
    .direction(scrollable::Direction::Horizontal(
        scrollable::Scrollbar::new(),
    ))
    .into();

    column![toolbar, content]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}
