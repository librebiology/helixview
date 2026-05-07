//! Mutual Information analysis view.
//! Two sub-views: top-pairs table and NxN heatmap canvas.

use std::cell::Cell;

use iced::{
    mouse,
    widget::{button, canvas, column, container, row, scrollable, slider, text},
    Alignment, Background, Border, Color, Element, Length, Point, Rectangle, Size,
};
use helixview_analysis::{MiResult, MiPair};

use crate::app::Message;
use crate::theme::palette;

// ── Public entry point ────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MiTab { TopPairs, Covariation, Heatmap }

pub fn mutual_info_view<'a>(
    result:              Option<&'a MiResult>,
    min_obs:             usize,
    top_n:               usize,
    stem_min_wc_types:   u8,
    stem_min_wc_frac:    f64,
    active_tab:          &MiTab,
    running:             bool,
) -> Element<'a, Message> {
    // ── Toolbar ───────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 10])
        .on_press(Message::CloseMutualInfo);

    let obs_ctrl = row![
        text("Min obs:").size(12).color(palette::TEXT_DIM),
        slider(2u32..=100, min_obs as u32, |v| Message::MiMinObs(v as usize)).width(90),
        text(format!("{min_obs}")).size(12).color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let top_ctrl = row![
        text("Top N:").size(12).color(palette::TEXT_DIM),
        slider(10u32..=500, top_n as u32, |v| Message::MiTopN(v as usize)).width(90),
        text(format!("{top_n}")).size(12).color(palette::TEXT),
    ]
    .spacing(6)
    .align_y(Alignment::Center);

    let wc_types_ctrl = row![
        text("Min WC types:").size(12).color(palette::TEXT_DIM),
        slider(1u32..=6, stem_min_wc_types as u32, |v| Message::MiStemMinWcTypes(v as u8)).width(70),
        text(format!("{stem_min_wc_types}")).size(12).color(palette::TEXT),
    ]
    .spacing(4)
    .align_y(Alignment::Center);

    let wc_frac_ctrl = row![
        text("Min WC %:").size(12).color(palette::TEXT_DIM),
        slider(0u32..=100, (stem_min_wc_frac * 100.0) as u32, |v| Message::MiStemMinWcFrac(v as f64 / 100.0)).width(70),
        text(format!("{:.0}%", stem_min_wc_frac * 100.0)).size(12).color(palette::TEXT),
    ]
    .spacing(4)
    .align_y(Alignment::Center);

    let run_btn = if running {
        button(text("Computing…").size(12)).padding([4, 12])
    } else {
        button(text("Run MI").size(12)).padding([4, 12])
            .on_press(Message::RunMutualInfo)
    };

    let tab_pairs = {
        let is_active = *active_tab == MiTab::TopPairs;
        let btn = button(text("Top Pairs").size(12)).padding([3, 10]);
        if is_active { btn.style(button::primary) } else { btn.style(button::secondary) }
            .on_press(Message::MiSetTab(MiTab::TopPairs))
    };
    let tab_cov = {
        let is_active = *active_tab == MiTab::Covariation;
        let btn = button(text("Covariation").size(12)).padding([3, 10]);
        if is_active { btn.style(button::primary) } else { btn.style(button::secondary) }
            .on_press(Message::MiSetTab(MiTab::Covariation))
    };
    let tab_heat = {
        let is_active = *active_tab == MiTab::Heatmap;
        let btn = button(text("Heatmap").size(12)).padding([3, 10]);
        if is_active { btn.style(button::primary) } else { btn.style(button::secondary) }
            .on_press(Message::MiSetTab(MiTab::Heatmap))
    };

    let status_txt = match &result {
        None => "No results — press Run MI".to_string(),
        Some(r) => format!(
            "{} cols · {} top pairs · max MI {:.3}",
            r.n_cols, r.top_pairs.len(), r.max_mi()
        ),
    };

    let toolbar_row1 = container(
        row![
            back_btn,
            text("Mutual Information").size(14).color(palette::TEXT),
            iced::widget::horizontal_space(),
            text(status_txt).size(11).color(palette::TEXT_DIM),
            run_btn,
        ]
        .spacing(10)
        .align_y(Alignment::Center)
        .padding([4, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        border: Border { color: crate::theme::palette::BORDER, width: 1.0, ..Default::default() },
        ..Default::default()
    })
    .width(Length::Fill)
    .height(Length::Shrink);

    let toolbar_row2 = container(
        row![
            tab_pairs,
            tab_cov,
            tab_heat,
            obs_ctrl,
            top_ctrl,
            wc_types_ctrl,
            wc_frac_ctrl,
        ]
        .spacing(10)
        .align_y(Alignment::Center)
        .padding([3, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgb(0.95, 0.96, 0.98))),
        border: Border { color: crate::theme::palette::BORDER, width: 1.0, ..Default::default() },
        ..Default::default()
    })
    .width(Length::Fill)
    .height(Length::Shrink);

    let content: Element<Message> = match result {
        None => container(
            text("Press \"Run MI\" to compute mutual information for all column pairs.\n\
                  Large alignments (>500 columns) may take a few seconds.")
                .size(13)
                .color(palette::TEXT_DIM),
        )
        .padding([24, 24])
        .width(Length::Fill)
        .height(Length::Fill)
        .into(),

        Some(r) => match active_tab {
            MiTab::TopPairs    => top_pairs_table(&r.top_pairs, false, stem_min_wc_types, stem_min_wc_frac),
            MiTab::Covariation => top_pairs_table(&r.top_covarying, true, stem_min_wc_types, stem_min_wc_frac),
            MiTab::Heatmap     => heatmap_view(r),
        },
    };

    column![toolbar_row1, toolbar_row2, content]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Top-pairs table ───────────────────────────────────────────────────────────

/// `show_cov` — show covariation columns and WC highlighting (Covariation tab).
fn top_pairs_table<'a>(pairs: &'a [MiPair], show_cov: bool, stem_min_wc_types: u8, stem_min_wc_frac: f64) -> Element<'a, Message> {
    if pairs.is_empty() {
        return container(
            text("No correlated pairs found.\n\
                  • Try lowering Min obs (currently set in the toolbar above)\n\
                  • Check the status bar — it shows seqs · cols · pairs after computation\n\
                  • Fully conserved alignments produce zero MI by design")
                .size(13)
                .color(palette::TEXT_DIM),
        )
        .padding([24, 24])
        .width(Length::Fill)
        .height(Length::Fill)
        .into();
    }

    fn col<'a>(s: impl Into<String>, w: f32) -> Element<'a, Message> {
        text(s.into()).size(13).color(palette::TEXT).width(Length::Fixed(w)).into()
    }
    fn col_dim<'a>(s: impl Into<String>, w: f32) -> Element<'a, Message> {
        text(s.into()).size(13).color(palette::TEXT_DIM).width(Length::Fixed(w)).into()
    }

    let mut hdr_cells: Vec<Element<Message>> = vec![
        col_dim("#",       40.0),
        col_dim("Col A",   65.0),
        col_dim("Col B",   65.0),
        col_dim("MI",      85.0),
        col_dim("N obs",   65.0),
    ];
    if show_cov {
        hdr_cells.push(col_dim("WC types", 80.0));
        hdr_cells.push(col_dim("WC %",     65.0));
        hdr_cells.push(col_dim("Cov",      70.0));
    }
    let header: Element<Message> = container(
        row(hdr_cells).spacing(0).padding([5, 12]).align_y(Alignment::Center).width(Length::Fill)
    )
    .width(Length::Fill)
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgb(0.22, 0.24, 0.30))),
        ..Default::default()
    })
    .into();

    let mut rows: Vec<Element<Message>> = vec![header];

    for (rank, pair) in pairs.iter().enumerate() {
        let is_stem = show_cov
            && pair.wc_types >= stem_min_wc_types
            && pair.wc_frac >= stem_min_wc_frac;
        let bg = if is_stem {
            Color::from_rgb(0.87, 0.96, 0.87)
        } else if rank % 2 == 0 {
            Color::WHITE
        } else {
            Color::from_rgb(0.96, 0.97, 0.99)
        };

        let jump_btn = button(text("Go").size(11))
            .padding([2, 8])
            .on_press(Message::JumpToColumn(pair.col_a));

        let mut cells: Vec<Element<Message>> = vec![
            col(format!("{}", rank + 1),            40.0),
            col(format!("{}", pair.col_a + 1),      65.0),
            col(format!("{}", pair.col_b + 1),      65.0),
            col(format!("{:.4}", pair.mi),           85.0),
            col(format!("{}", pair.n_obs),           65.0),
        ];
        if show_cov {
            cells.push(col(wc_types_label(pair.wc_types),              80.0));
            cells.push(col(format!("{:.0}%", pair.wc_frac * 100.0),    65.0));
            cells.push(col(format!("{:.3}", pair.cov_score),            70.0));
        }

        let row_el = container(
            row(cells)
                .push(iced::widget::horizontal_space())
                .push(jump_btn)
                .spacing(0)
                .padding([4, 12])
                .align_y(Alignment::Center)
                .width(Length::Fill),
        )
        .width(Length::Fill)
        .style(move |_| container::Style {
            background: Some(Background::Color(bg)),
            ..Default::default()
        });

        rows.push(row_el.into());
    }

    scrollable(column(rows).spacing(0).width(Length::Fill))
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

/// Human-readable WC type count with asterisk for strong signal.
fn wc_types_label(n: u8) -> String {
    match n {
        0 => "-".to_string(),
        1 => "1".to_string(),
        2 => "2 **".to_string(),
        3 => "3 ***".to_string(),
        4..=6 => format!("{n} ****"),
        _ => format!("{n}"),
    }
}


// ── Heatmap canvas ────────────────────────────────────────────────────────────

fn heatmap_view<'a>(result: &'a MiResult) -> Element<'a, Message> {
    // Pre-extract what the canvas needs so it doesn't borrow `result`.
    let n_cols = result.n_cols;
    let max_mi = result.max_mi();
    let values = result.values.clone();

    canvas(MiHeatmap { n_cols, max_mi, values })
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

struct MiHeatmap {
    n_cols: usize,
    max_mi: f64,
    values: Vec<f64>,
}

pub struct MiHeatmapState {
    cache:    canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for MiHeatmapState {
    fn default() -> Self {
        Self { cache: canvas::Cache::default(), last_key: Cell::new(u64::MAX) }
    }
}

impl canvas::Program<Message> for MiHeatmap {
    type State = MiHeatmapState;

    fn draw(
        &self,
        state:    &MiHeatmapState,
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        let key = (self.n_cols as u64)
            .wrapping_add(self.values.len() as u64)
            .wrapping_add(bounds.width.to_bits() as u64)
            .wrapping_add(bounds.height.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let n_cols  = self.n_cols;
        let max_mi  = self.max_mi;
        let values  = &self.values;
        let geom = state.cache.draw(renderer, bounds.size(), |frame| {
            draw_heatmap(frame, bounds.size(), n_cols, max_mi, values);
        });
        vec![geom]
    }
}

fn draw_heatmap(
    frame:  &mut iced::widget::canvas::Frame<iced::Renderer>,
    size:   Size,
    n:      usize,
    max_mi: f64,
    values: &[f64],
) {
    frame.fill_rectangle(Point::ORIGIN, size, Color::WHITE);

    if n == 0 { return; }
    let max_mi = max_mi.max(1e-9);
    let cell_w = size.width  / n as f32;
    let cell_h = size.height / n as f32;

    let get = |i: usize, j: usize| -> f64 {
        if i == j { return 0.0; }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        let idx = a * n - a * (a + 1) / 2 + (b - a - 1);
        values.get(idx).copied().unwrap_or(0.0)
    };

    // Draw cells — only upper triangle (mirror it).
    for i in 0..n {
        for j in (i+1)..n {
            let mi = get(i, j);
            if mi < 1e-9 { continue; }

            let t = (mi / max_mi) as f32;
            let color = mi_color(t);

            let x = i as f32 * cell_w;
            let y = j as f32 * cell_h;
            frame.fill_rectangle(Point::new(x, y), Size::new(cell_w.max(1.0), cell_h.max(1.0)), color);

            // Mirror.
            let xm = j as f32 * cell_w;
            let ym = i as f32 * cell_h;
            frame.fill_rectangle(Point::new(xm, ym), Size::new(cell_w.max(1.0), cell_h.max(1.0)), color);
        }
    }

    // Diagonal (self = 0).
    for i in 0..n {
        let x = i as f32 * cell_w;
        let y = i as f32 * cell_h;
        frame.fill_rectangle(
            Point::new(x, y),
            Size::new(cell_w.max(1.0), cell_h.max(1.0)),
            Color::from_rgb(0.92, 0.92, 0.92),
        );
    }

    // Axis labels (if cells are large enough).
    if cell_w >= 12.0 {
        for i in (0..n).step_by((n / 10).max(1)) {
            let x = i as f32 * cell_w + cell_w / 2.0;
            frame.fill_text(canvas::Text {
                content: format!("{}", i + 1),
                position: Point::new(x, size.height - 12.0),
                size: iced::Pixels(9.0),
                color: Color::from_rgb(0.3, 0.3, 0.35),
                horizontal_alignment: iced::alignment::Horizontal::Center,
                vertical_alignment:   iced::alignment::Vertical::Top,
                ..Default::default()
            });
        }
    }
}

/// Map a 0..1 fraction to a white→blue→red color scale (like Matplotlib "hot" rotated).
fn mi_color(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    if t < 0.5 {
        let s = t * 2.0;
        Color::from_rgb(s * 0.18, s * 0.45, 1.0 - s * 0.20)
    } else {
        let s = (t - 0.5) * 2.0;
        Color::from_rgb(0.18 + s * 0.72, 0.45 - s * 0.35, 0.80 - s * 0.70)
    }
}
