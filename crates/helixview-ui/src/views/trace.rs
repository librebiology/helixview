//! ABI chromatogram trace viewer.
//!
//! Renders a four-channel Sanger sequencing trace with:
//! - Coloured polylines for each detection channel (A=green, C=blue, G=black, T=red)
//! - Base-call letters drawn above their peak positions
//! - Horizontal scroll (0.0–1.0) and zoom (1×–16×) controls
//! - Peak-height quantification panel for detecting heterozygous/edited sites

use std::cell::Cell;
use std::sync::Arc;

use iced::{
    alignment,
    mouse,
    widget::{button, canvas, column, container, row, scrollable, slider, text},
    Alignment as IAlign, Background, Color, Element, Font, Length, Pixels, Point,
    Rectangle, Size,
};
use iced::widget::canvas::{Frame, Path, Text as CanvasText};

use helixview_formats::abi::AbiTrace;
use crate::app::Message;
use crate::theme::palette;

// ── Base channel colors ───────────────────────────────────────────────────────

fn base_color(base: u8) -> Color {
    match base.to_ascii_uppercase() {
        b'A' => Color::from_rgb(0.07, 0.70, 0.12), // green
        b'C' => Color::from_rgb(0.08, 0.30, 0.88), // blue
        b'T' => Color::from_rgb(0.86, 0.12, 0.12), // red
        b'G' => Color::from_rgb(0.10, 0.10, 0.10), // near-black
        _    => Color::from_rgb(0.60, 0.60, 0.60), // grey fallback
    }
}

// ── Public entry-point ────────────────────────────────────────────────────────

/// Build the full trace viewer: toolbar + chromatogram canvas + quant panel.
pub fn trace_view<'a>(
    trace:  &Arc<AbiTrace>,
    scroll: f32,
    zoom:   f32,
    ch_a:   u8,
    ch_b:   u8,
) -> Element<'a, Message> {
    // ── Toolbar ──────────────────────────────────────────────────────────────
    let back_btn = button(text("< Back").size(13))
        .padding([4, 10])
        .on_press(Message::CloseTrace);

    let info_lbl = text(format!(
        "ABI Trace  ·  {} bases  ·  {} samples",
        trace.bases.len(),
        trace.num_samples,
    ))
    .size(13)
    .color(palette::TEXT);

    // Channel labels in order
    let channel_key = trace.channel_bases.iter().enumerate().fold(
        row![].spacing(6).align_y(IAlign::Center),
        |r, (i, &b)| {
            let dot = container(text("■").size(12).color(base_color(b)))
                .padding([0, 2]);
            let lbl = text(format!(
                "ch{}: {}",
                i + 1,
                b as char
            ))
            .size(11)
            .color(palette::TEXT_DIM);
            r.push(dot).push(lbl)
        },
    );

    let zoom_lbl  = text(format!("{:.1}×", zoom)).size(11).color(palette::TEXT_DIM);
    let zoom_sl   = slider(1.0f32..=16.0, zoom, Message::TraceZoom).width(80);
    let zoom_row  = row![text("Zoom:").size(11).color(palette::TEXT_DIM), zoom_sl, zoom_lbl]
        .spacing(4).align_y(IAlign::Center);

    let scroll_sl = slider(0.0f32..=1.0, scroll, Message::TraceScroll).width(120);
    let scroll_row = row![
        text("Scroll:").size(11).color(palette::TEXT_DIM),
        scroll_sl,
    ]
    .spacing(4).align_y(IAlign::Center);

    let toolbar = container(
        row![back_btn, info_lbl,
             iced::widget::horizontal_space(),
             channel_key,
             zoom_row,
             scroll_row,
        ]
        .spacing(12)
        .align_y(IAlign::Center)
        .padding([5, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        ..Default::default()
    })
    .width(Length::Fill);

    // ── Chromatogram canvas ───────────────────────────────────────────────────
    let chrom = ChromatogramCanvas {
        trace: Arc::clone(trace),
        scroll,
        zoom,
    };
    let chrom_canvas = canvas(chrom)
        .width(Length::Fill)
        .height(Length::Fill);

    // ── Quantification panel ─────────────────────────────────────────────────
    let quant = quant_panel(trace, ch_a, ch_b);

    column![toolbar, chrom_canvas, quant]
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Quantification panel ──────────────────────────────────────────────────────

fn quant_panel<'a>(trace: &Arc<AbiTrace>, ch_a: u8, ch_b: u8) -> Element<'a, Message> {
    let bases_for: [u8; 4] = [b'A', b'C', b'G', b'T'];

    // ── Channel selector row ──────────────────────────────────────────────────
    let make_cycle_btns = |current: u8, mk_msg: fn(u8) -> Message| {
        let btns = bases_for.iter().map(|&b| {
            let label = (b as char).to_string();
            let btn = button(text(label).size(11).color(base_color(b)))
                .padding([2, 7]);
            if b == current {
                btn.style(|_theme, _state| button::Style {
                    background: Some(Background::Color(Color::from_rgb(0.85, 0.92, 1.0))),
                    border: iced::Border { radius: 3.0.into(), width: 1.0, color: Color::from_rgb(0.4, 0.6, 0.9) },
                    ..Default::default()
                })
            } else {
                btn.on_press(mk_msg(b))
            }
        });
        let mut r = row![].spacing(2).align_y(IAlign::Center);
        for b in btns { r = r.push(b); }
        r
    };

    let ch_a_sel = row![
        text("Channel A:").size(11).color(palette::TEXT_DIM),
        make_cycle_btns(ch_a, Message::TraceQuantChA),
    ].spacing(4).align_y(IAlign::Center);

    let ch_b_sel = row![
        text("Channel B:").size(11).color(palette::TEXT_DIM),
        make_cycle_btns(ch_b, Message::TraceQuantChB),
    ].spacing(4).align_y(IAlign::Center);

    // ── Per-base-type summary ─────────────────────────────────────────────────
    let idx_a = trace.channel_for_base(ch_a);
    let idx_b = trace.channel_for_base(ch_b);

    // Compute means per called-base type.
    struct BaseStat { count: u32, sum_a: i64, sum_b: i64 }
    let mut stats: std::collections::HashMap<u8, BaseStat> = std::collections::HashMap::new();

    for (i, &called) in trace.bases.iter().enumerate() {
        let ha = idx_a.map(|c| trace.peak_height(i, c)).unwrap_or(0) as i64;
        let hb = idx_b.map(|c| trace.peak_height(i, c)).unwrap_or(0) as i64;
        let e  = stats.entry(called.to_ascii_uppercase()).or_insert(BaseStat { count: 0, sum_a: 0, sum_b: 0 });
        e.count += 1;
        e.sum_a += ha;
        e.sum_b += hb;
    }

    // Build summary table rows.
    let header = container(
        row![
            cell("Base".to_string(),  50),
            cell("Count".to_string(), 55),
            cell("h_A".to_string(),   60),
            cell("h_B".to_string(),   60),
            cell("P(%)".to_string(),  60),
        ]
        .spacing(0),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgb(0.92, 0.94, 0.97))),
        ..Default::default()
    })
    .padding([2, 6]);

    let mut sorted_bases: Vec<u8> = stats.keys().copied().collect();
    sorted_bases.sort();

    let rows: Vec<Element<Message>> = sorted_bases.iter().map(|&b| {
        let st = &stats[&b];
        let mean_a = if st.count > 0 { st.sum_a / st.count as i64 } else { 0 };
        let mean_b = if st.count > 0 { st.sum_b / st.count as i64 } else { 0 };
        let pct = if mean_a + mean_b > 0 {
            format!("{:.1}%", mean_a as f64 / (mean_a + mean_b) as f64 * 100.0)
        } else {
            "—".to_string()
        };
        container(
            row![
                cell_colored((b as char).to_string(), 50, base_color(b)),
                cell(st.count.to_string(),  55),
                cell(mean_a.to_string(),    60),
                cell(mean_b.to_string(),    60),
                cell(pct,                   60),
            ]
            .spacing(0),
        )
        .padding([1, 6])
        .into()
    }).collect();

    let mut table_col = column![header].spacing(0);
    for r in rows { table_col = table_col.push(r); }

    // Formula note
    let formula_note = text(format!(
        "P(%) = h_{} / (h_{} + h_{}) × 100   (peak heights at called bases)",
        ch_a as char, ch_a as char, ch_b as char,
    ))
    .size(10)
    .color(palette::TEXT_DIM);

    let panel = container(
        column![
            row![ch_a_sel, iced::widget::horizontal_space(), ch_b_sel]
                .spacing(16)
                .align_y(IAlign::Center),
            table_col,
            formula_note,
        ]
        .spacing(4)
        .padding([6, 10]),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(Color::from_rgb(0.97, 0.97, 0.99))),
        border: iced::Border {
            color: palette::BORDER,
            width: 1.0,
            radius: 0.0.into(),
        },
        ..Default::default()
    })
    .width(Length::Fill);

    scrollable(panel)
        .height(Length::Fixed(160.0))
        .into()
}

fn cell<'a>(s: String, w: u16) -> Element<'a, Message> {
    container(text(s).size(11))
        .width(Length::Fixed(w as f32))
        .padding([1, 4])
        .into()
}

fn cell_colored<'a>(s: String, w: u16, color: Color) -> Element<'a, Message> {
    container(text(s).size(11).color(color))
        .width(Length::Fixed(w as f32))
        .padding([1, 4])
        .into()
}

// ── Chromatogram canvas ───────────────────────────────────────────────────────

struct ChromatogramCanvas {
    trace:  Arc<AbiTrace>,
    scroll: f32,
    zoom:   f32,
}

pub struct ChromatogramState {
    cache:    canvas::Cache,
    last_key: Cell<u64>,
}

impl Default for ChromatogramState {
    fn default() -> Self {
        Self { cache: canvas::Cache::default(), last_key: Cell::new(u64::MAX) }
    }
}

impl canvas::Program<Message> for ChromatogramCanvas {
    type State = ChromatogramState;

    fn draw(
        &self,
        state:    &ChromatogramState,
        renderer: &iced::Renderer,
        _theme:   &iced::Theme,
        bounds:   Rectangle,
        _cursor:  mouse::Cursor,
    ) -> Vec<canvas::Geometry<iced::Renderer>> {
        // Cache invalidation: any change to parameters or window size triggers redraw.
        let key = (Arc::as_ptr(&self.trace) as u64)
            .wrapping_add(self.scroll.to_bits() as u64)
            .wrapping_add(self.zoom.to_bits() as u64)
            .wrapping_add(bounds.width.to_bits()  as u64)
            .wrapping_add(bounds.height.to_bits() as u64);

        if key != state.last_key.get() {
            state.cache.clear();
            state.last_key.set(key);
        }

        let geom = state.cache.draw(renderer, bounds.size(), |frame| {
            draw_chromatogram(frame, bounds.size(), &self.trace, self.scroll, self.zoom);
        });

        vec![geom]
    }
}

// ── Drawing ───────────────────────────────────────────────────────────────────

const LABEL_AREA_H: f32 = 18.0;  // height reserved above trace for base labels
const AXIS_PAD_L:   f32 = 6.0;
const AXIS_PAD_B:   f32 = 4.0;

fn draw_chromatogram(
    frame:  &mut Frame<iced::Renderer>,
    size:   Size,
    trace:  &AbiTrace,
    scroll: f32,
    zoom:   f32,
) {
    // White background
    frame.fill_rectangle(Point::ORIGIN, size, Color::WHITE);

    if trace.num_samples == 0 {
        draw_placeholder(frame, size, "No trace data");
        return;
    }

    let total     = trace.num_samples;
    let zoom      = zoom.max(1.0_f32).min(64.0_f32);
    // Number of samples visible in the current window.
    let visible   = ((total as f32) / zoom).ceil() as usize;
    let visible   = visible.min(total).max(1);
    // Start sample derived from scroll fraction.
    let max_start = total.saturating_sub(visible);
    let start     = ((scroll.clamp(0.0, 1.0) * max_start as f32) as usize).min(max_start);
    let end       = (start + visible).min(total);

    // Drawing area
    let draw_w = (size.width  - AXIS_PAD_L).max(1.0);
    let draw_h = (size.height - LABEL_AREA_H - AXIS_PAD_B).max(1.0);
    let max_val = trace.max_value().max(1) as f32;

    let px_per_sample = draw_w / (end - start) as f32;

    // ── Draw each channel as a polyline ───────────────────────────────────────
    for (ch_idx, channel) in trace.channels.iter().enumerate() {
        if channel.is_empty() { continue; }
        let base  = trace.channel_bases[ch_idx];
        let color = base_color(base);

        // Collect visible (x, y) pairs.
        let points: Vec<Point> = (start..end)
            .filter_map(|s| {
                let val  = *channel.get(s)? as f32;
                let x    = AXIS_PAD_L + (s - start) as f32 * px_per_sample;
                let y    = LABEL_AREA_H + draw_h * (1.0 - (val.max(0.0) / max_val));
                Some(Point::new(x, y))
            })
            .collect();

        if points.len() < 2 { continue; }

        let path = Path::new(|b| {
            b.move_to(points[0]);
            for &pt in &points[1..] {
                b.line_to(pt);
            }
        });

        frame.stroke(
            &path,
            canvas::Stroke::default()
                .with_color(color)
                .with_width(1.5),
        );
    }

    // ── Base call labels at peak positions ────────────────────────────────────
    for (i, &loc) in trace.peak_locs.iter().enumerate() {
        let s = loc as usize;
        if s < start || s >= end { continue; }
        let Some(&base) = trace.bases.get(i) else { continue };
        let color = base_color(base);
        let x = AXIS_PAD_L + (s - start) as f32 * px_per_sample;

        // Vertical tick mark at peak
        frame.fill_rectangle(
            Point::new(x, LABEL_AREA_H),
            Size::new(1.0, 4.0),
            Color { a: 0.35, ..color },
        );

        // Base letter
        frame.fill_text(CanvasText {
            content:  (base as char).to_string(),
            position: Point::new(x, LABEL_AREA_H * 0.5),
            color,
            size:      Pixels(10.0),
            font:      Font::MONOSPACE,
            horizontal_alignment: alignment::Horizontal::Center,
            vertical_alignment:   alignment::Vertical::Center,
            ..CanvasText::default()
        });
    }

    // ── Bottom axis line ──────────────────────────────────────────────────────
    frame.fill_rectangle(
        Point::new(AXIS_PAD_L, LABEL_AREA_H + draw_h),
        Size::new(draw_w, 1.0),
        palette::BORDER,
    );

    // ── Sample-range annotation (bottom-right) ────────────────────────────────
    frame.fill_text(CanvasText {
        content:  format!("{start}–{end}"),
        position: Point::new(size.width - 4.0, size.height - 2.0),
        color:    palette::TEXT_DIM,
        size:     Pixels(9.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Right,
        vertical_alignment:   alignment::Vertical::Bottom,
        ..CanvasText::default()
    });
}

fn draw_placeholder(frame: &mut Frame<iced::Renderer>, size: Size, msg: &str) {
    frame.fill_text(CanvasText {
        content:  msg.to_string(),
        position: Point::new(size.width * 0.5, size.height * 0.5),
        color:    palette::TEXT_DIM,
        size:     Pixels(14.0),
        font:     Font::MONOSPACE,
        horizontal_alignment: alignment::Horizontal::Center,
        vertical_alignment:   alignment::Vertical::Center,
        ..CanvasText::default()
    });
}
