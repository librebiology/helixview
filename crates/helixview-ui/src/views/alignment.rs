use std::sync::Arc;

use helixview_core::Alignment as SeqAlignment;
use iced::widget::horizontal_space;
use iced::{
    widget::{button, canvas, column, container, row, text, text_input},
    Alignment, Background, Border, Color, Element, Length,
};

use crate::app::{ColorScheme, Message, SearchMode, ToolbarTab, ViewState};
use crate::color_table::ColorTable;
use crate::theme::palette;
use crate::views::analysis_panel;
use crate::widgets::grid::{CELL_H, CELL_W, FONT_SIZE};
use crate::widgets::{
    AlignmentGrid, EntropyStrip, PairingArcsStrip, PairingPair, RestrMapStrip, ARCS_H, RESTR_H,
    STRIP_H,
};
use helixview_external::ExternalAligner;

/// Full alignment editor: toolbar + canvas grid + horizontal scrollbar + status bar.
pub fn alignment_view<'a>(
    aln: &'a Arc<SeqAlignment>,
    view: &'a ViewState,
    status: &'a str,
    title_menu: Option<usize>,
    rename_input: &'a str,
    show_restr_map: bool,
    color_table: &'a Arc<ColorTable>,
    computing: Option<&'a str>,
    pairing_pairs: Vec<PairingPair>,
    pairing_max: f64,
    show_pairing_arcs: bool,
    pairing_threshold: f64,
    mi_available: bool,
    mi_stale: bool,
) -> Element<'a, Message> {
    // ── Toolbar: two-row ribbon ───────────────────────────────────────────────
    use crate::theme::buttons as btn_style;

    // Row 1: always-visible essentials
    let open_btn = button(text("Open").size(12))
        .padding([4, 12])
        .style(btn_style::primary)
        .on_press(Message::OpenFileDialog);
    let save_btn = button(text("Save").size(12))
        .padding([4, 12])
        .style(btn_style::save)
        .on_press(Message::SaveFileDialog);
    let proj_open_btn = button(text("Open Project").size(12))
        .padding([3, 8])
        .style(btn_style::secondary)
        .on_press(Message::OpenProjectDialog);
    let proj_save_btn = button(text("Save Project").size(12))
        .padding([3, 8])
        .style(btn_style::secondary)
        .on_press(Message::SaveProjectDialog);
    let fetch_btn = button(text("Fetch…").size(12))
        .padding([3, 8])
        .style(btn_style::secondary)
        .on_press(Message::ToggleFetchBar);
    let export_btn = button(text("Export…").size(12))
        .padding([3, 8])
        .style(btn_style::secondary)
        .on_press(Message::ExportFastaDialog);
    let find_btn = button(text("Find  ⌕").size(12))
        .padding([3, 8])
        .style(btn_style::ghost)
        .on_press(Message::ToggleSearch);

    let zoom_out_btn = button(text("−").size(14))
        .padding([2, 7])
        .style(btn_style::ghost)
        .on_press(Message::ZoomOut);
    let zoom_in_btn = button(text("+").size(14))
        .padding([2, 7])
        .style(btn_style::ghost)
        .on_press(Message::ZoomIn);
    let zoom_label = text(format!("{:.0}%", view.zoom * 100.0))
        .size(11)
        .color(palette::TEXT_DIM);
    let zoom_row = row![zoom_out_btn, zoom_label, zoom_in_btn]
        .spacing(1)
        .align_y(Alignment::Center);

    // Tab selector buttons — active tab is green, inactive is muted
    let tab_btn = |label: &str, tab: ToolbarTab| -> iced::widget::Button<'static, Message> {
        let active = view.toolbar_tab == tab;
        button(text(label.to_string()).size(12))
            .padding([4, 12])
            .style(if active {
                btn_style::tab_active
            } else {
                btn_style::tab_inactive
            })
            .on_press(Message::SetToolbarTab(tab))
    };

    let doc_info = row![
        text(&aln.name).size(13).color(palette::TEXT),
        text(format!("  {}×{}", aln.seq_count(), aln.col_count()))
            .size(11)
            .color(palette::TEXT_DIM),
    ]
    .align_y(Alignment::Center);

    let toolbar_row1 = container(
        row![
            open_btn,
            save_btn,
            vsep(),
            proj_open_btn,
            proj_save_btn,
            vsep(),
            fetch_btn,
            export_btn,
            vsep(),
            doc_info,
            horizontal_space(),
            find_btn,
            vsep(),
            zoom_row,
            vsep(),
            button(text("Prefs").size(12))
                .padding([3, 8])
                .style(btn_style::ghost)
                .on_press(Message::OpenPrefs),
            button(text("?").size(12))
                .padding([3, 8])
                .style(btn_style::ghost)
                .on_press(Message::ShowHelp),
        ]
        .spacing(5)
        .align_y(Alignment::Center)
        .padding([4, 10])
        .width(Length::Fill),
    )
    .style(|_| container::Style {
        background: Some(Background::Color(palette::HEADER_BG)),
        border: Border {
            color: palette::BORDER,
            width: 1.0,
            ..Default::default()
        },
        ..Default::default()
    })
    .width(Length::Fill);

    // Row 2: tab-dependent buttons
    let toolbar_row2 = {
        // Shorthand helpers that apply a consistent style per category.
        let sec = |label: &str, msg: Message| -> iced::widget::Button<'static, Message> {
            button(text(label.to_string()).size(12))
                .padding([3, 9])
                .style(btn_style::secondary)
                .on_press(msg)
        };
        let ana = |label: &str, msg: Message| -> iced::widget::Button<'static, Message> {
            button(text(label.to_string()).size(12))
                .padding([3, 9])
                .style(btn_style::analysis)
                .on_press(msg)
        };

        let inner: Element<Message> = match view.toolbar_tab {
            ToolbarTab::Edit => {
                let features_label = if view.show_features {
                    "Features ●"
                } else {
                    "Features ○"
                };
                let restr_label = if show_restr_map {
                    "RE Map ●"
                } else {
                    "RE Map ○"
                };
                let feat_style = if view.show_features {
                    btn_style::toggle_on
                } else {
                    btn_style::toggle_off
                };
                let re_style = if show_restr_map {
                    btn_style::toggle_on
                } else {
                    btn_style::toggle_off
                };
                row![
                    sec("Minimize", Message::MinimizeAlignment),
                    sec("Del Gaps", Message::DeleteGapColumns),
                    sec("Sort A–Z", Message::SortByName),
                    sec("Rev-Comp", Message::ReverseComplementSelected),
                    sec("Translate", Message::TranslateSelected),
                    sec("Consensus", Message::AddConsensus),
                    vsep(),
                    sec("Annotate", Message::ToggleAnnotateBar),
                    button(text(features_label).size(12))
                        .padding([3, 9])
                        .style(feat_style)
                        .on_press(Message::ToggleFeatures),
                    button(text(restr_label).size(12))
                        .padding([3, 9])
                        .style(re_style)
                        .on_press(Message::ToggleRestrMap),
                ]
                .spacing(4)
                .align_y(Alignment::Center)
                .into()
            }

            ToolbarTab::Analysis => {
                let protein_label = if view.protein_view_frame.is_some() {
                    "< DNA"
                } else {
                    "Protein"
                };
                let protein_msg = if view.protein_view_frame.is_some() {
                    Message::ExitProteinView
                } else {
                    Message::SetProteinViewFrame(0)
                };
                row![
                    ana("Pairwise", Message::RunPairwise),
                    ana("ORFs", Message::RunOrfFinder),
                    ana("6-Frame", Message::RunSixFrame),
                    ana(protein_label, protein_msg),
                    vsep(),
                    ana("BLAST", Message::ToggleBlastPanel),
                    ana("ABI Trace", Message::OpenTraceDialog),
                    ana("Phylo Tree", Message::BuildTreeFromAlignment),
                    ana("Plasmid Map", Message::ShowPlasmid),
                    vsep(),
                    ana("Dot Plot", Message::ShowDotPlot),
                    ana("Col Stats", Message::ShowColSummary),
                    ana("Conserved", Message::ShowConservation),
                    ana("Identity", Message::ShowIdentityMatrix),
                    ana("Composition", Message::ShowComposition),
                    ana("Mutual Info", Message::ShowMutualInfo),
                    ana("Hydrophob.", Message::ShowHydrophobicity),
                    ana("Oligo Tm", Message::ShowOligoTm),
                    ana("Taxonomy", Message::ShowTaxonomy),
                ]
                .spacing(4)
                .align_y(Alignment::Center)
                .into()
            }

            ToolbarTab::View => {
                let analysis_label = if view.show_analysis {
                    "Side <"
                } else {
                    "Side >"
                };
                let arc_style = if show_pairing_arcs {
                    btn_style::toggle_on
                } else {
                    btn_style::toggle_off
                };
                let arc_label = if show_pairing_arcs {
                    "Cov Arcs ●"
                } else {
                    "Cov Arcs ○"
                };
                row![
                    scheme_btn("Residue", ColorScheme::ResidueType, &view.color_scheme),
                    scheme_btn("Identity", ColorScheme::Identity, &view.color_scheme),
                    scheme_btn("Plain", ColorScheme::Plain, &view.color_scheme),
                    scheme_btn("Strength", ColorScheme::Strength, &view.color_scheme),
                    vsep(),
                    sec("Colors", Message::OpenColorEditor),
                    sec(analysis_label, Message::ToggleAnalysis),
                    sec("Export SVG", Message::ExportSvgDialog),
                    sec("Shaded Fig", Message::ShowShadedExportDialog),
                    sec("Export PDF", Message::ExportPdfDialog),
                    sec("Export Text", Message::ShowTextExport),
                    vsep(),
                    button(text(arc_label).size(12))
                        .padding([3, 9])
                        .style(arc_style)
                        .on_press(Message::TogglePairingArcs),
                    iced::widget::slider(1u32..=20, pairing_threshold as u32, |v| {
                        Message::PairingArcThreshold(v as f64)
                    },)
                    .width(70),
                    text(format!("≥{:.0}", pairing_threshold))
                        .size(10)
                        .color(palette::TEXT_DIM),
                    vsep(),
                    sec("Accessories", Message::ShowAccessories),
                    sec("Create BLAST DB", Message::CreateBlastDbDialog),
                    vsep(),
                    sec(
                        "ClustalW",
                        Message::RunExternalAlign(ExternalAligner::ClustalW)
                    ),
                    sec(
                        "ClustalΩ",
                        Message::RunExternalAlign(ExternalAligner::ClustalOmega)
                    ),
                    sec("MUSCLE", Message::RunExternalAlign(ExternalAligner::Muscle)),
                    sec("MAFFT", Message::RunExternalAlign(ExternalAligner::Mafft)),
                ]
                .spacing(4)
                .align_y(Alignment::Center)
                .into()
            }
        };

        container(
            row![
                tab_btn("Edit", ToolbarTab::Edit),
                tab_btn("Analysis", ToolbarTab::Analysis),
                tab_btn("View", ToolbarTab::View),
                vsep(),
                inner,
            ]
            .padding([4, 10])
            .spacing(5)
            .width(Length::Fill)
            .align_y(Alignment::Center),
        )
        .style(|_| container::Style {
            background: Some(Background::Color(palette::HEADER_BG)),
            border: Border {
                color: palette::BORDER,
                width: 1.0,
                ..Default::default()
            },
            ..Default::default()
        })
        .width(Length::Fill)
    };

    let toolbar = column![toolbar_row1, toolbar_row2].width(Length::Fill);

    // ── Search bar (shown when view.show_search) ─────────────────────────────
    let search_bar: Option<Element<Message>> = if view.show_search {
        let hit_label = if view.pattern_matches.is_empty() {
            if view.search_query.is_empty() {
                String::new()
            } else {
                "0 hits".to_string()
            }
        } else {
            format!(
                "{}/{}",
                view.pattern_match_idx + 1,
                view.pattern_matches.len()
            )
        };

        let placeholder = match view.search_mode {
            SearchMode::Sequence => "IUPAC / amino acid pattern…",
            SearchMode::Name => "sequence name…",
        };

        let submit_msg = match view.search_mode {
            SearchMode::Sequence => Message::SearchPattern(view.search_query.clone()),
            SearchMode::Name => Message::SearchNames(view.search_query.clone()),
        };

        // Mode toggle: "Seq" / "Name"
        let seq_mode_btn = {
            let active = view.search_mode == SearchMode::Sequence;
            button(text("Seq").size(11))
                .padding([2, 6])
                .style(if active {
                    button::primary
                } else {
                    button::secondary
                })
                .on_press(Message::SetSearchMode(SearchMode::Sequence))
        };
        let name_mode_btn = {
            let active = view.search_mode == SearchMode::Name;
            button(text("Name").size(11))
                .padding([2, 6])
                .style(if active {
                    button::primary
                } else {
                    button::secondary
                })
                .on_press(Message::SetSearchMode(SearchMode::Name))
        };

        let prev_btn = button(text("<").size(11))
            .padding([2, 6])
            .on_press(Message::PatternMatchJump(-1));
        let next_btn = button(text(">").size(11))
            .padding([2, 6])
            .on_press(Message::PatternMatchJump(1));

        Some(
            container(
                row![
                    text("Find:").size(12).color(palette::TEXT_DIM),
                    seq_mode_btn,
                    name_mode_btn,
                    text_input(placeholder, &view.search_query)
                        .on_input(Message::SearchInput)
                        .on_submit(submit_msg)
                        .size(12)
                        .width(Length::Fill),
                    text(hit_label).size(11).color(palette::TEXT_DIM),
                    prev_btn,
                    next_btn,
                    button(text("✕").size(11))
                        .padding([2, 6])
                        .on_press(Message::SearchClose),
                ]
                .spacing(6)
                .align_y(Alignment::Center)
                .padding([3, 10])
                .width(Length::Fill),
            )
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(1.0, 0.98, 0.88))),
                border: Border {
                    color: Color::from_rgb(0.90, 0.80, 0.40),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .width(Length::Fill)
            .into(),
        )
    } else {
        None
    };

    // ── Alignment canvas ─────────────────────────────────────────────────────
    let zoom = view.zoom;

    // Version stamp: changes when any visually relevant field changes.
    // Hover position intentionally excluded so mouse-move never forces a redraw.
    let grid_version: u64 = {
        let ptr = Arc::as_ptr(aln) as u64;
        let cptr = Arc::as_ptr(color_table) as u64;
        let sc = view.scroll_col as u64;
        let sr = view.scroll_row as u64;
        let ec = view
            .edit_cell
            .map(|(r, c)| (r as u64) << 32 | c as u64)
            .unwrap_or(u64::MAX);
        ptr.wrapping_add(sc.wrapping_mul(0x9e3779b9_7f4a7c15))
            .wrapping_add(sr.wrapping_mul(0x6c62272e_07bb0142))
            .wrapping_add(zoom.to_bits() as u64)
            .wrapping_add(view.selected.len() as u64 * 0x517cc1b7_27220a95)
            .wrapping_add(view.selected_cols.len() as u64 * 0x43499fdf_f23952dd)
            .wrapping_add(view.pattern_matches.len() as u64 * 0x0b5a4bcf_b9b4e4a3)
            .wrapping_add(ec)
            .wrapping_add(view.show_features as u64)
            .wrapping_add(cptr)
    };

    let grid = AlignmentGrid {
        aln: Arc::clone(aln),
        scroll_col: view.scroll_col,
        scroll_row: view.scroll_row,
        selected: view.selected.clone(),
        scheme: view.color_scheme.clone(),
        ctrl_held: view.ctrl_held,
        shift_held: view.shift_held,
        edit_cell: view.edit_cell,
        cell_w: CELL_W * zoom,
        cell_h: CELL_H * zoom,
        font_size: FONT_SIZE * zoom,
        selected_cols: view.selected_cols.iter().copied().collect(),
        show_features: view.show_features,
        pattern_matches: view.pattern_matches.clone(),
        pattern_match_idx: view.pattern_match_idx,
        version: grid_version,
        color_table: Arc::clone(color_table),
    };

    let grid_canvas = canvas(grid).width(Length::Fill).height(Length::Fill);

    // ── Covariation pairing arcs strip (optional) ────────────────────────────
    // Four states: off | no MI (run hint) | stale MI (re-run hint) | show arcs.
    let pairing_strip: Option<Element<Message>> = if show_pairing_arcs && !mi_available {
        drop(pairing_pairs);
        Some(
            container(
                text("Run MI (Analysis → Mutual Information) to show covariation arcs")
                    .size(10)
                    .color(Color::from_rgb(0.55, 0.45, 0.65)),
            )
            .width(Length::Fill)
            .height(ARCS_H)
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(0.97, 0.97, 1.0))),
                border: iced::Border {
                    color: Color::from_rgb(0.80, 0.78, 0.90),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .align_y(iced::alignment::Vertical::Center)
            .padding([0, 8])
            .into(),
        )
    } else if show_pairing_arcs && mi_stale {
        drop(pairing_pairs);
        Some(
            container(
                text(
                    "⚠ MI results are stale — alignment changed since last run. Re-run to update.",
                )
                .size(10)
                .color(Color::from_rgb(0.55, 0.38, 0.05)),
            )
            .width(Length::Fill)
            .height(ARCS_H)
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(1.0, 0.97, 0.88))),
                border: iced::Border {
                    color: Color::from_rgb(0.80, 0.65, 0.30),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .align_y(iced::alignment::Vertical::Center)
            .padding([0, 8])
            .into(),
        )
    } else if show_pairing_arcs && !pairing_pairs.is_empty() {
        Some(
            canvas(PairingArcsStrip {
                pairs: pairing_pairs,
                scroll_col: view.scroll_col,
                cell_w: CELL_W * zoom,
                max_score: pairing_max,
            })
            .width(Length::Fill)
            .height(ARCS_H)
            .into(),
        )
    } else {
        drop(pairing_pairs);
        None
    };

    // ── Entropy ruler ────────────────────────────────────────────────────────
    let entropy_strip = canvas(EntropyStrip {
        aln: Arc::clone(aln),
        scroll_col: view.scroll_col,
        cell_w: CELL_W * zoom,
    })
    .width(Length::Fill)
    .height(STRIP_H);

    // ── Restriction map strip (optional) ─────────────────────────────────────
    let restr_seq_idx = view.selected.iter().next().copied().unwrap_or(0);
    let restr_strip: Option<iced::widget::Canvas<RestrMapStrip, Message>> = if show_restr_map {
        Some(
            canvas(RestrMapStrip {
                aln: Arc::clone(aln),
                scroll_col: view.scroll_col,
                cell_w: CELL_W * zoom,
                seq_idx: restr_seq_idx,
            })
            .width(Length::Fill)
            .height(RESTR_H),
        )
    } else {
        None
    };

    // ── Horizontal scrollbar ─────────────────────────────────────────────────
    let total_cols = aln.col_count();
    // Approx visible columns at current zoom (~100 at 100% zoom, matches rest of codebase).
    let approx_vis = ((100.0 / view.zoom).round() as usize)
        .max(10)
        .min(total_cols);
    let max_col = total_cols.saturating_sub(1) as u32;
    // Cap scrolling so at least approx_vis columns remain visible.
    let scroll_max = total_cols.saturating_sub(approx_vis) as u32;
    let cur_col = view.scroll_col.min(max_col as usize) as u32;

    let h_scroll: Element<Message> = if max_col > 0 {
        use crate::theme::buttons as bs;
        let total = total_cols;
        // Proportion of total that's currently visible.
        let thumb_w_pct = (approx_vis as f32 / total as f32).min(1.0);
        let thumb_off_pct = view.scroll_col as f32 / (scroll_max.max(1)) as f32;

        // ── Scroll buttons ──
        let btn_start = button(text("|<").size(11))
            .padding([3, 7])
            .style(bs::ghost)
            .on_press(Message::ScrollColSet(0));
        let btn_back = button(text("<").size(12))
            .padding([3, 7])
            .style(bs::ghost)
            .on_press(Message::ScrollColSet(cur_col.saturating_sub(20)));
        let btn_fwd = button(text(">").size(12))
            .padding([3, 7])
            .style(bs::ghost)
            .on_press(Message::ScrollColSet((cur_col + 20).min(scroll_max)));
        let btn_end = button(text(">|").size(11))
            .padding([3, 7])
            .style(bs::ghost)
            .on_press(Message::ScrollColSet(scroll_max));

        // ── Inline progress bar ──
        // Two side-by-side containers whose widths represent the scrolled position.
        let track = container(
            row![
                container(iced::widget::Space::with_height(0))
                    .width(Length::FillPortion((thumb_off_pct * 1000.0) as u16 + 1))
                    .height(Length::Fixed(10.0))
                    .style(|_| container::Style {
                        background: Some(Background::Color(palette::BORDER)),
                        ..Default::default()
                    }),
                container(iced::widget::Space::with_height(0))
                    .width(Length::FillPortion((thumb_w_pct * 1000.0) as u16 + 1))
                    .height(Length::Fixed(10.0))
                    .style(|_| container::Style {
                        background: Some(Background::Color(palette::ANALYSIS)),
                        border: Border {
                            radius: 3.0.into(),
                            ..Default::default()
                        },
                        ..Default::default()
                    }),
                container(iced::widget::Space::with_height(0))
                    .width(Length::Fill)
                    .height(Length::Fixed(10.0))
                    .style(|_| container::Style {
                        background: Some(Background::Color(palette::BORDER)),
                        ..Default::default()
                    }),
            ]
            .width(Length::Fill),
        )
        .width(Length::Fill)
        .height(Length::Fixed(10.0))
        .style(|_| container::Style {
            background: Some(Background::Color(palette::BORDER)),
            border: Border {
                radius: 3.0.into(),
                color: palette::BORDER_STRONG,
                width: 1.0,
            },
            ..Default::default()
        });

        container(
            row![
                btn_start,
                btn_back,
                track,
                text(format!(" col {} / {} ", view.scroll_col + 1, total))
                    .size(12)
                    .color(palette::TEXT),
                btn_fwd,
                btn_end,
            ]
            .spacing(6)
            .align_y(Alignment::Center)
            .padding([5, 10])
            .width(Length::Fill),
        )
        .style(|_| container::Style {
            background: Some(Background::Color(palette::BG_PANEL)),
            border: Border {
                color: palette::BORDER_STRONG,
                width: 1.0,
                radius: 0.0.into(),
            },
            ..Default::default()
        })
        .width(Length::Fill)
        .into()
    } else {
        iced::widget::Space::with_height(0).into()
    };

    // ── Status bar ───────────────────────────────────────────────────────────
    let selected_info = if view.selected.is_empty() {
        String::new()
    } else {
        format!("  ·  {} selected", view.selected.len())
    };

    let hover_info = match (view.hover_col, view.hover_row) {
        (Some(c), Some(r)) => {
            let residue_ch = aln
                .sequences
                .get(r)
                .and_then(|s| s.residues.get(c))
                .map(|&b| format!(" ({})", b as char))
                .unwrap_or_default();
            // Skip per-column frequency computation for large alignments — iterating
            // every sequence on every mouse-move event causes severe UI stutter.
            let freq_info = if aln.seq_count() <= 500 {
                let mut counts = [0u32; 256];
                let mut total = 0u32;
                for seq in &aln.sequences {
                    if let Some(&b) = seq.residues.get(c) {
                        let b = b.to_ascii_uppercase();
                        if !matches!(b, b'-' | b'~' | b'.') {
                            counts[b as usize] += 1;
                            total += 1;
                        }
                    }
                }
                if total > 0 {
                    let mut sorted: Vec<(u8, u32)> = counts
                        .iter()
                        .enumerate()
                        .filter(|(_, &n)| n > 0)
                        .map(|(i, &n)| (i as u8, n))
                        .collect();
                    sorted.sort_unstable_by(|a, b| b.1.cmp(&a.1));
                    let top: String = sorted
                        .iter()
                        .take(3)
                        .map(|(b, n)| {
                            format!("{}:{:.0}%", *b as char, *n as f32 / total as f32 * 100.0)
                        })
                        .collect::<Vec<_>>()
                        .join(" ");
                    format!("  [{}]", top)
                } else {
                    String::new()
                }
            } else {
                String::new()
            };
            format!(
                "  ·  col {} row {}{}{}",
                c + 1,
                r + 1,
                residue_ch,
                freq_info
            )
        }
        _ => String::new(),
    };

    // Show current view position as fractions
    let pos_info = format!(
        "  rows {}-{}/{} cols {}-{}/{}",
        view.scroll_row + 1,
        (view.scroll_row + 30).min(aln.seq_count()),
        aln.seq_count(),
        view.scroll_col + 1,
        (view.scroll_col + 100).min(aln.col_count()),
        aln.col_count(),
    );

    let edit_hint = "  ·  ^F find  ^Z undo  ^S save";

    let status_str = if let Some(c) = computing {
        format!("⏳ {c}")
    } else {
        format!("{status}{pos_info}{selected_info}{hover_info}{edit_hint}")
    };
    let status_bar = container(text(status_str).size(12).color(palette::TEXT_DIM))
        .padding([3, 10])
        .width(Length::Fill)
        .style(|_| container::Style {
            background: Some(Background::Color(palette::BG_PANEL)),
            border: Border {
                color: palette::BORDER,
                width: 1.0,
                radius: 0.0.into(),
            },
            ..Default::default()
        });

    // ── Protein view banner ───────────────────────────────────────────────────
    let protein_banner: Option<Element<Message>> = if let Some(frame) = view.protein_view_frame {
        let label: i8 = if frame < 3 {
            (frame as i8) + 1
        } else {
            -((frame as i8) - 2)
        };
        let frame_btns = row![
            button(text("F+1").size(10))
                .padding([1, 5])
                .on_press(Message::SetProteinViewFrame(0)),
            button(text("F+2").size(10))
                .padding([1, 5])
                .on_press(Message::SetProteinViewFrame(1)),
            button(text("F+3").size(10))
                .padding([1, 5])
                .on_press(Message::SetProteinViewFrame(2)),
            button(text("F-1").size(10))
                .padding([1, 5])
                .on_press(Message::SetProteinViewFrame(3)),
            button(text("F-2").size(10))
                .padding([1, 5])
                .on_press(Message::SetProteinViewFrame(4)),
            button(text("F-3").size(10))
                .padding([1, 5])
                .on_press(Message::SetProteinViewFrame(5)),
        ]
        .spacing(2);
        Some(
            container(
                row![
                    text(format!("Protein View — Frame {:+}  (read-only)", label))
                        .size(12)
                        .color(palette::TEXT),
                    iced::widget::horizontal_space(),
                    frame_btns,
                    button(text("< DNA View").size(12))
                        .padding([2, 8])
                        .on_press(Message::ExitProteinView),
                ]
                .spacing(8)
                .align_y(Alignment::Center)
                .padding([3, 10])
                .width(Length::Fill),
            )
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(1.0, 0.95, 0.80))),
                border: Border {
                    color: Color::from_rgb(0.90, 0.70, 0.20),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .width(Length::Fill)
            .into(),
        )
    } else {
        None
    };

    // ── BLAST submit panel (shown when view.show_blast_panel) ─────────────────
    let blast_panel: Option<Element<Message>> = if view.show_blast_panel {
        Some(
            container(
                row![
                    text("BLAST:").size(12).color(palette::TEXT),
                    text("Program:").size(11).color(palette::TEXT_DIM),
                    button(
                        text(if view.blast_program == "blastn" {
                            "blastn ●"
                        } else {
                            "blastn"
                        })
                        .size(11)
                    )
                    .padding([2, 6])
                    .on_press(Message::BlastProgramChanged("blastn".into())),
                    button(
                        text(if view.blast_program == "blastp" {
                            "blastp ●"
                        } else {
                            "blastp"
                        })
                        .size(11)
                    )
                    .padding([2, 6])
                    .on_press(Message::BlastProgramChanged("blastp".into())),
                    button(
                        text(if view.blast_program == "blastx" {
                            "blastx ●"
                        } else {
                            "blastx"
                        })
                        .size(11)
                    )
                    .padding([2, 6])
                    .on_press(Message::BlastProgramChanged("blastx".into())),
                    text("DB:").size(11).color(palette::TEXT_DIM),
                    button(
                        text(if view.blast_database == "nt" {
                            "nt ●"
                        } else {
                            "nt"
                        })
                        .size(11)
                    )
                    .padding([2, 6])
                    .on_press(Message::BlastDatabaseChanged("nt".into())),
                    button(
                        text(if view.blast_database == "nr" {
                            "nr ●"
                        } else {
                            "nr"
                        })
                        .size(11)
                    )
                    .padding([2, 6])
                    .on_press(Message::BlastDatabaseChanged("nr".into())),
                    button(
                        text(if view.blast_database == "refseq_rna" {
                            "refseq_rna ●"
                        } else {
                            "refseq_rna"
                        })
                        .size(11)
                    )
                    .padding([2, 6])
                    .on_press(Message::BlastDatabaseChanged("refseq_rna".into())),
                    iced::widget::horizontal_space(),
                    button(text("Submit selected sequence").size(12))
                        .padding([2, 10])
                        .on_press(Message::SubmitBlast),
                    button(text("✕").size(11))
                        .padding([2, 6])
                        .on_press(Message::ToggleBlastPanel),
                ]
                .spacing(4)
                .align_y(Alignment::Center)
                .padding([3, 10])
                .width(Length::Fill),
            )
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(0.95, 0.88, 1.0))),
                border: Border {
                    color: Color::from_rgb(0.70, 0.40, 0.90),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .width(Length::Fill)
            .into(),
        )
    } else {
        None
    };

    // ── Fetch bar (shown when view.show_fetch_bar) ────────────────────────────
    let fetch_bar: Option<Element<Message>> = if view.show_fetch_bar {
        Some(
            container(
                row![
                    text("Fetch:").size(12).color(palette::TEXT),
                    text_input(
                        "Accession — NCBI (NM_001256799) or UniProt (P02945), comma-separated…",
                        &view.fetch_input
                    )
                    .on_input(Message::FetchInput)
                    .on_submit(Message::FetchAccession)
                    .size(12)
                    .width(Length::Fill),
                    button(text("Fetch").size(12))
                        .padding([2, 10])
                        .on_press(Message::FetchAccession),
                    button(text("✕").size(11))
                        .padding([2, 6])
                        .on_press(Message::ToggleFetchBar),
                ]
                .spacing(6)
                .align_y(Alignment::Center)
                .padding([3, 10])
                .width(Length::Fill),
            )
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(0.88, 0.95, 1.0))),
                border: Border {
                    color: Color::from_rgb(0.40, 0.70, 0.90),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .width(Length::Fill)
            .into(),
        )
    } else {
        None
    };

    // ── Annotate bar (shown when view.show_annotate_bar) ─────────────────────
    let annotate_bar: Option<Element<Message>> = if view.show_annotate_bar {
        let col_range = if view.selected_cols.is_empty() {
            "no columns selected".to_string()
        } else {
            let lo = view.selected_cols.iter().next().unwrap() + 1;
            let hi = view.selected_cols.iter().next_back().unwrap() + 1;
            format!("cols {lo}–{hi}")
        };
        Some(
            container(
                row![
                    text("Annotate:").size(12).color(palette::TEXT),
                    text_input("Feature name…", &view.annotate_name)
                        .on_input(Message::AnnotateName)
                        .on_submit(Message::ApplyAnnotation)
                        .size(12)
                        .width(Length::Fill),
                    text("Type:").size(12).color(palette::TEXT_DIM),
                    text_input("misc_feature", &view.annotate_type)
                        .on_input(Message::AnnotateType)
                        .size(12)
                        .width(200),
                    text(col_range).size(11).color(palette::TEXT_DIM),
                    button(text("Apply").size(12))
                        .padding([2, 10])
                        .on_press(Message::ApplyAnnotation),
                    button(text("✕").size(11))
                        .padding([2, 6])
                        .on_press(Message::ToggleAnnotateBar),
                ]
                .spacing(6)
                .align_y(Alignment::Center)
                .padding([3, 10])
                .width(Length::Fill),
            )
            .style(|_| container::Style {
                background: Some(Background::Color(Color::from_rgb(0.92, 1.0, 0.90))),
                border: Border {
                    color: Color::from_rgb(0.40, 0.80, 0.40),
                    width: 1.0,
                    radius: 0.0.into(),
                },
                ..Default::default()
            })
            .width(Length::Fill)
            .into(),
        )
    } else {
        None
    };

    // ── Build the grid column (search bar + grid + strips) ───────────────────
    let grid_col: Element<Message> = {
        let mut parts: Vec<Element<Message>> = Vec::new();
        if let Some(pb) = protein_banner {
            parts.push(pb);
        }
        if let Some(bp) = blast_panel {
            parts.push(bp);
        }
        if let Some(sb) = search_bar {
            parts.push(sb);
        }
        if let Some(fb) = fetch_bar {
            parts.push(fb);
        }
        if let Some(ab) = annotate_bar {
            parts.push(ab);
        }
        parts.push(grid_canvas.into());
        if let Some(rs) = restr_strip {
            parts.push(rs.into());
        }
        if let Some(ps) = pairing_strip {
            parts.push(ps.into());
        }
        parts.push(entropy_strip.into());
        column(parts)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    };

    // ── Main content area (grid + optional analysis panel) ───────────────────
    let main_area: Element<Message> = if view.show_analysis {
        let panel = analysis_panel(aln);
        row![grid_col, panel]
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    } else {
        grid_col
    };

    // ── Title context menu (shown when right-clicking a sequence title) ──────
    let ctx_menu: Option<Element<Message>> = title_menu.and_then(|idx| {
        let seq = aln.sequences.get(idx)?;
        let lock_label = if seq.locked { "Unlock" } else { "Lock" };
        let name_for_rename = rename_input.to_owned();
        let panel = container(
            column![
                row![
                    text(format!("Sequence: {}", seq.name))
                        .size(12)
                        .color(palette::TEXT),
                    horizontal_space(),
                    button(text("✕").size(11))
                        .padding([2, 6])
                        .on_press(Message::TitleMenuClose),
                ]
                .spacing(6)
                .align_y(Alignment::Center),
                row![
                    text("Rename:").size(11).color(palette::TEXT_DIM),
                    text_input("new name…", rename_input)
                        .on_input(Message::RenameInput)
                        .on_submit(Message::RenameSequence {
                            idx,
                            name: name_for_rename
                        })
                        .size(11)
                        .width(Length::Fill),
                    button(text("OK").size(11))
                        .padding([2, 8])
                        .on_press(Message::RenameSequence {
                            idx,
                            name: rename_input.to_owned()
                        }),
                ]
                .spacing(6)
                .align_y(Alignment::Center),
                row![
                    button(text(lock_label).size(11))
                        .padding([3, 8])
                        .on_press(Message::ToggleLockSequence { idx }),
                    button(text("Duplicate").size(11))
                        .padding([3, 8])
                        .on_press(Message::DuplicateSequence { idx }),
                    button(text("Delete").size(11))
                        .padding([3, 8])
                        .on_press(Message::DeleteSequence { idx }),
                    button(text("Run ORFs").size(11))
                        .padding([3, 8])
                        .on_press(Message::RunOrfFinder),
                    button(text("Rev-Comp").size(11))
                        .padding([3, 8])
                        .on_press(Message::ReverseComplementSelected),
                ]
                .spacing(4),
            ]
            .spacing(6)
            .padding([6, 10])
            .width(Length::Fill),
        )
        .style(|_| container::Style {
            background: Some(Background::Color(Color::from_rgb(0.98, 0.98, 0.96))),
            border: Border {
                color: Color::from_rgb(0.70, 0.70, 0.50),
                width: 1.0,
                radius: 0.0.into(),
            },
            ..Default::default()
        })
        .width(Length::Fill);

        Some(panel.into())
    });

    // ── Layout ───────────────────────────────────────────────────────────────
    let mut layout_parts: Vec<Element<Message>> = vec![toolbar.into()];
    if let Some(menu) = ctx_menu {
        layout_parts.push(menu);
    }
    layout_parts.push(main_area);
    layout_parts.push(h_scroll);
    layout_parts.push(status_bar.into());

    column(layout_parts)
        .width(Length::Fill)
        .height(Length::Fill)
        .into()
}

// ── Helper: color scheme toggle button ───────────────────────────────────────

fn scheme_btn<'a>(
    label: &'a str,
    scheme: ColorScheme,
    current: &ColorScheme,
) -> Element<'a, Message> {
    use crate::theme::buttons as bs;
    let active = scheme == *current;
    let msg = Message::SetColorScheme(scheme);
    button(text(label).size(12))
        .padding([3, 9])
        .on_press(msg)
        .style(if active {
            bs::toggle_on
        } else {
            bs::toggle_off
        })
        .into()
}

/// Thin vertical separator for use inside toolbar rows.
fn vsep<'a>() -> Element<'a, Message> {
    container(iced::widget::vertical_space())
        .height(Length::Fixed(16.0))
        .width(Length::Fixed(1.0))
        .style(|_| container::Style {
            background: Some(Background::Color(Color::from_rgb(0.75, 0.76, 0.80))),
            ..Default::default()
        })
        .into()
}
