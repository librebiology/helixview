//! Project file (.hvproj) — saves and restores all open tabs.
//!
//! Format: JSON (via serde_json) containing a list of tab snapshots.
//! Each tab stores the serialized Alignment, view settings, and flags.
//! On restore every tab is reopened in a new slot; the active tab index is preserved.

use std::collections::BTreeSet;
use std::path::Path;
use std::sync::Arc;

use serde::{Deserialize, Serialize};
use helixview_core::Alignment as SeqAlignment;

use helixview_core::History as EditHistory;
use crate::app::{ColorScheme, HelixViewApp, TabRecord, ToolbarTab, ViewState};

// ── Serializable view state ───────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SavedViewState {
    scroll_col:         usize,
    scroll_row:         usize,
    selected:           Vec<usize>,
    color_scheme:       String,
    zoom:               f32,
    show_features:      bool,
    show_search:        bool,
    search_query:       String,
    show_blast_panel:   bool,
    blast_program:      String,
    blast_database:     String,
    toolbar_tab:        String,
}

impl From<&ViewState> for SavedViewState {
    fn from(v: &ViewState) -> Self {
        Self {
            scroll_col:       v.scroll_col,
            scroll_row:       v.scroll_row,
            selected:         v.selected.iter().copied().collect(),
            color_scheme:     match v.color_scheme {
                ColorScheme::ResidueType => "ResidueType",
                ColorScheme::Identity    => "Identity",
                ColorScheme::Strength    => "Strength",
                ColorScheme::Plain       => "Plain",
            }.to_string(),
            zoom:             v.zoom,
            show_features:    v.show_features,
            show_search:      v.show_search,
            search_query:     v.search_query.clone(),
            show_blast_panel: v.show_blast_panel,
            blast_program:    v.blast_program.clone(),
            blast_database:   v.blast_database.clone(),
            toolbar_tab:      match v.toolbar_tab {
                ToolbarTab::Edit     => "Edit",
                ToolbarTab::Analysis => "Analysis",
                ToolbarTab::View     => "View",
            }.to_string(),
        }
    }
}

impl From<SavedViewState> for ViewState {
    fn from(s: SavedViewState) -> Self {
        let mut v = ViewState::default();
        v.scroll_col  = s.scroll_col;
        v.scroll_row  = s.scroll_row;
        v.selected    = s.selected.into_iter().collect();
        v.color_scheme = match s.color_scheme.as_str() {
            "Identity" => ColorScheme::Identity,
            "Strength" => ColorScheme::Strength,
            "Plain"    => ColorScheme::Plain,
            _          => ColorScheme::ResidueType,
        };
        v.zoom            = s.zoom;
        v.show_features   = s.show_features;
        v.show_search     = s.show_search;
        v.search_query    = s.search_query;
        v.show_blast_panel = s.show_blast_panel;
        v.blast_program   = s.blast_program;
        v.blast_database  = s.blast_database;
        v.toolbar_tab     = match s.toolbar_tab.as_str() {
            "Analysis" => ToolbarTab::Analysis,
            "View"     => ToolbarTab::View,
            _          => ToolbarTab::Edit,
        };
        v
    }
}

// ── Per-tab snapshot ──────────────────────────────────────────────────────────

#[derive(Debug, Serialize, Deserialize)]
struct SavedTab {
    document:           Option<SeqAlignment>,
    view:               SavedViewState,
    show_restr_map:     bool,
    protein_view_frame: usize,
}

// ── Project file ──────────────────────────────────────────────────────────────

#[derive(Debug, Serialize, Deserialize)]
struct ProjectFile {
    version:    u32,
    active_tab: usize,
    tabs:       Vec<SavedTab>,
}

// ── Public API ────────────────────────────────────────────────────────────────

pub fn save(app: &HelixViewApp, path: &Path) -> Result<(), String> {
    // Build a snapshot of all tabs including the live active tab.
    let mut tabs: Vec<SavedTab> = Vec::new();

    for (i, rec) in app.tabs.iter().enumerate() {
        if i == app.active_tab {
            // Active tab: read from live fields, not the stale snapshot.
            tabs.push(SavedTab {
                document:           app.document.as_deref().cloned(),
                view:               SavedViewState::from(&app.view),
                show_restr_map:     app.show_restr_map,
                protein_view_frame: app.protein_view_frame,
            });
        } else {
            tabs.push(SavedTab {
                document:           rec.document.as_deref().cloned(),
                view:               SavedViewState::from(&rec.view),
                show_restr_map:     rec.show_restr_map,
                protein_view_frame: rec.protein_view_frame,
            });
        }
    }

    let proj = ProjectFile {
        version:    1,
        active_tab: app.active_tab,
        tabs,
    };

    let json = serde_json::to_string_pretty(&proj).map_err(|e| e.to_string())?;
    std::fs::write(path, json).map_err(|e| e.to_string())
}

/// Serialise the current session to a JSON string without writing to disk.
/// Used by `auto_save_session` to separate serialisation (main thread) from
/// the file write (background thread).
pub fn to_json(app: &HelixViewApp) -> Result<String, String> {
    let mut tabs: Vec<SavedTab> = Vec::new();
    for (i, rec) in app.tabs.iter().enumerate() {
        if i == app.active_tab {
            tabs.push(SavedTab {
                document:           app.document.as_deref().cloned(),
                view:               SavedViewState::from(&app.view),
                show_restr_map:     app.show_restr_map,
                protein_view_frame: app.protein_view_frame,
            });
        } else {
            tabs.push(SavedTab {
                document:           rec.document.as_deref().cloned(),
                view:               SavedViewState::from(&rec.view),
                show_restr_map:     rec.show_restr_map,
                protein_view_frame: rec.protein_view_frame,
            });
        }
    }
    let proj = ProjectFile { version: 1, active_tab: app.active_tab, tabs };
    serde_json::to_string_pretty(&proj).map_err(|e| e.to_string())
}

/// Restore a saved project into `app`, replacing all current tabs.
pub fn load(app: &mut HelixViewApp, path: &Path) -> Result<(), String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let proj: ProjectFile = serde_json::from_str(&text).map_err(|e| e.to_string())?;

    if proj.tabs.is_empty() {
        return Err("Project file contains no tabs".into());
    }

    // Build TabRecords for every tab except the active one.
    let active = proj.active_tab.min(proj.tabs.len() - 1);
    let mut new_tabs: Vec<TabRecord> = proj.tabs.iter().map(|t| TabRecord {
        document:           t.document.clone().map(Arc::new),
        view:               ViewState::from(t.view.clone()),
        history:            EditHistory::new(app.prefs.undo_depth),
        show_restr_map:     t.show_restr_map,
        protein_view_doc:   None,
        protein_view_frame: t.protein_view_frame,
    }).collect();

    // Pop the active slot and apply it to live fields.
    let live = new_tabs.remove(active);
    // Re-insert a placeholder so indices remain correct.
    new_tabs.insert(active, TabRecord {
        document:           live.document.clone(),
        view:               live.view.clone(),
        history:            EditHistory::new(app.prefs.undo_depth),
        show_restr_map:     live.show_restr_map,
        protein_view_doc:   None,
        protein_view_frame: live.protein_view_frame,
    });

    app.tabs       = new_tabs;
    app.active_tab = active;
    app.document   = live.document;
    app.view       = live.view;
    app.history    = EditHistory::new(app.prefs.undo_depth);
    app.show_restr_map     = live.show_restr_map;
    app.protein_view_frame = live.protein_view_frame;
    app.protein_view_doc   = None;

    // Clear transient state.
    app.pairwise_result  = None;
    app.orf_results      = None;
    app.six_frame_result = None;
    app.blast            = crate::app::BlastState::Idle;
    app.trace_data       = None;
    app.tree_data        = None;
    app.tree_selected.clear();

    Ok(())
}
