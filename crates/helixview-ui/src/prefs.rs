//! Persistent user preferences, stored at:
//!   Linux/macOS: ~/.config/helixview/preferences.toml
//!   Windows:     %APPDATA%\helixview\preferences.toml

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

// ── Preferences struct ────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Preferences {
    /// Default zoom level for new documents (0.4–3.0, default 1.0).
    pub default_zoom: f32,
    /// Undo history depth.
    pub undo_depth: usize,
    /// Show feature stripes by default.
    pub show_features: bool,
    /// Show restriction map strip by default.
    pub show_restr_map: bool,
    /// Default color scheme name ("ResidueType", "Identity", "Strength", "Plain").
    pub color_scheme: String,
    /// Identity shading threshold (0.0–1.0).
    pub identity_threshold: f32,
    /// Font size for the residue grid (8–20).
    pub font_size: f32,
    /// Number of recently opened files to remember.
    pub recent_files_max: usize,
    /// Recently opened file paths (most recent first).
    pub recent_files: Vec<PathBuf>,
}

impl Default for Preferences {
    fn default() -> Self {
        Self {
            default_zoom:       1.0,
            undo_depth:         100,
            show_features:      true,
            show_restr_map:     false,
            color_scheme:       "ResidueType".to_string(),
            identity_threshold: 0.5,
            font_size:          11.0,
            recent_files_max:   10,
            recent_files:       Vec::new(),
        }
    }
}

// ── File path ─────────────────────────────────────────────────────────────────

pub fn prefs_path() -> Option<PathBuf> {
    dirs::config_dir().map(|d| d.join("helixview").join("preferences.toml"))
}

// ── Load / save ───────────────────────────────────────────────────────────────

pub fn load() -> Preferences {
    let path = match prefs_path() {
        Some(p) => p,
        None    => return Preferences::default(),
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t)  => t,
        Err(_) => return Preferences::default(),
    };
    toml::from_str(&text).unwrap_or_default()
}

pub fn save(prefs: &Preferences) -> Result<(), String> {
    let path = prefs_path().ok_or("Cannot determine config directory")?;
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    }
    let text = toml::to_string_pretty(prefs).map_err(|e| e.to_string())?;
    std::fs::write(&path, text).map_err(|e| e.to_string())
}

/// Path to auto-saved last session project file.
pub fn session_path() -> Option<PathBuf> {
    dirs::config_dir().map(|d| d.join("helixview").join("last_session.hvproj"))
}

/// Push a path onto the recent-files list, deduplicating and capping at max.
pub fn push_recent(prefs: &mut Preferences, path: PathBuf) {
    prefs.recent_files.retain(|p| p != &path);
    prefs.recent_files.insert(0, path);
    prefs.recent_files.truncate(prefs.recent_files_max);
}
