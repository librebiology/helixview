//! Accessory app persistence — load/save from ~/.config/helixview/accessories.toml

use crate::app::AccessoryApp;
use std::path::PathBuf;

fn accessories_path() -> Option<PathBuf> {
    dirs::config_dir().map(|d| d.join("helixview").join("accessories.toml"))
}

pub fn load_accessories() -> Vec<AccessoryApp> {
    let path = match accessories_path() {
        Some(p) => p,
        None => return Vec::new(),
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(_) => return Vec::new(),
    };
    #[derive(serde::Deserialize)]
    struct File {
        accessories: Vec<AccessoryApp>,
    }
    toml::from_str::<File>(&text)
        .map(|f| f.accessories)
        .unwrap_or_default()
}

pub fn save_accessories(list: &[AccessoryApp]) -> Result<(), String> {
    let path = accessories_path().ok_or("Cannot determine config directory")?;
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    }
    #[derive(serde::Serialize)]
    struct File<'a> {
        accessories: &'a [AccessoryApp],
    }
    let text = toml::to_string_pretty(&File { accessories: list }).map_err(|e| e.to_string())?;
    std::fs::write(&path, text).map_err(|e| e.to_string())
}
