/// Shared palette constants and button style functions for HelixView.
/// Light-mode only — biologists expect publication-quality light interfaces.

pub mod palette {
    #![allow(dead_code)]
    use iced::Color;

    // ── Backgrounds ───────────────────────────────────────────────────────────
    pub const BG: Color = Color::from_rgb(0.97, 0.97, 0.99);
    pub const BG_PANEL: Color = Color::from_rgb(0.93, 0.94, 0.97);
    pub const BG_SEQ: Color = Color::WHITE;
    pub const BG_SEQ_ALT: Color = Color::from_rgb(0.96, 0.97, 1.00);
    pub const HEADER_BG: Color = Color::from_rgb(0.91, 0.93, 0.98);
    pub const TOOLBAR_BG: Color = Color::from_rgb(0.26, 0.28, 0.34);

    // ── Text ──────────────────────────────────────────────────────────────────
    pub const TEXT: Color = Color::from_rgb(0.10, 0.10, 0.14);
    pub const TEXT_DIM: Color = Color::from_rgb(0.46, 0.48, 0.55);
    pub const TEXT_ON_DARK: Color = Color::from_rgb(0.92, 0.93, 0.96);

    // ── Borders ───────────────────────────────────────────────────────────────
    pub const BORDER: Color = Color::from_rgb(0.78, 0.80, 0.86);
    pub const BORDER_STRONG: Color = Color::from_rgb(0.60, 0.63, 0.72);

    // ── Semantic accent colours ───────────────────────────────────────────────
    /// Primary brand green (biology / life-science associations).
    pub const ACCENT: Color = Color::from_rgb(0.10, 0.52, 0.28);
    pub const ACCENT_HOVER: Color = Color::from_rgb(0.12, 0.62, 0.33);
    pub const ACCENT_PRESS: Color = Color::from_rgb(0.07, 0.40, 0.21);

    /// Analysis / data-science blue.
    pub const ANALYSIS: Color = Color::from_rgb(0.09, 0.38, 0.72);
    pub const ANALYSIS_HOVER: Color = Color::from_rgb(0.11, 0.46, 0.85);

    /// Warning amber.
    pub const WARN: Color = Color::from_rgb(0.80, 0.50, 0.05);

    /// Destructive red.
    pub const DANGER: Color = Color::from_rgb(0.72, 0.12, 0.12);
    pub const DANGER_HOVER: Color = Color::from_rgb(0.84, 0.16, 0.16);
}

// ── Button style functions ────────────────────────────────────────────────────
// Each function matches the signature  `fn(&iced::Theme, iced::widget::button::Status)`
// so it can be passed directly to `.style(...)`.

pub mod buttons {
    use super::palette;
    use iced::{widget::button, Background, Border, Color, Shadow, Vector};

    fn shadow_sm() -> Shadow {
        Shadow {
            color: Color::from_rgba(0.0, 0.0, 0.0, 0.12),
            offset: Vector::new(0.0, 1.5),
            blur_radius: 3.0,
        }
    }

    // ── Primary (green) — most important call-to-action ───────────────────────
    pub fn primary(_t: &iced::Theme, s: button::Status) -> button::Style {
        let bg = match s {
            button::Status::Hovered => palette::ACCENT_HOVER,
            button::Status::Pressed => palette::ACCENT_PRESS,
            button::Status::Disabled => Color::from_rgba(0.10, 0.52, 0.28, 0.45),
            _ => palette::ACCENT,
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: Color::WHITE,
            border: Border {
                radius: 5.0.into(),
                ..Default::default()
            },
            shadow: shadow_sm(),
        }
    }

    // ── Save — slightly subdued green ─────────────────────────────────────────
    pub fn save(_t: &iced::Theme, s: button::Status) -> button::Style {
        let base = Color::from_rgb(0.20, 0.46, 0.30);
        let bg = match s {
            button::Status::Hovered => Color::from_rgb(0.24, 0.56, 0.36),
            button::Status::Pressed => Color::from_rgb(0.15, 0.37, 0.23),
            _ => base,
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: Color::WHITE,
            border: Border {
                radius: 5.0.into(),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    // ── Analysis (blue) — computation and results ─────────────────────────────
    pub fn analysis(_t: &iced::Theme, s: button::Status) -> button::Style {
        let bg = match s {
            button::Status::Hovered => palette::ANALYSIS_HOVER,
            button::Status::Pressed => Color::from_rgb(0.06, 0.28, 0.58),
            button::Status::Disabled => Color::from_rgba(0.09, 0.38, 0.72, 0.45),
            _ => palette::ANALYSIS,
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: Color::WHITE,
            border: Border {
                radius: 5.0.into(),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    // ── Secondary — default utility actions ───────────────────────────────────
    pub fn secondary(_t: &iced::Theme, s: button::Status) -> button::Style {
        let (bg, text) = match s {
            button::Status::Hovered => (
                Color::from_rgb(0.83, 0.85, 0.93),
                Color::from_rgb(0.08, 0.08, 0.18),
            ),
            button::Status::Pressed => (
                Color::from_rgb(0.74, 0.76, 0.86),
                Color::from_rgb(0.05, 0.05, 0.14),
            ),
            button::Status::Disabled => (
                Color::from_rgba(0.88, 0.90, 0.95, 0.5),
                Color::from_rgba(0.30, 0.30, 0.38, 0.5),
            ),
            _ => (
                Color::from_rgb(0.89, 0.91, 0.96),
                Color::from_rgb(0.16, 0.16, 0.26),
            ),
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: text,
            border: Border {
                radius: 5.0.into(),
                color: Color::from_rgb(0.72, 0.75, 0.84),
                width: 1.0,
            },
            ..Default::default()
        }
    }

    // ── Ghost — minimal, text-only feel for low-priority actions ─────────────
    pub fn ghost(_t: &iced::Theme, s: button::Status) -> button::Style {
        let (bg, text) = match s {
            button::Status::Hovered => (
                Color::from_rgba(0.0, 0.0, 0.0, 0.07),
                Color::from_rgb(0.10, 0.10, 0.20),
            ),
            button::Status::Pressed => (
                Color::from_rgba(0.0, 0.0, 0.0, 0.12),
                Color::from_rgb(0.05, 0.05, 0.14),
            ),
            _ => (Color::TRANSPARENT, Color::from_rgb(0.25, 0.28, 0.42)),
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: text,
            border: Border {
                radius: 4.0.into(),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    // ── Danger (red) ─────────────────────────────────────────────────────────
    pub fn danger(_t: &iced::Theme, s: button::Status) -> button::Style {
        let bg = match s {
            button::Status::Hovered => palette::DANGER_HOVER,
            button::Status::Pressed => Color::from_rgb(0.58, 0.09, 0.09),
            _ => palette::DANGER,
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: Color::WHITE,
            border: Border {
                radius: 5.0.into(),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    // ── Tab (active) ─────────────────────────────────────────────────────────
    pub fn tab_active(_t: &iced::Theme, s: button::Status) -> button::Style {
        let bg = match s {
            button::Status::Hovered => palette::ACCENT_HOVER,
            _ => palette::ACCENT,
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: Color::WHITE,
            border: Border {
                radius: 5.0.into(),
                ..Default::default()
            },
            shadow: Shadow {
                color: Color::from_rgba(0.0, 0.0, 0.0, 0.15),
                offset: Vector::new(0.0, 1.5),
                blur_radius: 3.0,
            },
        }
    }

    // ── Tab (inactive) ────────────────────────────────────────────────────────
    pub fn tab_inactive(_t: &iced::Theme, s: button::Status) -> button::Style {
        let (bg, text) = match s {
            button::Status::Hovered => (
                Color::from_rgb(0.82, 0.84, 0.92),
                Color::from_rgb(0.08, 0.08, 0.22),
            ),
            button::Status::Pressed => (
                Color::from_rgb(0.76, 0.78, 0.88),
                Color::from_rgb(0.05, 0.05, 0.18),
            ),
            _ => (
                Color::from_rgb(0.90, 0.91, 0.96),
                Color::from_rgb(0.26, 0.28, 0.40),
            ),
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: text,
            border: Border {
                radius: 5.0.into(),
                color: Color::from_rgb(0.74, 0.76, 0.86),
                width: 1.0,
            },
            ..Default::default()
        }
    }

    // ── Toggle (on/off indicator) — for Feat●/RE● type buttons ───────────────
    pub fn toggle_on(_t: &iced::Theme, s: button::Status) -> button::Style {
        let bg = match s {
            button::Status::Hovered => Color::from_rgb(0.10, 0.44, 0.55),
            button::Status::Pressed => Color::from_rgb(0.07, 0.34, 0.43),
            _ => Color::from_rgb(0.09, 0.40, 0.50),
        };
        button::Style {
            background: Some(Background::Color(bg)),
            text_color: Color::WHITE,
            border: Border {
                radius: 5.0.into(),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    pub fn toggle_off(_t: &iced::Theme, s: button::Status) -> button::Style {
        secondary(_t, s)
    }
}
