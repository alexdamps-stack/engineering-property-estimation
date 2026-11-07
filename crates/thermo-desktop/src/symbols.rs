//! Technical-drawing look: CAD colour scheme and process-flow-diagram symbols.
//!
//! Symbols are monochrome line work modelled on ISO 10628-2 graphical
//! symbols (vessel with dished heads, heat exchanger with heating/cooling
//! medium, shell-and-tube exchanger, control valve with diaphragm actuator,
//! centrifugal pump, turbo compressor and expander, static mixer). Arrowheads
//! are closed and filled with a 15° half angle as in ISO 128. The dark
//! scheme mirrors a CAD model space; the light scheme a plotted sheet.

use std::f32::consts::PI;

use eframe::egui::{Color32, Painter, Pos2, Rect, Shape, Stroke, Vec2, Visuals, pos2, vec2};
use flowsheet_core::UnitKind;

/// Colours used on the canvas, derived from the active light/dark visuals.
#[derive(Debug, Clone, Copy)]
pub struct Theme {
    /// Drawing background (model space or paper).
    pub background: Color32,
    /// Minor grid lines.
    pub grid_minor: Color32,
    /// Major grid lines (every fifth minor line).
    pub grid_major: Color32,
    /// Line work and lettering.
    pub line: Color32,
    /// De-emphasised line work (underspecified units, hints).
    pub muted: Color32,
    /// Selection highlight and grips.
    pub select: Color32,
    /// Snap marker on free nozzles.
    pub snap: Color32,
    /// Warnings in the surrounding panels.
    pub warning: Color32,
}

impl Theme {
    /// Theme for the given egui visuals.
    pub fn from_visuals(visuals: &Visuals) -> Self {
        if visuals.dark_mode {
            Self {
                background: Color32::from_rgb(33, 40, 48),
                grid_minor: Color32::from_rgb(40, 48, 57),
                grid_major: Color32::from_rgb(54, 64, 76),
                line: Color32::from_rgb(240, 240, 240),
                muted: Color32::from_rgb(128, 138, 150),
                select: Color32::from_rgb(30, 144, 255),
                snap: Color32::from_rgb(255, 200, 0),
                warning: Color32::from_rgb(210, 153, 34),
            }
        } else {
            Self {
                background: Color32::WHITE,
                grid_minor: Color32::from_rgb(240, 242, 245),
                grid_major: Color32::from_rgb(222, 226, 232),
                line: Color32::BLACK,
                muted: Color32::from_rgb(140, 140, 140),
                select: Color32::from_rgb(0, 102, 214),
                snap: Color32::from_rgb(214, 130, 0),
                warning: Color32::from_rgb(154, 103, 0),
            }
        }
    }
}

/// Closed, filled arrowhead with a 15° half angle (ISO 128), tip at `tip`
/// pointing along `dir`.
pub fn arrowhead(painter: &Painter, tip: Pos2, dir: Vec2, length: f32, color: Color32) {
    let d = dir.normalized();
    let n = d.rot90() * length * (15.0_f32).to_radians().tan();
    let base = tip - d * length;
    painter.add(Shape::convex_polygon(
        vec![tip, base + n, base - n],
        color,
        Stroke::NONE,
    ));
}

/// Draws the symbol of `kind` in the screen rectangle `r`. `fill` masks
/// grid and lines behind the symbol (use transparent for palette icons).
pub fn draw_unit(painter: &Painter, kind: UnitKind, r: Rect, stroke: Stroke, fill: Color32) {
    let c = r.center();
    match kind {
        UnitKind::Feed | UnitKind::Product => {
            // Flow boundary flag pointing in the flow direction.
            let tip = r.height() * 0.5;
            let pts = vec![
                r.left_top(),
                pos2(r.right() - tip, r.top()),
                pos2(r.right(), c.y),
                pos2(r.right() - tip, r.bottom()),
                r.left_bottom(),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        UnitKind::Mixer => {
            // Static mixer: pipe section with alternating mixing elements.
            painter.rect(r, 0, fill, stroke, eframe::egui::StrokeKind::Middle);
            let n = 4;
            let zig: Vec<Pos2> = (0..=n)
                .map(|i| {
                    let x = r.left() + r.width() * i as f32 / n as f32;
                    pos2(x, if i % 2 == 0 { r.top() } else { r.bottom() })
                })
                .collect();
            painter.line(zig, Stroke::new(stroke.width * 0.7, stroke.color));
        }
        UnitKind::Heater | UnitKind::Cooler => {
            // Heat exchanger: shell with the heating/cooling medium as a
            // zig-zag passing through; the medium arrow shows heating
            // (enters top) or cooling (enters bottom).
            let radius = r.width().min(r.height()) / 2.0;
            painter.circle(c, radius, fill, stroke);
            let ext = radius * 0.45;
            let w = radius * 0.4;
            let ys = [
                radius + ext,
                radius * 0.6,
                radius * 0.3,
                0.0,
                -radius * 0.3,
                -radius * 0.6,
                -radius - ext,
            ];
            let xs = [0.0, 0.0, -w, w, -w, 0.0, 0.0];
            let zig: Vec<Pos2> = xs.iter().zip(ys).map(|(&x, y)| c + vec2(x, y)).collect();
            painter.line(zig, stroke);
            let head = (radius * 0.5).max(5.0);
            if kind == UnitKind::Heater {
                arrowhead(
                    painter,
                    c - vec2(0.0, radius + ext * 0.1),
                    vec2(0.0, 1.0),
                    head,
                    stroke.color,
                );
            } else {
                arrowhead(
                    painter,
                    c + vec2(0.0, radius + ext * 0.1),
                    vec2(0.0, -1.0),
                    head,
                    stroke.color,
                );
            }
        }
        UnitKind::Valve => {
            // Control valve: globe body (bow-tie) with diaphragm actuator.
            let h = r.height() * 0.25;
            let body_c = pos2(c.x, c.y);
            let left = vec![
                pos2(r.left(), body_c.y - h),
                body_c,
                pos2(r.left(), body_c.y + h),
            ];
            let right = vec![
                pos2(r.right(), body_c.y - h),
                body_c,
                pos2(r.right(), body_c.y + h),
            ];
            painter.add(Shape::convex_polygon(left, fill, stroke));
            painter.add(Shape::convex_polygon(right, fill, stroke));
            let dome_r = r.width() * 0.22;
            let flat_y = r.top() + dome_r + stroke.width;
            painter.line_segment([body_c, pos2(c.x, flat_y)], stroke);
            let dome: Vec<Pos2> = (0..=16)
                .map(|i| {
                    let a = PI * i as f32 / 16.0;
                    pos2(c.x - dome_r * a.cos(), flat_y - dome_r * a.sin())
                })
                .collect();
            painter.add(Shape::convex_polygon(dome, fill, stroke));
        }
        UnitKind::Splitter => {
            // Branch: one inlet forking to the two outlet nozzles.
            painter.rect(r, 0, fill, stroke, eframe::egui::StrokeKind::Middle);
            let thin = Stroke::new(stroke.width * 0.7, stroke.color);
            let y1 = r.top() + r.height() / 3.0;
            let y2 = r.top() + r.height() * 2.0 / 3.0;
            let fork = pos2(c.x, c.y);
            painter.line_segment([pos2(r.left(), c.y), fork], thin);
            painter.line(
                vec![pos2(r.right(), y1), pos2(c.x + r.width() * 0.2, y1), fork],
                thin,
            );
            painter.line(
                vec![pos2(r.right(), y2), pos2(c.x + r.width() * 0.2, y2), fork],
                thin,
            );
        }
        UnitKind::HeatExchanger => {
            // Shell-and-tube exchanger: shell (hot side, upper nozzles) with
            // the tube bundle (cold side, lower nozzles) drawn as a zig-zag.
            let radius = (r.height() * 0.12).min(255.0) as u8;
            painter.rect(r, radius, fill, stroke, eframe::egui::StrokeKind::Middle);
            let y1 = r.top() + r.height() / 3.0;
            let y2 = r.top() + r.height() * 2.0 / 3.0;
            let peak = (y1 + y2) / 2.0 - r.height() * 0.05;
            let inset = r.width() * 0.1;
            let n = 6;
            let mut tubes = vec![pos2(r.left(), y2), pos2(r.left() + inset, y2)];
            for i in 1..n {
                let x = r.left() + inset + (r.width() - 2.0 * inset) * i as f32 / n as f32;
                tubes.push(pos2(x, if i % 2 == 1 { peak } else { y2 }));
            }
            tubes.push(pos2(r.right() - inset, y2));
            tubes.push(pos2(r.right(), y2));
            painter.line(tubes, Stroke::new(stroke.width * 0.8, stroke.color));
        }
        UnitKind::Pump => {
            // Centrifugal pump: casing with axial suction, tangential
            // discharge at the top, and a base.
            let radius = r.width() * crate::view::PUMP_RADIUS / crate::view::block_size(kind).x;
            let foot = vec![
                c + vec2(-radius * 0.7, radius * 0.7),
                pos2(c.x - radius, r.bottom()),
                pos2(c.x + radius, r.bottom()),
                c + vec2(radius * 0.7, radius * 0.7),
            ];
            painter.line(foot, stroke);
            painter.line_segment([pos2(r.left(), c.y), pos2(c.x - radius, c.y)], stroke);
            painter.circle(c, radius, fill, stroke);
            painter.line_segment(
                [pos2(c.x, c.y - radius), pos2(r.right(), c.y - radius)],
                stroke,
            );
        }
        UnitKind::Compressor | UnitKind::Expander => {
            // Turbomachine: trapezoid narrowing in the flow direction for a
            // compressor and widening for an expander (turbine).
            let small = r.height() * 0.3;
            let (left_h, right_h) = if kind == UnitKind::Compressor {
                (r.height() / 2.0, small)
            } else {
                (small, r.height() / 2.0)
            };
            let pts = vec![
                pos2(r.left(), c.y - left_h),
                pos2(r.right(), c.y - right_h),
                pos2(r.right(), c.y + right_h),
                pos2(r.left(), c.y + left_h),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        UnitKind::Separator => {
            // Vertical vessel with 2:1 semi-elliptical heads.
            let a = r.width() / 2.0;
            let b = a * 0.5;
            let top_tan = r.top() + b;
            let bot_tan = r.bottom() - b;
            let mut outline = Vec::with_capacity(40);
            for i in 0..=16 {
                let t = PI * i as f32 / 16.0;
                outline.push(pos2(c.x - a * t.cos(), top_tan - b * t.sin()));
            }
            for i in 0..=16 {
                let t = PI * i as f32 / 16.0;
                outline.push(pos2(c.x + a * t.cos(), bot_tan + b * t.sin()));
            }
            painter.add(Shape::convex_polygon(outline, fill, stroke));
            let thin = Stroke::new(stroke.width * 0.6, stroke.color);
            // Head-to-shell tangent lines.
            painter.line_segment([pos2(r.left(), top_tan), pos2(r.right(), top_tan)], thin);
            painter.line_segment([pos2(r.left(), bot_tan), pos2(r.right(), bot_tan)], thin);
            // Normal liquid level with level mark.
            let level = bot_tan - (bot_tan - top_tan) * 0.25;
            let inset = a * 0.25;
            painter.extend(Shape::dashed_line(
                &[
                    pos2(r.left() + inset, level),
                    pos2(r.right() - inset, level),
                ],
                thin,
                4.0,
                3.0,
            ));
            let m = a * 0.18;
            let mx = c.x + a * 0.35;
            painter.add(Shape::convex_polygon(
                vec![
                    pos2(mx, level),
                    pos2(mx - m, level - m * 1.2),
                    pos2(mx + m, level - m * 1.2),
                ],
                Color32::TRANSPARENT,
                thin,
            ));
        }
    }
}

/// Scales a symbol to fit an icon area, keeping its aspect ratio.
pub fn icon_rect(area: Rect, kind: UnitKind) -> Rect {
    let size = crate::view::block_size(kind);
    let scale = (area.width() / size.x).min(area.height() / size.y);
    Rect::from_center_size(area.center(), size * scale)
}
