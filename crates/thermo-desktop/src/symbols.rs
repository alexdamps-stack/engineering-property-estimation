//! Technical-drawing look: CAD colour scheme and process-flow-diagram symbols.
//!
//! Symbols are monochrome line work redrawn from the ISO 10628-2:2012 symbol
//! sheets by Steven Baltakatei Sandoval (Wikimedia Commons, CC BY-SA 4.0),
//! as extracted in `docs/design-reference/iso-10628`. Each symbol names the
//! asset ID it follows and keeps its proportions; nozzles sit on the
//! outline. The splitter and the feed/product flags are project symbols,
//! because the sheets have no equivalent. Arrowheads are closed and filled with a 15° half angle as in ISO 128. The dark
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
            // s4-x2673-in-line-static-mixer: square body with one diagonal
            // mixing element.
            painter.rect(r, 0, fill, stroke, eframe::egui::StrokeKind::Middle);
            let at = |fx: f32, fy: f32| pos2(r.left() + r.width() * fx, r.top() + r.height() * fy);
            painter.line(
                vec![
                    at(1.0 / 6.0, 5.0 / 12.0),
                    at(1.0 / 6.0, 1.0 / 12.0),
                    at(5.0 / 6.0, 11.0 / 12.0),
                    at(5.0 / 6.0, 7.0 / 12.0),
                ],
                stroke,
            );
        }
        UnitKind::Heater | UnitKind::Cooler => {
            // s2-x8079-heat-exchanger-general-condenser: circle crossed by
            // the process line with a Z-shaped bend, heating or cooling
            // medium nozzles at top and bottom. Project addition: an arrow
            // on the medium shows heat in (heater, top) or out (cooler,
            // bottom).
            let radius = r.width() / 2.0;
            let half = radius / 2.0;
            painter.circle(c, radius, fill, stroke);
            painter.line(
                vec![
                    c + vec2(-radius, 0.0),
                    c + vec2(-half, 0.0),
                    c + vec2(0.0, -half),
                    c + vec2(0.0, half),
                    c + vec2(half, 0.0),
                    c + vec2(radius, 0.0),
                ],
                stroke,
            );
            painter.line_segment([pos2(c.x, r.top()), pos2(c.x, c.y - radius)], stroke);
            painter.line_segment([pos2(c.x, r.bottom()), pos2(c.x, c.y + radius)], stroke);
            let head = (radius * 0.35).max(4.0);
            if kind == UnitKind::Heater {
                arrowhead(
                    painter,
                    pos2(c.x, c.y - radius),
                    vec2(0.0, 1.0),
                    head,
                    stroke.color,
                );
            } else {
                arrowhead(
                    painter,
                    pos2(c.x, r.bottom()),
                    vec2(0.0, 1.0),
                    head,
                    stroke.color,
                );
            }
        }
        UnitKind::Valve => {
            // s6-2101-valve-general: two triangles meeting at the seat.
            let left = vec![r.left_top(), c, r.left_bottom()];
            let right = vec![r.right_top(), c, r.right_bottom()];
            painter.add(Shape::convex_polygon(left, fill, stroke));
            painter.add(Shape::convex_polygon(right, fill, stroke));
        }
        UnitKind::Splitter => {
            // Branch: one inlet forking to the two outlet nozzles. ISO
            // 10628-2 has no unit symbol for a stream split; this is a
            // project symbol.
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
            // s2-2511-heat-exchanger-with-straight-tubes-fixed-tube-plates:
            // shell with a channel head and tube plate at each end and the
            // straight tube bundle between the plates.
            painter.rect(r, 0, fill, stroke, eframe::egui::StrokeKind::Middle);
            let plate = r.width() / 7.0;
            let (x1, x2) = (r.left() + plate, r.right() - plate);
            painter.line_segment([pos2(x1, r.top()), pos2(x1, r.bottom())], stroke);
            painter.line_segment([pos2(x2, r.top()), pos2(x2, r.bottom())], stroke);
            for i in 1..=3 {
                let y = r.top() + r.height() * i as f32 / 4.0;
                painter.line_segment([pos2(x1, y), pos2(x2, y)], stroke);
            }
        }
        UnitKind::Pump => {
            // s4-2301-pump-liquid-type-general: circle with two lines from
            // the top and bottom meeting at the discharge, so the apex
            // points in the flow direction.
            let radius = r.width().min(r.height()) / 2.0;
            painter.circle(c, radius, fill, stroke);
            painter.line(
                vec![
                    pos2(c.x, c.y - radius),
                    pos2(c.x + radius, c.y),
                    pos2(c.x, c.y + radius),
                ],
                stroke,
            );
        }
        UnitKind::Compressor => {
            // s4-2302-compressor-vacuum-pump-general: circle with two chords
            // converging in the flow direction (end points at 118 and 26
            // degrees either side of the flow axis, as in the reference).
            let radius = r.width().min(r.height()) / 2.0;
            painter.circle(c, radius, fill, stroke);
            let on_circle = |deg: f32| {
                let a = deg.to_radians();
                c + vec2(a.cos(), -a.sin()) * radius
            };
            let (wide, narrow) = (117.95_f32, 26.05_f32);
            painter.line_segment([on_circle(wide), on_circle(narrow)], stroke);
            painter.line_segment([on_circle(-wide), on_circle(-narrow)], stroke);
        }
        UnitKind::Expander => {
            // s5-2571-turbine-general: trapezoid widening in the flow
            // direction; the inlet side is one third as high as the outlet
            // side.
            let narrow = r.height() / 6.0;
            let pts = vec![
                pos2(r.left(), c.y - narrow),
                r.right_top(),
                r.right_bottom(),
                pos2(r.left(), c.y + narrow),
            ];
            painter.add(Shape::convex_polygon(pts, fill, stroke));
        }
        UnitKind::Separator => {
            // s1-2062-tank-vessel-with-dished-ends: vertical vessel with
            // shallow dished heads (depth one eighth of the diameter).
            let a = r.width() / 2.0;
            let b = r.width() / 8.0;
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
        }
    }
}

/// Scales a symbol to fit an icon area, keeping its aspect ratio.
pub fn icon_rect(area: Rect, kind: UnitKind) -> Rect {
    let size = crate::view::block_size(kind);
    let scale = (area.width() / size.x).min(area.height() / size.y);
    Rect::from_center_size(area.center(), size * scale)
}
