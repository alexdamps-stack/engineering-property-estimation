//! Interactive flowsheet canvas drawn as a technical process flow diagram:
//! orthogonal stream lines with ISO 128 arrowheads, numbered stream tags,
//! monochrome ISO 10628-style symbols, CAD lettering, a model-space grid,
//! and CAD-like grips and snap markers.

use eframe::egui::{
    Align2, Color32, FontId, Key, PointerButton, Pos2, Rect, Response, Sense, Shape, Stroke, Ui,
    Vec2, pos2, vec2,
};
use flowsheet_core::{PortDirection, PortRef, StreamId, UnitId, UnitKind, UnitStatus};

use crate::editor::{Drag, Editor, Selection};
use crate::symbols::{Theme, arrowhead, draw_unit};
use crate::view::{
    SNAP, block_rect, distance_to_polyline, place_tag, port_exit, port_pos, snap, stream_route,
};

const PORT_HIT_RADIUS: f32 = 9.0;
const STREAM_HIT_DISTANCE: f32 = 5.0;
const MAJOR_EVERY: i32 = 5;
const GRIP: f32 = 7.0;
const BADGE_RADIUS: f32 = 6.5;

/// Centre of the status badge of a symbol drawn in screen rectangle `r`.
fn badge_center(r: Rect, scale: f32) -> Pos2 {
    r.right_top() + vec2(4.0, -4.0) * scale
}

/// Shows the canvas filling the available space of `ui`.
pub fn show(ui: &mut Ui, ed: &mut Editor) {
    let (response, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let canvas = response.rect;
    ed.canvas_rect = canvas;
    let painter = painter.with_clip_rect(canvas);
    let theme = Theme::from_visuals(ui.visuals());

    if std::mem::take(&mut ed.fit_requested) {
        ed.fit_to_content(canvas);
    }

    handle_view_input(ui, &response, ed, canvas);
    handle_edit_input(ui, &response, ed, canvas);

    let z = ed.viewport.zoom;
    let scale = z.clamp(0.6, 1.6);
    let process_w = (1.6 * z).clamp(1.0, 3.0);
    let equipment_w = (1.3 * z).clamp(1.0, 2.5);
    let font = |size: f32| FontId::monospace((size * z).clamp(8.0, size * 1.5));

    draw_grid(&painter, ed, canvas, &theme);

    let hover = response.hover_pos();
    let idle = ed.drag.is_none();
    let hovered_port = hover.and_then(|p| port_at(ed, canvas, p));
    let hovered_unit = hover
        .filter(|_| idle && hovered_port.is_none())
        .and_then(|p| unit_at(ed, canvas, p));
    let hovered_stream = hover
        .filter(|_| idle && hovered_port.is_none() && hovered_unit.is_none())
        .and_then(|p| stream_at(ed, canvas, p));

    // Process lines.
    let routes: Vec<(StreamId, String, Vec<Pos2>)> = ed
        .fs
        .streams()
        .filter_map(|s| Some((s.id, s.name.clone(), stream_screen_route(ed, canvas, s.id)?)))
        .collect();
    let stream_color = |id: StreamId| {
        if ed.selection == Some(Selection::Stream(id)) {
            theme.select
        } else {
            theme.line
        }
    };
    for (id, _, route) in &routes {
        let color = stream_color(*id);
        let width = if hovered_stream == Some(*id) {
            process_w * 1.8
        } else {
            process_w
        };
        painter.line(route.clone(), Stroke::new(width, color));
        let n = route.len();
        arrowhead(
            &painter,
            route[n - 1],
            route[n - 1] - route[n - 2],
            11.0 * scale,
            color,
        );
        if ed.selection == Some(Selection::Stream(*id)) {
            for &p in route {
                grip(&painter, p, &theme);
            }
        }
    }

    // Stream numbers: beside their line and clear of equipment, other
    // numbers, and every stream line.
    if z >= 0.45 {
        let lines: Vec<&[Pos2]> = routes.iter().map(|(_, _, r)| r.as_slice()).collect();
        // Equipment symbols and their status badges are obstacles.
        let mut taken: Vec<Rect> = Vec::new();
        for u in ed.fs.units() {
            let Some(&center) = ed.layout.get(&u.id) else {
                continue;
            };
            let world = block_rect(u.kind, center);
            let r = Rect::from_min_max(
                ed.viewport.to_screen(canvas, world.min),
                ed.viewport.to_screen(canvas, world.max),
            );
            taken.push(r.expand(2.0));
            if ed.fs.unit_status(u.id) == Some(UnitStatus::Underspecified) {
                let d = 2.0 * BADGE_RADIUS * scale + 4.0;
                taken.push(Rect::from_center_size(
                    badge_center(r, scale),
                    Vec2::splat(d),
                ));
            }
        }
        for (id, name, route) in &routes {
            let color = stream_color(*id);
            let galley = painter.layout_no_wrap(name.clone(), font(11.0), color);
            let half = tag_half_size(galley.size());
            if let Some(at) = place_tag(route, half, 3.0, 16.0 * scale, &taken, &lines) {
                stream_tag(&painter, at, half, galley, color, &theme);
                taken.push(Rect::from_center_size(at, half * 2.0).expand(1.0));
            }
        }
    }

    // Rubber-band line while connecting.
    if let (Some(Drag::Connecting { from }), Some(pointer)) = (ed.drag, hover)
        && let Some(route) = pending_route(ed, canvas, from, pointer)
    {
        painter.extend(Shape::dashed_line(
            &route,
            Stroke::new(process_w, theme.select),
            8.0,
            5.0,
        ));
    }

    // Equipment.
    for unit in ed.fs.units() {
        let Some(&center) = ed.layout.get(&unit.id) else {
            continue;
        };
        let world = block_rect(unit.kind, center);
        let r = Rect::from_min_max(
            ed.viewport.to_screen(canvas, world.min),
            ed.viewport.to_screen(canvas, world.max),
        );
        let selected = ed.selection == Some(Selection::Unit(unit.id));
        // Proposal section 4.6: underspecified units get a grey outline and
        // a `!` badge; ready units are drawn normally.
        let underspecified = ed.fs.unit_status(unit.id) == Some(UnitStatus::Underspecified);
        let color = if selected {
            theme.select
        } else if underspecified {
            theme.muted
        } else {
            theme.line
        };
        let width = if hovered_unit == Some(unit.id) {
            equipment_w * 1.8
        } else {
            equipment_w
        };
        draw_unit(
            &painter,
            unit.kind,
            r,
            Stroke::new(width, color),
            theme.background,
        );

        if z >= 0.35 {
            let (pos, align) = match unit.kind {
                // Bottom nozzle is taken by the liquid line.
                UnitKind::Separator => (pos2(r.left() - 8.0, r.top() + 4.0), Align2::RIGHT_TOP),
                // Clear the heating/cooling medium stub below the shell.
                UnitKind::Heater | UnitKind::Cooler => (
                    pos2(r.right() + 2.0 * scale, r.bottom() + 2.0 * scale),
                    Align2::LEFT_TOP,
                ),
                _ => (
                    pos2(r.center().x, r.bottom() + 5.0 * scale),
                    Align2::CENTER_TOP,
                ),
            };
            painter.text(pos, align, &unit.name, font(12.0), theme.line);
        }
        if underspecified {
            let at = badge_center(r, scale);
            let rad = BADGE_RADIUS * scale;
            painter.circle(at, rad, theme.background, Stroke::new(1.0, theme.muted));
            painter.text(at, Align2::CENTER_CENTER, "!", font(10.0), theme.muted);
        }
        if selected {
            painter.extend(Shape::dashed_line(
                &[
                    r.left_top(),
                    r.right_top(),
                    r.right_bottom(),
                    r.left_bottom(),
                    r.left_top(),
                ],
                Stroke::new(1.0, theme.select),
                5.0,
                4.0,
            ));
            for p in [
                r.left_top(),
                r.right_top(),
                r.right_bottom(),
                r.left_bottom(),
                r.center(),
            ] {
                grip(&painter, p, &theme);
            }
        }

        // Nozzles: free ones carry a snap marker; connected ones are plain
        // line ends, as on a drawing.
        for (index, spec) in unit.kind.ports().iter().enumerate() {
            let port = PortRef::new(unit.id, index);
            let p = ed
                .viewport
                .to_screen(canvas, port_pos(unit.kind, center, index));
            let free = ed.fs.stream_at(port).is_none();
            let target =
                matches!(ed.drag, Some(Drag::Connecting { from }) if ed.can_connect(from, port));
            let hovered = hovered_port == Some(port);
            if free || hovered {
                let size = if hovered || target { 10.0 } else { 7.0 };
                let marker = Rect::from_center_size(p, Vec2::splat(size));
                let fill = if target {
                    theme.snap
                } else {
                    Color32::TRANSPARENT
                };
                painter.rect(
                    marker,
                    0,
                    fill,
                    Stroke::new(1.5, theme.snap),
                    eframe::egui::StrokeKind::Middle,
                );
            }
            if hovered {
                let exit = port_exit(unit.kind, index);
                let (align, offset) = match (spec.direction, exit.y) {
                    (_, y) if y < 0.0 => (Align2::LEFT_BOTTOM, vec2(8.0, -6.0)),
                    (_, y) if y > 0.0 => (Align2::LEFT_TOP, vec2(8.0, 6.0)),
                    (PortDirection::Inlet, _) => (Align2::RIGHT_BOTTOM, vec2(-8.0, -6.0)),
                    (PortDirection::Outlet, _) => (Align2::LEFT_BOTTOM, vec2(8.0, -6.0)),
                };
                painter.text(
                    p + offset,
                    align,
                    spec.name,
                    FontId::monospace(11.0),
                    theme.snap,
                );
            }
        }
    }

    if ed.fs.units().next().is_none() {
        painter.text(
            canvas.center(),
            Align2::CENTER_CENTER,
            "Place a unit from the palette, or load the demonstration flowsheet.",
            FontId::monospace(13.0),
            theme.muted,
        );
    }
}

/// Model-space grid with a major line every fifth snap increment.
fn draw_grid(painter: &eframe::egui::Painter, ed: &Editor, canvas: Rect, theme: &Theme) {
    painter.rect_filled(canvas, 0, theme.background);
    let minor = SNAP * ed.viewport.zoom;
    let major = minor * MAJOR_EVERY as f32;
    let (step, every) = if minor >= 7.0 {
        (minor, MAJOR_EVERY)
    } else if major >= 7.0 {
        (major, 1)
    } else {
        return;
    };
    let origin = ed.viewport.to_screen(canvas, Pos2::ZERO);
    let first_x = ((canvas.left() - origin.x) / step).floor() as i32;
    let first_y = ((canvas.top() - origin.y) / step).floor() as i32;
    let pick = |i: i32| {
        let c = if i.rem_euclid(every) == 0 {
            theme.grid_major
        } else {
            theme.grid_minor
        };
        Stroke::new(1.0, c)
    };
    let mut i = first_x;
    loop {
        let x = origin.x + i as f32 * step;
        if x > canvas.right() {
            break;
        }
        painter.vline(x, canvas.y_range(), pick(i));
        i += 1;
    }
    let mut j = first_y;
    loop {
        let y = origin.y + j as f32 * step;
        if y > canvas.bottom() {
            break;
        }
        painter.hline(canvas.x_range(), y, pick(j));
        j += 1;
    }
}

/// Half extents of the diamond enclosing a stream number of text size `text`.
fn tag_half_size(text: Vec2) -> Vec2 {
    let hw = text.x / 2.0 + text.y;
    vec2(hw, hw / 2.0 + 1.0)
}

/// Draws a stream number in a diamond centred at `at`.
fn stream_tag(
    painter: &eframe::egui::Painter,
    at: Pos2,
    half: Vec2,
    galley: std::sync::Arc<eframe::egui::Galley>,
    color: Color32,
    theme: &Theme,
) {
    let size = galley.size();
    let diamond = vec![
        at - vec2(half.x, 0.0),
        at - vec2(0.0, half.y),
        at + vec2(half.x, 0.0),
        at + vec2(0.0, half.y),
    ];
    painter.add(Shape::convex_polygon(
        diamond,
        theme.background,
        Stroke::new(1.0, color),
    ));
    painter.galley(at - size / 2.0, galley, color);
}

/// CAD grip: small filled square.
fn grip(painter: &eframe::egui::Painter, at: Pos2, theme: &Theme) {
    painter.rect(
        Rect::from_center_size(at, Vec2::splat(GRIP)),
        0,
        theme.select,
        Stroke::new(1.0, theme.background),
        eframe::egui::StrokeKind::Middle,
    );
}

fn handle_view_input(ui: &Ui, response: &Response, ed: &mut Editor, canvas: Rect) {
    if !response.hovered() {
        return;
    }
    let (zoom, scroll, pointer) =
        ui.input(|i| (i.zoom_delta(), i.smooth_scroll_delta, i.pointer.hover_pos()));
    if zoom != 1.0 {
        let anchor = pointer.unwrap_or(canvas.center());
        ed.viewport.zoom_about(canvas, anchor, zoom);
    }
    if scroll != Vec2::ZERO {
        ed.viewport.pan += scroll;
    }
}

fn handle_edit_input(ui: &Ui, response: &Response, ed: &mut Editor, canvas: Rect) {
    let press_origin = ui.input(|i| i.pointer.press_origin());

    if response.drag_started_by(PointerButton::Middle) {
        ed.drag = Some(Drag::Panning);
    } else if response.drag_started_by(PointerButton::Primary)
        && let Some(origin) = press_origin
    {
        ed.drag = Some(start_drag(ed, canvas, origin));
    }

    if response.dragged() {
        match ed.drag {
            Some(Drag::Panning) => ed.viewport.pan += response.drag_delta(),
            Some(Drag::Moving { unit, grab }) => {
                if let Some(pointer) = response.interact_pointer_pos() {
                    let world = ed.viewport.to_world(canvas, pointer);
                    // Snap while dragging so lines stay square on the grid.
                    ed.layout.insert(unit, snap(world - grab));
                }
            }
            _ => {}
        }
    }

    if response.drag_stopped()
        && let Some(Drag::Connecting { from }) = ed.drag.take()
        && let Some(to) = response
            .interact_pointer_pos()
            .and_then(|p| port_at(ed, canvas, p))
    {
        ed.connect(from, to);
    }

    if response.clicked()
        && let Some(p) = response.interact_pointer_pos()
    {
        ed.selection = if let Some(unit) = unit_at(ed, canvas, p) {
            Some(Selection::Unit(unit))
        } else {
            stream_at(ed, canvas, p).map(Selection::Stream)
        };
    }

    if !ui.ctx().egui_wants_keyboard_input() {
        let (delete, escape) = ui.input(|i| {
            (
                i.key_pressed(Key::Delete) || i.key_pressed(Key::Backspace),
                i.key_pressed(Key::Escape),
            )
        });
        if delete {
            ed.delete_selection();
        }
        if escape {
            ed.drag = None;
            ed.selection = None;
        }
    }
}

fn start_drag(ed: &mut Editor, canvas: Rect, origin: Pos2) -> Drag {
    if let Some(port) = port_at(ed, canvas, origin)
        && ed.fs.stream_at(port).is_none()
    {
        return Drag::Connecting { from: port };
    }
    if let Some(unit) = unit_at(ed, canvas, origin) {
        ed.selection = Some(Selection::Unit(unit));
        let center = ed.layout.get(&unit).copied().unwrap_or_default();
        let grab = ed.viewport.to_world(canvas, origin) - center;
        return Drag::Moving { unit, grab };
    }
    Drag::Panning
}

fn port_world(ed: &Editor, port: PortRef) -> Option<(Pos2, Vec2)> {
    let unit = ed.fs.unit(port.unit)?;
    let center = *ed.layout.get(&port.unit)?;
    Some((
        port_pos(unit.kind, center, port.port),
        port_exit(unit.kind, port.port),
    ))
}

fn to_screen(ed: &Editor, canvas: Rect, route: Vec<Pos2>) -> Vec<Pos2> {
    route
        .into_iter()
        .map(|p| ed.viewport.to_screen(canvas, p))
        .collect()
}

fn stream_screen_route(ed: &Editor, canvas: Rect, id: StreamId) -> Option<Vec<Pos2>> {
    let stream = ed.fs.stream(id)?;
    let (a, exit) = port_world(ed, stream.from)?;
    let (b, _) = port_world(ed, stream.to)?;
    Some(to_screen(ed, canvas, stream_route(a, exit, b)))
}

/// Route from the nozzle being connected to the pointer.
fn pending_route(ed: &Editor, canvas: Rect, from: PortRef, pointer: Pos2) -> Option<Vec<Pos2>> {
    let (p, exit) = port_world(ed, from)?;
    let cursor = ed.viewport.to_world(canvas, pointer);
    let route = match ed.fs.port_spec(from)?.direction {
        PortDirection::Outlet => stream_route(p, exit, cursor),
        PortDirection::Inlet => stream_route(cursor, vec2(1.0, 0.0), p),
    };
    Some(to_screen(ed, canvas, route))
}

fn port_at(ed: &Editor, canvas: Rect, p: Pos2) -> Option<PortRef> {
    ed.fs.units().find_map(|unit| {
        (0..unit.kind.ports().len()).find_map(|i| {
            let port = PortRef::new(unit.id, i);
            let (pos, _) = port_world(ed, port)?;
            (ed.viewport.to_screen(canvas, pos).distance(p) <= PORT_HIT_RADIUS).then_some(port)
        })
    })
}

fn unit_at(ed: &Editor, canvas: Rect, p: Pos2) -> Option<UnitId> {
    let world = ed.viewport.to_world(canvas, p);
    let tolerance = 3.0 / ed.viewport.zoom;
    // Topmost (last drawn) first.
    let units: Vec<_> = ed.fs.units().collect();
    units.into_iter().rev().find_map(|unit| {
        let center = *ed.layout.get(&unit.id)?;
        block_rect(unit.kind, center)
            .expand(tolerance)
            .contains(world)
            .then_some(unit.id)
    })
}

fn stream_at(ed: &Editor, canvas: Rect, p: Pos2) -> Option<StreamId> {
    ed.fs.streams().find_map(|s| {
        let route = stream_screen_route(ed, canvas, s.id)?;
        (distance_to_polyline(p, &route) <= STREAM_HIT_DISTANCE).then_some(s.id)
    })
}
