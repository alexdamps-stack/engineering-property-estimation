//! Canvas geometry: world/screen transform, block shapes, ports, and
//! orthogonal stream routing.
//!
//! World coordinates are unscaled drawing units (1 unit ≈ 1 px at 100 %);
//! the [`Viewport`] maps them to screen points with a uniform scale, so
//! horizontal and vertical segments stay axis-aligned. Everything here is
//! pure geometry so it can be unit-tested without a running GUI.

use eframe::egui::{Pos2, Rect, Vec2, pos2, vec2};
use flowsheet_core::{PortDirection, UnitKind};

/// Minimum and maximum zoom factors.
pub const ZOOM_RANGE: (f32, f32) = (0.25, 4.0);

/// Snap increment for block positions and route bends (world units).
pub const SNAP: f32 = 10.0;

/// Straight lead-out of a stream from a nozzle before its first bend.
const STUB: f32 = 20.0;

/// Pan and zoom state of the flowsheet canvas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Viewport {
    /// Screen-space offset of the world origin from the canvas' top-left corner.
    pub pan: Vec2,
    /// Screen points per world unit.
    pub zoom: f32,
}

impl Default for Viewport {
    fn default() -> Self {
        Self {
            pan: vec2(40.0, 40.0),
            zoom: 1.0,
        }
    }
}

impl Viewport {
    /// Maps a world position to a screen position.
    pub fn to_screen(self, canvas: Rect, world: Pos2) -> Pos2 {
        canvas.min + self.pan + world.to_vec2() * self.zoom
    }

    /// Maps a screen position to a world position.
    pub fn to_world(self, canvas: Rect, screen: Pos2) -> Pos2 {
        ((screen - canvas.min - self.pan) / self.zoom).to_pos2()
    }

    /// Multiplies the zoom by `factor`, keeping the world point under
    /// `anchor` (screen) fixed.
    pub fn zoom_about(&mut self, canvas: Rect, anchor: Pos2, factor: f32) {
        let world = self.to_world(canvas, anchor);
        self.zoom = (self.zoom * factor).clamp(ZOOM_RANGE.0, ZOOM_RANGE.1);
        self.pan = anchor - canvas.min - world.to_vec2() * self.zoom;
    }

    /// Fits a world-space rectangle into the canvas with a margin.
    pub fn fit(&mut self, canvas: Rect, world: Rect) {
        let margin = 60.0;
        let avail = (canvas.size() - Vec2::splat(2.0 * margin)).max(Vec2::splat(1.0));
        let size = world.size().max(Vec2::splat(1.0));
        self.zoom = (avail.x / size.x)
            .min(avail.y / size.y)
            .clamp(ZOOM_RANGE.0, 1.5);
        self.pan = canvas.size() / 2.0 - world.center().to_vec2() * self.zoom;
    }
}

/// Rounds a world position to the snap grid.
pub fn snap(p: Pos2) -> Pos2 {
    pos2((p.x / SNAP).round() * SNAP, (p.y / SNAP).round() * SNAP)
}

/// World-space size of a unit symbol.
pub fn block_size(kind: UnitKind) -> Vec2 {
    match kind {
        UnitKind::Feed | UnitKind::Product => vec2(70.0, 30.0),
        UnitKind::Mixer => vec2(60.0, 60.0),
        UnitKind::Heater | UnitKind::Cooler => vec2(50.0, 50.0),
        UnitKind::Valve => vec2(40.0, 40.0),
        UnitKind::Separator => vec2(50.0, 130.0),
    }
}

/// World-space rectangle of a unit symbol centred at `center`.
pub fn block_rect(kind: UnitKind, center: Pos2) -> Rect {
    Rect::from_center_size(center, block_size(kind))
}

/// Direction in which a stream leaves (outlet) or arrives at (inlet,
/// reversed) a nozzle. Inlets are entered from the left; separator vapor
/// leaves through the top head and liquid through the bottom head, as on a
/// process flow diagram.
pub fn port_exit(kind: UnitKind, index: usize) -> Vec2 {
    match (kind, kind.ports().get(index).map(|p| p.direction)) {
        (UnitKind::Separator, Some(PortDirection::Outlet)) if index == 1 => vec2(0.0, -1.0),
        (UnitKind::Separator, Some(PortDirection::Outlet)) => vec2(0.0, 1.0),
        (_, Some(PortDirection::Inlet)) => vec2(-1.0, 0.0),
        _ => vec2(1.0, 0.0),
    }
}

/// World-space position of nozzle `index` of a unit centred at `center`.
pub fn port_pos(kind: UnitKind, center: Pos2, index: usize) -> Pos2 {
    let rect = block_rect(kind, center);
    let ports = kind.ports();
    let Some(spec) = ports.get(index) else {
        return center;
    };
    let exit = port_exit(kind, index);
    if exit.y < 0.0 {
        return rect.center_top();
    }
    if exit.y > 0.0 {
        return rect.center_bottom();
    }
    let same_side: Vec<usize> = ports
        .iter()
        .enumerate()
        .filter(|&(i, p)| p.direction == spec.direction && port_exit(kind, i) == exit)
        .map(|(i, _)| i)
        .collect();
    let rank = same_side.iter().position(|&i| i == index).unwrap_or(0);
    let n = same_side.len() as f32;
    let x = if exit.x < 0.0 {
        rect.left()
    } else {
        rect.right()
    };
    pos2(
        x,
        rect.top() + rect.height() * (rank as f32 + 1.0) / (n + 1.0),
    )
}

/// Orthogonal route of a stream from a nozzle at `from` leaving in direction
/// `exit` to an inlet nozzle at `to` that is entered from the left.
///
/// Every segment is horizontal or vertical; bends lie on the snap grid where
/// possible. The route always ends with a horizontal segment pointing in +x
/// so the arrowhead enters the inlet correctly.
pub fn stream_route(from: Pos2, exit: Vec2, to: Pos2) -> Vec<Pos2> {
    let lead_in = to - vec2(STUB, 0.0);
    let mut pts = vec![from];
    if exit.y != 0.0 {
        let lead_out = from + exit * STUB;
        let same_way = (to.y - lead_out.y) * exit.y >= 0.0;
        if same_way && to.x >= from.x + STUB {
            pts.push(pos2(from.x, to.y));
        } else {
            pts.push(lead_out);
            pts.push(pos2(lead_in.x, lead_out.y));
            pts.push(pos2(lead_in.x, to.y));
        }
    } else {
        let lead_out = from + vec2(STUB, 0.0);
        if lead_in.x >= lead_out.x {
            let mid = ((from.x + to.x) / 2.0 / SNAP).round() * SNAP;
            let mid = mid.clamp(lead_out.x, lead_in.x);
            pts.push(pos2(mid, from.y));
            pts.push(pos2(mid, to.y));
        } else {
            // Backward run (for example a recycle): go round below both nozzles,
            // or between them when they are far enough apart vertically.
            let dy = (to.y - from.y).abs();
            let y = if dy >= 4.0 * STUB {
                ((from.y + to.y) / 2.0 / SNAP).round() * SNAP
            } else {
                from.y.max(to.y) + 4.0 * STUB
            };
            pts.extend([lead_out, pos2(lead_out.x, y), pos2(lead_in.x, y), lead_in]);
        }
    }
    pts.push(to);
    simplify(pts)
}

/// Removes repeated points and merges collinear segments.
fn simplify(pts: Vec<Pos2>) -> Vec<Pos2> {
    let mut out: Vec<Pos2> = Vec::with_capacity(pts.len());
    for p in pts {
        if out.last().is_some_and(|q| q.distance(p) < 1e-3) {
            continue;
        }
        if out.len() >= 2 {
            let a = out[out.len() - 2];
            let b = out[out.len() - 1];
            let collinear = (a.x == b.x && b.x == p.x) || (a.y == b.y && b.y == p.y);
            if collinear {
                out.pop();
            }
        }
        out.push(p);
    }
    out
}

/// Anchor for a stream number: the middle of the longest horizontal run,
/// falling back to the longest run. The last `arrow_clearance` of the route
/// is kept free for the arrowhead, and earlier runs win ties so tags sit
/// near the source, as is usual on a PFD.
///
/// Returns the anchor on the line and whether its run is horizontal.
pub fn route_label_anchor(route: &[Pos2], arrow_clearance: f32) -> Option<(Pos2, bool)> {
    let last = route.len().checked_sub(2)?;
    let usable = |i: usize, w: &[Pos2]| {
        let (a, mut b) = (w[0], w[1]);
        if i == last {
            let len = a.distance(b);
            b = a + (b - a) * ((len - arrow_clearance).max(0.0) / len.max(1e-6));
        }
        (a, b)
    };
    let best = |horizontal_only: bool| {
        let mut best: Option<(f32, Pos2)> = None;
        for (i, w) in route.windows(2).enumerate() {
            if horizontal_only && w[0].y != w[1].y {
                continue;
            }
            let (a, b) = usable(i, w);
            let len = a.distance(b);
            if best.is_none_or(|(l, _)| len > l + 0.5) {
                best = Some((len, a + (b - a) / 2.0));
            }
        }
        best.map(|(_, p)| p)
    };
    best(true)
        .map(|p| (p, true))
        .or_else(|| best(false).map(|p| (p, false)))
}

/// Centre of a stream tag with half extents `half` placed beside the line at
/// `anchor`, leaving `gap` between tag and line: above a horizontal run,
/// right of a vertical run. The tag never covers the line it labels.
pub fn tag_center(anchor: Pos2, horizontal: bool, half: Vec2, gap: f32) -> Pos2 {
    if horizontal {
        anchor - vec2(0.0, half.y + gap)
    } else {
        anchor + vec2(half.x + gap, 0.0)
    }
}

/// Shortest distance from `p` to a polyline.
pub fn distance_to_polyline(p: Pos2, line: &[Pos2]) -> f32 {
    line.windows(2)
        .map(|w| distance_to_segment(p, w[0], w[1]))
        .fold(f32::INFINITY, f32::min)
}

fn distance_to_segment(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let len2 = ab.length_sq();
    if len2 == 0.0 {
        return p.distance(a);
    }
    let t = ((p - a).dot(ab) / len2).clamp(0.0, 1.0);
    p.distance(a + ab * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn canvas() -> Rect {
        Rect::from_min_size(pos2(100.0, 50.0), vec2(800.0, 600.0))
    }

    fn assert_orthogonal_into_inlet(route: &[Pos2], from: Pos2, to: Pos2) {
        assert_eq!(route.first().copied(), Some(from));
        assert_eq!(route.last().copied(), Some(to));
        for w in route.windows(2) {
            assert!(
                w[0].x == w[1].x || w[0].y == w[1].y,
                "diagonal segment {w:?} in {route:?}"
            );
        }
        let n = route.len();
        let last = route[n - 1] - route[n - 2];
        assert!(
            last.x > 0.0 && last.y == 0.0,
            "must enter from the left: {route:?}"
        );
    }

    #[test]
    fn screen_world_round_trip() {
        let v = Viewport {
            pan: vec2(13.0, -7.0),
            zoom: 1.7,
        };
        let w = pos2(42.0, -3.5);
        let back = v.to_world(canvas(), v.to_screen(canvas(), w));
        assert!(back.distance(w) < 1e-4);
    }

    #[test]
    fn zoom_keeps_anchor_fixed_and_respects_limits() {
        let mut v = Viewport::default();
        let anchor = pos2(400.0, 300.0);
        let before = v.to_world(canvas(), anchor);
        v.zoom_about(canvas(), anchor, 1.5);
        assert!(v.to_world(canvas(), anchor).distance(before) < 1e-3);
        v.zoom_about(canvas(), anchor, 1000.0);
        assert_eq!(v.zoom, ZOOM_RANGE.1);
        v.zoom_about(canvas(), anchor, 1e-6);
        assert_eq!(v.zoom, ZOOM_RANGE.0);
    }

    #[test]
    fn fit_centres_the_world_rect() {
        let mut v = Viewport::default();
        let world = Rect::from_min_max(pos2(-200.0, 0.0), pos2(600.0, 300.0));
        v.fit(canvas(), world);
        let c = v.to_screen(canvas(), world.center());
        assert!(c.distance(canvas().center()) < 1e-3);
    }

    #[test]
    fn nozzles_lie_on_the_symbol_boundary_facing_their_exit() {
        for kind in UnitKind::ALL {
            let rect = block_rect(kind, Pos2::ZERO);
            for i in 0..kind.ports().len() {
                let p = port_pos(kind, Pos2::ZERO, i);
                let exit = port_exit(kind, i);
                let on_edge = match (exit.x, exit.y) {
                    (x, _) if x < 0.0 => p.x == rect.left(),
                    (x, _) if x > 0.0 => p.x == rect.right(),
                    (_, y) if y < 0.0 => p.y == rect.top(),
                    _ => p.y == rect.bottom(),
                };
                assert!(on_edge, "{kind} nozzle {i} at {p:?}");
            }
        }
    }

    #[test]
    fn separator_vapor_leaves_top_and_liquid_leaves_bottom() {
        let vapor = port_pos(UnitKind::Separator, Pos2::ZERO, 1);
        let liquid = port_pos(UnitKind::Separator, Pos2::ZERO, 2);
        assert_eq!(port_exit(UnitKind::Separator, 1), vec2(0.0, -1.0));
        assert_eq!(port_exit(UnitKind::Separator, 2), vec2(0.0, 1.0));
        assert!(vapor.y < liquid.y);
        assert_eq!(vapor.x, 0.0);
        let a = port_pos(UnitKind::Mixer, Pos2::ZERO, 0);
        let b = port_pos(UnitKind::Mixer, Pos2::ZERO, 1);
        assert!(a.y < b.y);
    }

    #[test]
    fn straight_and_forward_routes() {
        let right = vec2(1.0, 0.0);
        let (a, b) = (pos2(0.0, 0.0), pos2(200.0, 0.0));
        assert_eq!(stream_route(a, right, b), vec![a, b]);
        let (a, b) = (pos2(0.0, 0.0), pos2(200.0, 80.0));
        let r = stream_route(a, right, b);
        assert_orthogonal_into_inlet(&r, a, b);
        assert_eq!(r.len(), 4);
        assert_eq!(r[1], pos2(100.0, 0.0));
    }

    #[test]
    fn backward_routes_stay_orthogonal() {
        let right = vec2(1.0, 0.0);
        for (a, b) in [
            (pos2(300.0, 0.0), pos2(0.0, 0.0)),
            (pos2(300.0, 0.0), pos2(0.0, 200.0)),
            (pos2(10.0, 0.0), pos2(15.0, 5.0)),
        ] {
            let r = stream_route(a, right, b);
            assert_orthogonal_into_inlet(&r, a, b);
        }
    }

    #[test]
    fn vertical_nozzle_routes() {
        let up = vec2(0.0, -1.0);
        let down = vec2(0.0, 1.0);
        // Vapor up and over to a product above-right: one bend.
        let (a, b) = (pos2(0.0, 0.0), pos2(150.0, -60.0));
        let r = stream_route(a, up, b);
        assert_orthogonal_into_inlet(&r, a, b);
        assert_eq!(r, vec![a, pos2(0.0, -60.0), b]);
        // Target on the "wrong" side or to the left: still orthogonal.
        for b in [pos2(150.0, 60.0), pos2(-150.0, -60.0), pos2(-150.0, 60.0)] {
            assert_orthogonal_into_inlet(&stream_route(a, up, b), a, b);
            assert_orthogonal_into_inlet(&stream_route(a, down, b), a, b);
        }
    }

    #[test]
    fn label_sits_on_the_longest_horizontal_run() {
        let r = vec![
            pos2(0.0, 0.0),
            pos2(20.0, 0.0),
            pos2(20.0, 100.0),
            pos2(120.0, 100.0),
        ];
        assert_eq!(route_label_anchor(&r, 0.0), Some((pos2(70.0, 100.0), true)));
        // The arrowhead zone is excluded, and earlier runs win ties.
        assert_eq!(
            route_label_anchor(&r, 20.0),
            Some((pos2(60.0, 100.0), true))
        );
        let tie = [
            pos2(0.0, 0.0),
            pos2(100.0, 0.0),
            pos2(100.0, 50.0),
            pos2(220.0, 50.0),
        ];
        assert_eq!(
            route_label_anchor(&tie, 20.0),
            Some((pos2(50.0, 0.0), true))
        );
        let vertical = [pos2(0.0, 0.0), pos2(0.0, 100.0)];
        assert_eq!(
            route_label_anchor(&vertical, 0.0),
            Some((pos2(0.0, 50.0), false))
        );
        assert!(distance_to_polyline(pos2(70.0, 100.0), &r) < 1e-3);
        assert!(distance_to_polyline(pos2(70.0, 40.0), &r) > 40.0);
    }

    #[test]
    fn tags_sit_beside_the_line_without_covering_it() {
        let half = vec2(18.0, 10.0);
        let gap = 3.0;
        // Horizontal run: tag above, its bottom corner `gap` above the line.
        let c = tag_center(pos2(50.0, 0.0), true, half, gap);
        assert_eq!(c, pos2(50.0, -13.0));
        assert!(c.y + half.y < 0.0);
        // Vertical run: tag to the right, its left corner `gap` right of the line.
        let c = tag_center(pos2(0.0, 50.0), false, half, gap);
        assert_eq!(c, pos2(21.0, 50.0));
        assert!(c.x - half.x > 0.0);
    }
}
