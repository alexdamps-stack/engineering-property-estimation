//! Editor state: the flowsheet model plus GUI-only layout, viewport, and selection.

use std::collections::HashMap;

use eframe::egui::{Pos2, Rect, Vec2, pos2};
use flowsheet_core::{Flowsheet, PortRef, StreamId, UnitId, UnitKind};

use crate::units::DisplayUnits;
use crate::view::{Viewport, block_rect};

/// The selected canvas element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Selection {
    Unit(UnitId),
    Stream(StreamId),
}

/// An in-progress pointer drag on the canvas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Drag {
    /// Panning the view.
    Panning,
    /// Moving a unit; `grab` is the pointer offset from the unit centre (world).
    Moving { unit: UnitId, grab: Vec2 },
    /// Drawing a new stream from a free port.
    Connecting { from: PortRef },
}

/// Complete editor state.
pub struct Editor {
    pub fs: Flowsheet,
    /// World-space centre of each unit block. Layout is GUI-only data.
    pub layout: HashMap<UnitId, Pos2>,
    pub viewport: Viewport,
    pub selection: Option<Selection>,
    pub drag: Option<Drag>,
    pub units: DisplayUnits,
    /// Fit the view to the content on the next canvas frame.
    pub fit_requested: bool,
    /// Last user-facing status message (e.g. a rejected connection).
    pub status: Option<String>,
    /// Screen rectangle of the canvas in the last frame (used by UI tests to
    /// convert drawing coordinates to pointer positions).
    pub canvas_rect: Rect,
}

impl Default for Editor {
    fn default() -> Self {
        Self {
            fs: Flowsheet::new(),
            layout: HashMap::new(),
            viewport: Viewport::default(),
            selection: None,
            drag: None,
            units: DisplayUnits::default(),
            fit_requested: false,
            status: None,
            canvas_rect: Rect::NOTHING,
        }
    }
}

impl Editor {
    /// Replaces the flowsheet with an empty one, keeping display preferences.
    pub fn clear(&mut self) {
        *self = Self {
            units: self.units,
            ..Self::default()
        };
    }

    /// Adds a unit at `center` (world) and selects it.
    pub fn add_unit_at(&mut self, kind: UnitKind, center: Pos2) -> UnitId {
        let id = self.fs.add_unit(kind);
        self.layout.insert(id, center);
        self.selection = Some(Selection::Unit(id));
        id
    }

    /// Adds a unit near `hint` (world), offset so it does not cover an existing block.
    pub fn add_unit_near(&mut self, kind: UnitKind, hint: Pos2) -> UnitId {
        let mut center = pos2(
            (hint.x / 10.0).round() * 10.0,
            (hint.y / 10.0).round() * 10.0,
        );
        while self.layout.values().any(|p| p.distance(center) < 30.0) {
            center += Vec2::splat(30.0);
        }
        self.add_unit_at(kind, center)
    }

    /// Whether a new stream between the two ports would be accepted, in either order.
    pub fn can_connect(&self, a: PortRef, b: PortRef) -> bool {
        self.oriented(a, b).is_some_and(|(from, to)| {
            from.unit != to.unit
                && self.fs.stream_at(from).is_none()
                && self.fs.stream_at(to).is_none()
        })
    }

    /// Orders two ports as (outlet, inlet), if they have opposite directions.
    fn oriented(&self, a: PortRef, b: PortRef) -> Option<(PortRef, PortRef)> {
        use flowsheet_core::PortDirection::{Inlet, Outlet};
        match (
            self.fs.port_spec(a)?.direction,
            self.fs.port_spec(b)?.direction,
        ) {
            (Outlet, Inlet) => Some((a, b)),
            (Inlet, Outlet) => Some((b, a)),
            _ => None,
        }
    }

    /// Connects two ports dragged in either direction; reports rejections in
    /// the status line and selects the new stream on success.
    pub fn connect(&mut self, a: PortRef, b: PortRef) -> Option<StreamId> {
        if a == b {
            return None;
        }
        let Some((from, to)) = self.oriented(a, b) else {
            self.status = Some("Connect an outlet to an inlet.".into());
            return None;
        };
        match self.fs.connect(from, to) {
            Ok(id) => {
                self.status = None;
                self.selection = Some(Selection::Stream(id));
                Some(id)
            }
            Err(e) => {
                self.status = Some(format!("Cannot connect: {e}."));
                None
            }
        }
    }

    /// Deletes the selected unit (with its streams) or stream.
    pub fn delete_selection(&mut self) {
        match self.selection.take() {
            Some(Selection::Unit(id)) => {
                self.fs.remove_unit(id);
                self.layout.remove(&id);
            }
            Some(Selection::Stream(id)) => {
                self.fs.disconnect(id);
            }
            None => {}
        }
    }

    /// World-space bounding box of all blocks.
    pub fn content_bounds(&self) -> Option<Rect> {
        self.fs
            .units()
            .filter_map(|u| Some(block_rect(u.kind, *self.layout.get(&u.id)?)))
            .reduce(|a, b| a.union(b))
    }

    /// Fits the viewport to the content (or resets it when empty).
    pub fn fit_to_content(&mut self, canvas: Rect) {
        match self.content_bounds() {
            Some(bounds) => self.viewport.fit(canvas, bounds.expand(30.0)),
            None => self.viewport = Viewport::default(),
        }
    }

    /// Replaces the flowsheet with the topology of the planned demonstration:
    /// mix two gas streams, cool, reduce pressure, and separate the phases
    /// (development plan §3; proposal section 4.1). Feed conditions, cooler
    /// outlet temperature, and valve outlet pressure are left unset because
    /// the demonstration process conditions are specified under A2.
    pub fn load_demonstration(&mut self) {
        self.clear();
        let feed1 = self.add_unit_at(UnitKind::Feed, pos2(-100.0, 0.0));
        let feed2 = self.add_unit_at(UnitKind::Feed, pos2(-100.0, 120.0));
        let mix = self.add_unit_at(UnitKind::Mixer, pos2(160.0, 60.0));
        let cool = self.add_unit_at(UnitKind::Cooler, pos2(300.0, 60.0));
        let valve = self.add_unit_at(UnitKind::Valve, pos2(440.0, 60.0));
        let sep = self.add_unit_at(UnitKind::Separator, pos2(580.0, 60.0));
        let vap = self.add_unit_at(UnitKind::Product, pos2(740.0, -60.0));
        let liq = self.add_unit_at(UnitKind::Product, pos2(740.0, 180.0));
        self.fs.rename_unit(vap, "Vapor");
        self.fs.rename_unit(liq, "Liquid");
        let links = [
            ((feed1, 0), (mix, 0)),
            ((feed2, 0), (mix, 1)),
            ((mix, 2), (cool, 0)),
            ((cool, 1), (valve, 0)),
            ((valve, 1), (sep, 0)),
            ((sep, 1), (vap, 0)),
            ((sep, 2), (liq, 0)),
        ];
        for ((a, pa), (b, pb)) in links {
            // The fixed topology above is valid by construction.
            let _ = self.fs.connect(PortRef::new(a, pa), PortRef::new(b, pb));
        }
        self.selection = None;
        self.fit_requested = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flowsheet_core::DiagnosticKind;

    #[test]
    fn connect_accepts_either_drag_direction() {
        let mut ed = Editor::default();
        let feed = ed.add_unit_at(UnitKind::Feed, pos2(-100.0, 0.0));
        let valve = ed.add_unit_at(UnitKind::Valve, pos2(200.0, 0.0));
        // Dragged from the inlet back to the outlet.
        let id = ed
            .connect(PortRef::new(valve, 0), PortRef::new(feed, 0))
            .expect("connected");
        let s = ed.fs.stream(id).unwrap();
        assert_eq!(s.from, PortRef::new(feed, 0));
        assert_eq!(s.to, PortRef::new(valve, 0));
        assert_eq!(ed.selection, Some(Selection::Stream(id)));
    }

    #[test]
    fn rejected_connection_sets_status() {
        let mut ed = Editor::default();
        let a = ed.add_unit_at(UnitKind::Feed, pos2(-100.0, 0.0));
        let b = ed.add_unit_at(UnitKind::Feed, pos2(200.0, 0.0));
        assert!(!ed.can_connect(PortRef::new(a, 0), PortRef::new(b, 0)));
        assert!(ed.connect(PortRef::new(a, 0), PortRef::new(b, 0)).is_none());
        assert!(ed.status.is_some());
        assert_eq!(ed.fs.streams().count(), 0);
    }

    #[test]
    fn deleting_a_unit_removes_layout_and_streams() {
        let mut ed = Editor::default();
        let feed = ed.add_unit_at(UnitKind::Feed, pos2(-100.0, 0.0));
        let product = ed.add_unit_at(UnitKind::Product, pos2(200.0, 0.0));
        ed.connect(PortRef::new(feed, 0), PortRef::new(product, 0));
        ed.selection = Some(Selection::Unit(feed));
        ed.delete_selection();
        assert!(ed.fs.unit(feed).is_none());
        assert!(!ed.layout.contains_key(&feed));
        assert_eq!(ed.fs.streams().count(), 0);
        assert_eq!(
            ed.fs.unit_status(product),
            Some(flowsheet_core::UnitStatus::Underspecified)
        );
    }

    #[test]
    fn add_unit_near_avoids_stacking() {
        let mut ed = Editor::default();
        let a = ed.add_unit_near(UnitKind::Heater, pos2(100.0, 100.0));
        let b = ed.add_unit_near(UnitKind::Heater, pos2(100.0, 100.0));
        assert_ne!(ed.layout[&a], ed.layout[&b]);
    }

    #[test]
    fn demonstration_is_fully_connected_with_specs_still_open() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        assert_eq!(ed.fs.units().count(), 8);
        assert_eq!(ed.fs.streams().count(), 7);
        let diagnostics = ed.fs.diagnostics();
        assert!(!diagnostics.is_empty());
        assert!(
            diagnostics
                .iter()
                .all(|d| matches!(d.kind, DiagnosticKind::MissingSpecification(..)))
        );
        assert!(ed.content_bounds().is_some());
    }
}
