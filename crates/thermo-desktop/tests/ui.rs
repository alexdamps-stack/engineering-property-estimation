//! Headless UI tests for the flowsheet editor.
//!
//! The whole application runs inside `egui_kittest`: no window is opened and
//! the real mouse and keyboard are never touched. Interaction tests inject
//! pointer and key events; snapshot tests render the UI off-screen with wgpu
//! and compare it with the PNG files in `tests/snapshots/`.
//!
//! Update snapshots after an intended visual change with
//! `UPDATE_SNAPSHOTS=1 cargo test -p thermo-desktop --test ui` and review the
//! changed PNGs before committing. On CI (`CI` set) the snapshots are rendered
//! but not compared pixel by pixel, because the baselines are recorded on a
//! developer machine and software renderers on CI differ slightly.

use eframe::egui::{self, Event, Key, Modifiers, PointerButton, Pos2, Vec2, vec2};
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;
use flowsheet_core::{PortRef, PressureSpec, StreamId, UnitId, UnitKind, UnitSpec};
use thermo_desktop::FlowsheetApp;
use thermo_desktop::editor::Selection;
use thermo_desktop::view::{port_exit, port_pos, stream_route};

const WINDOW: Vec2 = vec2(1280.0, 800.0);

fn harness(app: FlowsheetApp, theme: egui::Theme) -> Harness<'static, FlowsheetApp> {
    let mut h = Harness::builder()
        .with_size(WINDOW)
        .with_theme(theme)
        .build_eframe(move |_cc| app);
    h.run();
    h
}

fn unit_of_kind(h: &Harness<FlowsheetApp>, kind: UnitKind) -> UnitId {
    h.state()
        .editor()
        .fs
        .units()
        .find(|u| u.kind == kind)
        .map(|u| u.id)
        .expect("unit of kind")
}

fn to_screen(h: &Harness<FlowsheetApp>, world: Pos2) -> Pos2 {
    let ed = h.state().editor();
    ed.viewport.to_screen(ed.canvas_rect, world)
}

fn nozzle(h: &Harness<FlowsheetApp>, port: PortRef) -> Pos2 {
    let ed = h.state().editor();
    let kind = ed.fs.unit(port.unit).expect("unit").kind;
    to_screen(h, port_pos(kind, ed.layout[&port.unit], port.port))
}

fn press(h: &mut Harness<FlowsheetApp>, pos: Pos2, pressed: bool) {
    h.event(Event::PointerButton {
        pos,
        button: PointerButton::Primary,
        pressed,
        modifiers: Modifiers::default(),
    });
    h.step();
}

fn drag(h: &mut Harness<FlowsheetApp>, from: Pos2, to: Pos2) {
    h.event(Event::PointerMoved(from));
    h.step();
    press(h, from, true);
    for i in 1..=10 {
        h.event(Event::PointerMoved(from + (to - from) * (i as f32 / 10.0)));
        h.step();
    }
    press(h, to, false);
    h.run();
}

fn click_at(h: &mut Harness<FlowsheetApp>, pos: Pos2) {
    h.event(Event::PointerMoved(pos));
    h.step();
    press(h, pos, true);
    press(h, pos, false);
    h.run();
}

fn key(h: &mut Harness<FlowsheetApp>, key: Key) {
    for pressed in [true, false] {
        h.event(Event::Key {
            key,
            physical_key: None,
            pressed,
            repeat: false,
            modifiers: Modifiers::default(),
        });
    }
    h.run();
}

/// Serialises GPU rendering: concurrent wgpu renders from parallel test
/// threads crashed the Windows driver in a local run.
static RENDER: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Pixel comparison locally; render-only on CI (see module docs).
fn check_snapshot(h: &mut Harness<FlowsheetApp>, name: &str) {
    let _gpu = RENDER.lock().unwrap_or_else(|e| e.into_inner());
    if std::env::var_os("CI").is_some() {
        h.render().expect("headless render");
    } else {
        h.snapshot(name);
    }
}

fn two_units(a: UnitKind, b: UnitKind) -> FlowsheetApp {
    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    ed.add_unit_at(a, egui::pos2(0.0, 0.0));
    ed.add_unit_at(b, egui::pos2(200.0, 0.0));
    ed.selection = None;
    app
}

#[test]
fn palette_click_places_and_selects_a_unit() {
    let mut h = harness(FlowsheetApp::default(), egui::Theme::Dark);
    h.get_by_label("Valve").click();
    h.run();
    let ed = h.state().editor();
    let units: Vec<_> = ed.fs.units().collect();
    assert_eq!(units.len(), 1);
    assert_eq!(units[0].kind, UnitKind::Valve);
    assert_eq!(ed.selection, Some(Selection::Unit(units[0].id)));
}

#[test]
fn dragging_from_outlet_to_inlet_creates_a_stream() {
    let mut h = harness(
        two_units(UnitKind::Feed, UnitKind::Product),
        egui::Theme::Dark,
    );
    let feed = unit_of_kind(&h, UnitKind::Feed);
    let product = unit_of_kind(&h, UnitKind::Product);
    let (a, b) = (PortRef::new(feed, 0), PortRef::new(product, 0));
    let (pa, pb) = (nozzle(&h, a), nozzle(&h, b));
    drag(&mut h, pa, pb);
    let ed = h.state().editor();
    let s = ed.fs.stream_at(a).expect("stream created");
    assert_eq!((s.from, s.to), (a, b));
}

#[test]
fn dragging_from_inlet_back_to_outlet_also_connects() {
    let mut h = harness(
        two_units(UnitKind::Feed, UnitKind::Valve),
        egui::Theme::Dark,
    );
    let feed = unit_of_kind(&h, UnitKind::Feed);
    let valve = unit_of_kind(&h, UnitKind::Valve);
    let (out, inlet) = (PortRef::new(feed, 0), PortRef::new(valve, 0));
    let (pa, pb) = (nozzle(&h, inlet), nozzle(&h, out));
    drag(&mut h, pa, pb);
    let s = h.state().editor().fs.stream_at(inlet).cloned();
    assert_eq!(s.map(|s| (s.from, s.to)), Some((out, inlet)));
}

#[test]
fn dropping_on_empty_canvas_creates_nothing() {
    let mut h = harness(
        two_units(UnitKind::Feed, UnitKind::Product),
        egui::Theme::Dark,
    );
    let feed = unit_of_kind(&h, UnitKind::Feed);
    let start = nozzle(&h, PortRef::new(feed, 0));
    drag(&mut h, start, start + vec2(60.0, 150.0));
    assert_eq!(h.state().editor().fs.streams().count(), 0);
}

#[test]
fn moving_a_unit_keeps_it_on_the_snap_grid() {
    let mut h = harness(
        two_units(UnitKind::Mixer, UnitKind::Product),
        egui::Theme::Dark,
    );
    let mixer = unit_of_kind(&h, UnitKind::Mixer);
    let center = to_screen(&h, h.state().editor().layout[&mixer]);
    drag(&mut h, center, center + vec2(53.0, 27.0));
    let ed = h.state().editor();
    let p = ed.layout[&mixer];
    assert_ne!(p, egui::pos2(0.0, 0.0), "unit should have moved");
    assert_eq!((p.x % 10.0, p.y % 10.0), (0.0, 0.0), "snapped: {p:?}");
    assert_eq!(ed.selection, Some(Selection::Unit(mixer)));
}

#[test]
fn clicking_a_stream_selects_it_and_delete_removes_it() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    let ed = h.state().editor();
    let cooler = ed
        .fs
        .units()
        .find(|u| u.kind == UnitKind::Cooler)
        .unwrap()
        .id;
    let stream: StreamId = ed.fs.stream_at(PortRef::new(cooler, 0)).unwrap().id;
    let s = ed.fs.stream(stream).unwrap();
    let from_kind = ed.fs.unit(s.from.unit).unwrap().kind;
    let route = stream_route(
        port_pos(from_kind, ed.layout[&s.from.unit], s.from.port),
        port_exit(from_kind, s.from.port),
        port_pos(UnitKind::Cooler, ed.layout[&cooler], 0),
    );
    // A point a little way along the first run, clear of nozzles and tags.
    let on_line = route[0] + (route[1] - route[0]) * 0.25;
    let target = to_screen(&h, on_line);

    click_at(&mut h, target);
    assert_eq!(
        h.state().editor().selection,
        Some(Selection::Stream(stream))
    );
    key(&mut h, Key::Delete);
    let ed = h.state().editor();
    assert!(ed.fs.stream(stream).is_none());
    assert_eq!(ed.fs.streams().count(), 6);
}

#[test]
fn valve_specification_is_edited_from_the_properties_panel() {
    let mut app = FlowsheetApp::default();
    let valve = app
        .editor_mut()
        .add_unit_at(UnitKind::Valve, egui::pos2(0.0, 0.0));
    let mut h = harness(app, egui::Theme::Dark);

    h.get_by_label("Enter value").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(valve).unwrap().spec,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(101_325.0))
        }
    );

    h.get_by_label("Pressure drop").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(valve).unwrap().spec,
        UnitSpec::Valve {
            pressure: PressureSpec::PressureDrop(None)
        }
    );
}

#[test]
fn snapshot_demonstration_dark() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    check_snapshot(&mut h, "demonstration_dark");
}

#[test]
fn snapshot_demonstration_light() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Light);
    check_snapshot(&mut h, "demonstration_light");
}

#[test]
fn snapshot_selected_cooler_properties() {
    let mut app = FlowsheetApp::with_demonstration();
    let cooler = app
        .editor()
        .fs
        .units()
        .find(|u| u.kind == UnitKind::Cooler)
        .unwrap()
        .id;
    app.editor_mut().selection = Some(Selection::Unit(cooler));
    let mut h = harness(app, egui::Theme::Dark);
    check_snapshot(&mut h, "cooler_selected_dark");
}

#[test]
fn every_palette_entry_places_its_unit_kind() {
    let mut h = harness(FlowsheetApp::default(), egui::Theme::Dark);
    for kind in UnitKind::ALL {
        h.get_by_label(kind.label()).click();
        h.run();
        let ed = h.state().editor();
        let Some(Selection::Unit(id)) = ed.selection else {
            panic!("{kind}: nothing selected after placing");
        };
        assert_eq!(ed.fs.unit(id).map(|u| u.kind), Some(kind));
    }
    assert_eq!(h.state().editor().fs.units().count(), UnitKind::ALL.len());
}

/// Splitter, pump, compressor, expander, and two-stream heat exchanger
/// wired into one small flowsheet.
fn equipment_gallery() -> FlowsheetApp {
    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    let at = egui::pos2;
    let feed1 = ed.add_unit_at(UnitKind::Feed, at(-100.0, 0.0));
    let tee = ed.add_unit_at(UnitKind::Splitter, at(20.0, 0.0));
    let pump = ed.add_unit_at(UnitKind::Pump, at(140.0, -60.0));
    let p1 = ed.add_unit_at(UnitKind::Product, at(300.0, -80.0));
    let comp = ed.add_unit_at(UnitKind::Compressor, at(140.0, 60.0));
    let exp = ed.add_unit_at(UnitKind::Expander, at(260.0, 60.0));
    let hx = ed.add_unit_at(UnitKind::HeatExchanger, at(420.0, 70.0));
    let feed2 = ed.add_unit_at(UnitKind::Feed, at(260.0, 160.0));
    let p2 = ed.add_unit_at(UnitKind::Product, at(580.0, 60.0));
    let p3 = ed.add_unit_at(UnitKind::Product, at(580.0, 140.0));
    for ((a, pa), (b, pb)) in [
        ((feed1, 0), (tee, 0)),
        ((tee, 1), (pump, 0)),
        ((pump, 1), (p1, 0)),
        ((tee, 2), (comp, 0)),
        ((comp, 1), (exp, 0)),
        ((exp, 1), (hx, 0)),
        ((feed2, 0), (hx, 1)),
        ((hx, 2), (p2, 0)),
        ((hx, 3), (p3, 0)),
    ] {
        ed.connect(PortRef::new(a, pa), PortRef::new(b, pb))
            .expect("gallery connection");
    }
    ed.selection = None;
    ed.fit_requested = true;
    app
}

#[test]
fn equipment_gallery_is_fully_connected() {
    let app = equipment_gallery();
    let fs = &app.editor().fs;
    assert_eq!(fs.streams().count(), 9);
    assert!(
        !fs.diagnostics()
            .iter()
            .any(|d| matches!(d.kind, flowsheet_core::DiagnosticKind::UnconnectedPort(_)))
    );
}

#[test]
fn snapshot_equipment_gallery() {
    let mut h = harness(equipment_gallery(), egui::Theme::Dark);
    check_snapshot(&mut h, "equipment_gallery_dark");
}
