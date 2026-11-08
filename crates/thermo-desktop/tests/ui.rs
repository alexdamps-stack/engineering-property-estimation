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
use thermo_desktop::editor::{ItemStatus, Selection};
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
    // Step instead of `run`: a started solve keeps requesting repaints.
    h.step();
    h.step();
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
    let hx = ed.add_unit_at(UnitKind::HeatExchanger, at(440.0, 70.0));
    let feed2 = ed.add_unit_at(UnitKind::Feed, at(260.0, 160.0));
    let p2 = ed.add_unit_at(UnitKind::Product, at(660.0, 60.0));
    let p3 = ed.add_unit_at(UnitKind::Product, at(660.0, 140.0));
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

/// Solves synchronously and pins the completion time, so snapshots do not
/// depend on the clock.
fn solve_now(app: &mut FlowsheetApp) {
    let ed = app.editor_mut();
    ed.start_solve(None);
    ed.wait_for_solve();
    if let Some(r) = &mut ed.results {
        r.solved_at_unix_s = 1_794_150_131;
    }
}

/// Steps the UI until the background solve has finished.
fn finish_solve(h: &mut Harness<FlowsheetApp>) {
    let start = std::time::Instant::now();
    while h.state().editor().is_solving() {
        assert!(
            start.elapsed() < std::time::Duration::from_secs(60),
            "solve did not finish"
        );
        h.step();
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    h.step();
    h.step();
}

fn statuses(h: &mut Harness<FlowsheetApp>) -> Vec<(String, ItemStatus)> {
    let ed = h.state_mut().editor_mut();
    let fps = ed.fingerprints();
    let issues = ed.input_issues();
    ed.fs
        .units()
        .map(|u| (u.name.clone(), ed.unit_status(u.id, &fps, &issues)))
        .collect()
}

/// Proposal section 4.1, steps 5–7, driven through the interface: solve,
/// change the valve pressure (downstream turns stale), solve again, save,
/// reopen and find the results current.
#[test]
fn demonstration_is_solved_modified_saved_and_reopened() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    h.get_by_label("Solve (F5)").click();
    h.step();
    finish_solve(&mut h);
    assert!(
        statuses(&mut h)
            .iter()
            .all(|(_, s)| *s == ItemStatus::Solved),
        "{:?}",
        statuses(&mut h)
    );

    let valve = unit_of_kind(&h, UnitKind::Valve);
    h.state_mut().editor_mut().selection = Some(Selection::Unit(valve));
    h.state_mut().editor_mut().fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(2.0e6)),
        },
    );
    h.step();
    for (name, s) in statuses(&mut h) {
        let downstream = ["V-100", "D-100", "Vapor", "Liquid"].contains(&name.as_str());
        let expected = if downstream {
            ItemStatus::Stale
        } else {
            ItemStatus::Solved
        };
        assert_eq!(s, expected, "{name}");
    }

    key(&mut h, Key::F5);
    finish_solve(&mut h);
    assert!(
        statuses(&mut h)
            .iter()
            .all(|(_, s)| *s == ItemStatus::Solved)
    );

    let dir = std::env::temp_dir().join(format!("thermo-desktop-ui-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("demo.flowsheet.json");
    h.state_mut().editor_mut().save_to(&path).unwrap();
    let mut reopened = harness(FlowsheetApp::with_file(path.clone()), egui::Theme::Dark);
    assert_eq!(
        reopened.state().editor().path.as_deref(),
        Some(path.as_path())
    );
    assert!(!reopened.state().editor().is_modified());
    assert!(
        statuses(&mut reopened)
            .iter()
            .all(|(_, s)| *s == ItemStatus::Solved)
    );
    let saved = std::fs::read_to_string(&path).unwrap();
    assert!(saved.contains("\"format_version\": 1"));
    assert!(saved.contains("\"outlet_pressure_pa\": 2000000.0"));
}

#[test]
fn failed_solve_is_shown_at_the_unit_and_values_are_not_current() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    let valve = unit_of_kind(&h, UnitKind::Valve);
    {
        let ed = h.state_mut().editor_mut();
        ed.fs.set_spec(
            valve,
            UnitSpec::Valve {
                pressure: PressureSpec::OutletPressure(Some(6.0e6)),
            },
        );
        ed.selection = Some(Selection::Unit(valve));
    }
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    let s = statuses(&mut h);
    assert!(s.contains(&("V-100".into(), ItemStatus::Failed)), "{s:?}");
    assert!(s.contains(&("D-100".into(), ItemStatus::Stale)), "{s:?}");
    assert!(s.contains(&("M-100".into(), ItemStatus::Solved)), "{s:?}");
    // The properties panel names the cause.
    let message = h.state().editor().failure.as_ref().unwrap().message.clone();
    h.get_by_label(&format!("✖ {message}"));
}

#[test]
fn cancel_stops_a_running_solve_without_changing_results() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    h.state_mut().editor_mut().start_solve(None);
    h.step();
    assert!(h.state().editor().is_solving());
    h.get_by_label("Cancel").click();
    h.step();
    let ed = h.state().editor();
    assert!(!ed.is_solving());
    assert!(ed.results.is_none());
    assert!(ed.status.as_deref().unwrap().contains("cancelled"));
}

#[test]
fn new_project_dialog_creates_a_package_and_feeds_take_a_composition() {
    let mut h = harness(FlowsheetApp::default(), egui::Theme::Dark);
    h.get_by_label("New").click();
    h.step();
    {
        let dialog = h.state_mut().new_project_dialog().expect("dialog open");
        for (id, on) in &mut dialog.components {
            *on = id == "methane" || id == "propane";
        }
    }
    h.step();
    h.get_by_label("Create").click();
    h.step();
    h.run();
    let ed = h.state().editor();
    let package = ed.package.as_ref().expect("package created");
    assert_eq!(package.components, ["methane", "propane"]);
    assert!(package.interaction_source.contains("zero-kij"));

    h.get_by_label("Feed").click();
    h.run();
    let feed = unit_of_kind(&h, UnitKind::Feed);
    h.get_by_label("Enter composition").click();
    h.run();
    assert_eq!(h.state().editor().feeds[&feed], vec![0.5, 0.5]);
}

#[test]
fn stream_table_copies_values_with_units_and_status() {
    let mut app = FlowsheetApp::with_demonstration();
    app.editor_mut().start_solve(None);
    app.editor_mut().wait_for_solve();
    let tsv = thermo_desktop::panels::stream_table_tsv(app.editor_mut());
    let lines: Vec<&str> = tsv.lines().collect();
    assert!(lines[0].starts_with("Stream\tUnit\tA\tB\tmix"));
    assert!(lines[1].starts_with("Status\t\tsolved"));
    assert!(lines[2].starts_with("Temperature\t°C\t56.85"));
    assert!(tsv.contains("\nx methane\tmol/mol\t0.87999"), "{tsv}");
    assert!(!tsv.contains("NaN"));
}

#[test]
fn snapshot_solved_demonstration_with_stream_table() {
    let mut app = FlowsheetApp::with_demonstration();
    solve_now(&mut app);
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Stream table").click();
    h.run();
    check_snapshot(&mut h, "demonstration_solved_dark");
}

#[test]
fn snapshot_stale_and_failed_after_valve_edit() {
    let mut app = FlowsheetApp::with_demonstration();
    solve_now(&mut app);
    let ed = app.editor_mut();
    let valve = ed
        .fs
        .units()
        .find(|u| u.kind == UnitKind::Valve)
        .unwrap()
        .id;
    let sep = ed
        .fs
        .units()
        .find(|u| u.kind == UnitKind::Separator)
        .unwrap()
        .id;
    ed.fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(6.0e6)),
        },
    );
    ed.start_solve(None);
    ed.wait_for_solve();
    ed.selection = Some(Selection::Unit(sep));
    let mut h = harness(app, egui::Theme::Light);
    check_snapshot(&mut h, "failed_valve_light");
}

/// Framework evaluation (#52): CPU time per frame of the whole window with a
/// 100-unit flowsheet (50 feed → heater chains plus products, 150 units and
/// 100 streams in total). Not a correctness gate; run explicitly with
/// `cargo test --release -p thermo-desktop --test ui -- --ignored --nocapture`.
#[test]
#[ignore = "performance measurement, not a correctness gate"]
fn canvas_frame_time_with_100_units() {
    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    for i in 0..50 {
        let y = i as f32 * 90.0;
        let feed = ed.add_unit_at(UnitKind::Feed, egui::pos2(0.0, y));
        let heater = ed.add_unit_at(UnitKind::Heater, egui::pos2(200.0, y));
        let product = ed.add_unit_at(UnitKind::Product, egui::pos2(400.0, y));
        ed.connect(PortRef::new(feed, 0), PortRef::new(heater, 0));
        ed.connect(PortRef::new(heater, 1), PortRef::new(product, 0));
    }
    ed.selection = None;
    ed.fit_requested = true;
    let mut h = harness(app, egui::Theme::Dark);
    for _ in 0..5 {
        h.step();
    }
    let frames = 100;
    let start = std::time::Instant::now();
    for i in 0..frames {
        // Move the pointer so every frame does hover hit-testing.
        h.event(Event::PointerMoved(egui::pos2(300.0 + i as f32, 300.0)));
        h.step();
    }
    let per_frame = start.elapsed() / frames;
    println!(
        "units={} streams={} cpu_frame_time={per_frame:?}",
        h.state().editor().fs.units().count(),
        h.state().editor().fs.streams().count()
    );
}
