//! Headless UI tests for the flowsheet editor.
//!
//! The whole application runs inside `egui_kittest`: no window is opened and
//! the real mouse and keyboard are never touched. Interaction tests inject
//! pointer and key events; snapshot tests render the UI off-screen with wgpu
//! and compare it with the PNG files in `tests/snapshots/`.
//!
//! Update snapshots after an intended visual change with
//! `UPDATE_SNAPSHOTS=1 cargo test -p thermo-desktop --test ui` and review the
//! changed PNGs before committing. The snapshots are compared everywhere,
//! including CI (Mesa lavapipe on Linux; egui_kittest prefers a CPU adapter
//! on every platform), with a per-pixel
//! tolerance and a small allowed count of differing pixels that absorb
//! anti-aliasing differences between renderers (see [`PIXEL_THRESHOLD`] and
//! [`MAX_FAILED_PIXELS`]). A failing comparison writes `<name>.new.png` and
//! `<name>.diff.png` next to the baseline; CI uploads them as an artifact.
//! Snapshot content must not depend on solver round-off: values that are
//! round-off by nature, such as balance residuals, are displayed as a bound.

use eframe::egui::{self, Event, Key, Modifiers, PointerButton, Pos2, Vec2, vec2};
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;
use flowsheet_core::{
    PaletteCategory, PortRef, PressureSpec, StreamId, UnitId, UnitKind, UnitSpec,
};
use thermo_desktop::FlowsheetApp;
use thermo_desktop::editor::{ItemStatus, Selection};
use thermo_desktop::view::{port_exit, port_pos, stream_route};

const WINDOW: Vec2 = vec2(1280.0, 800.0);

#[test]
fn recycle_demo_solves_invalidates_reloads_and_reports_failure() {
    let mut h = harness(
        FlowsheetApp::with_recycle_demonstration(),
        egui::Theme::Dark,
    );
    assert!(h.state_mut().editor_mut().solve_blocker().is_none());
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    assert!(
        statuses(&mut h)
            .iter()
            .all(|(_, s)| *s == ItemStatus::Solved)
    );
    assert!(
        h.state()
            .editor()
            .results
            .as_ref()
            .unwrap()
            .recycle
            .as_ref()
            .unwrap()
            .max_scaled_residual
            <= 1e-10
    );
    let heater = unit_of_kind(&h, UnitKind::Heater);
    h.state_mut().editor_mut().fs.set_spec(
        heater,
        UnitSpec::HeatExchange {
            thermal: flowsheet_core::ThermalSpec::OutletTemperature(Some(360.0)),
            pressure_drop_pa: 0.0,
        },
    );
    let fps = h.state_mut().editor_mut().fingerprints();
    let ed = h.state().editor();
    for u in ed.fs.units() {
        assert_eq!(
            ed.unit_status(u.id, &fps, &[]),
            if u.kind == UnitKind::Feed {
                ItemStatus::Solved
            } else {
                ItemStatus::Stale
            }
        );
    }
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    assert!(
        statuses(&mut h)
            .iter()
            .all(|(_, s)| *s == ItemStatus::Solved)
    );
    let path = std::env::temp_dir().join(format!(
        "thermo-recycle-ui-{}.flowsheet.json",
        std::process::id()
    ));
    h.state_mut().editor_mut().save_to(&path).unwrap();
    h.state_mut().editor_mut().open_path(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert!(
        statuses(&mut h)
            .iter()
            .all(|(_, s)| *s == ItemStatus::Solved)
    );
    // A valid specification with an impossible enthalpy target must leave the
    // previous loop values stale and show the actual failing heater.
    h.state_mut().editor_mut().fs.set_spec(
        heater,
        UnitSpec::HeatExchange {
            thermal: flowsheet_core::ThermalSpec::Duty(Some(1e20)),
            pressure_drop_pa: 0.0,
        },
    );
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    assert_eq!(
        h.state().editor().failure.as_ref().unwrap().unit,
        Some(heater)
    );
    assert!(
        statuses(&mut h)
            .iter()
            .any(|(_, s)| *s == ItemStatus::Failed)
    );
    assert!(
        statuses(&mut h)
            .iter()
            .any(|(_, s)| *s == ItemStatus::Stale)
    );
}

#[test]
fn snapshot_solved_recycle_demonstration() {
    let mut app = FlowsheetApp::with_recycle_demonstration();
    solve_now(&mut app);
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Stream table").click();
    h.run();
    check_snapshot(&mut h, "recycle_solved_dark");
}

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

/// Per-pixel tolerance: the weighted squared YIQ distance from which a pixel
/// counts as differing. Renderer drift between the Windows baselines and
/// Mesa lavapipe stays below it except for a few glyph edge pixels, and all
/// drift is gone at 3.0.
const PIXEL_THRESHOLD: f32 = 0.6;

/// Pixels allowed to exceed [`PIXEL_THRESHOLD`]: 6–14 glyph edge pixels per
/// image were measured between the baselines and lavapipe, plus a small
/// margin. A real change such as a moved label or a changed value differs
/// by far more (the residual text change of #150 alone was 261 pixels).
const MAX_FAILED_PIXELS: usize = 20;

fn snapshot_options() -> egui_kittest::SnapshotOptions {
    egui_kittest::SnapshotOptions::new()
        .threshold(PIXEL_THRESHOLD)
        .max_failed_pixels(MAX_FAILED_PIXELS)
}

/// Renders and compares with the baseline (see module docs for tolerances).
fn check_snapshot(h: &mut Harness<FlowsheetApp>, name: &str) {
    let _gpu = RENDER.lock().unwrap_or_else(|e| e.into_inner());
    h.snapshot_options(name, &snapshot_options());
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
    // The grouped palette (A17) is taller than the window: scroll to the entry.
    h.get_by_label("Valve").scroll_to_me();
    h.run();
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
fn delete_needs_the_delete_key_over_the_canvas() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    let cooler = unit_of_kind(&h, UnitKind::Cooler);
    let center = to_screen(&h, h.state().editor().layout[&cooler]);
    click_at(&mut h, center);
    assert_eq!(h.state().editor().selection, Some(Selection::Unit(cooler)));
    // Backspace never deletes: it follows text edits in the properties panel.
    key(&mut h, Key::Backspace);
    assert!(h.state().editor().fs.unit(cooler).is_some());
    // Delete with the pointer over the properties panel does nothing either.
    h.event(Event::PointerMoved(egui::pos2(WINDOW.x - 20.0, 300.0)));
    h.step();
    key(&mut h, Key::Delete);
    assert!(h.state().editor().fs.unit(cooler).is_some());
    // Over the canvas, Delete removes the unit.
    h.event(Event::PointerMoved(center));
    h.step();
    key(&mut h, Key::Delete);
    assert!(h.state().editor().fs.unit(cooler).is_none());
}

#[test]
fn keys_for_an_open_dialog_do_not_edit_the_flowsheet_or_cancel_the_solve() {
    let mut h = harness(FlowsheetApp::with_demonstration(), egui::Theme::Dark);
    let cooler = unit_of_kind(&h, UnitKind::Cooler);
    let center = to_screen(&h, h.state().editor().layout[&cooler]);
    click_at(&mut h, center);
    // The demonstration is unsaved, so New asks about unsaved changes.
    h.get_by_label("New").click();
    h.step();
    h.state_mut().editor_mut().start_solve(None);
    h.step();
    assert!(h.state().editor().is_solving());
    h.event(Event::PointerMoved(center));
    h.step();
    key(&mut h, Key::Delete);
    key(&mut h, Key::Escape);
    let ed = h.state().editor();
    assert!(ed.fs.unit(cooler).is_some(), "Delete went to the canvas");
    assert_eq!(ed.selection, Some(Selection::Unit(cooler)));
    assert!(ed.is_solving(), "Esc for the dialog cancelled the solve");
    h.state_mut().editor_mut().cancel_solve();
}

#[test]
fn a_package_without_components_does_not_offer_a_composition() {
    let mut app = FlowsheetApp::default();
    app.editor_mut()
        .new_project(flowsheet_project::PackageSelection::zero_interactions(
            flowsheet_project::CubicVariantChoice::PengRobinsonBostonMathias,
            Vec::new(),
            "test",
        ));
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Feed").click();
    h.run();
    assert!(h.query_by_label("Enter composition").is_none());
    h.get_by_label_contains("no components");
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
fn component_splitter_recovery_is_edited_from_the_properties_panel() {
    use flowsheet_core::units::SimpleLogicalSpec;
    use flowsheet_core::units::simple_logical::ComponentSplitterSpec;
    use flowsheet_project::{CubicVariantChoice, PackageSelection};

    let mut app = FlowsheetApp::default();
    app.editor_mut()
        .new_project(PackageSelection::zero_interactions(
            CubicVariantChoice::PengRobinson,
            vec!["methane".into()],
            "single-component UI test; no interactions",
        ));
    let splitter = app
        .editor_mut()
        .add_unit_at(UnitKind::ComponentSplitter, egui::pos2(0.0, 0.0));
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("methane");
    h.get_by_label("Enter value").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(splitter).unwrap().spec,
        UnitSpec::SimpleLogical(SimpleLogicalSpec::ComponentSplitter(
            ComponentSplitterSpec {
                overhead_recovery: [("methane".into(), 0.5)].into(),
            }
        ))
    );
}

#[test]
fn logical_stream_specs_and_exchange_use_the_properties_panel() {
    use flowsheet_core::units::SimpleLogicalSpec;
    use flowsheet_core::units::simple_logical::{MultiplierSpec, SelectorSpec};
    use flowsheet_project::{CubicVariantChoice, PackageSelection};

    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    ed.new_project(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        vec!["methane".into()],
        "single-component UI test",
    ));
    let first = ed.add_unit_at(UnitKind::Feed, egui::pos2(0.0, 0.0));
    let second = ed.add_unit_at(UnitKind::Feed, egui::pos2(0.0, 170.0));
    let scale = ed.add_unit_at(UnitKind::Multiplier, egui::pos2(150.0, 0.0));
    let select = ed.add_unit_at(UnitKind::Selector, egui::pos2(320.0, 0.0));
    let product = ed.add_unit_at(UnitKind::Product, egui::pos2(480.0, 0.0));
    for (feed, flow) in [(first, 2.0), (second, 3.0)] {
        ed.fs.set_spec(
            feed,
            UnitSpec::Feed {
                temperature_k: Some(330.0),
                pressure_pa: Some(2e6),
                molar_flow_mol_s: Some(flow),
            },
        );
        ed.feeds.insert(feed, vec![1.0]);
    }
    for (from, out, to, inlet) in [
        (first, 0, scale, 0),
        (scale, 1, select, 1),
        (second, 0, select, 2),
        (select, 0, product, 0),
    ] {
        ed.fs
            .connect(PortRef::new(from, out), PortRef::new(to, inlet))
            .unwrap();
    }
    ed.selection = Some(Selection::Unit(scale));
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Enter value").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(scale).unwrap().spec,
        UnitSpec::SimpleLogical(SimpleLogicalSpec::Multiplier(MultiplierSpec {
            factor: Some(1.0)
        })),
    );
    h.state_mut().editor_mut().selection = Some(Selection::Unit(select));
    h.run();
    h.get_by_label("Selected inlet").click();
    h.run();
    h.get_by_label("Inlet 2").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(select).unwrap().spec,
        UnitSpec::SimpleLogical(SimpleLogicalSpec::Selector(SelectorSpec {
            inlets: 2,
            selected_inlet: Some(2),
        })),
    );
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    h.get_by_label("External material (+ into unit)");
    h.get_by_label("Enthalpy carried by material");
    let r = h.state().editor().results.as_ref().unwrap();
    let exchange = r.units[&select].material_exchange.as_ref().unwrap();
    assert_eq!(exchange.component_flows_mol_s, [-2.0]);
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
    let available = UnitKind::ALL
        .into_iter()
        .filter(|kind| kind.category() != PaletteCategory::Solids);
    for kind in available.clone() {
        // Only the kind's palette category is expanded, so it is in view.
        for category in PaletteCategory::ALL {
            h.state_mut()
                .set_palette_category_open(category, category == kind.category());
        }
        h.run();
        h.get_by_label(kind.label()).scroll_to_me();
        h.run();
        h.get_by_label(kind.label()).click();
        h.run();
        let ed = h.state().editor();
        let Some(Selection::Unit(id)) = ed.selection else {
            panic!("{kind}: nothing selected after placing");
        };
        assert_eq!(ed.fs.unit(id).map(|u| u.kind), Some(kind));
    }
    assert_eq!(h.state().editor().fs.units().count(), available.count());
}

#[test]
fn solids_are_absent_from_add_palette_and_existing_projects_still_open() {
    use flowsheet_project::{CubicVariantChoice, Project, demonstration};

    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    let existing = project.flowsheet.add_unit(UnitKind::Cyclone);
    project.layout.insert(existing, [160.0, 190.0]);
    let path = std::env::temp_dir().join(format!(
        "thermo-existing-solids-{}.json",
        std::process::id()
    ));
    project.save(&path).unwrap();
    let mut h = harness(FlowsheetApp::with_file(path.clone()), egui::Theme::Dark);
    for category in PaletteCategory::ALL {
        h.state_mut()
            .set_palette_category_open(category, category == PaletteCategory::Solids);
    }
    h.run();
    assert!(h.query_by_label("Solids").is_none());
    for kind in PaletteCategory::Solids.kinds() {
        assert!(h.query_by_label(kind.label()).is_none(), "{kind}");
    }
    let ed = h.state_mut().editor_mut();
    let unit = ed.fs.unit(existing).unwrap();
    assert_eq!(unit.kind, UnitKind::Cyclone);
    assert_eq!(
        unit.ports(),
        project.flowsheet.unit(existing).unwrap().ports()
    );
    ed.save_to(&path).unwrap();
    let reopened = Project::load(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert!(
        reopened
            .flowsheet
            .same_calculation_inputs(&project.flowsheet)
    );
}

#[test]
fn tank_thermal_mode_is_editable_and_duty_is_reported() {
    use flowsheet_core::SeparatorThermal;
    use flowsheet_core::units::{SimpleLogicalSpec, simple_logical::TankSpec};
    use flowsheet_project::{CubicVariantChoice, PackageSelection};

    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    ed.new_project(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        vec!["methane".into()],
        "single-component tank UI test",
    ));
    let feed = ed.add_unit_at(UnitKind::Feed, egui::pos2(0.0, 0.0));
    let product = ed.add_unit_at(UnitKind::Product, egui::pos2(300.0, 0.0));
    let tank = ed.add_unit_at(UnitKind::Tank, egui::pos2(150.0, 0.0));
    ed.fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(2e6),
            molar_flow_mol_s: Some(2.0),
        },
    );
    ed.feeds.insert(feed, vec![1.0]);
    ed.fs
        .connect(PortRef::new(feed, 0), PortRef::new(tank, 0))
        .unwrap();
    let outlet = ed
        .fs
        .connect(PortRef::new(tank, 1), PortRef::new(product, 0))
        .unwrap();
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Outlet T").click();
    h.run();
    h.get_by_label("Enter value").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(tank).unwrap().spec,
        UnitSpec::SimpleLogical(SimpleLogicalSpec::Tank(TankSpec {
            thermal: SeparatorThermal::OutletTemperature(Some(298.15)),
            pressure_drop_pa: 0.0,
        })),
    );
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    h.get_by_label("Duty (+ into unit)");
    let r = h.state().editor().results.as_ref().unwrap();
    assert_eq!(
        r.streams[&outlet].state.as_ref().unwrap().temperature_k,
        298.15
    );
    assert!(r.units[&tank].duty_w < 0.0);
}

#[test]
fn heat_transport_forms_solve_and_report_process_duty_and_fan_power() {
    use flowsheet_core::units::HeatTransportSpec;
    use flowsheet_core::units::heat_transport::{AirCoolerSpec, FiredHeaterSpec};
    use flowsheet_core::{EnergySourceDuty, ThermalSpec};
    use flowsheet_project::{CubicVariantChoice, PackageSelection};

    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    ed.new_project(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        vec!["methane".into()],
        "single-component heat-transport UI test",
    ));
    let feed = ed.add_unit_at(UnitKind::Feed, egui::pos2(30.0, 150.0));
    let source = ed.add_unit_at(UnitKind::EnergySource, egui::pos2(210.0, 50.0));
    let furnace = ed.add_unit_at(UnitKind::FiredHeater, egui::pos2(210.0, 150.0));
    let cooler = ed.add_unit_at(UnitKind::AirCooler, egui::pos2(360.0, 150.0));
    let product = ed.add_unit_at(UnitKind::Product, egui::pos2(510.0, 150.0));
    ed.fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(320.0),
            pressure_pa: Some(2e6),
            molar_flow_mol_s: Some(2.0),
        },
    );
    ed.feeds.insert(feed, vec![1.0]);
    ed.fs.set_spec(
        source,
        UnitSpec::EnergySource {
            duty: EnergySourceDuty::Specified(Some(1_000.0)),
        },
    );
    ed.fs.set_spec(
        cooler,
        UnitSpec::HeatTransport(HeatTransportSpec::AirCooler(AirCoolerSpec {
            ambient_temperature_k: Some(280.0),
            air_mass_flow_kg_s: Some(2.0),
            fan_power_w: 75.0,
            ..Default::default()
        })),
    );
    for ((a, ap), (b, bp)) in [
        ((feed, 0), (furnace, 0)),
        ((source, 0), (furnace, 2)),
        ((furnace, 1), (cooler, 0)),
        ((cooler, 1), (product, 0)),
    ] {
        ed.fs
            .connect(PortRef::new(a, ap), PortRef::new(b, bp))
            .unwrap();
    }
    ed.selection = Some(Selection::Unit(furnace));

    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Energy stream").click();
    h.run();
    assert_eq!(
        h.state().editor().fs.unit(furnace).unwrap().spec,
        UnitSpec::HeatTransport(HeatTransportSpec::FiredHeater(FiredHeaterSpec {
            thermal: ThermalSpec::EnergyStream,
            ..Default::default()
        })),
    );
    h.state_mut().editor_mut().selection = Some(Selection::Unit(cooler));
    h.run();
    h.get_by_label("Dry-air mass flow");
    h.get_by_label("Fan power (≥ 0)");
    h.get_by_label("Enter value").click();
    h.run();
    key(&mut h, Key::F5);
    finish_solve(&mut h);
    h.get_by_label("Process duty (+ into fluid)");
    h.get_by_label("Air discharge temperature");
    h.get_by_label("Fan power");
    let r = h.state().editor().results.as_ref().unwrap();
    assert!(r.units[&cooler].duty_w < 0.0);
    assert_eq!(r.units[&furnace].duty_w, 1_000.0);
    h.state_mut().editor_mut().selection = Some(Selection::Unit(furnace));
    h.run();
    h.get_by_label("Duty (+ into unit)");
}

#[test]
fn multi_stream_exchanger_solves_a_fixed_terminal_approach_from_properties() {
    use flowsheet_core::units::HeatTransportSpec;
    use flowsheet_core::units::heat_transport::{ExchangerPassThermal, MultiStreamExchangerSpec};
    use flowsheet_project::{CubicVariantChoice, PackageSelection};

    let mut app = FlowsheetApp::default();
    let ed = app.editor_mut();
    ed.new_project(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        vec!["methane".into()],
        "multi-stream exchanger UI test",
    ));
    let exchanger = ed.add_unit_at(UnitKind::MultiStreamExchanger, egui::pos2(220.0, 150.0));
    let mut spec = MultiStreamExchangerSpec::with_pass_count(2);
    spec.passes[0].thermal = ExchangerPassThermal::Balance;
    spec.passes[1].thermal = ExchangerPassThermal::Duty(Some(0.0));
    ed.fs.set_spec(
        exchanger,
        UnitSpec::HeatTransport(HeatTransportSpec::MultiStreamExchanger(spec.clone())),
    );
    let mut cold_outlet = None;
    for (i, (t, n, y)) in [(390.0, 3.0, 120.0), (300.0, 1.0, 220.0)]
        .into_iter()
        .enumerate()
    {
        let feed = ed.add_unit_at(UnitKind::Feed, egui::pos2(40.0, y));
        let product = ed.add_unit_at(UnitKind::Product, egui::pos2(420.0, y));
        ed.fs.set_spec(
            feed,
            UnitSpec::Feed {
                temperature_k: Some(t),
                pressure_pa: Some(2e6),
                molar_flow_mol_s: Some(n),
            },
        );
        ed.feeds.insert(feed, vec![1.0]);
        ed.fs
            .connect(PortRef::new(feed, 0), PortRef::new(exchanger, i * 2))
            .unwrap();
        let outlet = ed
            .fs
            .connect(PortRef::new(exchanger, i * 2 + 1), PortRef::new(product, 0))
            .unwrap();
        if i == 1 {
            cold_outlet = Some(outlet);
        }
    }
    ed.selection = Some(Selection::Unit(exchanger));

    let mut h = harness(app, egui::Theme::Dark);
    h.get_all_by_label("Fixed approach").nth(1).unwrap().click();
    h.run();
    h.get_by_label("Hot inlet − cold outlet");
    h.get_by_label("Enter value").click();
    h.run();
    spec.passes[1].thermal = ExchangerPassThermal::TerminalApproach {
        reference_pass: 1,
        approach_k: Some(0.0),
    };
    assert_eq!(
        h.state().editor().fs.unit(exchanger).unwrap().spec,
        UnitSpec::HeatTransport(HeatTransportSpec::MultiStreamExchanger(spec)),
    );

    key(&mut h, Key::F5);
    finish_solve(&mut h);
    h.get_by_label("Transferred duty (hot to cold)");
    h.get_by_label("Pass 1 duty (+ into fluid)");
    h.get_by_label("Pass 2 duty (+ into fluid)");
    h.get_by_label("Pass 2 terminal approach");
    let result = h.state().editor().results.as_ref().unwrap();
    assert!(result.units[&exchanger].transferred_duty_w.unwrap() > 0.0);
    assert_eq!(
        result.streams[&cold_outlet.unwrap()]
            .state
            .as_ref()
            .unwrap()
            .temperature_k,
        390.0
    );
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
    assert!(saved.contains(&format!(
        "\"format_version\": {}",
        flowsheet_project::FORMAT_VERSION
    )));
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

/// The demonstration feeds routed through a heat exchanger that moves 1 kW
/// from the 350 K feed (hot side) to the 330 K feed (cold side).
fn exchanger_app() -> (FlowsheetApp, UnitId) {
    let mut project =
        flowsheet_project::demonstration(flowsheet_project::CubicVariantChoice::PengRobinson);
    let fs = &mut project.flowsheet;
    let feed_at = |fs: &flowsheet_core::Flowsheet, t: f64| {
        fs.units()
            .find(|u| matches!(u.spec, UnitSpec::Feed { temperature_k: Some(x), .. } if x == t))
            .unwrap()
            .id
    };
    let (cold, hot) = (feed_at(fs, 330.0), feed_at(fs, 350.0));
    let others: Vec<_> = fs
        .units()
        .filter(|u| u.kind != UnitKind::Feed)
        .map(|u| u.id)
        .collect();
    for id in others {
        fs.remove_unit(id);
    }
    let ex = fs.add_unit(UnitKind::HeatExchanger);
    fs.set_spec(
        ex,
        UnitSpec::HeatExchanger {
            thermal: flowsheet_core::ExchangerSpec::Duty(Some(1000.0)),
            hot_pressure_drop_pa: 0.0,
            cold_pressure_drop_pa: 0.0,
        },
    );
    fs.connect(PortRef::new(hot, 0), PortRef::new(ex, 0))
        .unwrap();
    fs.connect(PortRef::new(cold, 0), PortRef::new(ex, 1))
        .unwrap();
    for port in [2, 3] {
        let product = fs.add_unit(UnitKind::Product);
        fs.connect(PortRef::new(ex, port), PortRef::new(product, 0))
            .unwrap();
    }
    let mut app = FlowsheetApp::default();
    app.editor_mut().load_project(project);
    (app, ex)
}

#[test]
fn exchanger_shows_transferred_duty_in_display_units() {
    let (mut app, ex) = exchanger_app();
    solve_now(&mut app);
    let stored = &app.editor().results.as_ref().unwrap().units[&ex];
    assert_eq!(
        (stored.duty_w, stored.transferred_duty_w),
        (0.0, Some(1000.0))
    );
    app.editor_mut().selection = Some(Selection::Unit(ex));
    let mut h = harness(app, egui::Theme::Dark);
    h.get_by_label("Transferred duty (hot to cold)");
    let shown = thermo_desktop::units::format_quantity(
        thermo_desktop::units::DisplayUnits::ProcessMetric,
        thermo_desktop::units::Quantity::Duty,
        1000.0,
    );
    assert!(shown.ends_with(" kW"), "{shown}");
    h.get_by_label(&shown);

    // A result without a stored transferred duty (read from a version-2
    // file) asks for a recalculation instead of showing a number.
    h.state_mut()
        .editor_mut()
        .results
        .as_mut()
        .unwrap()
        .units
        .get_mut(&ex)
        .unwrap()
        .transferred_duty_w = None;
    h.run();
    h.get_by_label("not stored — recalculate");
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
