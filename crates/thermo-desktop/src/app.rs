//! Application shell following the main-window layout of
//! `docs/flowsheet-and-ui.md`, section 3: toolbar, palette, canvas,
//! properties panel, bottom panel (messages and stream table), status bar.

use eframe::egui::{self, DragValue, Key, Panel, RichText, ScrollArea, Sense, Stroke, Ui, vec2};
use flowsheet_core::{
    DiagnosticKind, ExchangerSpec, MixerPressure, PortRef, PressureRise, PressureSpec,
    SeparatorThermal, SpecProblem, StreamId, ThermalSpec, UnitId, UnitKind, UnitSpec, UnitStatus,
};

use crate::canvas;
use crate::editor::{Editor, Selection};
use crate::symbols::{Theme, draw_unit, icon_rect};
use crate::units::{DisplayUnits, Quantity, starting_value};

const ENGINE_PENDING: &str =
    "The flowsheet engine is available. Desktop solving is pending integration.";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum BottomTab {
    #[default]
    Messages,
    StreamTable,
}

/// The thermo-desktop application.
#[derive(Default)]
pub struct FlowsheetApp {
    editor: Editor,
    bottom_tab: BottomTab,
}

impl FlowsheetApp {
    /// Starts with the demonstration flowsheet loaded.
    pub fn with_demonstration() -> Self {
        let mut app = Self::default();
        app.editor.load_demonstration();
        app
    }

    /// Editor state (for tests and embedding).
    pub fn editor(&self) -> &Editor {
        &self.editor
    }

    /// Mutable editor state (for tests and embedding).
    pub fn editor_mut(&mut self) -> &mut Editor {
        &mut self.editor
    }
}

impl eframe::App for FlowsheetApp {
    fn ui(&mut self, ui: &mut Ui, _frame: &mut eframe::Frame) {
        let ed = &mut self.editor;
        if !ui.ctx().egui_wants_keyboard_input() && ui.input(|i| i.key_pressed(Key::F5)) {
            ed.status = Some(ENGINE_PENDING.into());
        }
        Panel::top("toolbar").show(ui, |ui| toolbar(ui, ed));
        Panel::bottom("status_bar").show(ui, |ui| status_bar(ui, ed));
        Panel::bottom("bottom_panel")
            .resizable(true)
            .default_size(160.0)
            .show(ui, |ui| bottom_panel(ui, ed, &mut self.bottom_tab));
        Panel::left("palette")
            .resizable(false)
            .exact_size(175.0)
            .show(ui, |ui| palette(ui, ed));
        Panel::right("properties")
            .resizable(true)
            .default_size(300.0)
            .show(ui, |ui| properties(ui, ed));
        egui::CentralPanel::no_frame().show(ui, |ui| canvas::show(ui, ed));
    }
}

fn toolbar(ui: &mut Ui, ed: &mut Editor) {
    ui.horizontal(|ui| {
        if ui
            .button("New")
            .on_hover_text("Start an empty flowsheet")
            .clicked()
        {
            ed.clear();
        }
        if ui
            .button("Demonstration")
            .on_hover_text("Mix two gas streams, cool, reduce pressure, separate phases")
            .clicked()
        {
            ed.load_demonstration();
        }
        ui.separator();
        ui.add_enabled(false, egui::Button::new("Solve (F5)"))
            .on_disabled_hover_text(ENGINE_PENDING);
        ui.separator();
        if ui.button("Fit view").clicked() {
            ed.fit_requested = true;
        }
        ui.label(RichText::new(format!("{:.0} %", ed.viewport.zoom * 100.0)).weak());
        ui.separator();
        let label = ui.label("Units:");
        egui::ComboBox::from_id_salt("display_units")
            .selected_text(ed.units.label())
            .show_ui(ui, |ui| {
                for u in DisplayUnits::ALL {
                    ui.selectable_value(&mut ed.units, u, u.label());
                }
            })
            .response
            .labelled_by(label.id);
    });
}

fn status_bar(ui: &mut Ui, ed: &Editor) {
    let open = ed.fs.diagnostics().len();
    ui.horizontal(|ui| {
        ui.label(RichText::new("Package: none selected").weak())
            .on_hover_text("Component selection and PR / PR-BM packages arrive with A2–A3.");
        ui.separator();
        ui.label(RichText::new("Not solved").weak());
        ui.separator();
        ui.label(
            RichText::new(format!(
                "{} units, {} streams, {open} open item(s)",
                ed.fs.units().count(),
                ed.fs.streams().count()
            ))
            .weak(),
        );
        ui.separator();
        ui.label(RichText::new(format!("revision {}", ed.fs.revision())).weak());
    });
}

fn palette(ui: &mut Ui, ed: &mut Editor) {
    ScrollArea::vertical().show(ui, |ui| palette_items(ui, ed));
}

fn palette_items(ui: &mut Ui, ed: &mut Editor) {
    ui.add_space(6.0);
    ui.strong("Unit operations");
    ui.add_space(4.0);
    for kind in UnitKind::ALL {
        let (rect, response) =
            ui.allocate_exact_size(vec2(ui.available_width(), 36.0), Sense::click());
        let visuals = ui.style().interact(&response);
        if response.hovered() {
            ui.painter().rect_filled(rect, 4, visuals.weak_bg_fill);
        }
        let icon_area = egui::Rect::from_min_size(rect.min + vec2(8.0, 5.0), vec2(36.0, 26.0));
        draw_unit(
            ui.painter(),
            kind,
            icon_rect(icon_area, kind),
            Stroke::new(1.2, visuals.text_color()),
            egui::Color32::TRANSPARENT,
        );
        ui.painter().text(
            rect.left_center() + vec2(54.0, 0.0),
            egui::Align2::LEFT_CENTER,
            kind.label(),
            egui::FontId::proportional(14.0),
            visuals.text_color(),
        );
        let label = kind.label();
        response.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, label));
        if response
            .on_hover_text(format!("Place a {}", kind.label().to_lowercase()))
            .clicked()
        {
            let hint = ed
                .content_bounds()
                .map(|b| egui::pos2(b.right() + 100.0, b.center().y))
                .unwrap_or(egui::Pos2::ZERO);
            ed.add_unit_near(kind, hint);
        }
    }
    ui.add_space(12.0);
    ui.separator();
    ui.label(
        RichText::new(
            "Drag from an outlet to an inlet to connect.\n\
             Drag blocks to move them.\n\
             Drag the background or scroll to pan; Ctrl+scroll to zoom.\n\
             Del removes the selection.",
        )
        .small()
        .weak(),
    );
}

fn properties(ui: &mut Ui, ed: &mut Editor) {
    ScrollArea::vertical().show(ui, |ui| {
        ui.add_space(6.0);
        match ed.selection {
            Some(Selection::Unit(id)) if ed.fs.unit(id).is_some() => unit_properties(ui, ed, id),
            Some(Selection::Stream(id)) if ed.fs.stream(id).is_some() => {
                stream_properties(ui, ed, id)
            }
            _ => {
                ui.strong("Properties");
                ui.label(RichText::new("Select a unit or stream to edit it.").weak());
                ui.add_space(8.0);
                ui.label(
                    RichText::new(
                        "Components and the property package are chosen when a project is \
                         created; both arrive with A2–A3.",
                    )
                    .small()
                    .weak(),
                );
            }
        }
    });
}

fn status_line(ui: &mut Ui, status: Option<UnitStatus>) {
    let theme = Theme::from_visuals(ui.visuals());
    ui.horizontal(|ui| {
        ui.label("Status:");
        match status {
            Some(UnitStatus::Underspecified) => {
                ui.label(RichText::new("! Underspecified").color(theme.warning));
            }
            Some(UnitStatus::Ready) => {
                ui.label("Ready (not solved)");
            }
            None => {}
        }
    });
}

fn unit_properties(ui: &mut Ui, ed: &mut Editor, id: UnitId) {
    let Some(unit) = ed.fs.unit(id).cloned() else {
        return;
    };
    let theme = Theme::from_visuals(ui.visuals());
    let mut name = unit.name.clone();
    ui.horizontal(|ui| {
        let label = ui.label("Tag");
        if ui
            .add(egui::TextEdit::singleline(&mut name).desired_width(110.0))
            .labelled_by(label.id)
            .changed()
        {
            ed.fs.rename_unit(id, name);
        }
        ui.label(RichText::new(format!("({})", unit.kind.label())).weak());
    });
    status_line(ui, ed.fs.unit_status(id));

    ui.add_space(8.0);
    ui.strong("Specification");
    let mut spec = unit.spec;
    let u = ed.units;
    egui::Grid::new(("spec", id))
        .num_columns(2)
        .spacing([8.0, 6.0])
        .show(ui, |ui| spec_form(ui, &mut spec, u));
    ed.fs.set_spec(id, spec);
    for issue in spec.issues() {
        if issue.problem != SpecProblem::Missing {
            ui.colored_label(
                theme.warning,
                format!("! {} {}", issue.field, issue.problem),
            );
        }
    }
    if unit.kind == UnitKind::Feed {
        ui.label(
            RichText::new("Composition: pending component data (A2).")
                .small()
                .weak(),
        );
    }

    ui.add_space(8.0);
    ui.strong("Connections");
    egui::Grid::new(("ports", id))
        .num_columns(2)
        .show(ui, |ui| {
            for (i, port) in unit.kind.ports().iter().enumerate() {
                ui.label(port.name);
                match ed.fs.stream_at(PortRef::new(id, i)) {
                    Some(s) => {
                        let sid = s.id;
                        if ui.link(s.name.clone()).clicked() {
                            ed.selection = Some(Selection::Stream(sid));
                        }
                    }
                    None => {
                        ui.label(RichText::new("not connected").color(theme.warning));
                    }
                }
                ui.end_row();
            }
        });

    results_placeholder(ui);
    ui.add_space(12.0);
    if ui.button("Delete unit").clicked() {
        ed.delete_selection();
    }
}

/// Specification form; mutually exclusive choices use radio buttons
/// (proposal section 3).
fn spec_form(ui: &mut Ui, spec: &mut UnitSpec, u: DisplayUnits) {
    match spec {
        UnitSpec::Product => {
            ui.label(RichText::new("No specification.").weak());
            ui.end_row();
        }
        UnitSpec::Feed {
            temperature_k,
            pressure_pa,
            molar_flow_mol_s,
        } => {
            opt_field(ui, "Temperature", temperature_k, Quantity::Temperature, u);
            opt_field(ui, "Pressure", pressure_pa, Quantity::Pressure, u);
            opt_field(ui, "Molar flow", molar_flow_mol_s, Quantity::MolarFlow, u);
        }
        UnitSpec::Mixer { outlet_pressure } => {
            let current = matches!(outlet_pressure, MixerPressure::Specified(_)) as usize;
            match choice(ui, "Outlet p", current, &["Lowest inlet", "Specified"]) {
                Some(0) => *outlet_pressure = MixerPressure::LowestInlet,
                Some(_) => *outlet_pressure = MixerPressure::Specified(None),
                None => {}
            }
            if let MixerPressure::Specified(p) = outlet_pressure {
                opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u);
            }
            note_row(ui, "Adiabatic.");
        }
        UnitSpec::HeatExchange {
            thermal,
            pressure_drop_pa,
        } => {
            let current = matches!(thermal, ThermalSpec::Duty(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet T", "Duty"]) {
                Some(0) => *thermal = ThermalSpec::OutletTemperature(None),
                Some(_) => *thermal = ThermalSpec::Duty(None),
                None => {}
            }
            match thermal {
                ThermalSpec::OutletTemperature(t) => {
                    opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u)
                }
                ThermalSpec::Duty(q) => opt_field(ui, "Duty (+ into unit)", q, Quantity::Duty, u),
            }
            value_field(
                ui,
                "Pressure drop",
                pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
        }
        UnitSpec::Valve { pressure } => {
            let current = matches!(pressure, PressureSpec::PressureDrop(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet p", "Pressure drop"]) {
                Some(0) => *pressure = PressureSpec::OutletPressure(None),
                Some(_) => *pressure = PressureSpec::PressureDrop(None),
                None => {}
            }
            match pressure {
                PressureSpec::OutletPressure(p) => {
                    opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u)
                }
                PressureSpec::PressureDrop(dp) => {
                    opt_field(ui, "Pressure drop", dp, Quantity::PressureDrop, u)
                }
            }
            note_row(ui, "Isenthalpic.");
        }
        UnitSpec::Separator {
            thermal,
            pressure_drop_pa,
        } => {
            let current = match thermal {
                SeparatorThermal::Adiabatic => 0,
                SeparatorThermal::OutletTemperature(_) => 1,
                SeparatorThermal::Duty(_) => 2,
            };
            match choice(ui, "Thermal", current, &["Adiabatic", "Outlet T", "Duty"]) {
                Some(0) => *thermal = SeparatorThermal::Adiabatic,
                Some(1) => *thermal = SeparatorThermal::OutletTemperature(None),
                Some(_) => *thermal = SeparatorThermal::Duty(None),
                None => {}
            }
            match thermal {
                SeparatorThermal::Adiabatic => {}
                SeparatorThermal::OutletTemperature(t) => {
                    opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u)
                }
                SeparatorThermal::Duty(q) => {
                    opt_field(ui, "Duty (+ into unit)", q, Quantity::Duty, u)
                }
            }
            value_field(
                ui,
                "Pressure drop",
                pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
        }
        UnitSpec::Splitter { fraction_to_first } => {
            opt_field(
                ui,
                "Fraction to out 1",
                fraction_to_first,
                Quantity::Fraction,
                u,
            );
            note_row(ui, "Out 2 receives the rest; same state in both.");
        }
        UnitSpec::HeatExchanger {
            thermal,
            hot_pressure_drop_pa,
            cold_pressure_drop_pa,
        } => {
            let current = match thermal {
                ExchangerSpec::HotOutletTemperature(_) => 0,
                ExchangerSpec::ColdOutletTemperature(_) => 1,
                ExchangerSpec::Duty(_) => 2,
            };
            match choice(ui, "Specify", current, &["Hot out T", "Cold out T", "Duty"]) {
                Some(0) => *thermal = ExchangerSpec::HotOutletTemperature(None),
                Some(1) => *thermal = ExchangerSpec::ColdOutletTemperature(None),
                Some(_) => *thermal = ExchangerSpec::Duty(None),
                None => {}
            }
            match thermal {
                ExchangerSpec::HotOutletTemperature(t) => {
                    opt_field(ui, "Hot outlet temperature", t, Quantity::Temperature, u)
                }
                ExchangerSpec::ColdOutletTemperature(t) => {
                    opt_field(ui, "Cold outlet temperature", t, Quantity::Temperature, u)
                }
                ExchangerSpec::Duty(q) => opt_field(ui, "Duty (hot to cold)", q, Quantity::Duty, u),
            }
            value_field(
                ui,
                "Hot-side pressure drop",
                hot_pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
            value_field(
                ui,
                "Cold-side pressure drop",
                cold_pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
        }
        UnitSpec::Pressurize {
            pressure,
            efficiency,
        } => {
            let current = matches!(pressure, PressureRise::PressureIncrease(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet p", "Pressure increase"]) {
                Some(0) => *pressure = PressureRise::OutletPressure(None),
                Some(_) => *pressure = PressureRise::PressureIncrease(None),
                None => {}
            }
            match pressure {
                PressureRise::OutletPressure(p) => {
                    opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u)
                }
                PressureRise::PressureIncrease(dp) => {
                    opt_field(ui, "Pressure increase", dp, Quantity::PressureDrop, u)
                }
            }
            value_field(ui, "Efficiency", efficiency, Quantity::Fraction, u);
            note_row(ui, ROTATING_DEFERRED);
        }
        UnitSpec::Expand {
            pressure,
            efficiency,
        } => {
            let current = matches!(pressure, PressureSpec::PressureDrop(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet p", "Pressure drop"]) {
                Some(0) => *pressure = PressureSpec::OutletPressure(None),
                Some(_) => *pressure = PressureSpec::PressureDrop(None),
                None => {}
            }
            match pressure {
                PressureSpec::OutletPressure(p) => {
                    opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u)
                }
                PressureSpec::PressureDrop(dp) => {
                    opt_field(ui, "Pressure drop", dp, Quantity::PressureDrop, u)
                }
            }
            value_field(
                ui,
                "Isentropic efficiency",
                efficiency,
                Quantity::Fraction,
                u,
            );
            note_row(ui, ROTATING_DEFERRED);
        }
    }
}

/// Pumps, compressors, and expanders need a PS flash (proposal section 2.3).
const ROTATING_DEFERRED: &str = "Calculation deferred: needs a PS flash and efficiency model.";

/// A row of radio buttons; returns the newly chosen index if it changed.
fn choice(ui: &mut Ui, label: &str, current: usize, options: &[&str]) -> Option<usize> {
    ui.label(label);
    let mut picked = current;
    ui.horizontal_wrapped(|ui| {
        for (i, text) in options.iter().enumerate() {
            ui.radio_value(&mut picked, i, *text);
        }
    });
    ui.end_row();
    (picked != current).then_some(picked)
}

fn note_row(ui: &mut Ui, text: &str) {
    ui.label("");
    ui.label(RichText::new(text).small().weak());
    ui.end_row();
}

fn drag_value<'a>(shown: &'a mut f64, q: Quantity, u: DisplayUnits) -> DragValue<'a> {
    DragValue::new(shown)
        .speed(u.drag_speed(q))
        .range(u.display_min(q)..=u.display_max(q))
        .max_decimals(6)
        .suffix(format!(" {}", u.symbol(q)))
}

/// An optional SI value edited in display units. Unset values are
/// highlighted and offer an "Enter value" button.
fn opt_field(ui: &mut Ui, label: &str, value: &mut Option<f64>, q: Quantity, u: DisplayUnits) {
    ui.label(label);
    ui.horizontal(|ui| match value {
        Some(si) => {
            let mut shown = u.si_to_display(q, *si);
            if ui.add(drag_value(&mut shown, q, u)).changed() {
                *si = u.display_to_si(q, shown);
            }
            if ui.small_button("×").on_hover_text("Clear").clicked() {
                *value = None;
            }
        }
        None => {
            let warning = Theme::from_visuals(ui.visuals()).warning;
            ui.label(RichText::new("required").italics().color(warning));
            if ui
                .small_button("Enter value")
                .on_hover_text("Enter a value (starts from standard conditions)")
                .clicked()
            {
                *value = Some(starting_value(q));
            }
        }
    });
    ui.end_row();
}

/// A value that always has a default (for example a pressure drop of zero).
fn value_field(ui: &mut Ui, label: &str, si: &mut f64, q: Quantity, u: DisplayUnits) {
    ui.label(label);
    let mut shown = u.si_to_display(q, *si);
    if ui.add(drag_value(&mut shown, q, u)).changed() {
        *si = u.display_to_si(q, shown);
    }
    ui.end_row();
}

fn stream_properties(ui: &mut Ui, ed: &mut Editor, id: StreamId) {
    let Some(stream) = ed.fs.stream(id).cloned() else {
        return;
    };
    let mut name = stream.name.clone();
    ui.horizontal(|ui| {
        let label = ui.label("Stream");
        if ui
            .add(egui::TextEdit::singleline(&mut name).desired_width(110.0))
            .labelled_by(label.id)
            .changed()
        {
            ed.fs.rename_stream(id, name);
        }
        ui.label(RichText::new("(Material stream)").weak());
    });
    let (from_unit, from) = endpoint(ed, stream.from);
    let (to_unit, to) = endpoint(ed, stream.to);
    egui::Grid::new(("stream", id))
        .num_columns(2)
        .show(ui, |ui| {
            ui.label("From");
            if ui.link(from).clicked() {
                ed.selection = Some(Selection::Unit(from_unit));
            }
            ui.end_row();
            ui.label("To");
            if ui.link(to).clicked() {
                ed.selection = Some(Selection::Unit(to_unit));
            }
            ui.end_row();
        });
    results_placeholder(ui);
    ui.add_space(12.0);
    if ui.button("Delete stream").clicked() {
        ed.delete_selection();
    }
}

fn endpoint(ed: &Editor, port: PortRef) -> (UnitId, String) {
    let unit_name = ed
        .fs
        .unit(port.unit)
        .map(|u| u.name.as_str())
        .unwrap_or("?");
    let port_name = ed.fs.port_spec(port).map(|p| p.name).unwrap_or("?");
    (port.unit, format!("{unit_name} · {port_name}"))
}

fn results_placeholder(ui: &mut Ui) {
    ui.add_space(8.0);
    ui.strong("Results");
    egui::Frame::group(ui.style()).show(ui, |ui| {
        ui.set_width(ui.available_width());
        ui.label(RichText::new("Not solved — no values.").weak());
        ui.label(RichText::new(ENGINE_PENDING).small().weak());
    });
}

fn bottom_panel(ui: &mut Ui, ed: &mut Editor, tab: &mut BottomTab) {
    let items = ed.fs.diagnostics();
    ui.horizontal(|ui| {
        ui.selectable_value(
            tab,
            BottomTab::Messages,
            format!("Messages ({})", items.len()),
        );
        ui.selectable_value(tab, BottomTab::StreamTable, "Stream table");
        if let Some(status) = &ed.status {
            ui.separator();
            ui.colored_label(Theme::from_visuals(ui.visuals()).warning, status);
            if ui.small_button("×").clicked() {
                ed.status = None;
            }
        }
    });
    ui.separator();
    match tab {
        BottomTab::Messages => messages(ui, ed, items),
        BottomTab::StreamTable => stream_table(ui, ed),
    }
}

fn messages(ui: &mut Ui, ed: &mut Editor, items: Vec<flowsheet_core::Diagnostic>) {
    ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        if items.is_empty() {
            ui.label(
                RichText::new("No open items: topology and specifications are complete.").weak(),
            );
        }
        for d in items {
            let target = match d.kind {
                DiagnosticKind::UnconnectedPort(p) => Some(p.unit),
                DiagnosticKind::MissingSpecification(u, _)
                | DiagnosticKind::InvalidSpecification(u, _) => Some(u),
                DiagnosticKind::RecycleUnsupported(ref c) => c.first().copied(),
                DiagnosticKind::NoFeed => None,
            };
            let text = format!("! {}", d.message);
            match target {
                Some(unit) => {
                    if ui.link(text).clicked() {
                        ed.selection = Some(Selection::Unit(unit));
                    }
                }
                None => {
                    ui.label(text);
                }
            }
        }
    });
}

/// Streams as columns, properties as rows (proposal section 4.2). Values stay
/// empty until a solve exists.
fn stream_table(ui: &mut Ui, ed: &mut Editor) {
    let streams: Vec<_> = ed.fs.streams().cloned().collect();
    if streams.is_empty() {
        ui.label(RichText::new("No streams yet.").weak());
        return;
    }
    let u = ed.units;
    let rows = [
        ("Temperature", u.symbol(Quantity::Temperature)),
        ("Pressure", u.symbol(Quantity::Pressure)),
        ("Molar flow", u.symbol(Quantity::MolarFlow)),
        ("Vapor fraction", "–"),
    ];
    ScrollArea::both().auto_shrink(false).show(ui, |ui| {
        egui::Grid::new("stream_table")
            .striped(true)
            .min_col_width(70.0)
            .show(ui, |ui| {
                ui.label("");
                ui.label("");
                for s in &streams {
                    let selected = ed.selection == Some(Selection::Stream(s.id));
                    if ui
                        .selectable_label(selected, RichText::new(&s.name).strong())
                        .clicked()
                    {
                        ed.selection = Some(Selection::Stream(s.id));
                    }
                }
                ui.end_row();
                for (label, unit) in [("From", ""), ("To", "")] {
                    ui.label(label);
                    ui.label(unit);
                    for s in &streams {
                        let port = if label == "From" { s.from } else { s.to };
                        let name = ed
                            .fs
                            .unit(port.unit)
                            .map(|u| u.name.as_str())
                            .unwrap_or("?");
                        ui.label(RichText::new(name).weak());
                    }
                    ui.end_row();
                }
                for (label, unit) in rows {
                    ui.label(label);
                    ui.label(RichText::new(unit).weak());
                    for _ in &streams {
                        ui.label(RichText::new("—").weak());
                    }
                    ui.end_row();
                }
            });
    });
}
