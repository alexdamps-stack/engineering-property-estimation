//! Properties panel (specifications, feed composition, results), message
//! list, and stream table (proposal sections 3, 4.2 and 4.6).

use eframe::egui::{self, Color32, DragValue, RichText, ScrollArea, Ui};
use flowsheet_core::{
    DiagnosticKind, ExchangerSpec, MixerPressure, PortDirection, PortRef, PressureRise,
    PressureSpec, SeparatorThermal, SpecProblem, StreamId, ThermalSpec, UnitId, UnitKind, UnitSpec,
};
use flowsheet_project::{CubicVariantChoice, Fingerprints, StreamResult, format_utc};
use thermo_types::COMPOSITION_SUM_TOLERANCE;

use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;
use crate::units::{DisplayUnits, Quantity, format_quantity, format_value, starting_value};

/// Tabs of the bottom panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BottomTab {
    #[default]
    Messages,
    StreamTable,
}

/// Glyph, label and colour of a status.
pub fn status_style(status: ItemStatus, theme: &Theme) -> (&'static str, &'static str, Color32) {
    match status {
        ItemStatus::Underspecified => ("!", "Underspecified", theme.warning),
        ItemStatus::Ready => ("-", "Ready (not solved)", theme.muted),
        ItemStatus::Solved => ("✔", "Solved", theme.ok),
        ItemStatus::Warning => ("⚠", "Warning", theme.warning),
        ItemStatus::Failed => ("✖", "Failed", theme.error),
        ItemStatus::Stale => (
            "~",
            "Stale: inputs changed since the last solve",
            theme.muted,
        ),
    }
}

fn status_line(ui: &mut Ui, status: ItemStatus) {
    let theme = Theme::from_visuals(ui.visuals());
    let (glyph, label, color) = status_style(status, &theme);
    ui.horizontal(|ui| {
        ui.label("Status:");
        let text = RichText::new(format!("{glyph} {label}")).color(color);
        ui.label(if status == ItemStatus::Stale {
            text.italics()
        } else {
            text
        });
    });
}

/// A result value; stale values are grey italic and marked.
fn value_text(text: String, stale: bool) -> RichText {
    if stale {
        RichText::new(format!("{text} (stale)")).italics().weak()
    } else {
        RichText::new(text)
    }
}

/// Properties panel. Returns `true` if the user asked for a new project.
pub fn properties(ui: &mut Ui, ed: &mut Editor) -> bool {
    let mut new_project = false;
    ScrollArea::vertical().show(ui, |ui| {
        ui.add_space(6.0);
        match ed.selection {
            Some(Selection::Unit(id)) if ed.fs.unit(id).is_some() => unit_properties(ui, ed, id),
            Some(Selection::Stream(id)) if ed.fs.stream(id).is_some() => {
                stream_properties(ui, ed, id)
            }
            _ => new_project = project_properties(ui, ed),
        }
    });
    new_project
}

fn project_properties(ui: &mut Ui, ed: &mut Editor) -> bool {
    ui.strong("Project");
    ui.label(
        RichText::new(match &ed.path {
            Some(p) => p.display().to_string(),
            None => "Not saved yet.".into(),
        })
        .small()
        .weak(),
    );
    ui.add_space(6.0);
    let Some(package) = ed.package.clone() else {
        ui.label("No property package selected.");
        ui.label(
            RichText::new(
                "Create a project to choose components and the PR or PR-BM package, \
                 or load the demonstration.",
            )
            .small()
            .weak(),
        );
        return ui.button("New project…").clicked();
    };
    ui.strong("Property package");
    let mut variant = package.variant;
    ui.horizontal(|ui| {
        ui.label("Equation of state");
        for v in CubicVariantChoice::ALL {
            ui.radio_value(&mut variant, v, v.label());
        }
    });
    if variant != package.variant
        && let Some(p) = &mut ed.package
    {
        // A different package makes every result stale (fingerprints).
        p.variant = variant;
    }
    egui::Grid::new("package").num_columns(2).show(ui, |ui| {
        ui.label("Components");
        ui.label(package.components.join(", "));
        ui.end_row();
        ui.label("Interactions");
        ui.label("kij matrix, see below");
        ui.end_row();
        ui.label("PH interval");
        ui.label(format!(
            "{} – {} K",
            package.ph_minimum_temperature_k, package.ph_maximum_temperature_k
        ));
        ui.end_row();
        if let Some(prov) = ed.package_cache.get(&package).provenance.clone() {
            ui.label("Model");
            ui.label(format!("{} ({})", prov.model, prov.model_version));
            ui.end_row();
            ui.label("Data");
            ui.label(format!("{} ({})", prov.data, prov.data_version));
            ui.end_row();
            ui.label("Engine");
            ui.label(prov.engine);
            ui.end_row();
        }
    });
    ui.add(egui::Label::new(RichText::new(&package.interaction_source).small()).wrap());
    ui.add_space(8.0);
    ui.label(RichText::new("Select a unit or stream to edit it.").weak());
    false
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
    let fps = ed.fingerprints();
    let issues = ed.input_issues();
    let status = ed.unit_status(id, &fps, &issues);
    status_line(ui, status);

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
        feed_composition(ui, ed, id);
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

    unit_results(ui, ed, id, status, &fps);
    ui.add_space(12.0);
    if ui.button("Delete unit").clicked() {
        ed.delete_selection();
    }
}

fn feed_composition(ui: &mut Ui, ed: &mut Editor, id: UnitId) {
    ui.add_space(8.0);
    ui.strong("Composition (mole fractions)");
    let Some(package) = ed.package.clone() else {
        ui.label(
            RichText::new("Create a project to choose components first (File > New project).")
                .small()
                .weak(),
        );
        return;
    };
    let theme = Theme::from_visuals(ui.visuals());
    let n = package.components.len();
    let Some(x) = ed.feeds.get_mut(&id) else {
        ui.horizontal(|ui| {
            ui.label(RichText::new("required").italics().color(theme.warning));
            if ui
                .small_button("Enter composition")
                .on_hover_text("Start from equal mole fractions")
                .clicked()
            {
                let mut equal = vec![1.0 / n as f64; n];
                let rest: f64 = equal[..n - 1].iter().sum();
                equal[n - 1] = 1.0 - rest;
                ed.feeds.insert(id, equal);
            }
        });
        return;
    };
    egui::Grid::new(("composition", id))
        .num_columns(2)
        .show(ui, |ui| {
            for (name, value) in package.components.iter().zip(x.iter_mut()) {
                let label = ui.label(name);
                ui.add(
                    DragValue::new(value)
                        .range(0.0..=1.0)
                        .speed(0.001)
                        .max_decimals(10),
                )
                .labelled_by(label.id);
                ui.end_row();
            }
        });
    let sum: f64 = x.iter().sum();
    ui.horizontal(|ui| {
        let off = (sum - 1.0).abs() > COMPOSITION_SUM_TOLERANCE;
        let text = format!("Sum {sum:.10}");
        if off {
            ui.colored_label(theme.warning, format!("! {text}"));
            if sum > 0.0 && ui.small_button("Normalize").clicked() {
                for v in x.iter_mut() {
                    *v /= sum;
                }
            }
        } else {
            ui.label(RichText::new(text).weak());
        }
    });
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
        .clamp_existing_to_range(false)
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

fn unit_results(ui: &mut Ui, ed: &mut Editor, id: UnitId, status: ItemStatus, fps: &Fingerprints) {
    ui.add_space(8.0);
    ui.strong("Results");
    let theme = Theme::from_visuals(ui.visuals());
    let u = ed.units;
    let kind = ed.fs.unit(id).map(|x| x.kind);
    egui::Frame::group(ui.style()).show(ui, |ui| {
        ui.set_width(ui.available_width());
        if status == ItemStatus::Failed
            && let Some(f) = &ed.failure
        {
            ui.colored_label(theme.error, format!("✖ {}", f.message));
            ui.label(
                RichText::new("No values are shown for a failed calculation.")
                    .small()
                    .weak(),
            );
            return;
        }
        let result = ed.results.as_ref().and_then(|r| r.units.get(&id));
        let Some(r) = result.filter(|_| {
            matches!(
                status,
                ItemStatus::Solved | ItemStatus::Warning | ItemStatus::Stale
            )
        }) else {
            ui.label(RichText::new("Not solved — no values.").weak());
            return;
        };
        let stale = status == ItemStatus::Stale;
        egui::Grid::new(("unit_results", id))
            .num_columns(2)
            .show(ui, |ui| {
                if matches!(
                    kind,
                    Some(UnitKind::Heater | UnitKind::Cooler | UnitKind::Separator)
                ) {
                    ui.label("Duty (+ into unit)");
                    ui.label(value_text(
                        format_quantity(u, Quantity::Duty, r.duty_w),
                        stale,
                    ));
                    ui.end_row();
                }
                // Outlet summary from the first outlet stream.
                if let Some(kind) = kind {
                    for (i, port) in kind.ports().iter().enumerate() {
                        if port.direction != PortDirection::Outlet {
                            continue;
                        }
                        let Some(s) = ed.fs.stream_at(PortRef::new(id, i)) else {
                            continue;
                        };
                        let Some(sr) = ed.results.as_ref().and_then(|r| r.streams.get(&s.id))
                        else {
                            continue;
                        };
                        let current = ed.stream_status(s.id, fps);
                        let st = current == Some(ItemStatus::Stale) || stale;
                        ui.label(format!("{} ({})", port.name, s.name));
                        ui.label(value_text(stream_summary(sr, u), st));
                        ui.end_row();
                    }
                }
                if let (Some(c), Some(e)) =
                    (r.max_scaled_component_residual, r.scaled_energy_residual)
                {
                    ui.label("Balance residuals");
                    ui.label(
                        value_text(format!("components {c:.1e}, energy {e:.1e}"), stale).small(),
                    );
                    ui.end_row();
                }
            });
        for w in &r.warnings {
            ui.colored_label(theme.warning, format!("⚠ {w}"));
        }
    });
}

fn stream_summary(s: &StreamResult, u: DisplayUnits) -> String {
    match &s.state {
        Some(state) => format!(
            "{}, {}, {}",
            format_quantity(u, Quantity::Temperature, state.temperature_k),
            format_quantity(u, Quantity::Pressure, state.pressure_pa),
            format_quantity(u, Quantity::MolarFlow, s.molar_flow_mol_s)
        ),
        None => "no flow (absent phase)".into(),
    }
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
    let fps = ed.fingerprints();
    let status = ed.stream_status(id, &fps);
    ui.add_space(8.0);
    ui.strong("Results");
    let theme = Theme::from_visuals(ui.visuals());
    let u = ed.units;
    let components = ed
        .package
        .as_ref()
        .map(|p| p.components.clone())
        .unwrap_or_default();
    egui::Frame::group(ui.style()).show(ui, |ui| {
        ui.set_width(ui.available_width());
        let (Some(status), Some(r)) =
            (status, ed.results.as_ref().and_then(|r| r.streams.get(&id)))
        else {
            ui.label(RichText::new("Not solved — no values.").weak());
            return;
        };
        status_line(ui, status);
        let stale = status == ItemStatus::Stale;
        let row = |ui: &mut Ui, label: &str, text: String| {
            ui.label(label);
            ui.label(value_text(text, stale));
            ui.end_row();
        };
        egui::Grid::new(("stream_results", id))
            .num_columns(2)
            .show(ui, |ui| {
                row(
                    ui,
                    "Molar flow",
                    format_quantity(u, Quantity::MolarFlow, r.molar_flow_mol_s),
                );
                row(
                    ui,
                    "Mass flow",
                    format_quantity(u, Quantity::MassFlow, r.mass_flow_kg_s),
                );
                if let Some(s) = &r.state {
                    row(
                        ui,
                        "Temperature",
                        format_quantity(u, Quantity::Temperature, s.temperature_k),
                    );
                    row(
                        ui,
                        "Pressure",
                        format_quantity(u, Quantity::Pressure, s.pressure_pa),
                    );
                    row(
                        ui,
                        "Vapor fraction",
                        s.vapor_fraction
                            .map_or_else(|| "— (supercritical)".to_string(), format_value),
                    );
                    row(
                        ui,
                        "Molar enthalpy",
                        format_quantity(u, Quantity::MolarEnthalpy, s.enthalpy_j_mol),
                    );
                    row(
                        ui,
                        "Molar entropy",
                        format_quantity(u, Quantity::MolarEntropy, s.entropy_j_mol_k),
                    );
                }
            });
        let Some(s) = &r.state else {
            ui.label(RichText::new("Absent phase: zero flow and no thermodynamic state.").weak());
            return;
        };
        ui.add_space(6.0);
        egui::Grid::new(("phases", id))
            .num_columns(s.phases.len() + 1)
            .striped(true)
            .show(ui, |ui| {
                ui.label("");
                ui.strong("Overall");
                for p in &s.phases {
                    ui.strong(&p.label);
                }
                ui.end_row();
                ui.label("Fraction");
                ui.label("");
                for p in &s.phases {
                    ui.label(value_text(format_value(p.fraction), stale));
                }
                ui.end_row();
                ui.label(format!("Density, {}", u.symbol(Quantity::Density)));
                ui.label("");
                for p in &s.phases {
                    ui.label(value_text(format_value(p.density_kg_m3), stale));
                }
                ui.end_row();
                ui.label("Z");
                ui.label("");
                for p in &s.phases {
                    ui.label(value_text(format_value(p.compressibility), stale));
                }
                ui.end_row();
                for (i, c) in components.iter().enumerate() {
                    ui.label(format!("x {c}"));
                    ui.label(value_text(
                        s.mole_fractions
                            .get(i)
                            .map_or_else(|| "—".into(), |v| format_value(*v)),
                        stale,
                    ));
                    for p in &s.phases {
                        ui.label(value_text(
                            p.mole_fractions
                                .get(i)
                                .map_or_else(|| "—".into(), |v| format_value(*v)),
                            stale,
                        ));
                    }
                    ui.end_row();
                }
            });
        for w in &r.warnings {
            ui.colored_label(theme.warning, format!("⚠ {w}"));
        }
        ui.collapsing("Diagnostics and provenance", |ui| {
            if let Some(results) = &ed.results {
                let p = &results.provenance;
                if let Some(r) = &results.recycle {
                    ui.label(format!(
                        "Recycle: {} iterations, max scaled residual {:.2e}",
                        r.iterations, r.max_scaled_residual
                    ));
                    for (id, residual) in &r.tears {
                        ui.label(format!("Tear S{id}: {residual:.2e}"));
                    }
                }
                ui.label(format!("Model {} ({})", p.model, p.model_version));
                ui.label(format!("Data {} ({})", p.data, p.data_version));
                ui.label(format!(
                    "{} · solved {}",
                    p.engine,
                    format_utc(results.solved_at_unix_s)
                ));
            }
            ui.label(format!(
                "Iterations {}, max |ln f residual| {:.1e}",
                s.iterations, s.max_ln_fugacity_residual
            ));
            for p in &s.phases {
                ui.label(format!("{} phase validity: {}", p.label, p.validity));
            }
            for d in &s.diagnostics {
                ui.label(RichText::new(d).small());
            }
        });
    });
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

/// One message in the list.
struct Message {
    unit: Option<UnitId>,
    status: ItemStatus,
    text: String,
}

fn messages_of(ed: &mut Editor) -> Vec<Message> {
    let mut out = Vec::new();
    if ed.results.as_ref().is_some_and(|r| r.recycle.is_some()) {
        let fps = ed.fingerprints();
        let current = ed.results.as_ref().is_some_and(|results| {
            results
                .units
                .keys()
                .all(|id| results.unit_currency(*id, &fps) == flowsheet_project::Currency::Current)
        });
        let r = ed
            .results
            .as_ref()
            .and_then(|r| r.recycle.as_ref())
            .expect("recycle result");
        out.push(Message {
            unit: None,
            status: if current {
                ItemStatus::Solved
            } else {
                ItemStatus::Stale
            },
            text: format!(
                "Recycle {}: {} iterations, max scaled residual {:.2e}; tears {:?}",
                if current { "converged" } else { "result stale" },
                r.iterations,
                r.max_scaled_residual,
                r.tears
                    .iter()
                    .map(|(id, _)| format!("S{id}"))
                    .collect::<Vec<_>>()
            ),
        });
    }
    if let Some(f) = &ed.failure {
        out.push(Message {
            unit: f.unit,
            status: ItemStatus::Failed,
            text: format!("Solve failed: {}", f.message),
        });
    }
    for d in ed.fs.recycle_diagnostics() {
        let unit = match d.kind {
            DiagnosticKind::UnconnectedPort(p) => Some(p.unit),
            DiagnosticKind::MissingSpecification(u, _)
            | DiagnosticKind::InvalidSpecification(u, _) => Some(u),
            DiagnosticKind::RecycleUnsupported(ref c) => c.first().copied(),
            DiagnosticKind::NoFeed => None,
        };
        out.push(Message {
            unit,
            status: ItemStatus::Underspecified,
            text: d.message,
        });
    }
    for i in ed.input_issues() {
        out.push(Message {
            unit: Some(i.unit),
            status: ItemStatus::Underspecified,
            text: i.message,
        });
    }
    let fps = ed.fingerprints();
    if let Some(r) = &ed.results {
        for (id, u) in &r.units {
            if ed.fs.unit(*id).is_none() || u.warnings.is_empty() {
                continue;
            }
            let current = r.unit_currency(*id, &fps) == flowsheet_project::Currency::Current;
            for w in &u.warnings {
                out.push(Message {
                    unit: Some(*id),
                    status: if current {
                        ItemStatus::Warning
                    } else {
                        ItemStatus::Stale
                    },
                    text: format!("{}: {w}", ed.fs.unit(*id).map_or("?", |u| u.name.as_str())),
                });
            }
        }
    }
    out
}

/// Bottom panel with messages and the stream table.
pub fn bottom_panel(ui: &mut Ui, ed: &mut Editor, tab: &mut BottomTab) {
    let items = messages_of(ed);
    ui.horizontal(|ui| {
        ui.selectable_value(
            tab,
            BottomTab::Messages,
            format!("Messages ({})", items.len()),
        );
        ui.selectable_value(tab, BottomTab::StreamTable, "Stream table");
        if *tab == BottomTab::StreamTable
            && ui
                .small_button("Copy")
                .on_hover_text("Copy the stream table as tab-separated text for spreadsheets")
                .clicked()
        {
            let tsv = stream_table_tsv(ed);
            ui.ctx().copy_text(tsv);
            ed.status = Some("Stream table copied as tab-separated text.".into());
        }
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

fn messages(ui: &mut Ui, ed: &mut Editor, items: Vec<Message>) {
    let theme = Theme::from_visuals(ui.visuals());
    ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        if items.is_empty() {
            ui.label(
                RichText::new("No open items: topology and specifications are complete.").weak(),
            );
        }
        for m in items {
            let (glyph, _, color) = status_style(m.status, &theme);
            let text = RichText::new(format!("{glyph} {}", m.text)).color(color);
            match m.unit {
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

/// Row definitions of the stream table: label, quantity (for the unit) and
/// value extraction in SI.
type Extract = fn(&StreamResult) -> Option<f64>;

fn table_rows() -> Vec<(&'static str, Option<Quantity>, Extract)> {
    vec![
        ("Temperature", Some(Quantity::Temperature), |r| {
            r.state.as_ref().map(|s| s.temperature_k)
        }),
        ("Pressure", Some(Quantity::Pressure), |r| {
            r.state.as_ref().map(|s| s.pressure_pa)
        }),
        ("Molar flow", Some(Quantity::MolarFlow), |r| {
            Some(r.molar_flow_mol_s)
        }),
        ("Mass flow", Some(Quantity::MassFlow), |r| {
            Some(r.mass_flow_kg_s)
        }),
        ("Vapor fraction", None, |r| {
            r.state.as_ref().and_then(|s| s.vapor_fraction)
        }),
        ("Molar enthalpy", Some(Quantity::MolarEnthalpy), |r| {
            r.state.as_ref().map(|s| s.enthalpy_j_mol)
        }),
        ("Molar entropy", Some(Quantity::MolarEntropy), |r| {
            r.state.as_ref().map(|s| s.entropy_j_mol_k)
        }),
    ]
}

struct TableColumn {
    name: String,
    status: Option<ItemStatus>,
    result: Option<StreamResult>,
}

fn table_columns(ed: &mut Editor) -> Vec<(StreamId, TableColumn)> {
    let fps = ed.fingerprints();
    let streams: Vec<_> = ed.fs.streams().map(|s| (s.id, s.name.clone())).collect();
    streams
        .into_iter()
        .map(|(id, name)| {
            let status = ed.stream_status(id, &fps);
            let result = status.and(
                ed.results
                    .as_ref()
                    .and_then(|r| r.streams.get(&id).cloned()),
            );
            (
                id,
                TableColumn {
                    name,
                    status,
                    result,
                },
            )
        })
        .collect()
}

/// Stream table as tab-separated text with display units, a status row and
/// one row per property; empty cells where no value exists.
pub fn stream_table_tsv(ed: &mut Editor) -> String {
    let u = ed.units;
    let cols = table_columns(ed);
    let components = ed
        .package
        .as_ref()
        .map(|p| p.components.clone())
        .unwrap_or_default();
    let mut out = String::from("Stream\tUnit");
    for (_, c) in &cols {
        out.push('\t');
        out.push_str(&c.name);
    }
    out.push_str("\nStatus\t");
    for (_, c) in &cols {
        out.push('\t');
        out.push_str(match c.status {
            Some(ItemStatus::Solved) => "solved",
            Some(ItemStatus::Warning) => "warning",
            Some(ItemStatus::Stale) => "stale",
            _ => "",
        });
    }
    for (label, q, get) in table_rows() {
        out.push('\n');
        out.push_str(label);
        out.push('\t');
        out.push_str(q.map_or("–", |q| u.symbol(q)));
        for (_, c) in &cols {
            out.push('\t');
            if let Some(v) = c.result.as_ref().and_then(get) {
                let shown = q.map_or(v, |q| u.si_to_display(q, v));
                if shown.is_finite() {
                    out.push_str(&shown.to_string());
                }
            }
        }
    }
    for (i, name) in components.iter().enumerate() {
        out.push_str(&format!("\nx {name}\tmol/mol"));
        for (_, c) in &cols {
            out.push('\t');
            if let Some(s) = c.result.as_ref().and_then(|r| r.state.as_ref())
                && let Some(value) = s.mole_fractions.get(i)
            {
                out.push_str(&value.to_string());
            }
        }
    }
    out.push('\n');
    out
}

/// Streams as columns, properties as rows (proposal section 4.2).
fn stream_table(ui: &mut Ui, ed: &mut Editor) {
    let cols = table_columns(ed);
    if cols.is_empty() {
        ui.label(RichText::new("No streams yet.").weak());
        return;
    }
    let theme = Theme::from_visuals(ui.visuals());
    let u = ed.units;
    let components = ed
        .package
        .as_ref()
        .map(|p| p.components.clone())
        .unwrap_or_default();
    ScrollArea::both().auto_shrink(false).show(ui, |ui| {
        egui::Grid::new("stream_table")
            .striped(true)
            .min_col_width(70.0)
            .show(ui, |ui| {
                ui.label("");
                ui.label("");
                for (id, c) in &cols {
                    let selected = ed.selection == Some(Selection::Stream(*id));
                    if ui
                        .selectable_label(selected, RichText::new(&c.name).strong())
                        .clicked()
                    {
                        ed.selection = Some(Selection::Stream(*id));
                    }
                }
                ui.end_row();
                for label in ["From", "To"] {
                    ui.label(label);
                    ui.label("");
                    for (id, _) in &cols {
                        let s = ed.fs.stream(*id).expect("listed stream");
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
                ui.label("Status");
                ui.label("");
                for (_, c) in &cols {
                    match c.status {
                        Some(s) => {
                            let (glyph, label, color) = status_style(s, &theme);
                            let short = if s == ItemStatus::Stale {
                                "stale"
                            } else {
                                label
                            };
                            ui.label(RichText::new(format!("{glyph} {short}")).color(color));
                        }
                        None => {
                            ui.label(RichText::new("—").weak());
                        }
                    }
                }
                ui.end_row();
                for (label, q, get) in table_rows() {
                    ui.label(label);
                    ui.label(RichText::new(q.map_or("–", |q| u.symbol(q))).weak());
                    for (_, c) in &cols {
                        let stale = c.status == Some(ItemStatus::Stale);
                        match c.result.as_ref().and_then(get) {
                            Some(v) => {
                                let text = format_value(q.map_or(v, |q| u.si_to_display(q, v)));
                                ui.label(if stale {
                                    RichText::new(text).italics().weak()
                                } else {
                                    RichText::new(text)
                                });
                            }
                            None => {
                                ui.label(RichText::new("—").weak());
                            }
                        }
                    }
                    ui.end_row();
                }
                for (i, name) in components.iter().enumerate() {
                    ui.label(format!("x {name}"));
                    ui.label(RichText::new("mol/mol").weak());
                    for (_, c) in &cols {
                        let stale = c.status == Some(ItemStatus::Stale);
                        match c.result.as_ref().and_then(|r| r.state.as_ref()) {
                            Some(s) => {
                                let text = s
                                    .mole_fractions
                                    .get(i)
                                    .map_or_else(|| "—".into(), |v| format_value(*v));
                                ui.label(if stale {
                                    RichText::new(text).italics().weak()
                                } else {
                                    RichText::new(text)
                                });
                            }
                            None => {
                                ui.label(RichText::new("—").weak());
                            }
                        }
                    }
                    ui.end_row();
                }
            });
    });
}
