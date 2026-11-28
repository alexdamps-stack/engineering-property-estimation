//! Unit properties view: tag, specification, feed composition, connections
//! and results.

use eframe::egui::{self, DragValue, RichText, Ui};
use flowsheet_core::units::HeatTransportSpec;
use flowsheet_core::units::heat_transport::{
    AirCoolerThermal, ExchangerPassDirection, ExchangerPassThermal,
};
use flowsheet_core::{
    ExchangerSpec, PortDirection, PortMedium, PortRef, SpecProblem, UnitId, UnitKind, UnitSpec,
};
use flowsheet_project::Fingerprints;
use thermo_types::COMPOSITION_SUM_TOLERANCE;

use super::presentation::stream_summary;
use super::reaction::chemistry_form;
use super::spec::spec_form;
use super::{status_line, value_text};
use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;
use crate::units::{DisplayUnits, Quantity, format_quantity, format_residual};

pub(super) fn unit_properties(ui: &mut Ui, ed: &mut Editor, id: UnitId) {
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
    let mut spec = unit.spec.clone();
    let u = ed.units;
    let components = ed
        .package
        .as_ref()
        .map_or(&[][..], |p| p.components.as_slice());
    egui::Grid::new(("spec", id))
        .num_columns(2)
        .spacing([8.0, 6.0])
        .show(ui, |ui| spec_form(ui, unit.kind, &mut spec, u, components));
    chemistry_form(ui, &mut spec, components, u);
    let issues = spec.issues();
    if let Err(e) = ed.fs.try_set_spec(id, spec) {
        // For example fewer column feeds while the last one is connected.
        ed.status = Some(format!("Specification not changed: {e}."));
    }
    for issue in issues {
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
            for (i, port) in unit.ports().iter().enumerate() {
                let label = match port.medium {
                    PortMedium::Material => port.label(),
                    PortMedium::Energy => format!("{port} (energy)"),
                };
                ui.label(label);
                match ed.fs.stream_at(PortRef::new(id, i)) {
                    Some(s) => {
                        let sid = s.id;
                        if ui.link(s.name.clone()).clicked() {
                            ed.selection = Some(Selection::Stream(sid));
                        }
                    }
                    None if port.optional => {
                        ui.label(RichText::new("optional, not connected").weak());
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
    if n == 0 {
        ui.colored_label(
            theme.warning,
            "! The property package has no components; create a project with at least one component.",
        );
        return;
    }
    let Some(x) = ed.feeds.get_mut(&id) else {
        ui.horizontal(|ui| {
            ui.label(RichText::new("required").italics().color(theme.warning));
            if ui
                .small_button("Enter composition")
                .on_hover_text("Start from equal mole fractions")
                .clicked()
            {
                ed.feeds.insert(id, equal_fractions(n));
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
                        .clamp_existing_to_range(false)
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

fn unit_results(ui: &mut Ui, ed: &mut Editor, id: UnitId, status: ItemStatus, fps: &Fingerprints) {
    ui.add_space(8.0);
    ui.strong("Results");
    let theme = Theme::from_visuals(ui.visuals());
    let u = ed.units;
    let unit = ed.fs.unit(id);
    let kind = unit.map(|x| x.kind);
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
                let stream = |port| {
                    let connected = ed.fs.stream_at(PortRef::new(id, port))?;
                    ed.results.as_ref()?.streams.get(&connected.id)
                };
                let temperature = |port| Some(stream(port)?.state.as_ref()?.temperature_k);
                if matches!(
                    kind,
                    Some(
                        UnitKind::Heater
                            | UnitKind::Cooler
                            | UnitKind::Separator
                            | UnitKind::ConversionReactor
                            | UnitKind::YieldReactor
                            | UnitKind::Tank
                            | UnitKind::AirCooler
                            | UnitKind::FiredHeater
                    )
                ) {
                    ui.label(if kind == Some(UnitKind::AirCooler) {
                        "Process duty (+ into fluid)"
                    } else {
                        "Duty (+ into unit)"
                    });
                    ui.label(value_text(
                        format_quantity(u, Quantity::Duty, r.duty_w),
                        stale,
                    ));
                    ui.end_row();
                }
                if let Some(reaction) = &r.reaction {
                    // Formation basis at 298.15 K: negative when exothermic.
                    ui.label("Standard heat of reaction")
                        .on_hover_text("Sum of generation × ideal-gas formation enthalpy at 298.15 K; negative when exothermic");
                    ui.label(value_text(
                        format_quantity(u, Quantity::Duty, reaction.standard_reaction_enthalpy_w),
                        stale,
                    ));
                    ui.end_row();
                    for (k, xi) in reaction.extents_mol_s.iter().enumerate() {
                        ui.label(format!("Extent, reaction {}", k + 1));
                        ui.label(value_text(
                            format_quantity(u, Quantity::MolarFlow, *xi),
                            stale,
                        ));
                        ui.end_row();
                    }
                }
                if let Some(m) = &r.rotating {
                    rotating_rows(ui, u, stale, m);
                }
                // These quantities combine the solved duty with the current specification.
                if !stale
                    && let Some(UnitSpec::HeatTransport(HeatTransportSpec::AirCooler(s))) =
                        unit.map(|x| &x.spec)
                    && let Some(t) = s.air_outlet_temperature_k(r.duty_w)
                {
                    ui.label("Air discharge temperature");
                    ui.label(value_text(
                        format_quantity(u, Quantity::Temperature, t),
                        false,
                    ));
                    ui.end_row();
                    if matches!(s.thermal, AirCoolerThermal::HotEndApproach(_) | AirCoolerThermal::ColdEndApproach(_))
                        && let Some(pre_fan) = s.pre_fan_air_outlet_temperature_k(r.duty_w)
                    {
                        ui.label("Pre-fan air outlet temperature");
                        ui.label(value_text(format_quantity(u, Quantity::Temperature, pre_fan), false));
                        ui.end_row();
                    }
                    ui.label("Fan power");
                    ui.label(value_text(
                        format_quantity(u, Quantity::Duty, s.fan_power_w),
                        false,
                    ));
                    ui.end_row();
                }
                if matches!(kind, Some(UnitKind::HeatExchanger | UnitKind::MultiStreamExchanger)) {
                    // An exchanger reports internal transfer and has zero external duty.
                    ui.label("Transferred duty (hot to cold)");
                    ui.label(match r.transferred_duty_w {
                        Some(q) => value_text(format_quantity(u, Quantity::Duty, q), stale),
                        None => RichText::new("not stored — recalculate").weak(),
                    });
                    ui.end_row();
                }
                if !stale {
                    let approach = match unit.map(|x| &x.spec) {
                        Some(UnitSpec::HeatExchanger { thermal: ExchangerSpec::HotEndApproach(_), .. }) => {
                            temperature(0).zip(temperature(3)).map(|(hot, cold)| ("Achieved hot-end approach", hot - cold))
                        }
                        Some(UnitSpec::HeatExchanger { thermal: ExchangerSpec::ColdEndApproach(_), .. }) => {
                            temperature(2).zip(temperature(1)).map(|(hot, cold)| ("Achieved cold-end approach", hot - cold))
                        }
                        Some(UnitSpec::HeatTransport(HeatTransportSpec::AirCooler(s))) => match s.thermal {
                            AirCoolerThermal::HotEndApproach(_) => temperature(0)
                                .zip(s.pre_fan_air_outlet_temperature_k(r.duty_w))
                                .map(|(hot, cold)| ("Achieved hot-end approach", hot - cold)),
                            AirCoolerThermal::ColdEndApproach(_) => temperature(1)
                                .zip(s.ambient_temperature_k)
                                .map(|(hot, cold)| ("Achieved cold-end approach", hot - cold)),
                            _ => None,
                        },
                        _ => None,
                    };
                    if let Some((label, delta)) = approach {
                        ui.label(label);
                        ui.label(value_text(format_quantity(u, Quantity::TemperatureDifference, delta), false));
                        ui.end_row();
                    }
                }
                if !stale
                    && let Some(UnitSpec::HeatTransport(HeatTransportSpec::MultiStreamExchanger(s))) =
                        unit.map(|x| &x.spec)
                {
                    for (i, pass) in s.passes.iter().enumerate() {
                        let (Some(inlet), Some(outlet)) = (stream(i * 2), stream(i * 2 + 1)) else {
                            continue;
                        };
                        let q = match (&inlet.state, &outlet.state) {
                            (Some(a), Some(b)) => {
                                inlet.molar_flow_mol_s * (b.enthalpy_j_mol - a.enthalpy_j_mol)
                            }
                            _ => 0.0,
                        };
                        ui.label(format!("Pass {} duty (+ into fluid)", i + 1));
                        ui.label(value_text(format_quantity(u, Quantity::Duty, q), false));
                        ui.end_row();
                        if let ExchangerPassThermal::TerminalApproach { reference_pass, .. } =
                            pass.thermal
                            && let Some(index) = reference_pass.checked_sub(1)
                            && let Some(reference) =
                                stream(usize::from(index) * 2).and_then(|s| s.state.as_ref())
                            && let Some(out) = &outlet.state
                        {
                            let approach = match pass.direction {
                                ExchangerPassDirection::Hot => out.temperature_k - reference.temperature_k,
                                ExchangerPassDirection::Cold => reference.temperature_k - out.temperature_k,
                            };
                            ui.label(format!("Pass {} terminal approach", i + 1));
                            ui.label(value_text(
                                format_quantity(u, Quantity::TemperatureDifference, approach),
                                false,
                            ));
                            ui.end_row();
                        }
                    }
                }
                if let Some(c) = &r.shortcut_column {
                    let number = |v: f64, digits: usize| value_text(format!("{v:.digits$}"), stale);
                    for (label, text) in [
                        ("Minimum stages (Fenske)", number(c.minimum_stages, 2)),
                        ("Minimum reflux ratio (Underwood)", number(c.minimum_reflux_ratio, 4)),
                        ("Reflux ratio L/D", number(c.reflux_ratio, 4)),
                        ("Theoretical stages (Gilliland)", number(c.theoretical_stages, 2)),
                        ("Feed stage from top (Kirkbride)", number(c.rectifying_stages + 1.0, 2)),
                        ("Key relative volatility", number(c.relative_volatility, 4)),
                        ("Feed condition q", number(c.feed_condition_q, 4)),
                        ("Gilliland X", number(c.gilliland_x, 4)),
                        (
                            "Condenser duty (+ removed)",
                            value_text(format_quantity(u, Quantity::Duty, c.condenser_duty_w), stale),
                        ),
                        (
                            "Reboiler duty (+ supplied)",
                            value_text(format_quantity(u, Quantity::Duty, c.reboiler_duty_w), stale),
                        ),
                        (
                            "Overhead vapour dew point",
                            value_text(
                                format_quantity(u, Quantity::Temperature, c.top_vapor_temperature_k),
                                stale,
                            ),
                        ),
                    ] {
                        ui.label(label);
                        ui.label(text);
                        ui.end_row();
                    }
                    // Key and non-key distribution (Fenske, total reflux).
                    if let Some(package) = &ed.package {
                        for (name, recovery) in package.components.iter().zip(&c.distillate_recoveries) {
                            if let Some(recovery) = recovery {
                                ui.label(format!("{name} to distillate"));
                                ui.label(value_text(
                                    format_quantity(u, Quantity::Fraction, *recovery),
                                    stale,
                                ));
                                ui.end_row();
                            }
                        }
                    }
                }
                if let Some(c) = &r.distillation_column {
                    for (label, q) in [
                        ("Condenser duty (+ removed)", c.condenser_duty_w),
                        ("Reboiler duty (+ supplied)", c.reboiler_duty_w),
                    ] {
                        ui.label(label);
                        ui.label(value_text(format_quantity(u, Quantity::Duty, q), stale));
                        ui.end_row();
                    }
                    for (label, q) in [
                        ("Distillate", c.distillate_mol_s),
                        ("Reflux", c.reflux_mol_s),
                        ("Bottoms", c.bottoms_mol_s),
                        ("Boilup", c.boilup_mol_s),
                    ] {
                        if let Some(q) = q {
                            ui.label(label);
                            ui.label(value_text(format_quantity(u, Quantity::MolarFlow, q), stale));
                            ui.end_row();
                        }
                    }
                    if let (Some(top), Some(bottom)) =
                        (c.stage_pressures_pa.first(), c.stage_pressures_pa.last())
                        && top != bottom
                    {
                        ui.label("Condenser / reboiler pressure");
                        ui.label(value_text(
                            format!(
                                "{} / {}",
                                format_quantity(u, Quantity::Pressure, *top),
                                format_quantity(u, Quantity::Pressure, *bottom)
                            ),
                            stale,
                        ));
                        ui.end_row();
                    }
                    stage_profile(ui, u, stale, &c.stage_temperatures_k, c.packed_height_m);
                    product_curve_rows(ui, u, stale, &c.product_curves);
                }
                if let Some(c) = &r.absorber {
                    match &c.rate_based {
                        None => {
                            stage_profile(ui, u, stale, &c.stage_temperatures_k, c.packed_height_m)
                        }
                        Some(rate) => rate_based_rows(ui, u, stale, c.packed_height_m, rate),
                    }
                }
                if let Some(c) = &r.liquid_extractor {
                    ui.label("Feed enters");
                    ui.label(value_text(
                        if c.feed_is_heavy { "top (heavy feed)" } else { "bottom (light feed)" }
                            .to_owned(),
                        stale,
                    ));
                    ui.end_row();
                    stage_profile(ui, u, stale, &c.stage_temperatures_k, c.packed_height_m);
                    for (j, (e, f)) in c
                        .extract_flows_mol_s
                        .iter()
                        .zip(&c.raffinate_flows_mol_s)
                        .enumerate()
                    {
                        ui.label(format!("Stage {} extract / raffinate", j + 1));
                        ui.label(value_text(
                            format!(
                                "{} / {}",
                                format_quantity(u, Quantity::MolarFlow, *e),
                                format_quantity(u, Quantity::MolarFlow, *f)
                            ),
                            stale,
                        ));
                        ui.end_row();
                    }
                }
                if let Some(c) = &r.three_phase_column {
                    for (label, q) in [
                        ("Condenser duty (+ removed)", c.condenser_duty_w),
                        ("Reboiler duty (+ supplied)", c.reboiler_duty_w),
                    ] {
                        ui.label(label);
                        ui.label(value_text(format_quantity(u, Quantity::Duty, q), stale));
                        ui.end_row();
                    }
                    for (label, q) in [
                        ("Decanter light liquid", c.light_liquid_mol_s),
                        ("Decanter heavy liquid", c.heavy_liquid_mol_s),
                        ("Reflux", c.liquid_flows_mol_s.first().copied().unwrap_or(0.0)),
                    ] {
                        ui.label(label);
                        ui.label(value_text(format_quantity(u, Quantity::MolarFlow, q), stale));
                        ui.end_row();
                    }
                    ui.label("Reflux liquid");
                    ui.label(value_text(
                        if c.light_refluxed { "light" } else { "heavy" }.to_owned(),
                        stale,
                    ));
                    ui.end_row();
                    stage_profile(ui, u, stale, &c.stage_temperatures_k, c.packed_height_m);
                    three_phase_stage_rows(ui, u, stale, c);
                }
                if let Some(x) = &r.material_exchange {
                    ui.label("External material (+ into unit)");
                    ui.label(value_text(
                        format_quantity(u, Quantity::MolarFlow, x.molar_flow_mol_s),
                        stale,
                    ));
                    ui.end_row();
                    if !stale && let Some(package) = &ed.package {
                        for (id, &flow) in package.components.iter().zip(&x.component_flows_mol_s) {
                            ui.label(id);
                            ui.label(value_text(
                                format_quantity(u, Quantity::MolarFlow, flow),
                                false,
                            ));
                            ui.end_row();
                        }
                    }
                    ui.label("Enthalpy carried by material");
                    ui.label(value_text(
                        format_quantity(u, Quantity::Duty, x.enthalpy_flow_w),
                        stale,
                    ));
                    ui.end_row();
                }
                // Energy streams at the unit (#171), positive in the stream
                // direction.
                let ports = ed.fs.unit(id).map(|x| x.ports()).unwrap_or_default();
                for (i, port) in ports.iter().enumerate() {
                    if port.medium != PortMedium::Energy {
                        continue;
                    }
                    let Some(s) = ed.fs.stream_at(PortRef::new(id, i)) else {
                        continue;
                    };
                    let Some(e) = ed
                        .results
                        .as_ref()
                        .and_then(|r| r.energy_streams.get(&s.id))
                    else {
                        continue;
                    };
                    let st = ed.stream_status(s.id, fps) == Some(ItemStatus::Stale) || stale;
                    let arrow = match port.direction {
                        PortDirection::Inlet => "in",
                        PortDirection::Outlet => "out",
                    };
                    ui.label(format!("{port} ({}, {arrow})", s.name));
                    ui.label(value_text(format_quantity(u, Quantity::Duty, e.duty_w), st));
                    ui.end_row();
                }
                // Outlet summary of every connected material outlet.
                if kind.is_some() {
                    for (i, port) in ports.iter().enumerate() {
                        if port.direction != PortDirection::Outlet
                            || port.medium != PortMedium::Material
                        {
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
                        ui.label(format!("{} ({})", port.label(), s.name));
                        ui.label(value_text(stream_summary(sr, u), st));
                        ui.end_row();
                    }
                }
                if let (Some(c), Some(e)) =
                    (r.max_scaled_component_residual, r.scaled_energy_residual)
                {
                    ui.label("Balance residuals");
                    ui.label(
                        value_text(
                            format!(
                                "components {}, energy {}",
                                format_residual(c),
                                format_residual(e)
                            ),
                            stale,
                        )
                        .small(),
                    );
                    ui.end_row();
                }
            });
        for w in &r.warnings {
            ui.colored_label(theme.warning, format!("⚠ {w}"));
        }
    });
}

/// TBP and ASTM D86 points of petroleum fractionator products (#190), one
/// row per product and curve: "T5 / T10 / T30 / T50 / T70 / T90 / T95".
fn product_curve_rows(
    ui: &mut Ui,
    units: DisplayUnits,
    stale: bool,
    curves: &[flowsheet_project::ProductCurveEstimate],
) {
    for curve in curves {
        for (name, points) in [("TBP", &curve.tbp_k), ("D86", &curve.d86_k)] {
            let percents: Vec<String> = points.iter().map(|p| format!("{}", p[0])).collect();
            ui.label(format!(
                "{} {name} {} vol%",
                curve.product,
                percents.join("/")
            ));
            let values: Vec<String> = points
                .iter()
                .map(|p| format_quantity(units, Quantity::Temperature, p[1]))
                .collect();
            ui.label(value_text(values.join(" / "), stale));
            ui.end_row();
        }
        if curve.excluded_mole_fraction > 0.0 {
            ui.label(format!("{} left out of the curves", curve.product));
            ui.label(value_text(
                format!(
                    "{} (non-pseudo-components)",
                    format_quantity(units, Quantity::Fraction, curve.excluded_mole_fraction)
                ),
                stale,
            ));
            ui.end_row();
        }
    }
}

fn stage_profile(
    ui: &mut Ui,
    units: DisplayUnits,
    stale: bool,
    stages: &[f64],
    height: Option<f64>,
) {
    if let Some(height) = height {
        ui.label("Packed height");
        ui.label(value_text(
            format_quantity(units, Quantity::Length, height),
            stale,
        ));
        ui.end_row();
    }
    for (j, &t) in stages.iter().enumerate() {
        ui.label(format!("Stage {} temperature", j + 1));
        ui.label(value_text(
            format_quantity(units, Quantity::Temperature, t),
            stale,
        ));
        ui.end_row();
    }
}

/// Hydraulics and segment profiles of a rate-based packed section (#241),
/// segments from the top.
fn rate_based_rows(
    ui: &mut Ui,
    units: DisplayUnits,
    stale: bool,
    height: Option<f64>,
    rate: &flowsheet_project::RateBasedSegmentsProfile,
) {
    let mut row = |label: String, text: String| {
        ui.label(label);
        ui.label(value_text(text, stale));
        ui.end_row();
    };
    if let Some(height) = height {
        row(
            "Packed height".into(),
            format_quantity(units, Quantity::Length, height),
        );
    }
    row("Packing".into(), rate.packing.clone());
    row(
        "Column diameter".into(),
        format_quantity(units, Quantity::Length, rate.diameter_m),
    );
    row(
        "Pressure drop (reported)".into(),
        format_quantity(units, Quantity::PressureDrop, rate.pressure_drop_pa),
    );
    row(
        "Max. fraction of flood".into(),
        format_quantity(units, Quantity::Fraction, rate.max_flooding_fraction),
    );
    if rate.max_neglected_transport_fraction > 0.0 {
        row(
            "Largest neglected transport fraction".into(),
            format!("{:.2e}", rate.max_neglected_transport_fraction),
        );
    }
    let boundary_pressures = &rate.estimated_boundary_pressures_pa;
    if let Some(&top) = boundary_pressures.first() {
        row(
            "Bed top estimated pressure".into(),
            format_quantity(units, Quantity::Pressure, top),
        );
    }
    for (j, t) in rate.liquid_temperatures_k.iter().enumerate() {
        let temperature = |t: f64| format_quantity(units, Quantity::Temperature, t);
        row(
            format!("Segment {} T vapor / interface / liquid", j + 1),
            format!(
                "{} / {} / {}",
                temperature(rate.vapor_temperatures_k[j]),
                temperature(rate.interface_temperatures_k[j]),
                temperature(*t)
            ),
        );
        if let Some(&pressure) = boundary_pressures.get(j + 1) {
            row(
                format!("Segment {} bottom estimated pressure", j + 1),
                format_quantity(units, Quantity::Pressure, pressure),
            );
        }
        row(
            format!("Segment {} flood / ΔP per m / a_e", j + 1),
            format!(
                "{} / {:.1} Pa/m / {:.1} m²/m³",
                format_quantity(units, Quantity::Fraction, rate.flooding_fractions[j]),
                rate.pressure_drops_pa_m[j],
                rate.interfacial_area_m2_m3[j]
            ),
        );
        row(
            format!("Segment {} transfer to liquid", j + 1),
            rate.fluxes_mol_s[j]
                .iter()
                .map(|n| format_quantity(units, Quantity::MolarFlow, *n))
                .collect::<Vec<_>>()
                .join(", "),
        );
    }
}

/// Phase counts and both liquids of every two-liquid stage of a three-phase
/// column (#188); stage 1 is the decanter.
fn three_phase_stage_rows(
    ui: &mut Ui,
    units: DisplayUnits,
    stale: bool,
    c: &flowsheet_project::ThreePhaseColumnProfile,
) {
    if !c.stage_phase_counts.is_empty() {
        ui.label("Phases per stage");
        ui.label(value_text(
            c.stage_phase_counts
                .iter()
                .map(u8::to_string)
                .collect::<Vec<_>>()
                .join(" "),
            stale,
        ));
        ui.end_row();
    }
    let fractions = |x: &[f64]| {
        x.iter()
            .map(|v| {
                if *v != 0.0 && v.abs() < 1e-3 {
                    format!("{v:.3e}")
                } else {
                    format!("{v:.4}")
                }
            })
            .collect::<Vec<_>>()
            .join(", ")
    };
    for (j, split) in c.stage_liquid_splits.iter().enumerate().skip(1) {
        let Some(s) = split else { continue };
        ui.label(format!("Stage {} light / heavy liquid", j + 1));
        ui.label(value_text(
            format!(
                "{} / {}",
                format_quantity(units, Quantity::MolarFlow, s.light_mol_s),
                format_quantity(units, Quantity::MolarFlow, s.heavy_mol_s)
            ),
            stale,
        ));
        ui.end_row();
        ui.label(format!("Stage {} light liquid x", j + 1));
        ui.label(value_text(fractions(&s.light_composition), stale));
        ui.end_row();
        ui.label(format!("Stage {} heavy liquid x", j + 1));
        ui.label(value_text(fractions(&s.heavy_composition), stale));
        ui.end_row();
    }
}

/// Shaft power, isentropic states and intercooling of a rotating machine.
fn rotating_rows(
    ui: &mut Ui,
    units: DisplayUnits,
    stale: bool,
    m: &flowsheet_project::RotatingPerformance,
) {
    let row = |ui: &mut Ui, label: String, q: Quantity, v: f64| {
        ui.label(label);
        ui.label(value_text(format_quantity(units, q, v), stale));
        ui.end_row();
    };
    row(
        ui,
        "Shaft power (+ into fluid)".into(),
        Quantity::Duty,
        m.power_w(),
    );
    let multistage = m.stages.len() > 1;
    for (k, s) in m.stages.iter().enumerate() {
        let prefix = if multistage {
            format!("Stage {} ", k + 1)
        } else {
            String::new()
        };
        if multistage {
            row(
                ui,
                format!("{prefix}discharge pressure"),
                Quantity::Pressure,
                s.outlet_pressure_pa,
            );
            row(ui, format!("{prefix}power"), Quantity::Duty, s.power_w);
        }
        row(
            ui,
            format!("{prefix}isentropic power"),
            Quantity::Duty,
            s.isentropic_power_w,
        );
        row(
            ui,
            format!("{prefix}isentropic outlet T"),
            Quantity::Temperature,
            s.isentropic_outlet_temperature_k,
        );
        if multistage {
            row(
                ui,
                format!("{prefix}discharge T"),
                Quantity::Temperature,
                s.outlet_temperature_k,
            );
        }
    }
    if multistage {
        row(
            ui,
            "Intercooler duty (+ removed)".into(),
            Quantity::Duty,
            m.intercooler_duty_w,
        );
    }
}

/// Equal mole fractions for `n > 0` components, summing to exactly one.
fn equal_fractions(n: usize) -> Vec<f64> {
    let mut equal = vec![1.0 / n as f64; n];
    if let Some((last, rest)) = equal.split_last_mut() {
        *last = 1.0 - rest.iter().sum::<f64>();
    }
    equal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_fractions_sum_to_one_and_handle_zero_components() {
        assert!(equal_fractions(0).is_empty());
        assert_eq!(equal_fractions(1), vec![1.0]);
        let x = equal_fractions(3);
        assert_eq!(x.len(), 3);
        assert_eq!(x.iter().sum::<f64>(), 1.0);
    }
}
