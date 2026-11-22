//! Unit properties view: tag, specification, feed composition, connections
//! and results.

use eframe::egui::{self, DragValue, RichText, Ui};
use flowsheet_core::{PortDirection, PortMedium, PortRef, SpecProblem, UnitId, UnitKind};
use flowsheet_project::Fingerprints;
use thermo_types::COMPOSITION_SUM_TOLERANCE;

use super::presentation::stream_summary;
use super::spec::spec_form;
use super::{status_line, value_text};
use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;
use crate::units::{Quantity, format_quantity, format_residual};

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
    egui::Grid::new(("spec", id))
        .num_columns(2)
        .spacing([8.0, 6.0])
        .show(ui, |ui| spec_form(ui, unit.kind, &mut spec, u));
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
                if kind == Some(UnitKind::HeatExchanger) {
                    // Internal heat moved from the hot to the cold side; the
                    // external duty of an exchanger is zero.
                    ui.label("Transferred duty (hot to cold)");
                    ui.label(match r.transferred_duty_w {
                        Some(q) => value_text(format_quantity(u, Quantity::Duty, q), stale),
                        None => RichText::new("not stored — recalculate").weak(),
                    });
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
