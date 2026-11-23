//! Stream properties view: connections, results, phases and provenance.

use eframe::egui::{self, RichText, Ui};
use flowsheet_core::{PortRef, StreamId, UnitId};
use flowsheet_project::format_utc;

use super::presentation::{StreamProperty, component_names, fraction_text, stream_result};
use super::{status_line, value_text};
use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;
use crate::units::{DisplayUnitsExt, Quantity, format_value};

pub(super) fn stream_properties(ui: &mut Ui, ed: &mut Editor, id: StreamId) {
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
    ui.add_space(8.0);
    ui.strong("Results");
    let theme = Theme::from_visuals(ui.visuals());
    let u = ed.units;
    let components = component_names(ed);
    egui::Frame::group(ui.style()).show(ui, |ui| {
        ui.set_width(ui.available_width());
        let Some((status, r)) = stream_result(ed, id, &fps) else {
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
                for property in StreamProperty::PANEL_FLOWS {
                    row(ui, property.label(), property.panel_text(r, u));
                }
                if r.state.is_some() {
                    for property in StreamProperty::PANEL_STATE {
                        row(ui, property.label(), property.panel_text(r, u));
                    }
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
                    ui.label(value_text(fraction_text(&s.mole_fractions, i), stale));
                    for p in &s.phases {
                        ui.label(value_text(fraction_text(&p.mole_fractions, i), stale));
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
