//! Project properties view: file, property package and provenance.

use eframe::egui::{self, RichText, Ui};
use flowsheet_project::CubicVariantChoice;

use crate::editor::Editor;

pub(super) fn project_properties(ui: &mut Ui, ed: &mut Editor) -> bool {
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
