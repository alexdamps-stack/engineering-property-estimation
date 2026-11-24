//! Project properties view: file, property package and provenance.

use eframe::egui::{self, RichText, Ui};
use flowsheet_project::{PackageChoice, PackageSelection};

use crate::app::DIALOG_INTERACTIONS;
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
                "Create a project to choose components and a property package, \
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
        let label = ui.label("Model");
        egui::ComboBox::from_id_salt("project_property_model")
            .selected_text(variant.label())
            .show_ui(ui, |ui| {
                for model in PackageChoice::ALL {
                    let compatible = model.uses_kij()
                        || PackageSelection::available_components_for(model)
                            .is_ok_and(|ids| package.components.iter().all(|id| ids.contains(id)));
                    ui.add_enabled_ui(compatible, |ui| {
                        ui.selectable_value(&mut variant, model, model.label());
                    });
                }
            })
            .response
            .labelled_by(label.id);
    });
    if variant != package.variant
        && let Some(p) = &mut ed.package
    {
        // A different package makes every result stale (fingerprints).
        *p = p.with_model(variant, DIALOG_INTERACTIONS);
    }
    egui::Grid::new("package")
        .num_columns(2)
        .max_col_width(220.0)
        .show(ui, |ui| {
            ui.label("Components");
            ui.add(egui::Label::new(package.components.join(", ")).wrap());
            ui.end_row();
            ui.label("Interactions");
            ui.label(match package.variant {
                PackageChoice::Gerg2008 => "Embedded GERG-2008 data",
                PackageChoice::EosCg2021 => "EOS-CG-2021 CO₂/H₂ coefficients",
                PackageChoice::NrtlIdealGas => "Embedded methanol–water NRTL data",
                _ => "kij matrix, see below",
            });
            ui.end_row();
            ui.label("PH interval");
            ui.label(format!(
                "{} – {} K",
                package.ph_minimum_temperature_k, package.ph_maximum_temperature_k
            ));
            ui.end_row();
            if let Some(prov) = ed.package_cache.get(&package).provenance.clone() {
                ui.label("Implementation");
                ui.add(egui::Label::new(format!("{} ({})", prov.model, prov.model_version)).wrap());
                ui.end_row();
                ui.label("Data");
                ui.add(egui::Label::new(format!("{} ({})", prov.data, prov.data_version)).wrap());
                ui.end_row();
                ui.label("Engine");
                ui.add(egui::Label::new(prov.engine).wrap());
                ui.end_row();
                if let Some(reference) = prov.caloric_reference {
                    ui.label("Caloric reference");
                    ui.add(egui::Label::new(reference).wrap());
                    ui.end_row();
                }
            }
        });
    if package.variant.uses_kij() {
        ui.add(egui::Label::new(RichText::new(&package.interaction_source).small()).wrap());
    }
    ui.add_space(8.0);
    ui.label(RichText::new("Select a unit or stream to edit it.").weak());
    false
}
