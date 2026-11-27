//! Specification forms of the columns group; owned by the columns work
//! package.

use eframe::egui::{ComboBox, Ui};
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::columns::{
    ColumnInternals, DistillationColumnSpec, MAX_COLUMN_PORTS, MAX_COLUMN_STAGES,
    MAX_DISTILLATION_INTERNAL_STAGES, ShortcutColumnSpec, ShortcutCondenser, ShortcutReflux,
};
use flowsheet_core::units::{ColumnSpec, GroupSpec};

use super::{choice, count_field, drawn_only_note, note_row, opt_field};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`ColumnSpec`]. Feeds and side draws are numbered ports; a
/// connected port cannot be removed (disconnect it first).
pub(super) fn form(ui: &mut Ui, spec: &mut ColumnSpec, u: DisplayUnits, components: &[String]) {
    let kind = spec.kind();
    match spec {
        ColumnSpec::Shortcut(s) => shortcut_form(ui, s, u, components),
        ColumnSpec::Distillation(s) => distillation_form(ui, s, u),
        ColumnSpec::PetroleumFractionator(s) => {
            count_field(ui, "Side draws", &mut s.side_draws, 0..=MAX_COLUMN_PORTS);
        }
        ColumnSpec::Absorber(s) => {
            internals_form(ui, &mut s.internals, u);
            opt_field(
                ui,
                "Column pressure",
                &mut s.pressure_pa,
                Quantity::Pressure,
                u,
            );
            note_row(
                ui,
                "Adiabatic contact. Liquid enters stage 1 and gas enters the last stage.",
            );
        }
    }
    drawn_only_note(ui, kind);
}

fn distillation_form(ui: &mut Ui, s: &mut DistillationColumnSpec, units: DisplayUnits) {
    count_field(ui, "Feeds", &mut s.feeds, 1..=MAX_COLUMN_PORTS);
    count_field(
        ui,
        "Liquid side draws",
        &mut s.side_draws,
        0..=MAX_COLUMN_PORTS,
    );
    s.feed_stages.truncate(usize::from(s.feeds));
    s.liquid_draws.truncate(usize::from(s.side_draws));

    internals_form(ui, &mut s.internals, units);
    opt_field(
        ui,
        "Column pressure",
        &mut s.pressure_pa,
        Quantity::Pressure,
        units,
    );
    opt_field(
        ui,
        "Reflux ratio L/D",
        &mut s.reflux_ratio,
        Quantity::Factor,
        units,
    );
    opt_field(
        ui,
        "Feed to distillate",
        &mut s.distillate_fraction,
        Quantity::Fraction,
        units,
    );
    note_row(
        ui,
        "Stage 1 is the total condenser. The last stage is the partial reboiler.",
    );
    let stage_range = 2..=MAX_DISTILLATION_INTERNAL_STAGES + 2;
    for i in 0..usize::from(s.feeds.min(MAX_COLUMN_PORTS)) {
        let mut stage = s.feed_stages.get(i).copied().unwrap_or(2);
        if count_field(
            ui,
            &format!("Feed {} stage", i + 1),
            &mut stage,
            stage_range.clone(),
        ) {
            s.feed_stages.resize(s.feed_stages.len().max(i + 1), 2);
            s.feed_stages[i] = stage;
        }
    }
    for i in 0..usize::from(s.side_draws.min(MAX_COLUMN_PORTS)) {
        let original = s.liquid_draws.get(i).copied().unwrap_or_default();
        let mut draw = original;
        count_field(
            ui,
            &format!("Draw {} stage", i + 1),
            &mut draw.stage,
            stage_range.clone(),
        );
        ui.push_id(("liquid side draw", i), |ui| {
            opt_field(
                ui,
                &format!("Draw {} molar flow", i + 1),
                &mut draw.molar_flow_mol_s,
                Quantity::MolarFlow,
                units,
            );
        });
        if draw != original {
            s.liquid_draws
                .resize(s.liquid_draws.len().max(i + 1), Default::default());
            s.liquid_draws[i] = draw;
        }
    }
}

fn internals_form(ui: &mut Ui, internals: &mut ColumnInternals, units: DisplayUnits) {
    let packed = matches!(*internals, ColumnInternals::Packed(_)) as usize;
    match choice(
        ui,
        "Internals",
        packed,
        &["Ideal trays", "Packed bed (HETP)"],
    ) {
        Some(0) => *internals = ColumnInternals::Trays(5),
        Some(_) => *internals = ColumnInternals::Packed(PackedSectionSpec::default()),
        None => {}
    }
    match internals {
        ColumnInternals::Trays(count) => {
            count_field(ui, "Internal trays", count, 1..=MAX_COLUMN_STAGES);
        }
        ColumnInternals::Packed(spec) => packed_form(ui, spec, units),
    }
}

fn packed_form(ui: &mut Ui, spec: &mut PackedSectionSpec, units: DisplayUnits) {
    opt_field(ui, "HETP", &mut spec.hetp_m, Quantity::Length, units);
    let stages = matches!(spec.basis, PackedSectionBasis::TheoreticalStages(_)) as usize;
    match choice(
        ui,
        "Packing basis",
        stages,
        &["Height", "Theoretical stages"],
    ) {
        Some(0) => spec.basis = PackedSectionBasis::default(),
        Some(_) => spec.basis = PackedSectionBasis::TheoreticalStages(5),
        None => {}
    }
    match &mut spec.basis {
        PackedSectionBasis::TheoreticalStages(count) => {
            count_field(ui, "Packed stages", count, 1..=MAX_COLUMN_STAGES);
        }
        PackedSectionBasis::PackedHeight { height_m, rounding } => {
            opt_field(ui, "Packed height", height_m, Quantity::Length, units);
            let choices = [
                StageRounding::Nearest,
                StageRounding::Down,
                StageRounding::Up,
            ];
            let current = choices
                .iter()
                .position(|value| value == rounding)
                .expect("known rounding rule");
            if let Some(index) = choice(ui, "Stage rounding", current, &["Nearest", "Down", "Up"]) {
                *rounding = choices[index];
            }
        }
    }
}

fn key_field(ui: &mut Ui, label: &str, key: &mut Option<String>, components: &[String]) {
    let text = ui.label(label);
    let shown = match key.as_deref() {
        None => "Choose a component".to_owned(),
        Some(id) if components.iter().any(|c| c == id) => id.to_owned(),
        Some(id) => format!("{id} (outside package)"),
    };
    ComboBox::from_id_salt(label)
        .selected_text(shown)
        .show_ui(ui, |ui| {
            ui.selectable_value(key, None, "Choose a component");
            for id in components {
                ui.selectable_value(key, Some(id.clone()), id);
            }
        })
        .response
        .labelled_by(text.id);
    ui.end_row();
}

fn shortcut_form(ui: &mut Ui, s: &mut ShortcutColumnSpec, u: DisplayUnits, components: &[String]) {
    key_field(ui, "Light key", &mut s.light_key, components);
    key_field(ui, "Heavy key", &mut s.heavy_key, components);
    opt_field(
        ui,
        "Light key to distillate",
        &mut s.light_key_recovery,
        Quantity::Fraction,
        u,
    );
    opt_field(
        ui,
        "Heavy key to bottoms",
        &mut s.heavy_key_recovery,
        Quantity::Fraction,
        u,
    );
    opt_field(
        ui,
        "Column pressure",
        &mut s.pressure_pa,
        Quantity::Pressure,
        u,
    );
    let current = match s.condenser {
        ShortcutCondenser::Total => 0,
        ShortcutCondenser::Partial => 1,
    };
    match choice(
        ui,
        "Condenser",
        current,
        &["Total (liquid distillate)", "Partial (vapour distillate)"],
    ) {
        Some(0) => s.condenser = ShortcutCondenser::Total,
        Some(_) => s.condenser = ShortcutCondenser::Partial,
        None => {}
    }
    let current = match s.reflux {
        ShortcutReflux::MinimumMultiple(_) => 0,
        ShortcutReflux::Ratio(_) => 1,
        ShortcutReflux::TheoreticalStages(_) => 2,
    };
    match choice(
        ui,
        "Specify",
        current,
        &["Multiple of minimum", "Reflux ratio", "Theoretical stages"],
    ) {
        Some(0) => s.reflux = ShortcutReflux::MinimumMultiple(None),
        Some(1) => s.reflux = ShortcutReflux::Ratio(None),
        Some(_) => s.reflux = ShortcutReflux::TheoreticalStages(None),
        None => {}
    }
    match &mut s.reflux {
        ShortcutReflux::MinimumMultiple(m) => {
            opt_field(ui, "R / R_min (> 1)", m, Quantity::Factor, u)
        }
        ShortcutReflux::Ratio(r) => opt_field(ui, "Reflux ratio L/D", r, Quantity::Factor, u),
        ShortcutReflux::TheoreticalStages(n) => {
            opt_field(ui, "Theoretical stages (> N_min)", n, Quantity::Factor, u)
        }
    }
    note_row(
        ui,
        match s.condenser {
            ShortcutCondenser::Total => {
                "Stages include the reboiler; the total condenser is not a stage."
            }
            ShortcutCondenser::Partial => {
                "Stages include the reboiler and the partial condenser (stage 1)."
            }
        },
    );
    note_row(
        ui,
        "Fenske–Underwood–Gilliland estimate. Recoveries lie strictly between 0 and 1 and sum to more than 1.",
    );
    note_row(
        ui,
        "Volatilities come from the package. Keys must be adjacent in volatility.",
    );
    note_row(
        ui,
        "Condenser and reboiler at column pressure. Bottoms leave as saturated liquid; the distillate as saturated liquid (total) or vapour (partial).",
    );
    note_row(
        ui,
        "Feed components are limited to methane, ethane, propane, nitrogen and CO2 between 200 and 350 K.",
    );
}
