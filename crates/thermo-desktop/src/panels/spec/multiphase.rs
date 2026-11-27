//! Specification forms of the multiphase group; owned by the multiphase
//! work package.

use eframe::egui::Ui;
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::columns::{ColumnInternals, MAX_COLUMN_STAGES};
use flowsheet_core::units::multiphase::{
    DecanterReflux, ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorSpec,
    MAX_EXTRACTOR_STAGES, PhaseReflux, ThreePhaseColumnSpec,
};
use flowsheet_core::units::{GroupSpec, MultiphaseSpec};

use super::{choice, count_field, drawn_only_note, note_row, opt_field};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`MultiphaseSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut MultiphaseSpec, u: DisplayUnits) {
    let kind = spec.kind();
    match spec {
        MultiphaseSpec::LiquidExtractor(s) => extractor_form(ui, s, u),
        MultiphaseSpec::ThreePhaseColumn(s) => three_phase_column_form(ui, s, u),
        _ => {}
    }
    drawn_only_note(ui, kind);
}

fn three_phase_column_form(ui: &mut Ui, s: &mut ThreePhaseColumnSpec, u: DisplayUnits) {
    let packed = matches!(s.internals, ColumnInternals::Packed(_)) as usize;
    match choice(
        ui,
        "Internals",
        packed,
        &["Ideal trays", "Packed bed (HETP)"],
    ) {
        Some(0) => s.internals = ColumnInternals::Trays(5),
        Some(_) => s.internals = ColumnInternals::Packed(PackedSectionSpec::default()),
        None => {}
    }
    match &mut s.internals {
        ColumnInternals::Trays(count) => {
            count_field(ui, "Internal trays", count, 1..=MAX_COLUMN_STAGES);
        }
        ColumnInternals::Packed(spec) => packed_form(ui, spec, u),
    }
    let last = s.stage_count().map_or(MAX_COLUMN_STAGES + 2, |n| n as u16);
    count_field(ui, "Feed stage", &mut s.feed_stage, 2..=last);
    opt_field(
        ui,
        "Column pressure",
        &mut s.pressure_pa,
        Quantity::Pressure,
        u,
    );
    opt_field(
        ui,
        "Decanter temperature",
        &mut s.decanter_temperature_k,
        Quantity::Temperature,
        u,
    );
    let current = match s.reflux {
        DecanterReflux::Light(_) => 0,
        DecanterReflux::Heavy(_) => 1,
        DecanterReflux::Fractions { .. } => 2,
    };
    match choice(
        ui,
        "Reflux liquid",
        current,
        &["Light liquid", "Heavy liquid", "Fractions of both"],
    ) {
        Some(0) => s.reflux = DecanterReflux::Light(PhaseReflux::Total),
        Some(1) => s.reflux = DecanterReflux::Heavy(PhaseReflux::Total),
        Some(_) => {
            s.reflux = DecanterReflux::Fractions {
                light: None,
                heavy: None,
            }
        }
        None => {}
    }
    match &mut s.reflux {
        DecanterReflux::Light(amount) | DecanterReflux::Heavy(amount) => {
            let ratio = matches!(amount, PhaseReflux::Ratio(_)) as usize;
            match choice(
                ui,
                "Reflux amount",
                ratio,
                &["Whole liquid (no distillate)", "Reflux ratio"],
            ) {
                Some(0) => *amount = PhaseReflux::Total,
                Some(_) => *amount = PhaseReflux::Ratio(None),
                None => {}
            }
            if let PhaseReflux::Ratio(r) = amount {
                opt_field(ui, "Reflux ratio L/D", r, Quantity::Factor, u);
            }
        }
        DecanterReflux::Fractions { light, heavy } => {
            opt_field(ui, "Light liquid refluxed", light, Quantity::Fraction, u);
            opt_field(ui, "Heavy liquid refluxed", heavy, Quantity::Fraction, u);
        }
    }
    opt_field(
        ui,
        "Feed to bottoms",
        &mut s.bottoms_fraction,
        Quantity::Fraction,
        u,
    );
    note_row(
        ui,
        "Stage 1 is the total condenser with a decanter, the last stage the partial reboiler. Distillate: unrefluxed part of the reflux liquid (light liquid for fractions); decant: the other liquid. Only the decanter may hold two liquids; needs a VLLE/LLE package (PR, PR-BM, NRTL-IG).",
    );
}

fn extractor_form(ui: &mut Ui, s: &mut LiquidExtractorSpec, u: DisplayUnits) {
    let packed = matches!(s.internals, ExtractorInternals::Packed(_)) as usize;
    match choice(
        ui,
        "Internals",
        packed,
        &["Ideal stages", "Packed section (HETP)"],
    ) {
        Some(0) => s.internals = ExtractorInternals::Stages(5),
        Some(_) => s.internals = ExtractorInternals::Packed(PackedSectionSpec::default()),
        None => {}
    }
    match &mut s.internals {
        ExtractorInternals::Stages(count) => {
            count_field(ui, "Ideal stages", count, 1..=MAX_EXTRACTOR_STAGES);
        }
        ExtractorInternals::Packed(spec) => packed_form(ui, spec, u),
    }
    opt_field(
        ui,
        "Column pressure",
        &mut s.pressure_pa,
        Quantity::Pressure,
        u,
    );
    let adiabatic = matches!(s.thermal, ExtractorThermal::Adiabatic) as usize;
    match choice(ui, "Thermal mode", adiabatic, &["Isothermal", "Adiabatic"]) {
        Some(0) => {
            s.thermal = ExtractorThermal::Isothermal {
                temperature_k: None,
            }
        }
        Some(_) => s.thermal = ExtractorThermal::Adiabatic,
        None => {}
    }
    if let ExtractorThermal::Isothermal { temperature_k } = &mut s.thermal {
        opt_field(ui, "Temperature", temperature_k, Quantity::Temperature, u);
    }
    let phases = [FeedPhase::Detect, FeedPhase::Heavy, FeedPhase::Light];
    let current = phases
        .iter()
        .position(|p| *p == s.feed_phase)
        .expect("known feed phase");
    if let Some(index) = choice(ui, "Feed phase", current, &["Detect", "Heavy", "Light"]) {
        s.feed_phase = phases[index];
    }
    note_row(
        ui,
        "Countercurrent LLE stages from the feed end: the feed enters stage 1, the solvent the last stage. Needs a VLLE/LLE package (PR, PR-BM, NRTL-IG).",
    );
}

fn packed_form(ui: &mut Ui, spec: &mut PackedSectionSpec, u: DisplayUnits) {
    opt_field(ui, "HETP", &mut spec.hetp_m, Quantity::Length, u);
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
            count_field(ui, "Packed stages", count, 1..=MAX_EXTRACTOR_STAGES);
        }
        PackedSectionBasis::PackedHeight { height_m, rounding } => {
            opt_field(ui, "Packed height", height_m, Quantity::Length, u);
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
