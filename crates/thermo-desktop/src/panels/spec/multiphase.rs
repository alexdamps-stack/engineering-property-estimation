//! Specification forms of the multiphase group; owned by the multiphase
//! work package.

use eframe::egui::Ui;
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::multiphase::{
    ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorSpec, MAX_EXTRACTOR_STAGES,
};
use flowsheet_core::units::{GroupSpec, MultiphaseSpec};

use super::{choice, count_field, drawn_only_note, note_row, opt_field};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`MultiphaseSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut MultiphaseSpec, u: DisplayUnits) {
    let kind = spec.kind();
    if let MultiphaseSpec::LiquidExtractor(s) = spec {
        extractor_form(ui, s, u);
    }
    drawn_only_note(ui, kind);
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
