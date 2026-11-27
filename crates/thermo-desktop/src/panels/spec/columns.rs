//! Specification forms of the columns group; owned by the columns work
//! package.

use eframe::egui::{ComboBox, Ui};
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::columns::{
    BottomSpecification, ColumnInternals, ColumnPressureDrop, ComponentTarget,
    DistillationColumnSpec, DistillationCondenser, DrawRate, FractionatorBottom, MAX_COLUMN_PORTS,
    MAX_COLUMN_STAGES, MAX_DISTILLATION_INTERNAL_STAGES, PetroleumFractionatorSpec,
    ShortcutColumnSpec, ShortcutCondenser, ShortcutReflux, TopSpecification,
};
use flowsheet_core::units::{ColumnSpec, GroupSpec};

use super::{choice, count_field, drawn_only_note, note_row, opt_field, value_field};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`ColumnSpec`]. Feeds and side draws are numbered ports; a
/// connected port cannot be removed (disconnect it first).
pub(super) fn form(ui: &mut Ui, spec: &mut ColumnSpec, u: DisplayUnits, components: &[String]) {
    let kind = spec.kind();
    match spec {
        ColumnSpec::Shortcut(s) => shortcut_form(ui, s, u, components),
        ColumnSpec::Distillation(s) => distillation_form(ui, s, u, components),
        ColumnSpec::PetroleumFractionator(s) => petroleum_form(ui, s, u, components),
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

fn distillation_form(
    ui: &mut Ui,
    s: &mut DistillationColumnSpec,
    units: DisplayUnits,
    components: &[String],
) {
    count_field(ui, "Feeds", &mut s.feeds, 1..=MAX_COLUMN_PORTS);
    count_field(
        ui,
        "Liquid side draws",
        &mut s.side_draws,
        0..=MAX_COLUMN_PORTS,
    );
    s.feed_stages.truncate(usize::from(s.feeds));
    s.liquid_draws.truncate(usize::from(s.side_draws));

    count_field(
        ui,
        "Vapor side draws",
        &mut s.vapor_side_draws,
        0..=MAX_COLUMN_PORTS,
    );
    s.vapor_draws.truncate(usize::from(s.vapor_side_draws));

    internals_form(ui, &mut s.internals, units);
    condenser_form(ui, &mut s.condenser, units);
    pressure_form(
        ui,
        &mut s.pressure_pa,
        &mut s.condenser_pressure_drop_pa,
        &mut s.pressure_drop,
        units,
    );
    top_form(ui, &mut s.top, units, components);
    bottom_form(ui, &mut s.bottom, units, components);
    note_row(
        ui,
        "Stage 1 is the condenser. The last stage is the partial reboiler. \
         Specifications other than flows and the reflux ratio are met by iteration.",
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
    for i in 0..usize::from(s.vapor_side_draws.min(MAX_COLUMN_PORTS)) {
        let original = s.vapor_draws.get(i).copied().unwrap_or_default();
        let mut draw = original;
        count_field(
            ui,
            &format!("Vapor draw {} stage", i + 1),
            &mut draw.stage,
            stage_range.clone(),
        );
        ui.push_id(("vapor side draw", i), |ui| {
            opt_field(
                ui,
                &format!("Vapor draw {} molar flow", i + 1),
                &mut draw.molar_flow_mol_s,
                Quantity::MolarFlow,
                units,
            );
        });
        if draw != original {
            s.vapor_draws
                .resize(s.vapor_draws.len().max(i + 1), Default::default());
            s.vapor_draws[i] = draw;
        }
    }
}

/// Petroleum fractionator (#190): the rigorous column with pseudo-component
/// feeds, side products by molar or standard liquid volume flow, and a
/// reboiler or stripping steam.
fn petroleum_form(
    ui: &mut Ui,
    s: &mut PetroleumFractionatorSpec,
    units: DisplayUnits,
    components: &[String],
) {
    count_field(ui, "Feeds", &mut s.feeds, 1..=MAX_COLUMN_PORTS);
    count_field(ui, "Side products", &mut s.side_draws, 0..=MAX_COLUMN_PORTS);
    s.feed_stages.truncate(usize::from(s.feeds));
    s.products.truncate(usize::from(s.side_draws));
    internals_form(ui, &mut s.internals, units);
    condenser_form(ui, &mut s.condenser, units);
    pressure_form(
        ui,
        &mut s.pressure_pa,
        &mut s.condenser_pressure_drop_pa,
        &mut s.pressure_drop,
        units,
    );
    top_form(ui, &mut s.top, units, components);
    let steam = s.steam() as usize;
    match choice(
        ui,
        "Bottom",
        steam,
        &["Reboiler", "Stripping steam (no reboiler)"],
    ) {
        Some(0) if steam != 0 => s.bottom = FractionatorBottom::default(),
        Some(1) if steam != 1 => s.bottom = FractionatorBottom::StrippingSteam,
        _ => {}
    }
    match &mut s.bottom {
        FractionatorBottom::Reboiler(bottom) => bottom_form(ui, bottom, units, components),
        FractionatorBottom::StrippingSteam => note_row(
            ui,
            "Steam enters the adiabatic bottom stage; the top specification is the only one.",
        ),
    }
    note_row(
        ui,
        "Stage 1 is the condenser and the last stage the reboiler or steam-stripped bottom tray. \
         Feeds must carry characterised pseudo-components. Side strippers and pump-arounds \
         are not modelled.",
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
        let original = s.products.get(i).copied().unwrap_or_default();
        let mut product = original;
        count_field(
            ui,
            &format!("Product {} stage", i + 1),
            &mut product.stage,
            stage_range.clone(),
        );
        ui.push_id(("side product", i), |ui| {
            let volume = matches!(product.rate, DrawRate::StandardLiquidVolume(_)) as usize;
            match choice(
                ui,
                &format!("Product {} basis", i + 1),
                volume,
                &["Molar flow", "Standard liquid volume (60 °F)"],
            ) {
                Some(0) if volume != 0 => product.rate = DrawRate::Molar(None),
                Some(1) if volume != 1 => product.rate = DrawRate::StandardLiquidVolume(None),
                _ => {}
            }
            match &mut product.rate {
                DrawRate::Molar(v) => opt_field(
                    ui,
                    &format!("Product {} molar flow", i + 1),
                    v,
                    Quantity::MolarFlow,
                    units,
                ),
                DrawRate::StandardLiquidVolume(v) => opt_field(
                    ui,
                    &format!("Product {} volume flow", i + 1),
                    v,
                    Quantity::VolumeFlow,
                    units,
                ),
            }
        });
        if product != original {
            s.products
                .resize(s.products.len().max(i + 1), Default::default());
            s.products[i] = product;
        }
    }
}

fn condenser_form(ui: &mut Ui, condenser: &mut DistillationCondenser, units: DisplayUnits) {
    let current = match *condenser {
        DistillationCondenser::Total => 0,
        DistillationCondenser::PartialVapor => 1,
        DistillationCondenser::PartialVaporLiquid { .. } => 2,
    };
    match choice(
        ui,
        "Condenser",
        current,
        &["Total", "Partial (vapor)", "Partial (vapor and liquid)"],
    ) {
        Some(0) => *condenser = DistillationCondenser::Total,
        Some(1) => *condenser = DistillationCondenser::PartialVapor,
        Some(_) if current != 2 => {
            *condenser = DistillationCondenser::PartialVaporLiquid {
                liquid_distillate_mol_s: None,
            }
        }
        _ => {}
    }
    if let DistillationCondenser::PartialVaporLiquid {
        liquid_distillate_mol_s,
    } = condenser
    {
        opt_field(
            ui,
            "Liquid distillate",
            liquid_distillate_mol_s,
            Quantity::MolarFlow,
            units,
        );
    }
}

fn pressure_form(
    ui: &mut Ui,
    pressure_pa: &mut Option<f64>,
    condenser_drop_pa: &mut f64,
    pressure_drop: &mut ColumnPressureDrop,
    units: DisplayUnits,
) {
    opt_field(
        ui,
        "Condenser pressure",
        pressure_pa,
        Quantity::Pressure,
        units,
    );
    value_field(
        ui,
        "Condenser pressure drop",
        condenser_drop_pa,
        Quantity::PressureDrop,
        units,
    );
    let current = matches!(pressure_drop, ColumnPressureDrop::Column(_)) as usize;
    match choice(
        ui,
        "Stage pressure drop",
        current,
        &["Per stage", "Whole column"],
    ) {
        Some(0) if current != 0 => *pressure_drop = ColumnPressureDrop::PerStage(0.0),
        Some(1) if current != 1 => *pressure_drop = ColumnPressureDrop::Column(0.0),
        _ => {}
    }
    let (label, drop) = match pressure_drop {
        ColumnPressureDrop::PerStage(drop) => ("Drop per stage", drop),
        ColumnPressureDrop::Column(drop) => ("Drop first tray to reboiler", drop),
    };
    value_field(ui, label, drop, Quantity::PressureDrop, units);
}

fn target_form(
    ui: &mut Ui,
    label: &str,
    target: &mut ComponentTarget,
    units: DisplayUnits,
    components: &[String],
) {
    key_field(
        ui,
        &format!("{label} component"),
        &mut target.component,
        components,
    );
    opt_field(ui, label, &mut target.value, Quantity::Fraction, units);
}

fn top_form(ui: &mut Ui, top: &mut TopSpecification, u: DisplayUnits, components: &[String]) {
    use TopSpecification as T;
    let current = match top {
        T::RefluxRatio(_) => 0,
        T::RefluxRate(_) => 1,
        T::DistillateRate(_) => 2,
        T::CondenserDuty(_) => 3,
        T::DistillatePurity(_) => 4,
        T::DistillateRecovery(_) => 5,
    };
    if let Some(picked) = choice(
        ui,
        "Top specification",
        current,
        &[
            "Reflux ratio",
            "Reflux rate",
            "Distillate rate",
            "Condenser duty",
            "Distillate purity",
            "Distillate recovery",
        ],
    ) && picked != current
    {
        *top = match picked {
            0 => T::RefluxRatio(None),
            1 => T::RefluxRate(None),
            2 => T::DistillateRate(None),
            3 => T::CondenserDuty(None),
            4 => T::DistillatePurity(Default::default()),
            _ => T::DistillateRecovery(Default::default()),
        };
    }
    match top {
        T::RefluxRatio(v) => opt_field(ui, "Reflux ratio L/D", v, Quantity::Factor, u),
        T::RefluxRate(v) => opt_field(ui, "Reflux molar flow", v, Quantity::MolarFlow, u),
        T::DistillateRate(v) => opt_field(ui, "Distillate molar flow", v, Quantity::MolarFlow, u),
        T::CondenserDuty(v) => opt_field(ui, "Condenser duty (removed)", v, Quantity::Duty, u),
        T::DistillatePurity(t) => target_form(ui, "Distillate mole fraction", t, u, components),
        T::DistillateRecovery(t) => target_form(ui, "Recovery to distillate", t, u, components),
    }
}

fn bottom_form(
    ui: &mut Ui,
    bottom: &mut BottomSpecification,
    u: DisplayUnits,
    components: &[String],
) {
    use BottomSpecification as B;
    let current = match bottom {
        B::DistillateFraction(_) => 0,
        B::BottomsRate(_) => 1,
        B::ReboilerDuty(_) => 2,
        B::BoilupRatio(_) => 3,
        B::BottomsPurity(_) => 4,
        B::BottomsRecovery(_) => 5,
    };
    if let Some(picked) = choice(
        ui,
        "Bottom specification",
        current,
        &[
            "Distillate to feed",
            "Bottoms rate",
            "Reboiler duty",
            "Boilup ratio",
            "Bottoms purity",
            "Bottoms recovery",
        ],
    ) && picked != current
    {
        *bottom = match picked {
            0 => B::DistillateFraction(None),
            1 => B::BottomsRate(None),
            2 => B::ReboilerDuty(None),
            3 => B::BoilupRatio(None),
            4 => B::BottomsPurity(Default::default()),
            _ => B::BottomsRecovery(Default::default()),
        };
    }
    match bottom {
        B::DistillateFraction(v) => opt_field(ui, "Feed to distillate", v, Quantity::Fraction, u),
        B::BottomsRate(v) => opt_field(ui, "Bottoms molar flow", v, Quantity::MolarFlow, u),
        B::ReboilerDuty(v) => opt_field(ui, "Reboiler duty (supplied)", v, Quantity::Duty, u),
        B::BoilupRatio(v) => opt_field(ui, "Boilup ratio V/B", v, Quantity::Factor, u),
        B::BottomsPurity(t) => target_form(ui, "Bottoms mole fraction", t, u, components),
        B::BottomsRecovery(t) => target_form(ui, "Recovery to bottoms", t, u, components),
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
