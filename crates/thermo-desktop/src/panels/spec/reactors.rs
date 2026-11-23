use eframe::egui::Ui;
use flowsheet_core::{BatchThermal, FugacityBasis, PlugFlowThermal, UnitSpec};

use super::{choice, note_row, opt_field, value_field, vessel_thermal_fields};
use crate::units::{DisplayUnits, Quantity};

pub(super) fn form(ui: &mut Ui, spec: &mut UnitSpec, u: DisplayUnits) {
    if let UnitSpec::Cstr { volume_m3, .. } = spec {
        opt_field(ui, "Volume", volume_m3, Quantity::Volume, u);
    }
    if let UnitSpec::EquilibriumReactor { equilibrium, .. }
    | UnitSpec::GibbsReactor { equilibrium, .. } = spec
    {
        let current = (equilibrium.fugacity == FugacityBasis::IdealGas) as usize;
        match choice(ui, "Fugacity", current, &["Package", "Ideal gas"]) {
            Some(0) => equilibrium.fugacity = FugacityBasis::Package,
            Some(_) => equilibrium.fugacity = FugacityBasis::IdealGas,
            None => {}
        }
        value_field(
            ui,
            "Temperature approach",
            &mut equilibrium.temperature_approach_k,
            Quantity::TemperatureDifference,
            u,
        );
    }
    match spec {
        UnitSpec::ConversionReactor {
            thermal,
            pressure_drop_pa,
            ..
        }
        | UnitSpec::YieldReactor {
            thermal,
            pressure_drop_pa,
            ..
        }
        | UnitSpec::EquilibriumReactor {
            thermal,
            pressure_drop_pa,
            ..
        }
        | UnitSpec::GibbsReactor {
            thermal,
            pressure_drop_pa,
            ..
        }
        | UnitSpec::Cstr {
            thermal,
            pressure_drop_pa,
            ..
        } => vessel_thermal_fields(ui, thermal, pressure_drop_pa, u),
        UnitSpec::PlugFlow {
            length_m,
            diameter_m,
            thermal,
            pressure_drop_pa,
            ..
        } => {
            opt_field(ui, "Length", length_m, Quantity::Length, u);
            opt_field(ui, "Diameter", diameter_m, Quantity::Length, u);
            let current = match thermal {
                PlugFlowThermal::Isothermal(_) => 0,
                PlugFlowThermal::Adiabatic => 1,
                PlugFlowThermal::Coolant { .. } => 2,
            };
            match choice(
                ui,
                "Thermal",
                current,
                &["Isothermal", "Adiabatic", "Coolant"],
            ) {
                Some(0) => *thermal = PlugFlowThermal::Isothermal(None),
                Some(1) => *thermal = PlugFlowThermal::Adiabatic,
                Some(_) => {
                    *thermal = PlugFlowThermal::Coolant {
                        temperature_k: None,
                        u_w_m2_k: None,
                    }
                }
                None => {}
            }
            match thermal {
                PlugFlowThermal::Isothermal(t) => {
                    opt_field(ui, "Temperature", t, Quantity::Temperature, u)
                }
                PlugFlowThermal::Adiabatic => {}
                PlugFlowThermal::Coolant {
                    temperature_k,
                    u_w_m2_k,
                } => {
                    opt_field(
                        ui,
                        "Coolant temperature",
                        temperature_k,
                        Quantity::Temperature,
                        u,
                    );
                    opt_field(
                        ui,
                        "U (wall)",
                        u_w_m2_k,
                        Quantity::HeatTransferCoefficient,
                        u,
                    );
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
        UnitSpec::Batch {
            reaction_time_s,
            thermal,
            ..
        } => {
            opt_field(ui, "Reaction time", reaction_time_s, Quantity::Time, u);
            let current = matches!(thermal, BatchThermal::Adiabatic) as usize;
            match choice(ui, "Thermal", current, &["Isothermal", "Adiabatic"]) {
                Some(0) => *thermal = BatchThermal::Isothermal(None),
                Some(_) => *thermal = BatchThermal::Adiabatic,
                None => {}
            }
            if let BatchThermal::Isothermal(t) = thermal {
                opt_field(ui, "Temperature", t, Quantity::Temperature, u);
            }
            note_row(
                ui,
                "Constant pressure; outlet is the time-averaged product.",
            );
        }
        _ => unreachable!("expected a reactor specification"),
    }
}
