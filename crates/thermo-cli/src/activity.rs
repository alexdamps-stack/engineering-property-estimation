//! NRTL-IG (NRTL liquid + ideal-gas vapor) CLI cases. Interactions come only from
//! the embedded dataset; a missing pair or pure-liquid correlation is an error.
use crate::{Case, phase_label};
use serde_json::{Value, json};
use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_models::{ActivityModel, ActivityPhaseState, VaporModel};
use thermo_solvers::{
    ActivityFlashPhase, ActivityOptions, ActivitySaturation, activity_bubble_pressure,
    activity_bubble_temperature, activity_dew_pressure, activity_dew_temperature,
    activity_ph_flash, activity_tp_flash,
};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, Diagnostic, NormalizationPolicy, Pressure,
    Temperature, ValidityStatus,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn state_json(s: &ActivityPhaseState) -> Value {
    json!({
        "ln_fugacity_coefficients": s.ln_fugacity_coefficients,
        "ln_activity_coefficients": s.excess.as_ref().map(|e| e.ln_activity_coefficients.clone()),
        "excess_enthalpy_J_per_mol": s.excess.as_ref().map(|e| e.excess_enthalpy),
        "enthalpy_J_per_mol": s.enthalpy, "entropy_J_per_mol_K": s.entropy,
        "gibbs_energy_J_per_mol": s.gibbs_energy, "cp_J_per_mol_K": s.cp,
        "molar_volume_m3_per_mol": s.molar_volume, "mass_density_kg_per_m3": s.mass_density,
        "validity": validity(s.validity),
    })
}

fn validity(v: ValidityStatus) -> &'static str {
    match v {
        ValidityStatus::WithinDeclaredDomain => "within NRTL regression temperature interval",
        ValidityStatus::Extrapolated => {
            "NRTL parameters extrapolated outside their regression temperature interval"
        }
        ValidityStatus::Unassessed => "unassessed",
    }
}

fn phases_json(phases: &[ActivityFlashPhase]) -> Vec<Value> {
    phases
        .iter()
        .map(|ph| {
            let mut v = state_json(&ph.state);
            v["label"] = json!(phase_label(ph.phase));
            v["molar_phase_fraction"] = json!(ph.phase_fraction);
            v["mole_fractions"] = json!(ph.composition);
            v
        })
        .collect()
}

fn provenance(model: &ActivityModel) -> Value {
    let m = model.metadata();
    let c = ActivityModel::caloric_reference();
    json!({
        "model":{"id":m.identifier,"implementation":m.implementation_version,"formulation":m.formulation},
        "caloric_reference":{"id":c.identifier,"definition":c.definition},
    })
}

fn diagnostics(d: &[Diagnostic]) -> Vec<Value> {
    d.iter()
        .map(|d| json!({"code":d.code,"severity":format!("{:?}",d.severity),"message":d.message}))
        .collect()
}

fn saturation_json(
    case: &Case,
    kind: &str,
    s: &ActivitySaturation,
    model: &ActivityModel,
) -> Value {
    json!({
        "case": case.id, "calculation": kind, "components": case.components,
        "temperature_K": s.temperature_k, "pressure_Pa_absolute": s.pressure_pa,
        "liquid_mole_fractions": s.liquid, "vapor_mole_fractions": s.vapor,
        "convergence": {"iterations": s.iterations, "max_abs_ln_fugacity_residual": s.max_ln_fugacity_residual},
        "incipient_phase": if kind.starts_with("bubble") { "vapor" } else { "liquid" },
        "provenance": provenance(model),
    })
}

pub(crate) fn calculate(case: Case) -> Result<Value> {
    if case.kij.is_some() || case.interaction_source.is_some() {
        return Err(
            "NRTL-IG takes interactions from the embedded dataset; omit kij and interaction_source"
                .into(),
        );
    }
    if case.root.is_some() {
        return Err("root is not an input of NRTL-IG calculations; omit it".into());
    }
    let entries = case
        .components
        .iter()
        .zip(&case.mole_fractions)
        .map(|(id, x)| Ok((ComponentId::new(id)?, *x)))
        .collect::<std::result::Result<_, thermo_types::InputError>>()?;
    let composition = Composition::new(
        CompositionBasis::Mole,
        entries,
        NormalizationPolicy::RequireUnitSum,
    )?;
    let model = ActivityModel::from_databases(
        &ComponentDatabase::embedded()?,
        &LiquidDatabase::embedded()?,
        &composition,
        VaporModel::IdealGas,
    )?;
    let o = ActivityOptions::default();
    let calculation = case.calculation.as_deref().unwrap_or("homogeneous");
    let temperature = || -> Result<Temperature> {
        Ok(Temperature::new(
            case.temperature_k.ok_or("temperature_k is required")?,
        )?)
    };
    let pressure = || -> Result<Pressure> {
        Ok(Pressure::new(
            case.pressure_pa.ok_or("pressure_pa is required")?,
        )?)
    };
    let enthalpy_only_ph = |name: &str| -> Result<()> {
        if case.enthalpy_j_per_mol.is_some() {
            return Err(format!("enthalpy_j_per_mol is not an input of {name}").into());
        }
        Ok(())
    };
    match calculation {
        "bubble-pressure" | "dew-pressure" => {
            enthalpy_only_ph(calculation)?;
            if case.pressure_pa.is_some() {
                return Err(format!("pressure is an output of {calculation}; omit pressure_pa").into());
            }
            let t = temperature()?;
            let s = if calculation == "bubble-pressure" {
                activity_bubble_pressure(&model, t)?
            } else {
                activity_dew_pressure(&model, t, o)?
            };
            Ok(saturation_json(&case, calculation, &s, &model))
        }
        "bubble-temperature" | "dew-temperature" => {
            enthalpy_only_ph(calculation)?;
            if case.temperature_k.is_some() {
                return Err(format!("temperature is an output of {calculation}; omit temperature_k").into());
            }
            let p = pressure()?;
            let s = if calculation == "bubble-temperature" {
                activity_bubble_temperature(&model, p, o)?
            } else {
                activity_dew_temperature(&model, p, o)?
            };
            Ok(saturation_json(&case, calculation, &s, &model))
        }
        "tp-flash" => {
            enthalpy_only_ph(calculation)?;
            let r = activity_tp_flash(&model, temperature()?, pressure()?, o)?;
            Ok(json!({
                "case": case.id, "calculation": "tp-flash", "components": case.components,
                "feed_mole_fractions": model.fractions(),
                "temperature_K": r.temperature_k, "pressure_Pa_absolute": r.pressure_pa,
                "enthalpy_J_per_mol": r.enthalpy(), "phase_count": r.phases.len(),
                "molar_vapor_fraction": r.vapor_fraction(),
                "bubble_pressure_Pa": r.bubble_pressure_pa, "dew_pressure_Pa": r.dew_pressure_pa,
                "phases": phases_json(&r.phases),
                "convergence": {"successive_substitution_iterations": r.iterations,
                    "max_abs_ln_fugacity_residual": r.max_ln_fugacity_residual,
                    "max_abs_component_balance_residual": r.max_component_balance_residual},
                "liquid_stability_finite_trials": r.liquid_stability.as_ref().map(|s| s.stable_in_trials),
                "provenance": provenance(&model),
                "diagnostics": diagnostics(&r.diagnostics),
            }))
        }
        "ph-flash" => {
            if case.temperature_k.is_some() {
                return Err("temperature_k is an output of ph-flash; omit it".into());
            }
            let h = case.enthalpy_j_per_mol.ok_or("ph-flash requires enthalpy_j_per_mol")?;
            let r = activity_ph_flash(&model, pressure()?, h, o)?;
            Ok(json!({
                "case": case.id, "calculation": "ph-flash", "components": case.components,
                "feed_mole_fractions": model.fractions(),
                "temperature_K": r.temperature_k, "pressure_Pa_absolute": r.pressure_pa,
                "target_enthalpy_J_per_mol": r.target_enthalpy, "enthalpy_J_per_mol": r.enthalpy,
                "bubble_temperature_K": r.bubble_temperature_k, "dew_temperature_K": r.dew_temperature_k,
                "phase_count": r.phases.len(), "phases": phases_json(&r.phases),
                "convergence": {"temperature_iterations": r.iterations,
                    "scaled_energy_residual": r.scaled_energy_residual, "energy_scale_floor_J_per_mol": 1.0},
                "provenance": provenance(&model),
                "diagnostics": diagnostics(&r.diagnostics),
            }))
        }
        other => Err(format!(
            "NRTL-IG calculation must be tp-flash, ph-flash, bubble-pressure, dew-pressure, bubble-temperature or dew-temperature (got {other})"
        )
        .into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(json: &str) -> Case {
        serde_json::from_str(json).unwrap()
    }

    #[test]
    fn nrtl_examples_and_schema_failures() {
        let tp = include_str!("../../../examples/nrtl-tp-flash-methanol-water.json");
        let result = crate::calculate(case(tp)).unwrap();
        assert_eq!(result["phase_count"], 2);
        assert!(
            result["convergence"]["max_abs_ln_fugacity_residual"]
                .as_f64()
                .unwrap()
                <= 1e-10
        );
        let mut ph = case(tp);
        ph.calculation = Some("ph-flash".into());
        ph.temperature_k = None;
        ph.enthalpy_j_per_mol = result["enthalpy_J_per_mol"].as_f64();
        let out = crate::calculate(ph).unwrap();
        assert!((out["temperature_K"].as_f64().unwrap() - 355.0).abs() < 1e-7);
        let bubble = include_str!("../../../examples/nrtl-bubble-temperature-methanol-water.json");
        let out = crate::calculate(case(bubble)).unwrap();
        let t = out["temperature_K"].as_f64().unwrap();
        assert!(t > 337.0 && t < 374.0);
        for invalid in [
            tp.replace(
                "\"model\": \"NRTL-IG\"",
                "\"model\": \"NRTL-IG\", \"kij\": [[0,0],[0,0]]",
            ),
            tp.replace(
                "\"model\": \"NRTL-IG\"",
                "\"model\": \"NRTL-IG\", \"root\": \"minimum-gibbs\"",
            ),
            tp.replace("\"tp-flash\"", "\"homogeneous\""),
            tp.replace("\"methanol\"", "\"ethane\""),
            bubble.replace("\"pressure_pa\": 101325.0", "\"temperature_k\": 350.0"),
        ] {
            assert!(
                serde_json::from_str::<Case>(&invalid)
                    .map(crate::calculate)
                    .is_ok_and(|r| r.is_err()),
                "{invalid}"
            );
        }
    }
}
