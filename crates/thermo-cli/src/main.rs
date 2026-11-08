//! Reproducible cubic and native GERG properties, TP and PH calculations, and (A7)
//! NRTL gamma-phi bubble/dew points and TP/PH flashes.
mod activity;
mod saturation;

use serde::Deserialize;
use serde_json::{Value, json};
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, GergModel, HomogeneousState, RootSelection};
use thermo_solvers::{
    FlashOptions, GergDensityOptions, GergPhFlashOptions, PhFlashOptions, StabilityAssessment,
    gerg_density, gerg_ph_flash, gerg_tp_flash, ph_flash, tp_flash,
};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase, Pressure, Temperature,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    id: String,
    /// Homogeneous, TP/PH flash, cubic or NRTL saturation, or bounded cubic envelope.
    #[serde(default)]
    calculation: Option<String>,
    components: Vec<String>,
    mole_fractions: Vec<f64>,
    #[serde(default)]
    temperature_k: Option<f64>,
    /// Required only for PH flash, on the reported model caloric reference.
    #[serde(default)]
    enthalpy_j_per_mol: Option<f64>,
    /// Absolute pressure; an output (omit it) for bubble-pressure and dew-pressure.
    #[serde(default)]
    pressure_pa: Option<f64>,
    model: String,
    /// Required for homogeneous evaluation; rejected for a flash.
    #[serde(default)]
    root: Option<String>,
    /// Cubic models only: complete symmetric kij and its source. Rejected for NRTL-IG,
    /// whose directed NRTL interactions come only from the embedded dataset.
    #[serde(default)]
    kij: Option<Vec<Vec<f64>>>,
    #[serde(default)]
    interaction_source: Option<String>,
    #[serde(default)]
    saturation_temperature_bounds_k: Option<[f64; 2]>,
    #[serde(default)]
    saturation_pressure_bounds_pa: Option<[f64; 2]>,
    #[serde(default)]
    envelope_step_k: Option<f64>,
    #[serde(default)]
    envelope_maximum_points: Option<usize>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GergCase {
    id: String,
    calculation: String,
    model: String,
    components: Vec<String>,
    mole_fractions: Vec<f64>,
    #[serde(default)]
    temperature_k: Option<f64>,
    #[serde(default)]
    molar_density_mol_per_m3: Option<f64>,
    #[serde(default)]
    pressure_pa: Option<f64>,
    #[serde(default)]
    density_bracket_mol_per_m3: Option<[f64; 2]>,
    #[serde(default)]
    enthalpy_j_per_mol: Option<f64>,
    #[serde(default)]
    temperature_bracket_k: Option<[f64; 2]>,
}

fn calculate_gerg(case: GergCase) -> Result<Value, Box<dyn std::error::Error>> {
    if case.model != "GERG-2008" {
        return Err("GERG case requires model GERG-2008".into());
    }
    if case.components.len() != case.mole_fractions.len() {
        return Err("component and fraction lengths differ".into());
    }
    let entries = case
        .components
        .iter()
        .zip(&case.mole_fractions)
        .map(|(id, x)| Ok((ComponentId::new(id)?, *x)))
        .collect::<Result<_, thermo_types::InputError>>()?;
    let composition = Composition::new(
        CompositionBasis::Mole,
        entries,
        NormalizationPolicy::RequireUnitSum,
    )?;
    let model = GergModel::new(&composition)?;
    if matches!(case.calculation.as_str(), "tp-flash" | "ph-flash") {
        if case.molar_density_mol_per_m3.is_some() || case.density_bracket_mol_per_m3.is_some() {
            return Err(
                "GERG flash determines density and phase roots; omit density inputs".into(),
            );
        }
        let p = Pressure::new(case.pressure_pa.ok_or("GERG flash requires pressure_pa")?)?;
        let fractions: Vec<_> = composition.entries().iter().map(|(_, x)| *x).collect();
        if case.calculation == "tp-flash" {
            if case.enthalpy_j_per_mol.is_some() || case.temperature_bracket_k.is_some() {
                return Err("GERG TP flash requires temperature; omit energy inputs".into());
            }
            let t = Temperature::new(
                case.temperature_k
                    .ok_or("GERG TP flash requires temperature_k")?,
            )?;
            let r = gerg_tp_flash(&model, t, p, FlashOptions::default())?;
            return Ok(
                json!({"case":case.id,"calculation":"tp-flash","components":case.components,"feed_mole_fractions":fractions,
                "temperature_K":r.temperature_k,"pressure_Pa_absolute":r.pressure_pa,"phase_count":r.phases.len(),"molar_vapor_fraction":r.vapor_fraction(),
                "enthalpy_J_per_mol":r.enthalpy(),"phases":flash_phases_json(&r.phases),"provenance":provenance(&r.phases[0].state),
                "convergence":{"iterations":r.iterations,"max_abs_ln_fugacity_residual":r.max_ln_fugacity_residual,
                    "max_abs_component_balance_residual":r.max_component_balance_residual,"reduced_gibbs_change_per_mol_feed":r.reduced_gibbs_change},
                "feed_stability":stability_json(&r.feed_stability),"phase_stability":r.phase_stability.iter().map(stability_json).collect::<Vec<_>>(),
                "validity":"bounded GERG dry-gas fluid domain; finite-trial stability evidence; no solids or global stability proof",
                "diagnostics":r.diagnostics.iter().map(|d| json!({"code":d.code,"message":d.message})).collect::<Vec<_>>() }),
            );
        }
        if case.temperature_k.is_some() {
            return Err("temperature is an output of GERG PH flash; omit temperature_k".into());
        }
        let h = case
            .enthalpy_j_per_mol
            .ok_or("GERG PH flash requires enthalpy_j_per_mol on the native GERG reference")?;
        let mut options = GergPhFlashOptions::default();
        if let Some([lo, hi]) = case.temperature_bracket_k {
            options.minimum_temperature_k = lo;
            options.maximum_temperature_k = hi;
        }
        let r = gerg_ph_flash(&model, p, h, options)?;
        return Ok(
            json!({"case":case.id,"calculation":"ph-flash","components":case.components,"feed_mole_fractions":fractions,
            "temperature_K":r.temperature_k,"pressure_Pa_absolute":r.pressure_pa,"target_enthalpy_J_per_mol":h,"enthalpy_J_per_mol":r.enthalpy,
            "phase_count":r.phases.len(),"phases":flash_phases_json(&r.phases),"provenance":provenance(&r.phases[0].state),
            "convergence":{"temperature_iterations":r.iterations,"scaled_energy_residual":r.scaled_energy_residual,
                "max_abs_ln_fugacity_residual":r.tp.as_ref().map(|f| f.max_ln_fugacity_residual).or_else(||r.pure_saturation.as_ref().map(|s|s.fugacity_residual)),
                "pure_saturation_iterations":r.pure_saturation.as_ref().map(|s|s.iterations)},
            "feed_stability":r.tp.as_ref().map(|f|stability_json(&f.feed_stability)),
            "phase_stability":r.tp.as_ref().map(|f|f.phase_stability.iter().map(stability_json).collect::<Vec<_>>()),
            "validity":"bounded GERG dry-gas fluid domain; finite-trial stability evidence or pure saturation lever rule",
            "diagnostics":r.diagnostics.iter().map(|d| json!({"code":d.code,"message":d.message})).collect::<Vec<_>>() }),
        );
    }
    if case.enthalpy_j_per_mol.is_some() || case.temperature_bracket_k.is_some() {
        return Err("GERG homogeneous calculations reject PH energy inputs".into());
    }
    let t = Temperature::new(case.temperature_k.ok_or("temperature_k is required")?)?;
    let (rho, convergence) = match case.calculation.as_str() {
        "homogeneous-trho" => {
            if case.pressure_pa.is_some() || case.density_bracket_mol_per_m3.is_some() {
                return Err("homogeneous-trho requires density; pressure and density bracket must be omitted".into());
            }
            (
                case.molar_density_mol_per_m3
                    .ok_or("homogeneous-trho requires molar_density_mol_per_m3")?,
                Value::Null,
            )
        }
        "homogeneous-tp" => {
            if case.molar_density_mol_per_m3.is_some() {
                return Err("density is an output of homogeneous-tp; omit it".into());
            }
            let p = Pressure::new(
                case.pressure_pa
                    .ok_or("homogeneous-tp requires pressure_pa")?,
            )?;
            let [lo, hi] = case
                .density_bracket_mol_per_m3
                .ok_or("homogeneous-tp requires density_bracket_mol_per_m3")?;
            let r = gerg_density(&model, t, p, GergDensityOptions::new(lo, hi))?;
            (
                r.state.molar_density,
                json!({"iterations":r.iterations,"scaled_pressure_residual":r.scaled_pressure_residual,
                "target_pressure_Pa":r.target_pressure_pa,"pressure_tolerance":1e-10,"pressure_scale_floor_Pa":1.0,
                "initial_density_bracket_mol_per_m3":r.initial_density_bracket,"final_density_bracket_mol_per_m3":r.final_density_bracket}),
            )
        }
        _ => {
            return Err(
                "GERG-2008 supports homogeneous-trho, homogeneous-tp, tp-flash and ph-flash".into(),
            );
        }
    };
    let fugacity = model.fugacity(t, rho)?;
    let s = fugacity.state;
    Ok(json!({
        "case":case.id, "calculation":case.calculation, "components":case.components,
        "mole_fractions":composition.entries().iter().map(|(_,x)| *x).collect::<Vec<_>>(),
        "temperature_K":s.temperature_k, "molar_density_mol_per_m3":s.molar_density,
        "pressure_Pa_absolute":s.pressure_pa, "Z":s.compressibility, "mass_density_kg_per_m3":s.mass_density,
        "enthalpy_J_per_mol":s.enthalpy, "entropy_J_per_mol_K":s.entropy,
        "internal_energy_J_per_mol":s.internal_energy, "cp_J_per_mol_K":s.cp, "cv_J_per_mol_K":s.cv,
        "ln_fugacity_coefficients":fugacity.ln_fugacity_coefficients,
        "speed_of_sound_m_per_s":s.speed_of_sound, "dp_drho_Pa_per_mol_per_m3":s.dp_drho,"dp_dt_Pa_per_K":s.dp_dt,
        "model":{"id":"GERG-2008", "implementation":s.model_version},
        "data":{"id":"NIST AGA8 2.01", "revision":s.data_revision},
        "caloric_reference":{"id":"gerg-native-ideal", "definition":s.caloric_reference},
        "validity":"unassessed physical phase/domain; local mechanical and thermal checks passed",
        "convergence":convergence,
        "diagnostics":["Homogeneous state; an explicit density bracket selects a TP root. Root completeness, phase identity and equilibrium stability are unassessed."]
    }))
}

fn flash_phases_json(phases: &[thermo_solvers::FlashPhase]) -> Vec<Value> {
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

fn calculate_input(input: &str) -> Result<Value, Box<dyn std::error::Error>> {
    let header: Value = serde_json::from_str(input)?;
    if header["model"] == "GERG-2008" {
        calculate_gerg(serde_json::from_str(input)?)
    } else {
        calculate(serde_json::from_str(input)?)
    }
}

fn state_json(s: &HomogeneousState) -> Value {
    json!({
        "Z":s.compressibility,"molar_volume_m3_per_mol":s.molar_volume,"mass_density_kg_per_m3":s.mass_density,
        "ln_fugacity_coefficients":s.ln_fugacity_coefficients,
        "residual_enthalpy_J_per_mol":s.residual_enthalpy,"residual_entropy_J_per_mol_K":s.residual_entropy,
        "enthalpy_J_per_mol":s.enthalpy,"entropy_J_per_mol_K":s.entropy,"internal_energy_J_per_mol":s.internal_energy,
        "cp_J_per_mol_K":s.cp,"cv_J_per_mol_K":s.cv,"speed_of_sound_m_per_s":s.speed_of_sound,
    })
}

fn provenance(s: &HomogeneousState) -> Value {
    json!({
        "model":{"id":s.model.identifier,"implementation":s.model.implementation_version,"formulation":s.model.formulation},
        "data":{"id":s.data.identifier,"version":s.data.version,"source":s.data.source},
        "caloric_reference":{"id":s.caloric_reference.identifier,"definition":s.caloric_reference.definition},
    })
}

fn stability_json(s: &StabilityAssessment) -> Value {
    json!({
        "stable_in_finite_trials": s.stable_in_trials,
        "tm_threshold": s.options.tm_threshold,
        "trials": s.trials.iter().map(|t| json!({
            "initialization": format!("{:?}", t.initialization),
            "outcome": format!("{:?}", t.outcome),
            "modified_tangent_plane_distance": t.tangent_plane_distance,
            "composition": t.composition, "Z": t.compressibility, "iterations": t.iterations,
        })).collect::<Vec<_>>(),
    })
}

fn phase_label(phase: Phase) -> &'static str {
    match phase {
        Phase::Vapor => "vapor",
        Phase::Liquid => "liquid",
        Phase::Supercritical => "supercritical",
    }
}

fn calculate(case: Case) -> Result<Value, Box<dyn std::error::Error>> {
    if case.components.len() != case.mole_fractions.len() {
        return Err("component and fraction lengths differ".into());
    }
    let saturation_controls = case.saturation_temperature_bounds_k.is_some()
        || case.saturation_pressure_bounds_pa.is_some()
        || case.envelope_step_k.is_some()
        || case.envelope_maximum_points.is_some();
    let is_cubic_saturation = case.model != "NRTL-IG"
        && case
            .calculation
            .as_deref()
            .is_some_and(|c| c.starts_with("bubble-") || c.starts_with("dew-"));
    if saturation_controls && !is_cubic_saturation {
        return Err(
            "A13 saturation bounds and envelope controls require a cubic saturation calculation"
                .into(),
        );
    }
    let variant = match case.model.as_str() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        "NRTL-IG" => return activity::calculate(case),
        _ => return Err("model must be PR, PR-BM or NRTL-IG".into()),
    };
    let calculation = case.calculation.as_deref().unwrap_or("homogeneous");
    let flash = match calculation {
        "homogeneous" => false,
        "tp-flash" | "ph-flash" | "bubble-pressure" | "dew-pressure" | "bubble-temperature"
        | "dew-temperature" | "bubble-envelope" | "dew-envelope" => true,
        _ => return Err("unknown calculation; use homogeneous, tp-flash, ph-flash, bubble/dew-pressure, bubble/dew-temperature or bubble/dew-envelope".into()),
    };
    let root = match (flash, case.root.as_deref()) {
        (true, None) => None,
        (true, Some(_)) => return Err("root is not an input of a flash; omit it".into()),
        (false, Some("smallest-volume")) => Some(RootSelection::SmallestVolume),
        (false, Some("largest-volume")) => Some(RootSelection::LargestVolume),
        (false, Some("minimum-gibbs")) => Some(RootSelection::MinimumGibbs),
        (false, _) => {
            return Err("root must be smallest-volume, largest-volume or minimum-gibbs".into());
        }
    };
    let entries = case
        .components
        .iter()
        .zip(&case.mole_fractions)
        .map(|(id, x)| Ok((ComponentId::new(id)?, *x)))
        .collect::<Result<_, thermo_types::InputError>>()?;
    let composition = Composition::new(
        CompositionBasis::Mole,
        entries,
        NormalizationPolicy::RequireUnitSum,
    )?;
    let db = ComponentDatabase::embedded()?;
    let model = CubicModel::new(
        &db,
        &composition,
        variant,
        case.kij.clone().ok_or("cubic models require kij")?,
        case.interaction_source
            .clone()
            .ok_or("cubic models require interaction_source")?,
    )?;
    if calculation.starts_with("bubble-") || calculation.starts_with("dew-") {
        return saturation::calculate(&case, &model);
    }
    let p = Pressure::new(case.pressure_pa.ok_or("pressure_pa is required")?)?;
    let fractions: Vec<f64> = composition.entries().iter().map(|(_, x)| *x).collect();
    if calculation == "ph-flash" {
        if case.temperature_k.is_some() {
            return Err("temperature_k is an output of ph-flash; omit it".into());
        }
        let h = case
            .enthalpy_j_per_mol
            .ok_or("ph-flash requires enthalpy_j_per_mol")?;
        let r = ph_flash(&model, p, h, PhFlashOptions::default())?;
        return Ok(json!({
            "case":case.id, "calculation":"ph-flash", "components":case.components,
            "feed_mole_fractions":fractions, "temperature_K":r.temperature_k,
            "pressure_Pa_absolute":r.pressure_pa, "target_enthalpy_J_per_mol":r.target_enthalpy,
            "enthalpy_J_per_mol":r.enthalpy, "phase_count":r.phases.len(),
            "phases": r.phases.iter().map(|ph| {
                let mut v = state_json(&ph.state);
                v["label"] = json!(phase_label(ph.phase));
                v["molar_phase_fraction"] = json!(ph.phase_fraction);
                v["mole_fractions"] = json!(ph.composition);
                v
            }).collect::<Vec<_>>(),
            "convergence": {
                "temperature_iterations":r.iterations,
                "scaled_energy_residual":r.scaled_energy_residual,
                "energy_scale_floor_J_per_mol":1.0,
                "max_abs_ln_fugacity_residual":r.tp.as_ref().map(|f| f.max_ln_fugacity_residual).or_else(|| r.pure_saturation.as_ref().map(|s| s.fugacity_residual)),
                "max_abs_component_balance_residual":r.tp.as_ref().map(|f| f.max_component_balance_residual).unwrap_or(0.0),
                "pure_saturation_iterations":r.pure_saturation.as_ref().map(|s| s.iterations),
            },
            "feed_stability":r.tp.as_ref().map(|f| stability_json(&f.feed_stability)),
            "phase_stability":r.tp.as_ref().map(|f| f.phase_stability.iter().map(stability_json).collect::<Vec<_>>()),
            "provenance":provenance(&r.phases[0].state),
            "validity":"unassessed physical phase/domain; evaluated inside computational envelope",
            "diagnostics":r.diagnostics.iter().map(|d| json!({"code":d.code,"severity":format!("{:?}",d.severity),"message":d.message})).collect::<Vec<_>>(),
        }));
    }
    if case.enthalpy_j_per_mol.is_some() {
        return Err("enthalpy_j_per_mol is only an input of ph-flash".into());
    }
    let t = Temperature::new(case.temperature_k.ok_or("temperature_k is required")?)?;
    let Some(root) = root else {
        let r = tp_flash(&model, t, p, FlashOptions::default())?;
        return Ok(json!({
            "case":case.id, "calculation":"tp-flash", "components":case.components,
            "feed_mole_fractions":fractions,
            "temperature_K":r.temperature_k, "pressure_Pa_absolute":r.pressure_pa,
            "enthalpy_J_per_mol":r.enthalpy(),
            "phase_count":r.phases.len(), "molar_vapor_fraction":r.vapor_fraction(),
            "phases": r.phases.iter().map(|ph| {
                let mut v = state_json(&ph.state);
                v["label"] = json!(phase_label(ph.phase));
                v["molar_phase_fraction"] = json!(ph.phase_fraction);
                v["mole_fractions"] = json!(ph.composition);
                v
            }).collect::<Vec<_>>(),
            "convergence":{
                "split_iterations":r.iterations,
                // Retained for existing consumers; includes Newton iterations.
                "successive_substitution_iterations":r.iterations,
                "max_abs_ln_fugacity_residual":r.max_ln_fugacity_residual,
                "max_abs_component_balance_residual":r.max_component_balance_residual,
                "reduced_gibbs_change_per_mol_feed":r.reduced_gibbs_change},
            "feed_stability":stability_json(&r.feed_stability),
            "phase_stability":r.phase_stability.iter().map(stability_json).collect::<Vec<_>>(),
            "provenance":provenance(&r.phases[0].state),
            "validity":"unassessed physical phase/domain; evaluated inside computational envelope",
            "diagnostics":r.diagnostics.iter().map(|d| json!({
                "code":d.code,"severity":format!("{:?}",d.severity),"message":d.message})).collect::<Vec<_>>(),
        }));
    };
    let roots = model.compressibility_roots(t, p)?;
    let s = model.evaluate(t, p, root)?;
    let mut out = json!({
        "case":case.id, "calculation":"homogeneous", "components":case.components,
        "mole_fractions":fractions,
        "temperature_K":s.temperature_k, "pressure_Pa_absolute":s.pressure_pa,
        "root_request":case.root, "all_admissible_Z_roots":roots,
        "validity":"unassessed physical phase/domain; evaluated inside computational envelope",
        "diagnostics":["Fixed-composition homogeneous EOS root. Global phase stability and equilibrium phase count are not assessed; use calculation tp-flash."]
    });
    for source in [state_json(&s), provenance(&s)] {
        for (key, value) in source.as_object().into_iter().flatten() {
            out[key] = value.clone();
        }
    }
    Ok(out)
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err("usage: thermo-cli [case.json]".into());
    }
    if args.first().is_some_and(|a| a == "--help" || a == "-h") {
        println!(
            "Usage: thermo-cli [case.json]\nDefault: embedded dry natural-gas PR-BM homogeneous example.\nSchema/examples: examples/homogeneous-*.json, examples/tp-flash-*.json, examples/ph-flash-*.json, examples/gerg-*.json.\nCubic calculation: homogeneous (requires root), tp-flash (requires temperature_k) or ph-flash (requires enthalpy_j_per_mol; temperature is an output). Flashes support at most two phases.\nGERG-2008: homogeneous-trho requires temperature_k and molar_density_mol_per_m3; homogeneous-tp requires temperature_k, pressure_pa and density_bracket_mol_per_m3. GERG tp-flash and ph-flash support the dry five-component fluid domain, 230..350 K and 0.1..10 MPa. PH uses native GERG enthalpy, optionally temperature_bracket_k.
model NRTL-IG (NRTL liquid, ideal-gas vapor; embedded methanol/water data, no kij): tp-flash, ph-flash, bubble-pressure, dew-pressure (temperature_k), bubble-temperature, dew-temperature. See examples/nrtl-*.json."
        );
        return Ok(());
    }
    let input = if let Some(path) = args.first() {
        std::fs::read_to_string(path)?
    } else {
        include_str!("../../../examples/homogeneous-natural-gas.json").into()
    };
    println!(
        "{}",
        serde_json::to_string_pretty(&calculate_input(&input)?)?
    );
    Ok(())
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("thermo-cli: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gerg_tp_ph_flash_outputs_and_strict_schema() {
        let input = include_str!("../../../examples/gerg-tp-flash-binary.json");
        let tp = calculate_input(input).unwrap();
        assert_eq!(tp["phase_count"], 2);
        assert!((tp["molar_vapor_fraction"].as_f64().unwrap() - 0.40175558515).abs() < 2e-8);
        assert_eq!(tp["provenance"]["model"]["id"], "GERG-2008");
        assert_eq!(
            tp["provenance"]["caloric_reference"]["id"],
            "gerg-native-ideal"
        );
        let ph_input = include_str!("../../../examples/gerg-ph-flash-natural-gas.json");
        let ph = calculate_input(ph_input).unwrap();
        assert_eq!(ph["phase_count"], 1);
        assert!((300.0..350.0).contains(&ph["temperature_K"].as_f64().unwrap()));
        assert!(
            ph["convergence"]["scaled_energy_residual"]
                .as_f64()
                .unwrap()
                .abs()
                <= 1e-10
        );
        for invalid in [
            input.replace(
                "\"model\":",
                "\"density_bracket_mol_per_m3\": [1,20000], \"model\":",
            ),
            input.replace("temperature_k", "wrong_temperature"),
            input.replace("250.0", "220.0"),
            input.replace("propane", "water"),
            ph_input.replace("\"model\":", "\"temperature_k\": 320, \"model\":"),
            ph_input.replace("enthalpy_j_per_mol", "wrong_enthalpy"),
            ph_input.replace("temperature_bracket_k", "wrong_temperature_bracket"),
        ] {
            assert!(calculate_input(&invalid).is_err(), "{invalid}");
        }
    }
    #[test]
    fn gerg_tp_density_fugacity_and_input_contract() {
        let input = include_str!("../../../examples/gerg-tp-natural-gas.json");
        let result = calculate_input(input).unwrap();
        assert!((result["pressure_Pa_absolute"].as_f64().unwrap() - 1e7).abs() < 1e-3);
        assert_eq!(
            result["ln_fugacity_coefficients"].as_array().unwrap().len(),
            5
        );
        assert!(
            result["convergence"]["scaled_pressure_residual"]
                .as_f64()
                .unwrap()
                .abs()
                <= 1e-10
        );
        let mut trho: Value = serde_json::from_str(input).unwrap();
        trho["calculation"] = json!("homogeneous-trho");
        trho["molar_density_mol_per_m3"] = result["molar_density_mol_per_m3"].clone();
        trho.as_object_mut().unwrap().remove("pressure_pa");
        trho.as_object_mut()
            .unwrap()
            .remove("density_bracket_mol_per_m3");
        let recovered = calculate_input(&trho.to_string()).unwrap();
        assert_eq!(
            recovered["ln_fugacity_coefficients"],
            result["ln_fugacity_coefficients"]
        );
        for invalid in [
            input.replace("density_bracket_mol_per_m3", "density_bracket"),
            input.replace("[1000.0, 7000.0]", "[7000.0, 1000.0]"),
            input.replace(
                "\"model\":",
                "\"molar_density_mol_per_m3\": 1000, \"model\":",
            ),
            input.replace("10000000.0", "36000000.0"),
            input.replace("homogeneous-tp", "tp-flash"),
        ] {
            assert!(calculate_input(&invalid).is_err(), "{invalid}");
        }
    }
    #[test]
    fn gerg_schema_and_capability_boundaries() {
        let input = include_str!("../../../examples/gerg-natural-gas.json");
        let result = calculate_input(input).unwrap();
        assert!(result["Z"].as_f64().unwrap() > 0.9);
        assert_eq!(result["model"]["id"], "GERG-2008");
        assert_eq!(result["caloric_reference"]["id"], "gerg-native-ideal");
        for invalid in [
            input.replace("homogeneous-trho", "tp-flash"),
            input.replace("homogeneous-trho", "ph-flash"),
            input.replace("\"model\":", "\"pressure_pa\": 100000, \"model\":"),
            input.replace("molar_density_mol_per_m3", "density_kg_per_m3"),
            input.replace("methane", "methanol"),
            input.replace("300.0", "500.0"),
        ] {
            assert!(calculate_input(&invalid).is_err(), "{invalid}");
        }
    }
    #[test]
    fn checked_cli_schema_and_failure_outcomes() {
        let input = include_str!("../../../examples/homogeneous-binary.json");
        let result = calculate(serde_json::from_str(input).unwrap()).unwrap();
        assert!(result["Z"].as_f64().unwrap() > 0.0);
        let mut case: Case = serde_json::from_str(input).unwrap();
        case.mole_fractions.pop();
        assert!(calculate(case).is_err());
        let mut case: Case = serde_json::from_str(input).unwrap();
        case.model = "typo".into();
        assert!(calculate(case).is_err());
        let mut case: Case = serde_json::from_str(input).unwrap();
        case.pressure_pa = Some(-1.0);
        assert!(calculate(case).is_err());
        assert!(
            serde_json::from_str::<Case>(&input.replace("\"pressure_pa\"", "\"gauge_pressure\""))
                .is_err()
        );
        let mut case: Case = serde_json::from_str(input).unwrap();
        case.root = None;
        assert!(
            calculate(case).is_err(),
            "homogeneous evaluation requires a root"
        );
    }

    #[test]
    fn tp_flash_cases_report_split_and_evidence() {
        let input = include_str!("../../../examples/tp-flash-binary.json");
        let result = calculate(serde_json::from_str(input).unwrap()).unwrap();
        assert_eq!(result["phase_count"], 2);
        assert!(result["convergence"]["split_iterations"].as_u64().unwrap() > 0);
        assert_eq!(
            result["convergence"]["split_iterations"],
            result["convergence"]["successive_substitution_iterations"]
        );
        let beta = result["molar_vapor_fraction"].as_f64().unwrap();
        assert!(beta > 0.0 && beta < 1.0);
        assert!(
            result["convergence"]["max_abs_ln_fugacity_residual"]
                .as_f64()
                .unwrap()
                <= 1e-10
        );
        assert_eq!(result["feed_stability"]["stable_in_finite_trials"], false);
        assert_eq!(result["phases"][0]["label"], "vapor");
        let mut case: Case = serde_json::from_str(input).unwrap();
        case.root = Some("minimum-gibbs".into());
        assert!(calculate(case).is_err(), "flash rejects a root request");
        let mut case: Case = serde_json::from_str(input).unwrap();
        case.calculation = Some("unknown-flash".into());
        assert!(calculate(case).is_err());
        let gas = include_str!("../../../examples/tp-flash-natural-gas.json");
        let result = calculate(serde_json::from_str(gas).unwrap()).unwrap();
        assert_eq!(result["phase_count"], 2);
    }

    #[test]
    fn ph_cli_round_trip_and_schema_errors() {
        let input = include_str!("../../../examples/tp-flash-binary.json");
        let tp = calculate(serde_json::from_str(input).unwrap()).unwrap();
        let mut ph: Case = serde_json::from_str(input).unwrap();
        ph.calculation = Some("ph-flash".into());
        ph.temperature_k = None;
        ph.enthalpy_j_per_mol = tp["enthalpy_J_per_mol"].as_f64();
        let output = calculate(ph).unwrap();
        assert!((output["temperature_K"].as_f64().unwrap() - 250.0).abs() < 1e-7);
        assert_eq!(output["phase_count"], 2);
        assert!(
            output["convergence"]["scaled_energy_residual"]
                .as_f64()
                .unwrap()
                .abs()
                <= 1e-10
        );
        let example = include_str!("../../../examples/ph-flash-natural-gas.json");
        assert_eq!(
            calculate(serde_json::from_str(example).unwrap()).unwrap()["phase_count"],
            2
        );
        for invalid in [
            example.replace("enthalpy_j_per_mol", "wrong_enthalpy"),
            example.replace(
                "\"calculation\": \"ph-flash\",",
                "\"calculation\": \"ph-flash\", \"temperature_k\": 250.0,",
            ),
            example.replace(
                "\"calculation\": \"ph-flash\",",
                "\"calculation\": \"ph-flash\", \"root\": \"minimum-gibbs\",",
            ),
        ] {
            assert!(
                serde_json::from_str::<Case>(&invalid)
                    .map(calculate)
                    .map_or(true, |r| r.is_err())
            );
        }
    }
}
