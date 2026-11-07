//! Reproducible homogeneous cubic calculations and A4 two-phase TP flashes.
use serde::Deserialize;
use serde_json::{Value, json};
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, HomogeneousState, RootSelection};
use thermo_solvers::{FlashOptions, StabilityAssessment, tp_flash};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase, Pressure, Temperature,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    id: String,
    /// "homogeneous" (default) or "tp-flash".
    #[serde(default)]
    calculation: Option<String>,
    components: Vec<String>,
    mole_fractions: Vec<f64>,
    temperature_k: f64,
    pressure_pa: f64,
    model: String,
    /// Required for homogeneous evaluation; rejected for a flash.
    #[serde(default)]
    root: Option<String>,
    kij: Vec<Vec<f64>>,
    interaction_source: String,
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
    let variant = match case.model.as_str() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        _ => return Err("model must be PR or PR-BM".into()),
    };
    let flash = match case.calculation.as_deref() {
        None | Some("homogeneous") => false,
        Some("tp-flash") => true,
        Some(_) => return Err("calculation must be homogeneous or tp-flash".into()),
    };
    let root = match (flash, case.root.as_deref()) {
        (true, None) => None,
        (true, Some(_)) => return Err("root is not an input of tp-flash; omit it".into()),
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
        case.kij,
        case.interaction_source,
    )?;
    let t = Temperature::new(case.temperature_k)?;
    let p = Pressure::new(case.pressure_pa)?;
    let fractions: Vec<f64> = composition.entries().iter().map(|(_, x)| *x).collect();
    let Some(root) = root else {
        let r = tp_flash(&model, t, p, FlashOptions::default())?;
        return Ok(json!({
            "case":case.id, "calculation":"tp-flash", "components":case.components,
            "feed_mole_fractions":fractions,
            "temperature_K":r.temperature_k, "pressure_Pa_absolute":r.pressure_pa,
            "phase_count":r.phases.len(), "molar_vapor_fraction":r.vapor_fraction(),
            "phases": r.phases.iter().map(|ph| {
                let mut v = state_json(&ph.state);
                v["label"] = json!(phase_label(ph.phase));
                v["molar_phase_fraction"] = json!(ph.phase_fraction);
                v["mole_fractions"] = json!(ph.composition);
                v
            }).collect::<Vec<_>>(),
            "convergence":{
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
            "Usage: thermo-cli [case.json]\nDefault: embedded dry natural-gas PR-BM homogeneous example.\nSchema/examples: examples/homogeneous-*.json and examples/tp-flash-*.json.\ncalculation: homogeneous (default; requires root) or tp-flash (at most two phases). No PH flash."
        );
        return Ok(());
    }
    let input = if let Some(path) = args.first() {
        std::fs::read_to_string(path)?
    } else {
        include_str!("../../../examples/homogeneous-natural-gas.json").into()
    };
    let case: Case = serde_json::from_str(&input)?;
    println!("{}", serde_json::to_string_pretty(&calculate(case)?)?);
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
        case.pressure_pa = -1.0;
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
        case.calculation = Some("ph-flash".into());
        assert!(calculate(case).is_err());
        let gas = include_str!("../../../examples/tp-flash-natural-gas.json");
        let result = calculate(serde_json::from_str(gas).unwrap()).unwrap();
        assert_eq!(result["phase_count"], 2);
    }
}
