//! Reproducible homogeneous cubic calculations. Does not perform a flash.
use serde::Deserialize;
use serde_json::{Value, json};
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, RootSelection};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure, Temperature,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    id: String,
    components: Vec<String>,
    mole_fractions: Vec<f64>,
    temperature_k: f64,
    pressure_pa: f64,
    model: String,
    root: String,
    kij: Vec<Vec<f64>>,
    interaction_source: String,
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
    let root = match case.root.as_str() {
        "smallest-volume" => RootSelection::SmallestVolume,
        "largest-volume" => RootSelection::LargestVolume,
        "minimum-gibbs" => RootSelection::MinimumGibbs,
        _ => return Err("root must be smallest-volume, largest-volume or minimum-gibbs".into()),
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
    let roots = model.compressibility_roots(t, p)?;
    let s = model.evaluate(t, p, root)?;
    Ok(json!({
        "case":case.id, "components":case.components,
        "mole_fractions":composition.entries().iter().map(|(_,x)| *x).collect::<Vec<_>>(),
        "temperature_K":s.temperature_k, "pressure_Pa_absolute":s.pressure_pa,
        "root_request":case.root, "all_admissible_Z_roots":roots,
        "Z":s.compressibility,"molar_volume_m3_per_mol":s.molar_volume,"mass_density_kg_per_m3":s.mass_density,
        "ln_fugacity_coefficients":s.ln_fugacity_coefficients,
        "residual_enthalpy_J_per_mol":s.residual_enthalpy,"residual_entropy_J_per_mol_K":s.residual_entropy,
        "enthalpy_J_per_mol":s.enthalpy,"entropy_J_per_mol_K":s.entropy,"internal_energy_J_per_mol":s.internal_energy,
        "cp_J_per_mol_K":s.cp,"cv_J_per_mol_K":s.cv,"speed_of_sound_m_per_s":s.speed_of_sound,
        "model":{"id":s.model.identifier,"implementation":s.model.implementation_version,"formulation":s.model.formulation},
        "data":{"id":s.data.identifier,"version":s.data.version,"source":s.data.source},
        "caloric_reference":{"id":s.caloric_reference.identifier,"definition":s.caloric_reference.definition},
        "validity":"unassessed physical phase/domain; evaluated inside computational envelope",
        "diagnostics":["Fixed-composition homogeneous EOS root. Global phase stability and equilibrium phase count are not assessed."]
    }))
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err("usage: thermo-cli [case.json]".into());
    }
    if args.first().is_some_and(|a| a == "--help" || a == "-h") {
        println!(
            "Usage: thermo-cli [case.json]\nDefault: embedded dry natural-gas PR-BM homogeneous example.\nSchema/examples: examples/homogeneous-*.json. No TP/PH flash or phase stability."
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
    }
}
