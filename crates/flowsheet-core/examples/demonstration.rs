//! Headless A8 process. Run with optional PR/PR-BM and positive flow multiplier.
mod support;
use serde_json::json;
use thermo_models::CubicVariant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let variant = match args.next().as_deref() {
        None | Some("PR-BM") => CubicVariant::PengRobinsonBostonMathias,
        Some("PR") => CubicVariant::PengRobinson,
        _ => return Err("expected PR or PR-BM".into()),
    };
    let scale = args
        .next()
        .map(|v| v.parse::<f64>())
        .transpose()?
        .unwrap_or(1.0);
    if !scale.is_finite() || scale <= 0.0 || args.next().is_some() {
        return Err("expected a finite positive flow multiplier, with no extra arguments".into());
    }
    let (fs, package, feeds) = support::demonstration(variant, scale);
    let r = fs.solve(&package, &feeds)?;
    let streams: Vec<_> = fs.streams().map(|s| {
        let result = &r.streams[&s.id];
        let state = result.state.as_ref().map(|a| json!({
            "temperature_K":a.temperature_k,"pressure_Pa_absolute":a.pressure_pa,"mole_fractions":a.mole_fractions,
            "enthalpy_J_per_mol":a.enthalpy(),"entropy_J_per_mol_K":a.entropy(),"iterations":a.iterations,
            "convergence":{"max_ln_fugacity_residual":a.convergence.max_ln_fugacity_residual,"max_component_balance_residual":a.convergence.max_component_balance_residual,"scaled_energy_residual":a.convergence.scaled_energy_residual},
            "diagnostics":a.diagnostics.iter().map(|d|json!({"code":d.code,"message":d.message})).collect::<Vec<_>>(),
            "phases":a.phases.iter().map(|p|json!({"label":format!("{:?}",p.phase),"fraction":p.phase_fraction,"mole_fractions":p.composition,
                "flow_mol_per_s":result.molar_flow_mol_s*p.phase_fraction,"density_kg_per_m3":p.state.mass_density,"Z":p.state.compressibility,
                "enthalpy_J_per_mol":p.state.enthalpy,"entropy_J_per_mol_K":p.state.entropy,
                "model":{"id":p.state.model.identifier,"version":p.state.model.implementation_version,"formulation":p.state.model.formulation},
                "data":{"id":p.state.data.identifier,"version":p.state.data.version,"source":p.state.data.source},
                "caloric_reference":{"id":p.state.caloric_reference.identifier,"definition":p.state.caloric_reference.definition},
                "validity":format!("{:?}",p.state.validity)})).collect::<Vec<_>>()
        }));
        json!({"id":s.name,"status":"solved","flow_mol_per_s":result.molar_flow_mol_s,"mass_flow_kg_per_s":result.mass_flow_kg_s(),"state":state})
    }).collect();
    let units: Vec<_> = r.order.iter().map(|id| {
        let u = &r.units[id];
        json!({"tag":fs.unit(*id).unwrap().name,"status":if u.warnings.is_empty(){"solved"}else{"warning"},"duty_W":u.duty_w,"warnings":u.warnings,
            "balance":u.balance.as_ref().map(|b|json!({"component_residuals_mol_per_s":b.component_residuals_mol_s,"max_scaled_component_residual":b.max_scaled_component_residual,"energy_residual_W":b.energy_residual_w,"scaled_energy_residual":b.scaled_energy_residual}))})
    }).collect();
    println!(
        "{}",
        serde_json::to_string_pretty(
            &json!({"schema":"a8-headless-demonstration-v1","components":support::COMPONENTS,"flow_multiplier":scale,"units":units,"streams":streams})
        )?
    );
    Ok(())
}
