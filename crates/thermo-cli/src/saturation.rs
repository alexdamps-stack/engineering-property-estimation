use serde_json::{Value, json};
use thermo_models::CubicModel;
use thermo_solvers::{
    CubicSaturationKind, EnvelopeStop, SaturationOptions, SaturationPoint, saturation_pressure,
    saturation_temperature, trace_phase_envelope,
};
use thermo_types::{Pressure, Temperature};

fn point(p: &SaturationPoint) -> Value {
    json!({"temperature_K":p.temperature_k,"pressure_Pa_absolute":p.pressure_pa,
        "liquid_mole_fractions":p.liquid_composition,"vapor_mole_fractions":p.vapor_composition,
        "liquid":super::state_json(&p.liquid),"vapor":super::state_json(&p.vapor),
        "max_abs_ln_fugacity_residual":p.max_ln_fugacity_residual,"iterations":p.iterations,
        "phase_stability":p.phase_stability.iter().map(super::stability_json).collect::<Vec<_>>(),
        "provenance":super::provenance(&p.liquid)})
}
fn stop_json(stop: &EnvelopeStop) -> Value {
    match stop {
        EnvelopeStop::RequestedRangeComplete => json!({"kind":"RequestedRangeComplete"}),
        EnvelopeStop::NearCritical {
            temperature_k,
            diagnostic,
        } => json!({"kind":"NearCritical","temperature_K":temperature_k,"diagnostic":diagnostic}),
        EnvelopeStop::UnverifiedBoundary {
            temperature_k,
            diagnostic,
        } => {
            json!({"kind":"UnverifiedBoundary","temperature_K":temperature_k,"diagnostic":diagnostic})
        }
    }
}
pub fn calculate(
    case: &super::Case,
    model: &CubicModel,
) -> Result<Value, Box<dyn std::error::Error>> {
    if case.root.is_some() || case.enthalpy_j_per_mol.is_some() {
        return Err("A13 saturation rejects root and enthalpy inputs".into());
    }
    let mut o = SaturationOptions::default();
    if let Some([lo, hi]) = case.saturation_temperature_bounds_k {
        o.minimum_temperature_k = lo;
        o.maximum_temperature_k = hi;
    }
    if let Some([lo, hi]) = case.saturation_pressure_bounds_pa {
        o.minimum_pressure_pa = lo;
        o.maximum_pressure_pa = hi;
    }
    let calculation = case.calculation.as_deref().unwrap();
    let kind = if calculation.starts_with("bubble") {
        CubicSaturationKind::Bubble
    } else {
        CubicSaturationKind::Dew
    };
    let mut result = if calculation.ends_with("envelope") {
        if case.temperature_k.is_some() || case.pressure_pa.is_some() {
            return Err("envelope uses saturation bounds; omit fixed T/p".into());
        }
        let r = trace_phase_envelope(
            model,
            kind,
            case.envelope_step_k.unwrap_or(5.0),
            case.envelope_maximum_points.unwrap_or(100),
            o,
        )?;
        json!({"points":r.points.iter().map(point).collect::<Vec<_>>(),"supported_temperature_range_K":r.supported_temperature_range_k,
            "supported_pressure_range_Pa_absolute":r.supported_pressure_range_pa,"stop":stop_json(&r.stop)})
    } else {
        if case.envelope_step_k.is_some() || case.envelope_maximum_points.is_some() {
            return Err("envelope controls only apply to envelope calculations".into());
        }
        let p = if calculation.ends_with("pressure") {
            if case.pressure_pa.is_some() {
                return Err("pressure is an output; omit pressure_pa".into());
            }
            saturation_pressure(
                model,
                Temperature::new(case.temperature_k.ok_or("temperature_k required")?)?,
                kind,
                o,
            )?
        } else {
            if case.temperature_k.is_some() {
                return Err("temperature is an output; omit temperature_k".into());
            }
            saturation_temperature(
                model,
                Pressure::new(case.pressure_pa.ok_or("pressure_pa required")?)?,
                kind,
                o,
            )?
        };
        point(&p)
    };
    result["case"] = json!(case.id);
    result["calculation"] = json!(calculation);
    result["components"] = json!(case.components);
    result["specified_mole_fractions"] = json!(case.mole_fractions);
    result["requested_temperature_range_K"] =
        json!([o.minimum_temperature_k, o.maximum_temperature_k]);
    result["requested_pressure_range_Pa_absolute"] =
        json!([o.minimum_pressure_pa, o.maximum_pressure_pa]);
    result["validity"] = json!(
        "A13 bounded dry cubic VLE; finite TPD trials, no global stability proof; volume ratio >=1.5, mixture composition separation >=0.03; no retrograde/critical extrapolation or experimental accuracy claim"
    );
    Ok(result)
}
#[cfg(test)]
mod tests {
    use serde_json::json;
    #[test]
    fn saturation_input_conflicts_are_rejected() {
        let input = include_str!("../../../examples/a13-bubble-pressure-binary.json");
        let mut value: serde_json::Value = serde_json::from_str(input).unwrap();
        value["pressure_pa"] = json!(1e6);
        let case = serde_json::from_value(value.clone()).unwrap();
        assert!(
            super::super::calculate(case)
                .unwrap_err()
                .to_string()
                .contains("output")
        );
        value.as_object_mut().unwrap().remove("pressure_pa");
        value["envelope_step_k"] = json!(5.);
        let case = serde_json::from_value(value.clone()).unwrap();
        assert!(super::super::calculate(case).is_err());
        value["calculation"] = json!("tp-flash");
        let case = serde_json::from_value(value).unwrap();
        assert!(super::super::calculate(case).is_err());
    }
}
