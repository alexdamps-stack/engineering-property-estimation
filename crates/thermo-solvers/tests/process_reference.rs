//! Fixed A2 demonstration states independently solved with teqp/Cantera/SciPy.
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::{FlashOptions, PhFlashOptions, ph_flash, tp_flash};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure, Temperature,
};

#[test]
fn all_demonstration_streams_match_independent_references() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/process-v1.json")).unwrap();
    assert_eq!(fixture["schema"], "a2-process-reference-v1");
    assert_eq!(fixture["packages"]["teqp"], "0.23.2");
    let close = |key: &str, actual: f64, expected: f64| {
        let b = &fixture["tolerances"][key];
        assert!(actual.is_finite() && expected.is_finite());
        assert!(
            (actual - expected).abs()
                <= b[0].as_f64().unwrap() + b[1].as_f64().unwrap() * expected.abs(),
            "{key}: {actual} != {expected}"
        );
    };
    assert_eq!(fixture["cases"].as_array().unwrap().len(), 2);
    for case in fixture["cases"].as_array().unwrap() {
        let streams = case["streams"].as_array().unwrap();
        assert_eq!(streams.len(), 7);
        let mut enthalpies = Vec::new();
        let mut separator_flows = Vec::new();
        for stream in streams {
            let z: Vec<f64> = stream["z"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            let composition = Composition::new(
                CompositionBasis::Mole,
                stream["ids"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .zip(&z)
                    .map(|(id, x)| (ComponentId::new(id.as_str().unwrap()).unwrap(), *x))
                    .collect(),
                NormalizationPolicy::RequireUnitSum,
            )
            .unwrap();
            let variant = if case["variant"] == "PR" {
                CubicVariant::PengRobinson
            } else {
                CubicVariant::PengRobinsonBostonMathias
            };
            let model = CubicModel::new(
                &ComponentDatabase::embedded().unwrap(),
                &composition,
                variant,
                vec![vec![0.0; z.len()]; z.len()],
                "explicit A2 zero kij approximation",
            )
            .unwrap();
            let pressure = Pressure::new(stream["p"].as_f64().unwrap()).unwrap();
            let t = stream["T"].as_f64().unwrap();
            let expected_h = stream["target_h"].as_f64().unwrap();
            let result = tp_flash(
                &model,
                Temperature::new(t).unwrap(),
                pressure,
                FlashOptions::default(),
            )
            .unwrap_or_else(|e| panic!("{} {}: {e}", case["variant"], stream["id"]));
            assert_eq!(
                result.phases.len(),
                stream["phases"].as_u64().unwrap() as usize
            );
            close("enthalpy_J_per_mol", result.enthalpy(), expected_h);
            if result.phases.len() == 2 {
                close(
                    "phase_fraction",
                    result.phases[0].phase_fraction,
                    stream["vapor_fraction"].as_f64().unwrap(),
                );
                for (phase, name) in result.phases.iter().zip(["vapor", "liquid"]) {
                    for (i, value) in phase.composition.iter().enumerate() {
                        close("composition", *value, stream[name][i].as_f64().unwrap());
                    }
                }
            }
            if stream["id"] == "mix" || stream["id"] == "expand" {
                let recovered = ph_flash(
                    &model,
                    pressure,
                    expected_h,
                    PhFlashOptions {
                        minimum_temperature_k: 200.0,
                        maximum_temperature_k: 400.0,
                        ..Default::default()
                    },
                )
                .unwrap();
                close("temperature_K", recovered.temperature_k, t);
            }
            enthalpies.push(result.enthalpy());
            if stream["id"].as_str().unwrap().starts_with("separator_") {
                separator_flows.push(stream["flow_mol_per_s"].as_f64().unwrap());
            }
        }
        close(
            "enthalpy_J_per_mol",
            (100.0 * enthalpies[0] + 20.0 * enthalpies[1]) / 120.0,
            enthalpies[2],
        );
        close("enthalpy_J_per_mol", enthalpies[3], enthalpies[4]);
        close(
            "enthalpy_J_per_mol",
            (separator_flows[0] * enthalpies[5] + separator_flows[1] * enthalpies[6]) / 120.0,
            enthalpies[4],
        );
        for ((vapor, liquid), feed) in streams[5]["z"]
            .as_array()
            .unwrap()
            .iter()
            .zip(streams[6]["z"].as_array().unwrap())
            .zip(streams[4]["z"].as_array().unwrap())
        {
            let recovered = (separator_flows[0] * vapor.as_f64().unwrap()
                + separator_flows[1] * liquid.as_f64().unwrap())
                / 120.0;
            close("composition", recovered, feed.as_f64().unwrap());
        }
        close(
            "enthalpy_J_per_mol",
            case["cooler_duty_W"].as_f64().unwrap() / 120.0,
            enthalpies[3] - enthalpies[2],
        );
    }
}
