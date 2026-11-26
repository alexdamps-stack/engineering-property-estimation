//! Transport properties (#225) against the committed independent fixture
//! `validation/fixtures/transport-v1.json` (see
//! `validation/transport/generate-transport-references.py`).
use serde_json::Value;
use thermo_data::{ComponentDatabase, DATA_VERSION, TRANSPORT_DATA_VERSION};
use thermo_models::{CubicModel, CubicVariant, RootSelection, TransportModel, TransportRegime};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure,
    Temperature, ValidityStatus,
};

const FIXTURE: &str = include_str!("../../../validation/fixtures/transport-v1.json");

fn fixture() -> Value {
    serde_json::from_str(FIXTURE).unwrap()
}

fn budget(f: &Value, key: &str) -> (f64, f64) {
    let b = &f["tolerances"][key];
    (b[0].as_f64().unwrap(), b[1].as_f64().unwrap())
}

fn strings(v: &Value) -> Vec<&str> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap())
        .collect()
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_f64().unwrap())
        .collect()
}

const PROPERTIES: [&str; 3] = ["viscosity", "thermal_conductivity", "surface_tension"];

fn property(t: &thermo_models::TransportProperties, name: &str) -> Result<f64, CalculationError> {
    match name {
        "viscosity" => t.viscosity.clone(),
        "thermal_conductivity" => t.thermal_conductivity.clone(),
        _ => t.surface_tension.clone(),
    }
}

#[test]
fn fixture_inputs_are_the_embedded_data_versions() {
    let f = fixture();
    let hashes = f["input_sha256"].as_object().unwrap();
    assert!(hashes.contains_key("data/source-extracts/transport.provisional.json"));
    let extract = include_str!("../../../data/source-extracts/transport.provisional.json");
    assert!(extract.contains(&format!("\"version\": \"{TRANSPORT_DATA_VERSION}\"")));
    let model = TransportModel::new(&["methane"]).unwrap();
    assert_eq!(
        model.data().version,
        format!("{TRANSPORT_DATA_VERSION}+components-{DATA_VERSION}")
    );
}

/// Matched formulation: every value, and every expected unavailability.
#[test]
fn implementation_matches_independent_python_evaluation() {
    let f = fixture();
    let mut checked = 0;
    for case in f["implementation"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let ids = strings(&case["components"]);
        let x = floats(&case["mole_fractions"]);
        let model = TransportModel::new(&ids).unwrap();
        let t = Temperature::new(case["temperature_K"].as_f64().unwrap()).unwrap();
        let result = if case["phase"] == "gas" {
            model
                .gas(&x, t, case["molar_volume_m3_per_mol"].as_f64().unwrap())
                .unwrap()
        } else {
            model
                .liquid(&x, t, case["pressure_Pa"].as_f64().unwrap())
                .unwrap()
        };
        let present = x.iter().filter(|v| **v > 0.0).count();
        assert_eq!(
            result.validity,
            if present == 1 {
                ValidityStatus::WithinDeclaredDomain
            } else {
                ValidityStatus::Unassessed
            },
            "{id}"
        );
        for name in PROPERTIES {
            let actual = property(&result, name);
            match case["expected"][name].as_f64() {
                Some(expected) => {
                    let (atol, rtol) = budget(&f, &format!("implementation_{name}"));
                    let actual = actual.unwrap_or_else(|e| panic!("{id} {name}: {e}"));
                    assert!(
                        (actual - expected).abs() <= atol + rtol * expected.abs(),
                        "{id} {name}: {actual} vs {expected}"
                    );
                    checked += 1;
                }
                None if case["phase"] == "gas" && name == "surface_tension" => {
                    assert!(matches!(
                        actual,
                        Err(CalculationError::UnsupportedCapability(_))
                    ));
                }
                None => assert!(
                    matches!(actual, Err(CalculationError::OutsideValidityDomain(_))),
                    "{id} {name}: expected outside the domain, got {actual:?}"
                ),
            }
        }
    }
    assert!(checked >= 36, "{checked}");
}

/// Model comparison with CoolProp 7.2.0 reference correlations for pure
/// components, gas volumes from the engine's own PR-exact root.
#[test]
fn pure_components_agree_with_coolprop_within_declared_method_budgets() {
    let f = fixture();
    let db = ComponentDatabase::embedded().unwrap();
    let mut failures = Vec::new();
    let mut worst: Vec<(String, f64)> = Vec::new();
    let cases = f["coolprop"].as_array().unwrap();
    assert!(cases.len() >= 100);
    for case in cases {
        let id = case["id"].as_str().unwrap();
        let component = case["component"].as_str().unwrap();
        let t = case["temperature_K"].as_f64().unwrap();
        let p = case["pressure_Pa"].as_f64().unwrap();
        let polar = case["polar"].as_bool().unwrap();
        let model = TransportModel::new(&[component]).unwrap();
        let temperature = Temperature::new(t).unwrap();
        let (result, phase) = if case["phase"] == "gas" {
            let composition = Composition::new(
                CompositionBasis::Mole,
                vec![(ComponentId::new(component).unwrap(), 1.0)],
                NormalizationPolicy::RequireUnitSum,
            )
            .unwrap();
            let cubic = CubicModel::new(
                &db,
                &composition,
                CubicVariant::PengRobinson,
                vec![vec![0.0]],
                "pure component",
            )
            .unwrap();
            let state = cubic
                .evaluate(
                    temperature,
                    Pressure::new(p).unwrap(),
                    RootSelection::LargestVolume,
                )
                .unwrap();
            (
                model.gas(&[1.0], temperature, state.molar_volume).unwrap(),
                "gas",
            )
        } else {
            (model.liquid(&[1.0], temperature, p).unwrap(), "liquid")
        };
        for name in PROPERTIES {
            let Some(expected) = case["expected"].get(name) else {
                continue;
            };
            let key = match (phase, name) {
                ("gas", "viscosity") if polar => "coolprop_gas_viscosity_polar",
                ("gas", "viscosity") if result.regime == TransportRegime::DenseFluid => {
                    "coolprop_dense_fluid_viscosity_nonpolar"
                }
                ("gas", "viscosity") => "coolprop_dilute_gas_viscosity_nonpolar",
                ("gas", _) => "coolprop_gas_thermal_conductivity_nonpolar",
                (_, "viscosity") => "coolprop_liquid_viscosity",
                (_, "thermal_conductivity") => "coolprop_liquid_thermal_conductivity",
                _ => "coolprop_surface_tension",
            };
            let actual = property(&result, name);
            let Some(expected) = expected.as_f64() else {
                // Polar gas conductivity: declared unavailable.
                if !matches!(actual, Err(CalculationError::OutsideValidityDomain(_))) {
                    failures.push(format!("{id} {name}: expected unavailable, got {actual:?}"));
                }
                continue;
            };
            let (atol, rtol) = budget(&f, key);
            match actual {
                Ok(v) => {
                    let deviation = v / expected - 1.0;
                    if let Some(w) = worst.iter_mut().find(|(k, _)| k == key) {
                        w.1 = w.1.max(deviation.abs());
                    } else {
                        worst.push((key.into(), deviation.abs()));
                    }
                    if (v - expected).abs() > atol + rtol * expected.abs() {
                        failures.push(format!(
                            "{id} {name}: {v:.6e} vs CoolProp {expected:.6e} ({:+.2}%)",
                            100.0 * deviation
                        ));
                    }
                }
                Err(e) => failures.push(format!("{id} {name}: unavailable inside the domain: {e}")),
            }
        }
    }
    for (key, w) in &worst {
        println!("{key}: largest |deviation| {:.2}%", 100.0 * w);
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
