//! Engineering liquid density (#228) against the committed independent
//! fixture `validation/fixtures/liquid-density-v1.json`:
//! - implementation verification against `chemicals` 1.5.2 (COSTALD,
//!   Hankinson-Thomson mixing, COSTALD-Tait), budget from the fixture;
//! - experimental validation against CoolProp 7.2.0 reference equations of
//!   state, per-band budgets from the fixture, qualified points only;
//! - the qualification policy (refused / qualified / unassessed) and the
//!   declared domain failures.
use serde_json::Value;
use thermo_data::{COSTALD_DATA_VERSION, CostaldDatabase};
use thermo_models::{
    COSTALD_MAXIMUM_PRESSURE_PA, COSTALD_MIXTURE_QUALIFIED_COMPONENTS,
    COSTALD_QUALIFIED_COMPONENTS, COSTALD_REDUCED_TEMPERATURE_RANGE, COSTALD_TRACE_MOLE_FRACTION,
    CostaldLiquidDensity,
};
use thermo_types::{CalculationError, MOLAR_GAS_CONSTANT, ValidityStatus};

fn fixture() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/fixtures/liquid-density-v1.json"
    ))
    .unwrap()
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

fn model(case: &Value) -> CostaldLiquidDensity {
    CostaldLiquidDensity::embedded(&strings(&case["components"])).unwrap()
}

#[test]
fn fixture_matches_engine_inputs_and_policy() {
    let f = fixture();
    assert_eq!(f["schema"], "liquid-density-reference-v1");
    assert_eq!(
        f["gas_constant_J_per_mol_K"].as_f64(),
        Some(MOLAR_GAS_CONSTANT)
    );
    let domain = &f["declared_domain"];
    assert_eq!(
        floats(&domain["reduced_temperature"]),
        COSTALD_REDUCED_TEMPERATURE_RANGE
    );
    assert_eq!(
        domain["maximum_pressure_Pa"].as_f64(),
        Some(COSTALD_MAXIMUM_PRESSURE_PA)
    );
    let policy = &f["qualification_policy"];
    assert_eq!(
        strings(&policy["qualified_components"]),
        COSTALD_QUALIFIED_COMPONENTS
    );
    assert_eq!(
        strings(&policy["mixture_qualified_components"]),
        COSTALD_MIXTURE_QUALIFIED_COMPONENTS
    );
    assert_eq!(
        policy["trace_mole_fraction"].as_f64(),
        Some(COSTALD_TRACE_MOLE_FRACTION)
    );
    // The fixture was generated from the committed extract.
    let costald =
        include_bytes!("../../../data/source-extracts/liquid-density-costald.provisional.json");
    assert!(
        String::from_utf8_lossy(costald)
            .contains(&format!("\"version\": \"{COSTALD_DATA_VERSION}\""))
    );
    assert_eq!(
        CostaldDatabase::embedded().unwrap().metadata().version,
        COSTALD_DATA_VERSION
    );
}

#[test]
fn implementation_agrees_with_independent_chemicals_evaluation() {
    let f = fixture();
    let budget = &f["tolerances"]["implementation"];
    let rtol = budget["rtol"].as_f64().unwrap();
    let atol = budget["atol"].as_f64().unwrap();
    let close = |a: f64, e: f64| (a - e).abs() <= atol + rtol * e.abs();
    let cases = f["implementation"].as_array().unwrap();
    assert!(cases.len() >= 100, "fixture has {} cases", cases.len());
    let (mut values, mut refused) = (0, 0);
    for case in cases {
        let id = case["id"].as_str().unwrap();
        let x = floats(&case["mole_fractions"]);
        let t = case["temperature_K"].as_f64().unwrap();
        let p = case["pressure_Pa"].as_f64().unwrap();
        let result = model(case).evaluate(&x, t, p);
        match case["engine_scope"].as_str().unwrap() {
            "refused" => {
                refused += 1;
                assert!(
                    matches!(result, Err(CalculationError::OutsideValidityDomain(_))),
                    "{id}: {result:?}"
                );
            }
            scope => {
                values += 1;
                let r = result.unwrap_or_else(|e| panic!("{id}: {e}"));
                let expected_validity = if scope == "qualified" {
                    ValidityStatus::WithinDeclaredDomain
                } else {
                    ValidityStatus::Unassessed
                };
                assert_eq!(r.validity, expected_validity, "{id}");
                assert_eq!(
                    r.pressure_corrected,
                    case["pressure_corrected"].as_bool().unwrap(),
                    "{id}"
                );
                for (actual, key) in [
                    (
                        r.saturated_molar_volume,
                        "saturated_molar_volume_m3_per_mol",
                    ),
                    (
                        r.pseudo_saturation_pressure_pa,
                        "pseudo_saturation_pressure_Pa",
                    ),
                    (r.molar_volume, "molar_volume_m3_per_mol"),
                    (r.mass_density, "mass_density_kg_per_m3"),
                    (r.reduced_temperature, "reduced_temperature"),
                ] {
                    let expected = case[key].as_f64().unwrap();
                    assert!(
                        close(actual, expected),
                        "{id} {key}: {actual} vs {expected}"
                    );
                }
            }
        }
    }
    assert!(values > 0 && refused > 0);
}

#[test]
fn declared_domain_failures() {
    let f = fixture();
    for case in f["failures"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        assert_eq!(case["error"], "OutsideValidityDomain");
        let result = model(case).evaluate(
            &floats(&case["mole_fractions"]),
            case["temperature_K"].as_f64().unwrap(),
            case["pressure_Pa"].as_f64().unwrap(),
        );
        assert!(
            matches!(result, Err(CalculationError::OutsideValidityDomain(_))),
            "{id}: {result:?}"
        );
    }
}

#[test]
fn qualified_scope_agrees_with_reference_equations_of_state() {
    let f = fixture();
    let points = f["experimental"].as_array().unwrap();
    let mut qualified = 0;
    let mut kinds = std::collections::BTreeSet::new();
    for point in points {
        let id = point["id"].as_str().unwrap();
        let x = floats(&point["mole_fractions"]);
        let t = point["temperature_K"].as_f64().unwrap();
        let p = point["pressure_Pa"].as_f64().unwrap();
        let result = model(point).evaluate(&x, t, p);
        match point["engine_scope"].as_str().unwrap() {
            "refused" => assert!(
                matches!(result, Err(CalculationError::OutsideValidityDomain(_))),
                "{id}: {result:?}"
            ),
            "unassessed" => {
                assert_eq!(result.unwrap().validity, ValidityStatus::Unassessed, "{id}")
            }
            "qualified" => {
                qualified += 1;
                kinds.insert(point["kind"].as_str().unwrap().to_string());
                let r = result.unwrap_or_else(|e| panic!("{id}: {e}"));
                assert_eq!(r.validity, ValidityStatus::WithinDeclaredDomain, "{id}");
                let reference = point["reference_molar_density_mol_per_m3"]
                    .as_f64()
                    .unwrap();
                let deviation = 1.0 / r.molar_volume / reference - 1.0;
                let budget = point["budget"].as_f64().unwrap();
                assert!(
                    deviation.abs() <= budget,
                    "{id}: deviation {:.3} % exceeds {:.1} %",
                    100.0 * deviation,
                    100.0 * budget
                );
            }
            other => panic!("{id}: unknown scope {other}"),
        }
    }
    assert!(qualified >= 100, "{qualified} qualified points");
    assert_eq!(
        kinds.into_iter().collect::<Vec<_>>(),
        ["bubble", "compressed", "saturated"]
    );
}
