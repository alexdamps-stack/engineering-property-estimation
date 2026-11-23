//! Independent PS comparisons on the PR and PR-BM caloric reference.
mod common;

use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::{PsFlashOptions, ps_flash};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase,
    Pressure,
};

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn model(case: &Value, database: &ComponentDatabase) -> CubicModel {
    let composition = Composition::new(
        CompositionBasis::Mole,
        case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .zip(case["z"].as_array().unwrap())
            .map(|(id, fraction)| {
                (
                    ComponentId::new(id.as_str().unwrap()).unwrap(),
                    number(fraction),
                )
            })
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let variant = match case["variant"].as_str().unwrap() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        other => panic!("unknown PS reference variant {other}"),
    };
    CubicModel::new(
        database,
        &composition,
        variant,
        case["kij"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| row.as_array().unwrap().iter().map(number).collect())
            .collect(),
        "explicit matched PS reference interactions",
    )
    .unwrap()
}

fn fixture() -> Value {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/rotating/ps-v1.json")).unwrap();
    assert_eq!(fixture["schema"], "a14-ps-reference-v1");
    fixture
}

#[test]
fn pressure_entropy_flash_matches_independent_references() {
    let fixture = fixture();
    let cases = fixture["ps"].as_array().unwrap();
    assert_eq!(cases.len(), 28);
    common::assert_input_sha256(
        "ps-v1.json",
        "data/source-extracts/components.provisional.json",
        fixture["input_sha256"].as_str().unwrap(),
    );
    for (path, hash) in fixture["adapter_sha256"].as_object().unwrap() {
        common::assert_input_sha256(
            "ps-v1.json",
            &format!("validation/{path}"),
            hash.as_str().unwrap(),
        );
    }
    let database = ComponentDatabase::embedded().unwrap();
    let options = PsFlashOptions::default();
    for case in cases {
        let context = format!("{} {}", case["variant"], case["id"]);
        let compare = |key: &str, actual: f64, expected: f64| {
            let budget = &fixture["tolerances"][key];
            let limit = number(&budget[0]) + number(&budget[1]) * expected.abs();
            assert!(
                (actual - expected).abs() <= limit,
                "{context}: {key}: {actual:e} differs from {expected:e} by more than {limit:e}"
            );
        };
        let pressure = Pressure::new(number(&case["p"])).unwrap();
        let state = ps_flash(
            &model(case, &database),
            pressure,
            number(&case["target_s"]),
            options,
        )
        .unwrap_or_else(|error| panic!("{context}: {error}"));
        compare("temperature_K", state.temperature_k, number(&case["T"]));
        compare(
            "enthalpy_J_per_mol",
            state.enthalpy,
            number(&case["target_h"]),
        );
        compare(
            "entropy_J_per_mol_K",
            state.entropy,
            number(&case["target_s"]),
        );
        assert!(state.scaled_entropy_residual.abs() <= options.energy_tolerance);
        assert_eq!(
            state.phases.len() as u64,
            case["phases"].as_u64().unwrap(),
            "{context}"
        );
        if let [phase] = state.phases.as_slice() {
            compare("phase_fraction", phase.phase_fraction, 1.0);
            let expected = case["z"].as_array().unwrap();
            assert_eq!(phase.composition.len(), expected.len(), "{context}");
            for (&actual, expected) in phase.composition.iter().zip(expected) {
                compare("composition", actual, number(expected));
            }
            continue;
        }
        for phase in &state.phases {
            let (key, fraction) = if phase.phase == Phase::Vapor {
                ("vapor", number(&case["vapor_fraction"]))
            } else {
                ("liquid", 1.0 - number(&case["vapor_fraction"]))
            };
            compare("phase_fraction", phase.phase_fraction, fraction);
            let expected = case[key].as_array().unwrap();
            assert_eq!(phase.composition.len(), expected.len(), "{context}");
            for (&actual, expected) in phase.composition.iter().zip(expected) {
                compare("composition", actual, number(expected));
            }
        }
    }
}

#[test]
fn pressure_entropy_flash_rejects_invalid_and_unbracketed_targets() {
    let database = ComponentDatabase::embedded().unwrap();
    let model = model(&fixture()["ps"][0], &database);
    let pressure = Pressure::new(3e6).unwrap();
    let options = PsFlashOptions::default();
    assert!(matches!(
        ps_flash(&model, pressure, f64::NAN, options),
        Err(CalculationError::InvalidInput(_))
    ));
    assert!(matches!(
        ps_flash(&model, pressure, 1e6, options),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        ps_flash(
            &model,
            pressure,
            0.0,
            PsFlashOptions {
                max_iterations: 0,
                ..options
            }
        ),
        Err(CalculationError::NonConvergence(_))
    ));
}
