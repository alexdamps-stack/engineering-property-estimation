//! Mandatory offline Clapeyron checks. The live CI check supplies a freshly
//! generated candidate through THERMO_CLAPEYRON_REFERENCE without replacing Git fixtures.
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, HomogeneousState, RootSelection};
use thermo_solvers::{PhFlashOptions, ph_flash};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure};

fn fixture() -> Value {
    let input = match std::env::var_os("THERMO_CLAPEYRON_REFERENCE") {
        Some(path) => std::fs::read_to_string(path).expect("read live Clapeyron candidate"),
        None => include_str!("../../../validation/fixtures/clapeyron-v1.json").into(),
    };
    let f: Value = serde_json::from_str(&input).expect("valid Clapeyron JSON");
    assert_eq!(f["schema"], "clapeyron-reference-v1");
    assert_eq!(f["clapeyron_version"], "0.6.30");
    assert_eq!(f["julia_version"], "1.11.7");
    assert_eq!(f["ph"].as_array().unwrap().len(), 28);
    assert_eq!(f["homogeneous"].as_array().unwrap().len(), 60);
    assert_eq!(
        f["homogeneous"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|case| case["evaluation_precision_bits"] == 256)
            .count(),
        2
    );
    f
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn model(case: &Value, fractions: &str) -> CubicModel {
    let composition = Composition::new(
        CompositionBasis::Mole,
        case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .zip(numbers(&case[fractions]))
            .map(|(id, x)| (ComponentId::new(id.as_str().unwrap()).unwrap(), x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let variant = match case["variant"].as_str().unwrap() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        other => panic!("unknown Julia fixture variant {other}"),
    };
    CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition,
        variant,
        case["kij"]
            .as_array()
            .unwrap()
            .iter()
            .map(numbers)
            .collect(),
        "explicit matched Clapeyron fixture interactions",
    )
    .unwrap()
}

fn close(f: &Value, key: &str, actual: f64, expected: f64, context: &str) {
    let b = numbers(&f["tolerances"][key]);
    let limit = b[0] + b[1] * expected.abs();
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= limit,
        "{context}: {key}: Rust {actual:e}, Clapeyron {expected:e}, deviation {:e}, budget {limit:e}",
        (actual - expected).abs()
    );
}

fn compare_state(f: &Value, s: &HomogeneousState, expected: &Value, context: &str) {
    for (key, actual) in [
        ("Z", s.compressibility),
        ("molar_volume", s.molar_volume),
        ("mass_density", s.mass_density),
        ("h", s.enthalpy),
        ("s", s.entropy),
        ("h_res", s.residual_enthalpy),
        ("s_res", s.residual_entropy),
        ("cp", s.cp),
        ("cv", s.cv),
        ("sound", s.speed_of_sound),
    ] {
        close(f, key, actual, expected[key].as_f64().unwrap(), context);
    }
    let ln_phi = numbers(&expected["ln_phi"]);
    assert_eq!(ln_phi.len(), s.ln_fugacity_coefficients.len());
    for (actual, expected) in s.ln_fugacity_coefficients.iter().zip(ln_phi) {
        close(f, "ln_phi", *actual, expected, context);
    }
}

#[test]
fn homogeneous_properties_match_clapeyron() {
    let f = fixture();
    let mut max_density: f64 = 0.0;
    let mut max_enthalpy: f64 = 0.0;
    let mut max_cp: f64 = 0.0;
    for case in f["homogeneous"].as_array().unwrap() {
        let m = model(case, "x");
        let t = thermo_types::Temperature::new(case["T"].as_f64().unwrap()).unwrap();
        let p = Pressure::new(case["p"].as_f64().unwrap()).unwrap();
        let context = format!(
            "{} {:?}, T={}, p={}",
            case["variant"],
            case["ids"],
            t.value(),
            p.value()
        );
        assert!(!case["roots"].as_array().unwrap().is_empty());
        for root in case["roots"].as_array().unwrap() {
            let selection = match root["selection"].as_str().unwrap() {
                "smallest-volume" => RootSelection::SmallestVolume,
                "largest-volume" => RootSelection::LargestVolume,
                other => panic!("unknown selection {other}"),
            };
            let state = m.evaluate(t, p, selection).unwrap();
            compare_state(&f, &state, &root["state"], &context);
            max_density = max_density
                .max((state.mass_density - root["state"]["mass_density"].as_f64().unwrap()).abs());
            max_enthalpy =
                max_enthalpy.max((state.enthalpy - root["state"]["h"].as_f64().unwrap()).abs());
            max_cp = max_cp.max((state.cp - root["state"]["cp"].as_f64().unwrap()).abs());
        }
    }
    println!(
        "60 homogeneous Clapeyron comparisons: max |d density|={max_density:e} kg/m3, |dh|={max_enthalpy:e} J/mol, |d cp|={max_cp:e} J/(mol K)"
    );
}

#[test]
fn ph_equilibrium_matches_native_clapeyron_solver() {
    let f = fixture();
    let mut max_temperature: f64 = 0.0;
    let mut max_fraction: f64 = 0.0;
    let mut max_enthalpy: f64 = 0.0;
    for case in f["ph"].as_array().unwrap() {
        let m = model(case, "z");
        let p = Pressure::new(case["p"].as_f64().unwrap()).unwrap();
        let target = case["target_h"].as_f64().unwrap();
        let r = ph_flash(&m, p, target, PhFlashOptions::default()).unwrap();
        let context = format!("{} {}", case["variant"], case["id"]);
        let expected_t = case["T"].as_f64().unwrap();
        close(&f, "temperature", r.temperature_k, expected_t, &context);
        close(&f, "h", r.enthalpy, case["h"].as_f64().unwrap(), &context);
        max_enthalpy = max_enthalpy.max((r.enthalpy - case["h"].as_f64().unwrap()).abs());
        assert!(r.scaled_energy_residual.abs() <= 1e-10);
        // Saturation endpoints may retain an incipient phase (Rust) or contain
        // only the material phase (Clapeyron). Compare material phases using
        // the declared phase-fraction budget, then their complete properties.
        let threshold = f["tolerances"]["fraction"][0].as_f64().unwrap();
        close(
            &f,
            "fraction",
            case["phases"]
                .as_array()
                .unwrap()
                .iter()
                .map(|ph| ph["fraction"].as_f64().unwrap())
                .sum(),
            1.0,
            &context,
        );
        let phases: Vec<_> = case["phases"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|ph| ph["fraction"].as_f64().unwrap() > threshold)
            .collect();
        let actual_phases: Vec<_> = r
            .phases
            .iter()
            .filter(|ph| ph.phase_fraction > threshold)
            .collect();
        assert!(!phases.is_empty(), "{context}: no material phase");
        assert_eq!(
            actual_phases.len(),
            phases.len(),
            "{context}: material phase count"
        );
        max_temperature = max_temperature.max((r.temperature_k - expected_t).abs());
        for (actual, expected) in actual_phases.iter().zip(phases) {
            let beta = expected["fraction"].as_f64().unwrap();
            close(&f, "fraction", actual.phase_fraction, beta, &context);
            max_fraction = max_fraction.max((actual.phase_fraction - beta).abs());
            assert_eq!(
                actual.composition.len(),
                expected["composition"].as_array().unwrap().len(),
                "{context}: component count"
            );
            for (x, y) in actual
                .composition
                .iter()
                .zip(numbers(&expected["composition"]))
            {
                close(&f, "composition", *x, y, &context);
            }
            compare_state(&f, &actual.state, &expected["state"], &context);
        }
    }
    println!(
        "28 native Clapeyron PH comparisons: max |dT|={max_temperature:e} K, |d beta|={max_fraction:e}, |dh|={max_enthalpy:e} J/mol"
    );
}
