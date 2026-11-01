use property_packages::{HomogeneousRequest, calculate};
use serde_json::{Value, json};

fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.16e}, expected={expected:.16e}, absolute tolerance={absolute}, relative tolerance={relative}"
    );
}

fn base_request() -> Value {
    serde_json::from_str(include_str!("../../../examples/binary-hydrocarbon.json")).unwrap()
}

fn run(value: Value) -> property_packages::HomogeneousResult {
    calculate(serde_json::from_value(value).unwrap()).unwrap()
}

#[test]
fn matches_independent_teqp_properties_and_all_homogeneous_roots() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/cubic-reference.json")).unwrap();
    let tolerances = &fixture["tolerances"];
    for case in fixture["cases"].as_array().unwrap() {
        let request: HomogeneousRequest = serde_json::from_value(case["request"].clone()).unwrap();
        let result = calculate(request).unwrap_or_else(|error| panic!("{}: {error}", case["id"]));
        assert_eq!(
            result.data_version,
            fixture["data_version"].as_str().unwrap()
        );
        assert_eq!(result.phase_stability, "not_assessed");
        let actual = serde_json::to_value(&result.properties).unwrap();
        for (field, tolerance) in [
            ("z", "z_absolute"),
            ("molar_density_mol_m3", "density_absolute"),
            ("mass_density_kg_m3", "density_absolute"),
            ("residual_enthalpy_j_mol", "enthalpy_absolute_j_mol"),
            ("residual_entropy_j_mol_k", "entropy_absolute_j_mol_k"),
        ] {
            close(
                actual[field].as_f64().unwrap(),
                case["expected"][field].as_f64().unwrap(),
                tolerances[tolerance].as_f64().unwrap(),
                tolerances["relative"].as_f64().unwrap(),
            );
        }
        for (field, tolerance) in [
            ("roots_z", "z_absolute"),
            ("ln_fugacity_coefficients", "ln_phi_absolute"),
        ] {
            let expected = case["expected"][field].as_array().unwrap();
            let values = actual[field].as_array().unwrap();
            assert_eq!(values.len(), expected.len(), "{}: {field}", case["id"]);
            for (value, expected) in values.iter().zip(expected) {
                close(
                    value.as_f64().unwrap(),
                    expected.as_f64().unwrap(),
                    tolerances[tolerance].as_f64().unwrap(),
                    tolerances["relative"].as_f64().unwrap(),
                );
            }
        }
        assert!(result.properties.pressure_relative_residual < 1e-8);
    }
}

#[test]
fn component_permutation_and_mass_basis_preserve_properties() {
    let original = base_request();
    let result = run(original.clone());
    let mut reversed = original.clone();
    reversed["components"] = json!(["ethane", "methane"]);
    reversed["fractions"] = json!([0.3, 0.7]);
    let reverse_result = run(reversed);
    close(
        result.properties.z,
        reverse_result.properties.z,
        1e-13,
        1e-12,
    );
    for i in 0..2 {
        close(
            result.properties.ln_fugacity_coefficients[i],
            reverse_result.properties.ln_fugacity_coefficients[1 - i],
            1e-12,
            1e-12,
        );
    }
    let catalog = thermo_data::Catalog::bundled().unwrap();
    let m1 = catalog.component("methane").unwrap().molar_mass().value();
    let m2 = catalog.component("ethane").unwrap().molar_mass().value();
    let total = 0.7 * m1 + 0.3 * m2;
    let mut mass = original;
    mass["basis"] = json!("mass");
    mass["fractions"] = json!([0.7 * m1 / total, 0.3 * m2 / total]);
    let mass_result = run(mass);
    close(
        result.properties.mass_density_kg_m3,
        mass_result.properties.mass_density_kg_m3,
        1e-11,
        1e-12,
    );
}

#[test]
fn zero_component_matches_pure_limit_and_low_pressure_approaches_ideal() {
    let mut pure = base_request();
    pure["temperature_k"] = json!(500.0);
    pure["pressure_pa"] = json!(1000.0);
    pure["components"] = json!(["methane"]);
    pure["fractions"] = json!([1.0]);
    pure["binary_interactions"] = json!([[0.0]]);
    let result = run(pure.clone());
    let mut with_zero = pure;
    with_zero["components"] = json!(["methane", "ethane"]);
    with_zero["fractions"] = json!([1.0, 0.0]);
    with_zero["binary_interactions"] = json!([[0.0, 0.02], [0.02, 0.0]]);
    let expanded = run(with_zero);
    close(result.properties.z, expanded.properties.z, 1e-13, 1e-12);
    close(
        result.properties.residual_enthalpy_j_mol,
        expanded.properties.residual_enthalpy_j_mol,
        1e-10,
        1e-10,
    );
    assert!((result.properties.z - 1.0).abs() < 1e-4);
    assert!(result.properties.ln_fugacity_coefficients[0].abs() < 1e-4);
}

#[test]
fn gibbs_identity_and_temperature_derivative_are_consistent() {
    let mut req = base_request();
    req["temperature_k"] = json!(400.0);
    let center = run(req.clone());
    let r = thermo_types::GAS_CONSTANT;
    let g_over_rt: f64 = center
        .mole_fractions
        .iter()
        .zip(&center.properties.ln_fugacity_coefficients)
        .map(|(x, phi)| x * phi)
        .sum();
    close(
        g_over_rt,
        (center.properties.residual_enthalpy_j_mol
            - 400.0 * center.properties.residual_entropy_j_mol_k)
            / (r * 400.0),
        1e-12,
        1e-12,
    );
    for step in [0.01, 0.001] {
        let mut lower = req.clone();
        let mut upper = req.clone();
        lower["temperature_k"] = json!(400.0 - step);
        upper["temperature_k"] = json!(400.0 + step);
        let gibbs = |v: Value| {
            let result = run(v);
            result
                .mole_fractions
                .iter()
                .zip(&result.properties.ln_fugacity_coefficients)
                .map(|(x, phi)| x * phi)
                .sum::<f64>()
        };
        let derivative = (gibbs(upper) - gibbs(lower)) / (2.0 * step);
        close(
            derivative,
            -center.properties.residual_enthalpy_j_mol / (r * 400.0_f64.powi(2)),
            1e-10,
            1e-6,
        );
    }
}

#[test]
fn unsupported_and_malformed_requests_do_not_return_results() {
    for (field, value) in [
        ("pressure_pa", json!(0.0)),
        ("temperature_k", json!(2000.0)),
        ("fractions", json!([0.2, 0.3])),
        ("binary_interactions", json!([])),
        ("components", json!(["unknown", "ethane"])),
    ] {
        let mut request = base_request();
        request[field] = value;
        assert!(calculate(serde_json::from_value(request).unwrap()).is_err());
    }
    let mut unknown = base_request();
    unknown["pressure_bar"] = json!(30.0);
    assert!(serde_json::from_value::<HomogeneousRequest>(unknown).is_err());
    let mut missing = base_request();
    missing
        .as_object_mut()
        .unwrap()
        .remove("binary_interactions");
    assert!(serde_json::from_value::<HomogeneousRequest>(missing).is_err());
}
