//! A deterministic robustness sweep, not an independent accuracy comparison.
use property_packages::{HomogeneousRequest, calculate};
use serde_json::json;

#[test]
#[ignore = "broader deterministic sweep; run manually or in scheduled validation CI"]
fn homogeneous_state_grid() {
    let mut passed = 0;
    let mut failures = Vec::new();
    for model in ["PR", "PR-BM"] {
        for (components, fractions, matrix) in [
            (json!(["methane"]), json!([1.0]), json!([[0.0]])),
            (
                json!(["methane", "ethane"]),
                json!([0.7, 0.3]),
                json!([[0.0, 0.02], [0.02, 0.0]]),
            ),
        ] {
            for temperature in [150.0, 250.0, 300.0, 500.0, 700.0] {
                for pressure in [1e5, 1e6, 5e6, 1e7, 5e7] {
                    for root in ["smallest", "largest"] {
                        let input = json!({"model": model, "components": components, "fractions": fractions,
                            "basis": "mole", "binary_interactions": matrix, "temperature_k": temperature,
                            "pressure_pa": pressure, "root": root});
                        let request: HomogeneousRequest =
                            serde_json::from_value(input.clone()).unwrap();
                        match calculate(request) {
                            Ok(result) => {
                                assert!(result.properties.pressure_relative_residual <= 1e-8);
                                assert_eq!(result.phase_stability, "not_assessed");
                                passed += 1;
                            }
                            Err(error) => {
                                failures.push(json!({"input": input, "error": error.to_string()}))
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "{}",
        json!({"kind": "homogeneous_robustness_not_accuracy", "cases": 200,
        "passed": passed, "failed": failures.len(), "failures": failures})
    );
    assert!(
        failures.is_empty(),
        "all unsupported states and failures are listed in the sweep output"
    );
    assert_eq!(passed, 200);
}
