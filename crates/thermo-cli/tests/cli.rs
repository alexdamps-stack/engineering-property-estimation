use std::{path::Path, process::Command};

#[test]
fn calculates_example_and_emits_explicit_stability_status() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/natural-gas.json");
    let output = Command::new(env!("CARGO_BIN_EXE_thermo-cli"))
        .args(["calculate", "--input"])
        .arg(path)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["phase_stability"], "not_assessed");
    assert!(result["properties"]["mass_density_kg_m3"].as_f64().unwrap() > 0.0);
}

#[test]
fn failure_has_nonzero_exit_and_no_success_output() {
    for args in [
        vec!["calculate", "--input", "this-file-does-not-exist.json"],
        vec!["unsupported-command"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_thermo-cli"))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        let error: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
        assert!(error["error"].is_string());
    }
}
