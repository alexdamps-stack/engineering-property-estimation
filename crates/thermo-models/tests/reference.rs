use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, RootSelection};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure, Temperature,
};

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}
fn close(actual: f64, expected: f64, atol: f64, rtol: f64, context: &str) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= atol + rtol * expected.abs(),
        "{context}: actual={actual:.16e}, expected={expected:.16e}, error={:.4e}, budget={:.4e}",
        (actual - expected).abs(),
        atol + rtol * expected.abs()
    );
}

#[test]
fn independent_teqp_roots_fugacity_and_caloric_grid() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/cubic-v1.json")).unwrap();
    let db = ComponentDatabase::embedded().unwrap();
    assert_eq!(fixture["cases"].as_array().unwrap().len(), 60);
    for case in fixture["cases"].as_array().unwrap() {
        let entries = case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .zip(case["x"].as_array().unwrap())
            .map(|(id, x)| (ComponentId::new(id.as_str().unwrap()).unwrap(), number(x)))
            .collect();
        let comp = Composition::new(
            CompositionBasis::Mole,
            entries,
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let variant = if case["variant"] == "PR" {
            CubicVariant::PengRobinson
        } else {
            CubicVariant::PengRobinsonBostonMathias
        };
        let kij = case["kij"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| r.as_array().unwrap().iter().map(number).collect())
            .collect();
        let model = CubicModel::new(
            &db,
            &comp,
            variant,
            kij,
            "fixture explicit synthetic interactions",
        )
        .unwrap();
        let t = Temperature::new(number(&case["T"])).unwrap();
        let p = Pressure::new(number(&case["p"])).unwrap();
        let roots = model.compressibility_roots(t, p).unwrap();
        let expected = case["roots"].as_array().unwrap();
        assert_eq!(roots.len(), expected.len(), "{case}");
        for (i, root) in expected.iter().enumerate() {
            let context = format!(
                "{:?}, {:?}, T={t:?}, p={p:?}, root {i}",
                variant, case["ids"]
            );
            close(roots[i], number(&root["Z"]), 2e-10, 2e-9, &context);
            if root["mechanically_stable"] == false {
                continue;
            }
            let selection = if i == 0 {
                RootSelection::SmallestVolume
            } else {
                RootSelection::LargestVolume
            };
            let state = model.evaluate(t, p, selection).unwrap();
            for (key, actual) in [
                ("Z", state.compressibility),
                ("mass_density", state.mass_density),
                ("h_res", state.residual_enthalpy),
                ("s_res", state.residual_entropy),
                ("h", state.enthalpy),
                ("s", state.entropy),
                ("cp", state.cp),
                ("cv", state.cv),
                ("sound", state.speed_of_sound),
            ] {
                let tolerance = &fixture["tolerances"][key];
                close(
                    actual,
                    number(&root[key]),
                    number(&tolerance[0]),
                    number(&tolerance[1]),
                    &format!("{context} {key}"),
                );
            }
            for (actual, expected) in state
                .ln_fugacity_coefficients
                .iter()
                .zip(root["ln_phi"].as_array().unwrap())
            {
                close(
                    *actual,
                    number(expected),
                    2e-9,
                    2e-9,
                    &format!("{context} ln_phi"),
                );
            }
        }
        let state = model.evaluate(t, p, RootSelection::MinimumGibbs).unwrap();
        let gibbs = |r: &Value| -> f64 {
            r["ln_phi"]
                .as_array()
                .unwrap()
                .iter()
                .zip(comp.entries())
                .map(|(v, (_, x))| number(v) * x)
                .sum()
        };
        let best = expected
            .iter()
            .filter(|r| r["mechanically_stable"] == true)
            .min_by(|a, b| gibbs(a).total_cmp(&gibbs(b)))
            .unwrap();
        close(
            state.compressibility,
            number(&best["Z"]),
            2e-10,
            2e-9,
            "minimum fixed-composition Gibbs root",
        );
    }
}

#[test]
fn nasa7_against_cantera_offsets_and_both_intervals() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/cubic-v1.json")).unwrap();
    let db = ComponentDatabase::embedded().unwrap();
    for case in fixture["ideal_gas"].as_array().unwrap() {
        let c = db.get(case["id"].as_str().unwrap()).unwrap();
        let ideal = c
            .ideal_gas(Temperature::new(number(&case["T"])).unwrap())
            .unwrap();
        close(ideal.enthalpy, number(&case["h"]), 2e-6, 1e-11, "NASA7 h");
        close(ideal.entropy, number(&case["s"]), 2e-8, 1e-11, "NASA7 s");
        close(
            ideal.heat_capacity,
            number(&case["cp"]),
            2e-8,
            1e-11,
            "NASA7 cp",
        );
    }
}

#[test]
fn helmholtz_matches_high_precision_reference_through_dilute_limit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/fixtures/helmholtz-dilute-v1.json"
    ))
    .unwrap();
    assert_eq!(fixture["cases"].as_array().unwrap().len(), 6);
    let db = ComponentDatabase::embedded().unwrap();
    let composition = Composition::new(
        CompositionBasis::Mole,
        vec![(ComponentId::new("methane").unwrap(), 1.0)],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let scalar = |key: &str| case[key].as_str().unwrap().parse::<f64>().unwrap();
        let variant = if case["variant"] == "PR" {
            CubicVariant::PengRobinson
        } else {
            CubicVariant::PengRobinsonBostonMathias
        };
        let model = CubicModel::new(
            &db,
            &composition,
            variant,
            vec![vec![0.0]],
            "pure component",
        )
        .unwrap();
        let t = Temperature::new(scalar("temperature_K")).unwrap();
        let energy = model
            .residual_helmholtz(t, scalar("volume_m3"), scalar("moles"))
            .unwrap();
        // Relative agreement is essential: an absolute energy floor would hide
        // the wrong-sign dilute result. References were evaluated at 80 digits.
        close(energy, scalar("energy_J"), 0.0, 5e-14, &case.to_string());
    }
}
