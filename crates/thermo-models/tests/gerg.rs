use serde_json::Value;
use thermo_models::{GERG_GAS_CONSTANT, GergModel};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Temperature,
};

fn model(entries: &[(&str, f64)]) -> GergModel {
    GergModel::new(
        &Composition::new(
            CompositionBasis::Mole,
            entries
                .iter()
                .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
                .collect(),
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap(),
    )
    .unwrap()
}
fn near(name: &str, actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{name}: actual={actual:.16e}, reference={expected:.16e}"
    );
}

#[test]
fn independent_nist_properties_and_derivatives() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/gerg-v1.json")).unwrap();
    let atol = fixture["tolerance"]["absolute"].as_f64().unwrap();
    let rtol = fixture["tolerance"]["relative"].as_f64().unwrap();
    let mut admissible = 0;
    for case in fixture["cases"].as_array().unwrap() {
        let entries: Vec<_> = case["composition"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(id, x)| (id.as_str(), x.as_f64().unwrap()))
            .collect();
        let m = model(&entries);
        let t = Temperature::new(case["temperature_k"].as_f64().unwrap()).unwrap();
        let rho = case["molar_density"].as_f64().unwrap();
        let a = m.helmholtz(t, rho).unwrap();
        let e = &case["expected"];
        for (key, v) in [
            ("ar", a.ar),
            ("ar_d", a.ar_d),
            ("ar_dd", a.ar_dd),
            ("ar_t", a.ar_t),
            ("ar_tt", a.ar_tt),
            ("ar_dt", a.ar_dt),
            ("a0", a.a0),
            ("a0_t", a.a0_t),
            ("a0_tt", a.a0_tt),
            ("reducing_temperature", m.reducing_temperature()),
            ("reducing_density", m.reducing_density()),
        ] {
            near(
                &format!("{} {key}", case["id"]),
                v,
                e[key].as_f64().unwrap(),
                atol,
                rtol,
            );
        }
        let regular = e["dp_drho"].as_f64().unwrap() > 0.
            && e["cv"].as_f64().unwrap() > 0.
            && e["pressure_pa"].as_f64().unwrap() > 0.
            && e["pressure_pa"].as_f64().unwrap() <= 35e6;
        if !regular {
            assert!(matches!(
                m.evaluate(t, rho),
                Err(CalculationError::OutsideValidityDomain(_))
            ));
            continue;
        }
        let s = m.evaluate(t, rho).unwrap();
        admissible += 1;
        for (key, v) in [
            ("pressure_pa", s.pressure_pa),
            ("compressibility", s.compressibility),
            ("dp_drho", s.dp_drho),
            ("dp_dt", s.dp_dt),
            ("internal_energy", s.internal_energy),
            ("enthalpy", s.enthalpy),
            ("entropy", s.entropy),
            ("cv", s.cv),
            ("cp", s.cp),
            ("speed_of_sound", s.speed_of_sound),
            ("helmholtz_energy", s.helmholtz_energy),
            ("gibbs_energy", s.gibbs_energy),
        ] {
            near(
                &format!("{} {key}", case["id"]),
                v,
                e[key].as_f64().unwrap(),
                atol,
                rtol,
            );
        }
    }
    println!("273 derivative states; {admissible} locally admissible property states");
    assert!(admissible >= 250);
}

#[test]
fn identities_finite_differences_and_order_invariance() {
    let entries = [
        ("methane", 0.9),
        ("ethane", 0.04),
        ("propane", 0.02),
        ("nitrogen", 0.02),
        ("carbon-dioxide", 0.02),
    ];
    let m = model(&entries);
    let mut reversed = entries;
    reversed.reverse();
    let other = model(&reversed);
    for (t, rho) in [(250., 100.), (300., 1000.), (400., 4000.)] {
        let eval = |temp, d| m.evaluate(Temperature::new(temp).unwrap(), d).unwrap();
        let s = eval(t, rho);
        let permuted = other.evaluate(Temperature::new(t).unwrap(), rho).unwrap();
        near("permutation", s.enthalpy, permuted.enthalpy, 1e-10, 1e-12);
        near(
            "h=u+p/rho",
            s.enthalpy,
            s.internal_energy + s.pressure_pa / rho,
            1e-10,
            1e-12,
        );
        near(
            "g=h-Ts",
            s.gibbs_energy,
            s.enthalpy - t * s.entropy,
            1e-10,
            1e-12,
        );
        let dt = 0.001;
        let dr = rho * 1e-5;
        near(
            "dp/dT",
            s.dp_dt,
            (eval(t + dt, rho).pressure_pa - eval(t - dt, rho).pressure_pa) / (2. * dt),
            1e-6,
            2e-8,
        );
        near(
            "dp/drho",
            s.dp_drho,
            (eval(t, rho + dr).pressure_pa - eval(t, rho - dr).pressure_pa) / (2. * dr),
            1e-6,
            2e-8,
        );
        near(
            "cv=du/dT",
            s.cv,
            (eval(t + dt, rho).internal_energy - eval(t - dt, rho).internal_energy) / (2. * dt),
            1e-6,
            2e-8,
        );
        near(
            "ds/drho=-dp/dT/rho^2",
            -s.dp_dt / rho.powi(2),
            (eval(t, rho + dr).entropy - eval(t, rho - dr).entropy) / (2. * dr),
            1e-8,
            2e-8,
        );
    }
    let pure = model(&[("methane", 1.)]);
    let zeros = model(&[("argon", 0.), ("methane", 1.), ("water", 0.)]);
    let t = Temperature::new(300.).unwrap();
    near(
        "retained zero pure limit",
        pure.evaluate(t, 100.).unwrap().cp,
        zeros.evaluate(t, 100.).unwrap().cp,
        0.,
        0.,
    );
    let trace = model(&[("methane", 1.0 - 1e-12), ("water", 1e-12)]);
    near(
        "vanishing-component continuity",
        pure.evaluate(t, 100.).unwrap().enthalpy,
        trace.evaluate(t, 100.).unwrap().enthalpy,
        1e-7,
        1e-10,
    );
    let subnormal = model(&[("methane", 1.0), ("water", 1e-320)]);
    assert!(subnormal.evaluate(t, 0.001).unwrap().entropy.is_finite());
    near(
        "dilute gas",
        pure.evaluate(t, 0.001).unwrap().compressibility,
        1.,
        1e-6,
        0.,
    );
    near(
        "cp-cv dilute",
        pure.evaluate(t, 0.001).unwrap().cp - pure.evaluate(t, 0.001).unwrap().cv,
        GERG_GAS_CONSTANT,
        1e-5,
        0.,
    );
}

#[test]
fn explicit_unsupported_inputs_and_local_instability() {
    assert_eq!(GergModel::supported_components().len(), 21);
    for (id, x) in [("methanol", 1.), ("unknown", 0.)] {
        let c = Composition::new(
            CompositionBasis::Mole,
            vec![
                (ComponentId::new("methane").unwrap(), 1. - x),
                (ComponentId::new(id).unwrap(), x),
            ],
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        assert!(matches!(
            GergModel::new(&c),
            Err(CalculationError::MissingParameter(_))
        ));
    }
    let c = Composition::new(
        CompositionBasis::Mass,
        vec![(ComponentId::new("methane").unwrap(), 1.)],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    assert!(matches!(
        GergModel::new(&c),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    let m = model(&[("methane", 1.)]);
    assert!(matches!(
        m.evaluate(Temperature::new(300.).unwrap(), f64::NAN),
        Err(CalculationError::InvalidInput(_))
    ));
    for rho in [0., -1., f64::NAN, f64::INFINITY, 0.0001, 30001.] {
        assert!(m.evaluate(Temperature::new(300.).unwrap(), rho).is_err());
    }
    for t in [199., 451.] {
        assert!(m.evaluate(Temperature::new(t).unwrap(), 1.).is_err());
    }
    assert!(
        m.evaluate(Temperature::new(450.).unwrap(), 30000.).is_err(),
        "pressure limit"
    );
    let propane = model(&[("propane", 1.)]);
    assert!(
        propane
            .evaluate(Temperature::new(200.).unwrap(), 2000.)
            .is_err(),
        "unstable state"
    );
}

#[test]
fn analytic_fugacity_matches_teqp_ad_including_infinite_dilution() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/fixtures/gerg-fugacity-v1.json"
    ))
    .unwrap();
    assert_eq!(fixture["cases"].as_array().unwrap().len(), 45);
    let absolute = fixture["tolerance"][0].as_f64().unwrap();
    let relative = fixture["tolerance"][1].as_f64().unwrap();
    let mut maximum_error = 0.0_f64;
    for case in fixture["cases"].as_array().unwrap() {
        // Reference values use source JSON order, not serde_json's sorted map order.
        let original: Vec<String> = case["component_ids"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect();
        let entries: Vec<_> = original
            .iter()
            .map(|id| (id.as_str(), case["composition"][id].as_f64().unwrap()))
            .collect();
        let m = model(&entries);
        let r = m
            .fugacity(
                Temperature::new(case["temperature_k"].as_f64().unwrap()).unwrap(),
                case["molar_density"].as_f64().unwrap(),
            )
            .unwrap();
        assert_eq!(r.component_ids, original);
        for (i, actual) in r.ln_fugacity_coefficients.iter().enumerate() {
            maximum_error = maximum_error
                .max((*actual - case["ln_fugacity_coefficients"][i].as_f64().unwrap()).abs());
            near(
                &format!("{} {} ln(phi)", case["id"], original[i]),
                *actual,
                case["ln_fugacity_coefficients"][i].as_f64().unwrap(),
                absolute,
                relative,
            );
        }
        let weighted: f64 = entries
            .iter()
            .zip(&r.ln_fugacity_coefficients)
            .map(|((_, x), v)| x * v)
            .sum();
        let s = &r.state;
        near(
            "Euler residual Gibbs",
            weighted,
            s.derivatives.ar + s.compressibility - 1.0 - s.compressibility.ln(),
            2e-10,
            2e-10,
        );
    }
    println!("45 teqp fugacity states: maximum absolute ln(phi) error {maximum_error:e}");
}

#[test]
fn fugacity_composition_derivatives_match_total_helmholtz_differences() {
    use thermo_models::GERG_GAS_CONSTANT as R;
    let t = Temperature::new(300.).unwrap();
    let ids = ["methane", "ethane", "nitrogen"];
    let amounts = [0.85, 0.1, 0.05];
    let volume = 0.001;
    let entries: Vec<_> = ids.iter().zip(amounts).map(|(id, n)| (*id, n)).collect();
    let r = model(&entries).fugacity(t, 1000.).unwrap();
    let permuted = model(&[
        ("nitrogen", 0.05),
        ("methane", 0.85),
        ("water", 0.),
        ("ethane", 0.1),
    ])
    .fugacity(t, 1000.)
    .unwrap();
    for (i, id) in ids.iter().enumerate() {
        let j = permuted
            .component_ids
            .iter()
            .position(|name| name == id)
            .unwrap();
        near(
            "fugacity order and retained zeros",
            r.ln_fugacity_coefficients[i],
            permuted.ln_fugacity_coefficients[j],
            1e-11,
            1e-11,
        );
    }
    let energy = |n: [f64; 3]| {
        let total: f64 = n.iter().sum();
        let composition = Composition::new(
            CompositionBasis::Mole,
            ids.iter()
                .zip(n)
                .map(|(id, v)| (ComponentId::new(*id).unwrap(), v))
                .collect(),
            NormalizationPolicy::Normalize,
        )
        .unwrap();
        total
            * GergModel::new(&composition)
                .unwrap()
                .helmholtz(t, total / volume)
                .unwrap()
                .ar
    };
    for i in 0..3 {
        for step in [1e-5, 1e-6] {
            let mut plus = amounts;
            let mut minus = amounts;
            plus[i] += step;
            minus[i] -= step;
            let mu_res = (energy(plus) - energy(minus)) / (2. * step);
            near(
                "d(n ar)/dn at fixed T,V",
                mu_res,
                r.ln_fugacity_coefficients[i] + r.state.compressibility.ln(),
                2e-8,
                2e-8,
            );
        }
    }
    near(
        "pure ln(phi)",
        model(&[("methane", 1.)])
            .fugacity(t, 1000.)
            .unwrap()
            .ln_fugacity_coefficients[0],
        {
            let s = model(&[("methane", 1.)]).evaluate(t, 1000.).unwrap();
            s.gibbs_energy / (R * t.value()) - (s.derivatives.a0 + 1.0 + s.compressibility.ln())
        },
        1e-10,
        1e-10,
    );
}
