use serde_json::Value;
use thermo_models::GergModel;
use thermo_solvers::{GergDensityOptions, gerg_density};
use thermo_types::{
    CalculationError as Error, ComponentId, Composition, CompositionBasis, NormalizationPolicy,
    Pressure, Temperature, ValidityStatus,
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

#[test]
fn density_and_properties_match_nist_and_teqp_roots() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/fixtures/gerg-density-v1.json"
    ))
    .unwrap();
    let cases = fixture["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 42);
    let mut maximum_density_error = 0.0_f64;
    let mut maximum_iterations = 0;
    for case in cases {
        let m = model(
            &case["composition"]
                .as_object()
                .unwrap()
                .iter()
                .map(|(id, x)| (id.as_str(), x.as_f64().unwrap()))
                .collect::<Vec<_>>(),
        );
        let lo = case["density_bracket"][0].as_f64().unwrap();
        let hi = case["density_bracket"][1].as_f64().unwrap();
        let p = Pressure::new(case["pressure_pa"].as_f64().unwrap()).unwrap();
        let result = gerg_density(
            &m,
            Temperature::new(case["temperature_k"].as_f64().unwrap()).unwrap(),
            p,
            GergDensityOptions::new(lo, hi),
        )
        .unwrap();
        let s = &result.state;
        maximum_density_error = maximum_density_error
            .max((s.molar_density - case["expected"]["molar_density"].as_f64().unwrap()).abs());
        maximum_iterations = maximum_iterations.max(result.iterations);
        assert_eq!(s.validity, ValidityStatus::Unassessed);
        assert!(result.scaled_pressure_residual.abs() <= 1e-10);
        assert!(result.iterations > 1 && result.iterations <= 100);
        assert!(
            s.molar_density >= result.final_density_bracket[0]
                && s.molar_density <= result.final_density_bracket[1]
        );
        assert_eq!(result.initial_density_bracket, [lo, hi]);
        for (key, actual) in [
            ("molar_density", s.molar_density),
            ("pressure_pa", s.pressure_pa),
            ("compressibility", s.compressibility),
            ("enthalpy", s.enthalpy),
            ("entropy", s.entropy),
            ("cp", s.cp),
            ("cv", s.cv),
            ("speed_of_sound", s.speed_of_sound),
            ("dp_drho", s.dp_drho),
            ("dp_dt", s.dp_dt),
        ] {
            let expected = case["expected"][key].as_f64().unwrap();
            let tol = &fixture["tolerances"][key];
            assert!(
                (actual - expected).abs()
                    <= tol[0].as_f64().unwrap() + tol[1].as_f64().unwrap() * expected.abs(),
                "{} {key}: actual={actual}, NIST={expected}",
                case["id"]
            );
        }
        let independent = case["teqp_molar_density"].as_f64().unwrap();
        assert!(
            (s.molar_density - independent).abs() <= 2e-7 + 2e-9 * independent,
            "{} teqp root",
            case["id"]
        );
    }
    println!(
        "42 independent GERG TP roots: maximum density error {maximum_density_error:e} mol/m³; maximum iterations {maximum_iterations}"
    );
}

#[test]
fn explicit_brackets_select_distinct_homogeneous_roots() {
    let m = model(&[("propane", 1.)]);
    let t = Temperature::new(250.).unwrap();
    let p = Pressure::new(100000.).unwrap();
    let low = gerg_density(&m, t, p, GergDensityOptions::new(1., 100.)).unwrap();
    let high = gerg_density(&m, t, p, GergDensityOptions::new(12000., 13000.)).unwrap();
    assert!(low.state.molar_density < 100. && high.state.molar_density > 12000.);
    assert!(
        low.scaled_pressure_residual.abs() <= 1e-10 && high.scaled_pressure_residual.abs() <= 1e-10
    );
    assert_eq!(low.state.validity, ValidityStatus::Unassessed);
    assert_eq!(high.state.validity, ValidityStatus::Unassessed);
}

#[test]
fn endpoints_broad_brackets_and_component_invariance() {
    let m = model(&[("methane", 0.9), ("ethane", 0.1)]);
    let t = Temperature::new(300.).unwrap();
    let original = m.evaluate(t, 1000.).unwrap();
    let p = Pressure::new(original.pressure_pa).unwrap();
    for bracket in [[1000., 2000.], [500., 1000.]] {
        let r = gerg_density(&m, t, p, GergDensityOptions::new(bracket[0], bracket[1])).unwrap();
        assert_eq!(r.iterations, 0);
        assert_eq!(r.state.molar_density, 1000.);
    }
    let options = GergDensityOptions::new(0.001, 20000.);
    let r = gerg_density(&m, t, p, options).unwrap();
    assert!((r.state.molar_density - 1000.).abs() < 1e-6);
    let permuted = model(&[("ethane", 0.1), ("water", 0.), ("methane", 0.9)]);
    let other = gerg_density(&permuted, t, p, options).unwrap();
    assert!((r.state.molar_density - other.state.molar_density).abs() < 1e-6);
    assert!((r.state.enthalpy - other.state.enthalpy).abs() < 1e-5);
}

#[test]
fn options_unbracketed_states_instability_and_exhaustion_are_explicit() {
    let m = model(&[("methane", 1.)]);
    let t = Temperature::new(300.).unwrap();
    let p = Pressure::new(1e6).unwrap();
    let base = GergDensityOptions::new(1., 1000.);
    for options in [
        GergDensityOptions::new(f64::NAN, 1000.),
        GergDensityOptions::new(1., f64::INFINITY),
        GergDensityOptions::new(1000., 1.),
        GergDensityOptions::new(1., 1.),
        GergDensityOptions::new(0., 1000.),
        GergDensityOptions::new(1., 30001.),
        GergDensityOptions {
            max_iterations: 0,
            ..base
        },
        GergDensityOptions {
            pressure_tolerance: f64::NAN,
            ..base
        },
        GergDensityOptions {
            pressure_tolerance: 0.,
            ..base
        },
        GergDensityOptions {
            pressure_tolerance: 1e-6,
            ..base
        },
    ] {
        assert!(matches!(
            gerg_density(&m, t, p, options),
            Err(Error::UnsupportedCapability(_))
        ));
    }
    assert!(matches!(
        gerg_density(&m, t, p, GergDensityOptions::new(1., 2.)),
        Err(Error::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        gerg_density(&m, t, Pressure::new(36e6).unwrap(), base),
        Err(Error::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        gerg_density(&m, Temperature::new(451.).unwrap(), p, base),
        Err(Error::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        gerg_density(
            &m,
            t,
            p,
            GergDensityOptions {
                max_iterations: 1,
                ..base
            }
        ),
        Err(Error::NonConvergence(_))
    ));
    let propane = model(&[("propane", 1.)]);
    assert!(matches!(
        gerg_density(
            &propane,
            Temperature::new(250.).unwrap(),
            p,
            GergDensityOptions::new(1200., 2000.)
        ),
        Err(Error::OutsideValidityDomain(_))
    ));
}
