use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::{
    CubicSaturationKind as Kind, EnvelopeStop, SaturationOptions, saturation_pressure,
    saturation_temperature, trace_phase_envelope,
};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure, Temperature,
};
fn values(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}
fn inputs() -> Value {
    serde_json::from_str(include_str!("../../../validation/a13-inputs.json")).unwrap()
}
fn model(c: &Value) -> CubicModel {
    let ids = c["ids"].as_array().unwrap();
    let z = values(&c["z"]);
    let composition = Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(z)
            .map(|(id, x)| (ComponentId::new(id.as_str().unwrap()).unwrap(), x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition,
        if c["variant"] == "PR" {
            CubicVariant::PengRobinson
        } else {
            CubicVariant::PengRobinsonBostonMathias
        },
        c["kij"].as_array().unwrap().iter().map(values).collect(),
        "explicit A13 fixture interactions",
    )
    .unwrap()
}
fn close(f: &Value, key: &str, a: f64, b: f64, id: &str) {
    let limits = values(&f["tolerances"][key]);
    assert!(
        a.is_finite() && (a - b).abs() <= limits[0] + limits[1] * b.abs(),
        "{id} {key}: Rust={a:e} Clapeyron={b:e}"
    );
}
#[test]
fn independent_native_clapeyron_boundaries() {
    let text = std::env::var_os("THERMO_A13_REFERENCE")
        .map(|p| std::fs::read_to_string(p).unwrap())
        .unwrap_or_else(|| include_str!("../../../validation/fixtures/a13-v1.json").into());
    let f: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(f["schema"], "a13-clapeyron-v1");
    assert_eq!(f["julia_version"], "1.11.7");
    assert_eq!(f["clapeyron_version"], "0.6.30");
    let inputs = inputs();
    let cases = inputs["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 62);
    assert_eq!(f["cases"].as_array().unwrap().len(), 62);
    for (c, reference) in cases.iter().zip(f["cases"].as_array().unwrap()) {
        assert_eq!(c["id"], reference["id"]);
        let id = c["id"].as_str().unwrap();
        let m = model(c);
        let kind = if c["kind"] == "bubble" {
            Kind::Bubble
        } else {
            Kind::Dew
        };
        let fixed = c["fixed"].as_f64().unwrap();
        // At 225 K this mixture has more than one dew-pressure root in the
        // broad default interval. The independently tabulated lower-pressure
        // branch is explicitly bounded below 4 MPa; no implicit selection.
        let bounds = SaturationOptions {
            maximum_pressure_pa: if c["trace_group"] == "ng-dew" {
                4e6
            } else {
                1e7
            },
            ..Default::default()
        };
        let r = if c["unknown"] == "pressure" {
            saturation_pressure(&m, Temperature::new(fixed).unwrap(), kind, bounds)
        } else {
            saturation_temperature(&m, Pressure::new(fixed).unwrap(), kind, bounds)
        }
        .unwrap_or_else(|e| panic!("{id}: {e}"));
        close(
            &f,
            "temperature",
            r.temperature_k,
            reference["T"].as_f64().unwrap(),
            id,
        );
        close(
            &f,
            "pressure",
            r.pressure_pa,
            reference["p"].as_f64().unwrap(),
            id,
        );
        close(
            &f,
            "volume",
            r.liquid.molar_volume,
            reference["vl"].as_f64().unwrap(),
            id,
        );
        close(
            &f,
            "volume",
            r.vapor.molar_volume,
            reference["vv"].as_f64().unwrap(),
            id,
        );
        for (a, b) in r.liquid_composition.iter().zip(values(&reference["x"])) {
            close(&f, "composition", *a, b, id);
        }
        for (a, b) in r.vapor_composition.iter().zip(values(&reference["y"])) {
            close(&f, "composition", *a, b, id);
        }
        assert!(r.max_ln_fugacity_residual < 1e-9);
        println!("PASS {id}: T={} p={}", r.temperature_k, r.pressure_pa);
    }
}
#[test]
fn bounded_trace_retains_verified_points_and_stops_near_critical() {
    let c = &inputs()["cases"][12];
    assert_eq!(c["ids"][0], "propane");
    let m = model(c);
    let o = SaturationOptions {
        minimum_temperature_k: 300.,
        maximum_temperature_k: 360.,
        ..Default::default()
    };
    let trace = trace_phase_envelope(&m, Kind::Bubble, 5., 20, o).unwrap();
    assert_eq!(trace.points.len(), 11);
    assert_eq!(trace.supported_temperature_range_k, Some([300., 350.]));
    assert!(matches!(
        trace.stop,
        EnvelopeStop::NearCritical {
            temperature_k: 355.,
            ..
        }
    ));
    for p in trace.points {
        assert!(
            p.max_ln_fugacity_residual < 1e-9
                && p.vapor.molar_volume / p.liquid.molar_volume >= 1.5
        );
    }
    assert!(trace_phase_envelope(&m, Kind::Dew, 0., 10, o).is_err());
    assert!(trace_phase_envelope(&m, Kind::Dew, 1e-300, 10, o).is_err());
    assert!(trace_phase_envelope(&m, Kind::Dew, 5., 2, o).is_err());
    let bad = SaturationOptions {
        minimum_temperature_k: f64::NAN,
        ..o
    };
    assert!(saturation_pressure(&m, Temperature::new(300.).unwrap(), Kind::Bubble, bad).is_err());
    let bad = SaturationOptions {
        max_iterations: 0,
        ..o
    };
    assert!(saturation_temperature(&m, Pressure::new(1e6).unwrap(), Kind::Dew, bad).is_err());
    let empty = trace_phase_envelope(
        &m,
        Kind::Dew,
        1.,
        10,
        SaturationOptions {
            minimum_temperature_k: 355.,
            maximum_temperature_k: 360.,
            ..o
        },
    )
    .unwrap();
    assert!(empty.points.is_empty() && empty.supported_pressure_range_pa.is_none());
    assert!(matches!(empty.stop, EnvelopeStop::NearCritical { .. }));
}
#[test]
fn iteration_limit_bounds_and_inactive_components_are_explicit() {
    let c = &inputs()["cases"][0];
    let m = model(c);
    let limited = SaturationOptions {
        max_iterations: 1,
        ..Default::default()
    };
    assert!(
        saturation_pressure(&m, Temperature::new(250.).unwrap(), Kind::Bubble, limited).is_err()
    );
    assert!(
        saturation_pressure(
            &m,
            Temperature::new(199.).unwrap(),
            Kind::Dew,
            SaturationOptions::default()
        )
        .is_err()
    );
    assert!(
        saturation_temperature(
            &m,
            Pressure::new(1e9).unwrap(),
            Kind::Dew,
            SaturationOptions::default()
        )
        .is_err()
    );
    let m = m.with_fractions(&[0., 1.]).unwrap();
    let p = saturation_pressure(
        &m,
        Temperature::new(300.).unwrap(),
        Kind::Dew,
        SaturationOptions::default(),
    )
    .unwrap();
    assert_eq!(p.liquid_composition, vec![0., 1.]);
    assert_eq!(p.vapor_composition, vec![0., 1.]);
}

#[test]
fn sampled_mixture_envelopes_match_independent_native_points() {
    let text = std::env::var_os("THERMO_A13_REFERENCE")
        .map(|p| std::fs::read_to_string(p).unwrap())
        .unwrap_or_else(|| include_str!("../../../validation/fixtures/a13-v1.json").into());
    let f: Value = serde_json::from_str(&text).unwrap();
    let input = inputs();
    let all = input["cases"].as_array().unwrap();
    for (group, kind, lo, hi, step) in [
        ("binary-bubble", Kind::Bubble, 250., 290., 5.),
        ("propane-bubble", Kind::Bubble, 300., 350., 5.),
        ("ng-dew", Kind::Dew, 200., 225., 5.),
        ("ng-bubble", Kind::Bubble, 200., 206., 2.),
    ] {
        let cases: Vec<_> = all.iter().filter(|c| c["trace_group"] == group).collect();
        assert!(!cases.is_empty());
        let m = model(cases[0]);
        let o = SaturationOptions {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            maximum_pressure_pa: if group == "ng-dew" { 4e6 } else { 1e7 },
            ..Default::default()
        };
        let trace = trace_phase_envelope(&m, kind, step, 100, o).unwrap();
        assert!(
            matches!(trace.stop, EnvelopeStop::RequestedRangeComplete),
            "{group}: {:?}",
            trace.stop
        );
        assert_eq!(trace.points.len(), cases.len());
        for (p, c) in trace.points.iter().zip(cases) {
            let reference = f["cases"]
                .as_array()
                .unwrap()
                .iter()
                .find(|r| r["id"] == c["id"])
                .unwrap();
            close(
                &f,
                "pressure",
                p.pressure_pa,
                reference["p"].as_f64().unwrap(),
                group,
            );
        }
        if group == "ng-dew" {
            assert!(matches!(
                saturation_pressure(
                    &m,
                    Temperature::new(225.).unwrap(),
                    Kind::Dew,
                    SaturationOptions::default()
                ),
                Err(thermo_types::CalculationError::UnsupportedCapability(_))
            ));
        }
        if group == "ng-bubble" {
            let near = trace_phase_envelope(
                &m,
                kind,
                step,
                100,
                SaturationOptions {
                    maximum_temperature_k: 220.,
                    ..o
                },
            )
            .unwrap();
            assert!(
                matches!(near.stop, EnvelopeStop::NearCritical { .. }),
                "{:?}",
                near.stop
            );
            assert!(near.points.last().unwrap().temperature_k < 210.);
        }
    }
}

#[test]
fn component_order_and_zero_fraction_do_not_change_boundaries() {
    let mut c = inputs()["cases"][16].clone();
    let original = model(&c);
    let opts = SaturationOptions::default();
    let p = saturation_pressure(
        &original,
        Temperature::new(250.).unwrap(),
        Kind::Bubble,
        opts,
    )
    .unwrap();
    c["ids"] = serde_json::json!(["propane", "nitrogen", "methane"]);
    c["z"] = serde_json::json!([0.5, 0., 0.5]);
    c["kij"] = serde_json::json!([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]);
    let reordered = saturation_pressure(
        &model(&c),
        Temperature::new(250.).unwrap(),
        Kind::Bubble,
        opts,
    )
    .unwrap();
    assert!((p.pressure_pa - reordered.pressure_pa).abs() < 0.01);
    assert!((p.vapor_composition[0] - reordered.vapor_composition[2]).abs() < 1e-9);
    assert_eq!(reordered.vapor_composition[1], 0.);
    c["ids"] = serde_json::json!(["propane", "water", "methane"]);
    c["z"] = serde_json::json!([0.49, 0.01, 0.5]);
    assert!(matches!(
        saturation_pressure(
            &model(&c),
            Temperature::new(250.).unwrap(),
            Kind::Bubble,
            opts
        ),
        Err(thermo_types::CalculationError::UnsupportedCapability(_))
    ));
}

#[test]
fn near_unit_pure_fractions_are_canonicalized_in_both_saturation_paths() {
    for index in [12, 28, 0, 16] {
        let base = model(&inputs()["cases"][index]);
        let mut canonical = vec![0.; base.fractions().len()];
        let propane = canonical.len() - 1;
        canonical[propane] = 1.;
        let exact = base.with_fractions(&canonical).unwrap();
        for kind in [Kind::Bubble, Kind::Dew] {
            let expected_p = saturation_pressure(
                &exact,
                Temperature::new(300.).unwrap(),
                kind,
                SaturationOptions::default(),
            )
            .unwrap();
            let expected_t = saturation_temperature(
                &exact,
                Pressure::new(1e6).unwrap(),
                kind,
                SaturationOptions::default(),
            )
            .unwrap();
            for amount in [1. - 5e-11, 1. + 5e-11] {
                let mut fractions = canonical.clone();
                fractions[propane] = amount;
                let near = base.with_fractions(&fractions).unwrap();
                let p = saturation_pressure(
                    &near,
                    Temperature::new(300.).unwrap(),
                    kind,
                    SaturationOptions::default(),
                )
                .unwrap();
                let t = saturation_temperature(
                    &near,
                    Pressure::new(1e6).unwrap(),
                    kind,
                    SaturationOptions::default(),
                )
                .unwrap();
                assert_eq!(p.pressure_pa, expected_p.pressure_pa);
                assert_eq!(t.temperature_k, expected_t.temperature_k);
                for result in [p, t] {
                    assert_eq!(result.liquid_composition, canonical);
                    assert_eq!(result.vapor_composition, canonical);
                    assert!(result.max_ln_fugacity_residual < 1e-9);
                }
                let exhausted = SaturationOptions {
                    max_iterations: 1,
                    ..Default::default()
                };
                assert!(matches!(
                    saturation_temperature(&near, Pressure::new(1e6).unwrap(), kind, exhausted),
                    Err(thermo_types::CalculationError::NonConvergence(_))
                ));
            }
        }
    }
}
