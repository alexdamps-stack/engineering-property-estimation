//! A4 stability, TP flash and pure saturation against the offline independent
//! fixture `validation/fixtures/flash-v1.json`, plus consistency and failure checks.
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, RootSelection};
use thermo_solvers::*;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase,
    Pressure, Temperature,
};

fn fixture() -> Value {
    serde_json::from_str(include_str!("../../../validation/fixtures/flash-v1.json")).unwrap()
}

fn variant(name: &str) -> CubicVariant {
    match name {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        other => panic!("unknown variant {other}"),
    }
}

fn model(ids: &[&str], x: &[f64], v: CubicVariant, kij: Vec<Vec<f64>>) -> CubicModel {
    let c = Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(x)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::Normalize,
    )
    .unwrap();
    CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &c,
        v,
        kij,
        "explicit synthetic kij for A4 verification",
    )
    .unwrap()
}

fn zeros(n: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; n]; n]
}

fn tp(t: f64, p: f64) -> (Temperature, Pressure) {
    (Temperature::new(t).unwrap(), Pressure::new(p).unwrap())
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn budget(f: &Value, key: &str) -> (f64, f64) {
    let t = &f["tolerances"][key];
    (t[0].as_f64().unwrap(), t[1].as_f64().unwrap())
}

fn close(label: &str, actual: f64, expected: f64, (atol, rtol): (f64, f64)) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= atol + rtol * expected.abs(),
        "{label}: {actual:e} != {expected:e}, delta={:e}",
        (actual - expected).abs()
    );
}

fn case_model(case: &Value) -> CubicModel {
    let ids: Vec<&str> = case["ids"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    let kij: Vec<Vec<f64>> = case["kij"].as_array().unwrap().iter().map(floats).collect();
    model(
        &ids,
        &floats(&case["z"]),
        variant(case["variant"].as_str().unwrap()),
        kij,
    )
}

#[test]
fn fixture_provenance_matches_embedded_data() {
    let f = fixture();
    assert_eq!(f["schema"], "a4-flash-reference-v1");
    assert_eq!(f["packages"]["teqp"], "0.23.2");
    assert_eq!(f["flash"].as_array().unwrap().len(), 22);
    assert_eq!(f["saturation"].as_array().unwrap().len(), 30);
    // The generator records the SHA-256 of the source extract it read; reviewers
    // compare it with `validation/a4-checks.md` when regenerating.
    assert_eq!(f["input_sha256"].as_str().unwrap().len(), 64);
}

#[test]
fn tp_flash_matches_independent_references() {
    let f = fixture();
    let (beta_tol, x_tol, z_tol) = (
        budget(&f, "vapor_fraction"),
        budget(&f, "composition"),
        budget(&f, "Z"),
    );
    for case in f["flash"].as_array().unwrap() {
        let id = format!("{} {}", case["variant"], case["id"]);
        let m = case_model(case);
        let (t, p) = tp(case["T"].as_f64().unwrap(), case["p"].as_f64().unwrap());
        let result =
            tp_flash(&m, t, p, FlashOptions::default()).unwrap_or_else(|e| panic!("{id}: {e}"));
        let expected_phases = case["phases"].as_u64().unwrap() as usize;
        assert_eq!(result.phases.len(), expected_phases, "{id}: phase count");
        assert!(
            result
                .diagnostics
                .iter()
                .any(|d| d.code == "A4-STABILITY-FINITE-TRIALS")
        );
        if expected_phases == 1 {
            assert!(result.feed_stability.stable_in_trials);
            close(
                &id,
                result.phases[0].state.compressibility,
                case["Z"].as_f64().unwrap(),
                z_tol,
            );
            continue;
        }
        let tm = result.feed_stability.most_negative().unwrap();
        assert!(tm.tangent_plane_distance < -1e-8, "{id}: feed not unstable");
        let [vapor, liquid] = result.phases.as_slice() else {
            unreachable!()
        };
        assert_eq!((vapor.phase, liquid.phase), (Phase::Vapor, Phase::Liquid));
        close(
            &format!("{id} beta"),
            vapor.phase_fraction,
            case["vapor_fraction"].as_f64().unwrap(),
            beta_tol,
        );
        for (label, phase, key) in [("y", vapor, "vapor"), ("x", liquid, "liquid")] {
            for (a, e) in phase.composition.iter().zip(floats(&case[key])) {
                close(&format!("{id} {label}"), *a, e, x_tol);
            }
        }
        close(
            &format!("{id} Zv"),
            vapor.state.compressibility,
            case["Z_vapor"].as_f64().unwrap(),
            z_tol,
        );
        close(
            &format!("{id} Zl"),
            liquid.state.compressibility,
            case["Z_liquid"].as_f64().unwrap(),
            z_tol,
        );
        assert!(result.max_ln_fugacity_residual <= 1e-10, "{id}");
        assert!(result.max_component_balance_residual <= 1e-12, "{id}");
        assert!(result.reduced_gibbs_change < 0.0, "{id}");
        assert_eq!(result.phase_stability.len(), 2);
        assert!(result.phase_stability.iter().all(|s| s.stable_in_trials));
        assert_eq!(result.vapor_fraction(), Some(vapor.phase_fraction));
    }
}

#[test]
fn newton_acceleration_reduces_iterations_near_the_critical_region() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        for (fractions, pressures) in [
            ([0.75, 0.25], &[8.5e6, 8.9e6, 9.025e6][..]),
            (
                [0.8, 0.2],
                &[9.0e6, 9.1e6, 9.15e6, 9.2e6, 9.25e6, 9.3e6][..],
            ),
        ] {
            let m = model(&["methane", "propane"], &fractions, variant, zeros(2));
            for &pressure in pressures {
                let (t, p) = tp(250.0, pressure);
                let accelerated = tp_flash(&m, t, p, FlashOptions::default()).unwrap();
                let ss = tp_flash(
                    &m,
                    t,
                    p,
                    FlashOptions {
                        newton: false,
                        ..FlashOptions::default()
                    },
                )
                .unwrap();
                eprintln!(
                    "{variant:?} p={pressure}: Newton {}, eigenvalue SS {}",
                    accelerated.iterations, ss.iterations
                );
                assert_eq!(accelerated.phases.len(), 2);
                assert!(accelerated.iterations * 2 < ss.iterations);
                assert!(
                    accelerated
                        .diagnostics
                        .iter()
                        .any(|d| d.code == "A4-NEWTON-STEPS"
                            && !d.message.starts_with("Accepted 0 "))
                );
                close(
                    "beta",
                    accelerated.phases[0].phase_fraction,
                    ss.phases[0].phase_fraction,
                    (2e-8, 0.0),
                );
                for (a, b) in accelerated.phases.iter().zip(&ss.phases) {
                    for (x, y) in a.composition.iter().zip(&b.composition) {
                        close("phase composition", *x, *y, (2e-8, 1e-7));
                    }
                }
                assert!(accelerated.max_ln_fugacity_residual <= 1e-10);
                assert!(accelerated.max_component_balance_residual <= 1e-12);
            }
        }
    }
}

#[test]
fn split_phases_satisfy_fugacity_equality_from_property_evaluation() {
    let m = model(
        &["methane", "propane"],
        &[0.5, 0.5],
        CubicVariant::PengRobinsonBostonMathias,
        zeros(2),
    );
    let (t, p) = tp(250.0, 3e6);
    let r = tp_flash(&m, t, p, FlashOptions::default()).unwrap();
    let [v, l] = r.phases.as_slice() else {
        panic!()
    };
    for i in 0..2 {
        let lhs = v.composition[i].ln() + v.state.ln_fugacity_coefficients[i];
        let rhs = l.composition[i].ln() + l.state.ln_fugacity_coefficients[i];
        assert!((lhs - rhs).abs() <= 1e-10);
        let balance = v.phase_fraction * v.composition[i] + l.phase_fraction * l.composition[i];
        assert!((balance - m.fractions()[i]).abs() <= 1e-13);
    }
    assert!((v.phase_fraction + l.phase_fraction - 1.0).abs() <= 1e-15);
    assert!(v.state.molar_volume > l.state.molar_volume);
    // Phase enthalpies are reported with the shared caloric reference.
    assert_eq!(v.state.caloric_reference, l.state.caloric_reference);
}

#[test]
fn permutation_and_inactive_components_do_not_change_the_split() {
    let (t, p) = tp(250.0, 3e6);
    let base = tp_flash(
        &model(
            &["methane", "propane"],
            &[0.5, 0.5],
            CubicVariant::PengRobinson,
            zeros(2),
        ),
        t,
        p,
        FlashOptions::default(),
    )
    .unwrap();
    let swapped = tp_flash(
        &model(
            &["propane", "methane"],
            &[0.5, 0.5],
            CubicVariant::PengRobinson,
            zeros(2),
        ),
        t,
        p,
        FlashOptions::default(),
    )
    .unwrap();
    let trace = tp_flash(
        &model(
            &["methane", "nitrogen", "propane"],
            &[0.5, 0.0, 0.5],
            CubicVariant::PengRobinson,
            zeros(3),
        ),
        t,
        p,
        FlashOptions::default(),
    )
    .unwrap();
    let beta = base.phases[0].phase_fraction;
    assert!((swapped.phases[0].phase_fraction - beta).abs() <= 1e-10);
    assert!((trace.phases[0].phase_fraction - beta).abs() <= 1e-10);
    assert!((swapped.phases[0].composition[1] - base.phases[0].composition[0]).abs() <= 1e-10);
    for phase in &trace.phases {
        assert_eq!(
            phase.composition[1], 0.0,
            "inactive component must stay absent"
        );
    }
}

#[test]
fn phase_appearance_and_disappearance_along_an_isotherm() {
    // Methane/propane 50/50 at 250 K: vapor at low pressure, two phases, liquid at
    // high pressure. Vapor fraction must stay in [0, 1] and fall monotonically.
    let m = model(
        &["methane", "propane"],
        &[0.5, 0.5],
        CubicVariant::PengRobinson,
        zeros(2),
    );
    let mut previous = f64::INFINITY;
    let mut counts = Vec::new();
    for step in 0..=40 {
        let pressure = 2e5 + step as f64 * 2e5;
        let (t, p) = tp(250.0, pressure);
        let r = tp_flash(&m, t, p, FlashOptions::default())
            .unwrap_or_else(|e| panic!("p={pressure}: {e}"));
        let beta = r.vapor_fraction().expect("labelled vapor or liquid");
        assert!((0.0..=1.0).contains(&beta));
        assert!(
            beta <= previous + 1e-12,
            "vapor fraction increased at p={pressure}"
        );
        previous = beta;
        if counts.last() != Some(&r.phases.len()) {
            counts.push(r.phases.len());
        }
    }
    assert_eq!(
        counts,
        vec![1, 2, 1],
        "expected dew-point appearance and bubble-point disappearance"
    );
}

#[test]
fn single_phase_labels_and_pure_components() {
    let (t, p) = tp(330.0, 5e6);
    let gas = model(
        &["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
        &[0.88, 0.06, 0.03, 0.02, 0.01],
        CubicVariant::PengRobinsonBostonMathias,
        zeros(5),
    );
    let r = tp_flash(&gas, t, p, FlashOptions::default()).unwrap();
    assert_eq!(r.phases.len(), 1);
    assert_eq!(r.phases[0].phase, Phase::Supercritical);
    assert_eq!(r.vapor_fraction(), None);

    let (t, p) = tp(250.0, 5e6);
    let liquid = model(
        &["methane", "propane"],
        &[0.05, 0.95],
        CubicVariant::PengRobinson,
        zeros(2),
    );
    let r = tp_flash(&liquid, t, p, FlashOptions::default()).unwrap();
    assert_eq!(r.phases[0].phase, Phase::Liquid);

    let (t, p) = tp(300.0, 1e5);
    let pure = model(&["propane"], &[1.0], CubicVariant::PengRobinson, zeros(1));
    let r = tp_flash(&pure, t, p, FlashOptions::default()).unwrap();
    assert_eq!(r.phases.len(), 1);
    assert_eq!(r.phases[0].phase, Phase::Vapor);
    assert!(r.feed_stability.trials.is_empty());
}

#[test]
fn saturation_matches_superancillary_references() {
    let f = fixture();
    let (p_tol, z_tol) = (
        budget(&f, "saturation_pressure"),
        budget(&f, "saturation_Z"),
    );
    for case in f["saturation"].as_array().unwrap() {
        let id = format!(
            "{} {} Tr={}",
            case["variant"], case["component"], case["Tr"]
        );
        let m = model(
            &[case["component"].as_str().unwrap()],
            &[1.0],
            variant(case["variant"].as_str().unwrap()),
            zeros(1),
        );
        let t = Temperature::new(case["T"].as_f64().unwrap()).unwrap();
        let s = pure_saturation_pressure(&m, t, SATURATION_MAX_ITERATIONS)
            .unwrap_or_else(|e| panic!("{id}: {e}"));
        close(
            &format!("{id} p"),
            s.pressure_pa,
            case["p"].as_f64().unwrap(),
            p_tol,
        );
        close(
            &format!("{id} Zl"),
            s.liquid_compressibility,
            case["Z_liquid"].as_f64().unwrap(),
            z_tol,
        );
        close(
            &format!("{id} Zv"),
            s.vapor_compressibility,
            case["Z_vapor"].as_f64().unwrap(),
            z_tol,
        );
        assert!(s.fugacity_residual <= 1e-12, "{id}");
    }
}

#[test]
fn explicit_failures_and_iteration_limits() {
    let m = model(
        &["methane", "propane"],
        &[0.5, 0.5],
        CubicVariant::PengRobinson,
        zeros(2),
    );
    let (t, p) = tp(250.0, 3e6);
    let limited = FlashOptions {
        max_iterations: 1,
        ..FlashOptions::default()
    };
    assert!(matches!(
        tp_flash(&m, t, p, limited),
        Err(CalculationError::NonConvergence(_))
    ));
    let (t_stable, p_stable) = tp(330.0, 2e6);
    let no_stability_iterations = FlashOptions {
        stability: StabilityOptions {
            max_iterations: 1,
            ..StabilityOptions::default()
        },
        ..FlashOptions::default()
    };
    assert!(matches!(
        tp_flash(&m, t_stable, p_stable, no_stability_iterations),
        Err(CalculationError::NonConvergence(_))
    ));
    let (cold, _) = tp(150.0, 3e6);
    assert!(matches!(
        tp_flash(&m, cold, p, FlashOptions::default()),
        Err(CalculationError::OutsideValidityDomain(_))
    ));

    let methane = model(&["methane"], &[1.0], CubicVariant::PengRobinson, zeros(1));
    assert!(matches!(
        pure_saturation_pressure(&methane, Temperature::new(180.0).unwrap(), 100),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    let propane = model(&["propane"], &[1.0], CubicVariant::PengRobinson, zeros(1));
    assert!(matches!(
        pure_saturation_pressure(&propane, Temperature::new(360.0).unwrap(), 100),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        pure_saturation_pressure(&propane, Temperature::new(300.0).unwrap(), 0),
        Err(CalculationError::NonConvergence(_))
    ));
    assert!(matches!(
        pure_saturation_pressure(&m, Temperature::new(250.0).unwrap(), 100),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    assert!(m.with_fractions(&[0.6, 0.6]).is_err());
    assert!(m.with_fractions(&[1.0]).is_err());
}

#[test]
fn stability_trials_report_evidence_not_proof() {
    let m = model(
        &["methane", "propane"],
        &[0.5, 0.5],
        CubicVariant::PengRobinson,
        zeros(2),
    );
    let (t, p) = tp(330.0, 2e6);
    let s = assess_stability(&m, t, p, StabilityOptions::default()).unwrap();
    assert!(s.stable_in_trials);
    assert_eq!(
        s.trials.len(),
        4,
        "two Wilson trials plus one per active component"
    );
    assert!(s.trials.iter().all(|t| t.tangent_plane_distance >= -1e-8));
    let fug = m.ln_fugacity(t, p, RootSelection::MinimumGibbs).unwrap();
    assert!((s.feed_compressibility - fug.compressibility).abs() == 0.0);
    let no_pure = StabilityOptions {
        pure_component_trials: false,
        ..StabilityOptions::default()
    };
    assert_eq!(assess_stability(&m, t, p, no_pure).unwrap().trials.len(), 2);
}

#[test]
fn negative_tpd_survives_another_slow_trial() {
    for ids in [["methane", "propane"], ["propane", "methane"]] {
        let m = model(&ids, &[0.5, 0.5], CubicVariant::PengRobinson, zeros(2));
        for pressure in [7.0795e6, 7.08e6, 7.1e6] {
            let (t, p) = tp(295.0, pressure);
            let stability = assess_stability(&m, t, p, StabilityOptions::default()).unwrap();
            assert!(!stability.stable_in_trials);
            assert!(stability.most_negative().unwrap().tangent_plane_distance < -1e-8);
            let split = tp_flash(&m, t, p, FlashOptions::default()).unwrap();
            assert_eq!(split.phases.len(), 2);
            assert!(split.max_ln_fugacity_residual <= 1e-10);
            assert!(split.max_component_balance_residual < 1e-12);
            assert!(split.reduced_gibbs_change < 0.0);
            assert!(split.phase_stability.iter().all(|s| s.stable_in_trials));
        }
    }
}

#[test]
fn invalid_options_are_rejected_before_stable_or_pure_shortcuts() {
    let binary = model(
        &["methane", "propane"],
        &[0.5, 0.5],
        CubicVariant::PengRobinson,
        zeros(2),
    );
    let pure = model(&["propane"], &[1.0], CubicVariant::PengRobinson, zeros(1));
    for m in [&binary, &pure] {
        let (t, p) = tp(250.0, 3e6);
        let mut invalid = vec![StabilityOptions {
            max_iterations: 0,
            ..Default::default()
        }];
        for tolerance in [f64::NAN, f64::INFINITY, 0.0, -1e-10, 1e-6] {
            invalid.push(StabilityOptions {
                tolerance,
                ..Default::default()
            });
        }
        for tm_threshold in [f64::NAN, f64::NEG_INFINITY, 0.0, 1e-8, -1e-5] {
            invalid.push(StabilityOptions {
                tm_threshold,
                ..Default::default()
            });
        }
        for stability in invalid {
            assert!(matches!(
                assess_stability(m, t, p, stability),
                Err(CalculationError::UnsupportedCapability(_))
            ));
            let options = FlashOptions {
                stability,
                ..Default::default()
            };
            assert!(matches!(
                tp_flash(m, t, p, options),
                Err(CalculationError::UnsupportedCapability(_))
            ));
            assert!(matches!(
                ph_flash(
                    m,
                    p,
                    0.0,
                    PhFlashOptions {
                        tp: options,
                        ..Default::default()
                    }
                ),
                Err(CalculationError::UnsupportedCapability(_))
            ));
        }
        let mut invalid = vec![FlashOptions {
            max_iterations: 0,
            ..Default::default()
        }];
        for tolerance in [f64::NAN, f64::INFINITY, 0.0, -1e-10, 1e-6] {
            invalid.push(FlashOptions {
                tolerance,
                ..Default::default()
            });
        }
        for options in invalid {
            assert!(matches!(
                tp_flash(m, t, p, options),
                Err(CalculationError::UnsupportedCapability(_))
            ));
        }
    }
}
