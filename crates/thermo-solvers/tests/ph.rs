//! A5 recovery, coexistence, energy balance and explicit failure contracts.
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, RootSelection};
use thermo_solvers::*;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure,
    Temperature,
};

fn model(ids: &[&str], fractions: &[f64], variant: CubicVariant) -> CubicModel {
    let composition = Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(fractions)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition,
        variant,
        vec![vec![0.0; ids.len()]; ids.len()],
        "explicit zero kij approximation for A5 tests",
    )
    .unwrap()
}
fn t(t: f64) -> Temperature {
    Temperature::new(t).unwrap()
}
fn p(p: f64) -> Pressure {
    Pressure::new(p).unwrap()
}
fn variants() -> [CubicVariant; 2] {
    [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ]
}
fn close(actual: f64, expected: f64, atol: f64) {
    assert!(
        (actual - expected).abs() <= atol,
        "{actual} != {expected}, delta {}",
        actual - expected
    );
}

#[test]
fn unique_tp_ph_round_trips_including_zero_enthalpy() {
    for v in variants() {
        for (ids, x, temperature, pressure) in [
            (vec!["methane", "propane"], vec![0.5, 0.5], 250.0, 3e6),
            (vec!["methane", "propane"], vec![0.5, 0.5], 330.0, 2e6),
            (vec!["methane", "propane"], vec![0.05, 0.95], 250.0, 5e6),
            (vec!["methane"], vec![1.0], 400.0, 1e6),
            (vec!["methane"], vec![1.0], 400.0, 1e7),
            (vec!["water"], vec![1.0], 400.0, 1e6),
            (vec!["propane"], vec![1.0], 200.0, 1e5),
            (vec!["propane"], vec![1.0], 600.0, 1e5),
        ] {
            let m = model(&ids, &x, v);
            let tp = tp_flash(&m, t(temperature), p(pressure), FlashOptions::default()).unwrap();
            let ph = ph_flash(&m, p(pressure), tp.enthalpy(), PhFlashOptions::default()).unwrap();
            close(ph.temperature_k, temperature, 2e-7);
            assert!(ph.scaled_energy_residual.abs() <= 1e-10);
            assert_eq!(ph.phases.len(), tp.phases.len());
            for (a, b) in ph.phases.iter().zip(&tp.phases) {
                close(a.phase_fraction, b.phase_fraction, 1e-8);
                for (a, b) in a.composition.iter().zip(&b.composition) {
                    close(*a, *b, 1e-8);
                }
            }
        }
        let m = model(&["methane"], &[1.0], v);
        let zero = ph_flash(&m, p(1e5), 0.0, PhFlashOptions::default()).unwrap();
        assert!(
            zero.enthalpy.abs() <= 1e-10,
            "absolute energy floor near zero"
        );
    }
}

#[test]
fn pure_saturation_lever_rule_and_bubble_dew_endpoints() {
    for v in variants() {
        // Retain a zero component to test the pure active-component limit.
        let m = model(&["propane", "nitrogen"], &[1.0, 0.0], v);
        let sat = pure_saturation_pressure(&m, t(280.0), 200).unwrap();
        let pressure = p(sat.pressure_pa);
        let l = m
            .evaluate(t(280.0), pressure, RootSelection::SmallestVolume)
            .unwrap();
        let g = m
            .evaluate(t(280.0), pressure, RootSelection::LargestVolume)
            .unwrap();
        for beta in [0.0, 0.01, 0.5, 0.99, 1.0] {
            let h = (1.0 - beta) * l.enthalpy + beta * g.enthalpy;
            for (minimum_temperature_k, maximum_temperature_k) in
                [(200.0, 600.0), (280.0, 300.0), (250.0, 280.0)]
            {
                let options = PhFlashOptions {
                    minimum_temperature_k,
                    maximum_temperature_k,
                    ..Default::default()
                };
                let r = ph_flash(&m, pressure, h, options).unwrap();
                close(r.temperature_k, 280.0, 2e-8);
                assert!(r.pure_saturation.as_ref().unwrap().fugacity_residual <= 1e-12);
                close(r.phases[0].phase_fraction, beta, 2e-10);
                assert!(r.scaled_energy_residual.abs() <= 1e-10);
                assert_eq!(r.phases[0].composition, vec![1.0, 0.0]);
            }
        }
    }
}

#[test]
fn narrow_boiling_mixture_and_component_permutation() {
    for v in variants() {
        // 0.1% ethane in propane: finite but narrow coexistence temperature range.
        let m = model(&["ethane", "propane"], &[0.001, 0.999], v);
        let pure = model(&["propane"], &[1.0], v);
        let sat = pure_saturation_pressure(&pure, t(280.0), 200).unwrap();
        let pressure = p(sat.pressure_pa);
        // Search deterministically for an interior TP state, then require PH to
        // recover its phase split and energy without mistaking a small step for convergence.
        let tp = (0..1000)
            .map(|i| 280.0 - i as f64 * 0.0001)
            .find_map(|temperature| {
                let r = tp_flash(&m, t(temperature), pressure, FlashOptions::default()).unwrap();
                (r.phases.len() == 2
                    && r.phases[0].phase_fraction > 0.05
                    && r.phases[0].phase_fraction < 0.95)
                    .then_some(r)
            })
            .expect("narrow-boiling interior state");
        let ph = ph_flash(&m, pressure, tp.enthalpy(), PhFlashOptions::default()).unwrap();
        close(ph.temperature_k, tp.temperature_k, 2e-7);
        close(
            ph.phases[0].phase_fraction,
            tp.phases[0].phase_fraction,
            2e-7,
        );
        assert!(ph.scaled_energy_residual.abs() <= 1e-10);
        let permuted = model(&["propane", "ethane"], &[0.999, 0.001], v);
        let other = ph_flash(
            &permuted,
            pressure,
            tp.enthalpy(),
            PhFlashOptions::default(),
        )
        .unwrap();
        close(other.temperature_k, ph.temperature_k, 2e-8);
        close(
            other.phases[0].phase_fraction,
            ph.phases[0].phase_fraction,
            2e-7,
        );
    }
}

#[test]
fn isenthalpic_valve_preserves_feed_enthalpy() {
    for v in variants() {
        let m = model(
            &["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
            &[0.70, 0.15, 0.13, 0.01, 0.01],
            v,
        );
        let inlet = tp_flash(&m, t(240.0), p(5e6), FlashOptions::default()).unwrap();
        let outlet = ph_flash(&m, p(1.5e6), inlet.enthalpy(), PhFlashOptions::default()).unwrap();
        assert!(outlet.temperature_k < inlet.temperature_k);
        assert_eq!(outlet.phases.len(), 2);
        assert!(outlet.scaled_energy_residual.abs() <= 1e-10);
        let tp = outlet.tp.unwrap();
        assert!(tp.max_component_balance_residual <= 1e-8);
        assert!(tp.max_ln_fugacity_residual <= 1e-7);
    }
}

#[test]
fn failures_are_explicit_and_iterations_are_bounded() {
    let m = model(
        &["methane", "propane"],
        &[0.5, 0.5],
        CubicVariant::PengRobinson,
    );
    for h in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            ph_flash(&m, p(1e6), h, PhFlashOptions::default()),
            Err(CalculationError::InvalidInput(_))
        ));
    }
    for h in [-1e9, 1e9] {
        assert!(matches!(
            ph_flash(&m, p(1e6), h, PhFlashOptions::default()),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
    }
    for options in [
        PhFlashOptions {
            minimum_temperature_k: 150.0,
            ..Default::default()
        },
        PhFlashOptions {
            maximum_temperature_k: f64::NAN,
            ..Default::default()
        },
        PhFlashOptions {
            energy_tolerance: 0.0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            ph_flash(&m, p(1e6), 0.0, options),
            Err(CalculationError::UnsupportedCapability(_))
        ));
    }
    for max_iterations in [0, 1] {
        assert!(matches!(
            ph_flash(
                &m,
                p(1e6),
                0.0,
                PhFlashOptions {
                    max_iterations,
                    ..Default::default()
                }
            ),
            Err(CalculationError::NonConvergence(_))
        ));
    }
    let options = PhFlashOptions {
        tp: FlashOptions {
            stability: StabilityOptions {
                max_iterations: 0,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(
        matches!(
            ph_flash(&m, p(1e6), 0.0, options),
            Err(CalculationError::NonConvergence(_))
        ),
        "inner stability failure must propagate"
    );
    let pure = model(&["propane"], &[1.0], CubicVariant::PengRobinson);
    let sat = pure_saturation_pressure(&pure, t(280.0), 200).unwrap();
    assert!(matches!(
        ph_flash(
            &pure,
            p(sat.pressure_pa),
            -5000.0,
            PhFlashOptions {
                max_iterations: 1,
                ..Default::default()
            }
        ),
        Err(CalculationError::NonConvergence(_))
    ));
    assert!(matches!(
        ph_flash(&m, p(1e9), 0.0, PhFlashOptions::default()),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}

#[test]
fn ph_matches_independent_teqp_cantera_references() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/ph-v1.json")).unwrap();
    assert_eq!(fixture["schema"], "a5-ph-reference-v1");
    assert_eq!(fixture["packages"]["teqp"], "0.23.2");
    assert_eq!(fixture["packages"]["cantera"], "3.2.0");
    assert_eq!(fixture["cases"].as_array().unwrap().len(), 28);
    let budget = |key: &str, actual: f64, expected: f64| {
        let b = &fixture["tolerances"][key];
        close(
            actual,
            expected,
            b[0].as_f64().unwrap() + b[1].as_f64().unwrap() * expected.abs(),
        );
    };
    let mut max_t: f64 = 0.0;
    let mut max_beta: f64 = 0.0;
    let mut max_h: f64 = 0.0;
    for case in fixture["cases"].as_array().unwrap() {
        let ids: Vec<_> = case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .map(|id| id.as_str().unwrap())
            .collect();
        let z: Vec<_> = case["z"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_f64().unwrap())
            .collect();
        let variant = if case["variant"] == "PR" {
            CubicVariant::PengRobinson
        } else {
            CubicVariant::PengRobinsonBostonMathias
        };
        let m = model(&ids, &z, variant);
        let pressure = p(case["p"].as_f64().unwrap());
        let target = case["target_h"].as_f64().unwrap();
        let r = ph_flash(&m, pressure, target, PhFlashOptions::default())
            .unwrap_or_else(|e| panic!("{} {}: {e}", case["variant"], case["id"]));
        let expected_t = case["T"].as_f64().unwrap();
        budget("temperature_K", r.temperature_k, expected_t);
        budget("enthalpy_J_per_mol", r.enthalpy, target);
        assert_eq!(
            r.phases.len(),
            case["phases"].as_u64().unwrap() as usize,
            "{}",
            case["id"]
        );
        assert!(r.scaled_energy_residual.abs() <= 1e-10);
        max_t = max_t.max((r.temperature_k - expected_t).abs());
        max_h = max_h.max((r.enthalpy - target).abs());
        if r.phases.len() == 2 {
            let expected_beta = case["vapor_fraction"].as_f64().unwrap();
            budget("phase_fraction", r.phases[0].phase_fraction, expected_beta);
            max_beta = max_beta.max((r.phases[0].phase_fraction - expected_beta).abs());
            for (phase, key, hkey) in [
                (&r.phases[0], "vapor", "vapor_h"),
                (&r.phases[1], "liquid", "liquid_h"),
            ] {
                budget(
                    "enthalpy_J_per_mol",
                    phase.state.enthalpy,
                    case[hkey].as_f64().unwrap(),
                );
                for (i, fraction) in phase.composition.iter().enumerate() {
                    budget("composition", *fraction, case[key][i].as_f64().unwrap());
                }
            }
            for (i, feed) in z.iter().enumerate() {
                close(
                    r.phases
                        .iter()
                        .map(|ph| ph.phase_fraction * ph.composition[i])
                        .sum(),
                    *feed,
                    1e-8,
                );
            }
        }
        if let Some(tin) = case["inlet_T"].as_f64() {
            let inlet = tp_flash(
                &m,
                t(tin),
                p(case["inlet_p"].as_f64().unwrap()),
                FlashOptions::default(),
            )
            .unwrap();
            budget("enthalpy_J_per_mol", inlet.enthalpy(), target);
        }
    }
    println!(
        "28 independent PH states: max |dT|={max_t:e} K, |dh|={max_h:e} J/mol, |d beta|={max_beta:e}"
    );
}
