use serde_json::Value;
use thermo_models::GergModel;
use thermo_solvers::{FlashOptions, gerg_tp_flash};
use thermo_solvers::{GergPhFlashOptions, gerg_ph_flash};
use thermo_solvers::{StabilityOptions, gerg_assess_stability, gerg_outer_roots};
use thermo_types::CalculationError as Error;
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure, Temperature,
};

fn model(ids: &[&str], x: &[f64]) -> GergModel {
    GergModel::new(
        &Composition::new(
            CompositionBasis::Mole,
            ids.iter()
                .zip(x)
                .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
                .collect(),
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn regular_binary_and_natural_gas_flash() {
    for (ids, x, t, p, phases) in [
        (vec!["methane", "propane"], vec![0.5, 0.5], 250., 3e6, 2),
        (vec!["methane", "propane"], vec![0.5, 0.5], 330., 2e6, 1),
        (vec!["methane", "propane"], vec![0.05, 0.95], 250., 5e6, 1),
        (vec!["methane", "propane"], vec![0.5, 0.5], 290., 5e6, 2),
        (
            vec!["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
            vec![0.70, 0.15, 0.13, 0.01, 0.01],
            250.,
            3e6,
            2,
        ),
    ] {
        let m = model(&ids, &x);
        let r = gerg_tp_flash(
            &m,
            Temperature::new(t).unwrap(),
            Pressure::new(p).unwrap(),
            FlashOptions::default(),
        )
        .unwrap();
        println!(
            "{ids:?} {t} {p} phases={} beta={:?}",
            r.phases.len(),
            r.vapor_fraction()
        );
        assert_eq!(r.phases.len(), phases);
        assert!(r.max_ln_fugacity_residual <= 1e-10);
        assert!(r.max_component_balance_residual <= 1e-8);
    }
}

fn fixture() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/fixtures/gerg-equilibrium-v1.json"
    ))
    .unwrap()
}
fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}
fn case_model(c: &Value) -> GergModel {
    model(
        &c["ids"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>(),
        &numbers(&c["z"]),
    )
}
fn close(f: &Value, key: &str, actual: f64, expected: f64) {
    let budget = &f["tolerances"][key];
    assert!(
        (actual - expected).abs()
            <= budget[0].as_f64().unwrap() + budget[1].as_f64().unwrap() * expected.abs(),
        "{key}: {actual} != {expected}"
    );
}
fn compare_tp(f: &Value, c: &Value, r: &thermo_solvers::TpFlash) {
    assert_eq!(
        r.phases.len(),
        c["phases"].as_u64().unwrap() as usize,
        "{}",
        c["id"]
    );
    close(f, "h", r.enthalpy(), c["enthalpy"].as_f64().unwrap());
    assert!(r.max_ln_fugacity_residual <= 1e-10 && r.max_component_balance_residual <= 1e-8);
    assert_eq!(r.feed_stability.stable_in_trials, r.phases.len() == 1);
    for (i, phase) in r.phases.iter().enumerate() {
        assert_eq!(phase.state.model.identifier, "GERG-2008");
        assert_eq!(
            phase.state.caloric_reference.identifier,
            "gerg-native-ideal"
        );
        let props = if r.phases.len() == 1 {
            c
        } else if i == 0 {
            &c["vapor_properties"]
        } else {
            &c["liquid_properties"]
        };
        let composition = if r.phases.len() == 1 {
            &c["z"]
        } else if i == 0 {
            &c["vapor"]
        } else {
            &c["liquid"]
        };
        for (x, expected) in phase.composition.iter().zip(numbers(composition)) {
            close(f, "composition", *x, expected);
        }
        for (key, v) in [
            ("rho", 1. / phase.state.molar_volume),
            ("Z", phase.state.compressibility),
            ("h", phase.state.enthalpy),
            ("s", phase.state.entropy),
            ("cp", phase.state.cp),
        ] {
            close(f, key, v, props[key].as_f64().unwrap());
        }
    }
    if r.phases.len() == 2 {
        close(
            f,
            "vapor_fraction",
            r.phases[0].phase_fraction,
            c["vapor_fraction"].as_f64().unwrap(),
        );
        assert!(r.reduced_gibbs_change < 0.);
        assert!(r.phase_stability.iter().all(|s| s.stable_in_trials));
    }
}

#[test]
fn independent_tp_equilibrium_states() {
    let f = fixture();
    assert_eq!(f["tp"].as_array().unwrap().len(), 13);
    for c in f["tp"].as_array().unwrap() {
        println!("GERG TP reference {}", c["id"]);
        let r = gerg_tp_flash(
            &case_model(c),
            Temperature::new(c["T"].as_f64().unwrap()).unwrap(),
            Pressure::new(c["p"].as_f64().unwrap()).unwrap(),
            FlashOptions::default(),
        )
        .unwrap();
        compare_tp(&f, c, &r);
    }
}

#[test]
fn independent_ph_mixtures_pure_coexistence_and_valve() {
    let f = fixture();
    assert_eq!(f["ph"].as_array().unwrap().len(), 15);
    for c in f["ph"].as_array().unwrap() {
        println!("GERG PH reference {}", c["id"]);
        let bracket = numbers(&c["temperature_bracket"]);
        let r = gerg_ph_flash(
            &case_model(c),
            Pressure::new(c["p"].as_f64().unwrap()).unwrap(),
            c["target_h"].as_f64().unwrap(),
            GergPhFlashOptions {
                minimum_temperature_k: bracket[0],
                maximum_temperature_k: bracket[1],
                ..Default::default()
            },
        )
        .unwrap();
        close(&f, "temperature", r.temperature_k, c["T"].as_f64().unwrap());
        assert!(r.scaled_energy_residual.abs() <= 1e-10);
        if !c["saturation"].is_null() {
            assert!(r.pure_saturation.is_some() && r.tp.is_none());
            close(
                &f,
                "vapor_fraction",
                r.phases[0].phase_fraction,
                c["vapor_fraction"].as_f64().unwrap(),
            );
            for (i, phase) in r.phases.iter().enumerate() {
                let props = if i == 0 {
                    &c["saturation"]["vapor_properties"]
                } else {
                    &c["saturation"]["liquid_properties"]
                };
                close(&f, "h", phase.state.enthalpy, props["h"].as_f64().unwrap());
                close(
                    &f,
                    "rho",
                    1. / phase.state.molar_volume,
                    props["rho"].as_f64().unwrap(),
                );
            }
        } else {
            compare_tp(&f, &c["expected_tp"], r.tp.as_ref().unwrap());
        }
    }
}

#[test]
fn root_policy_permutation_zeros_and_failure_paths() {
    let m = model(&["methane", "propane"], &[0.5, 0.5]);
    let t = Temperature::new(250.).unwrap();
    let p = Pressure::new(3e6).unwrap();
    let roots = gerg_outer_roots(&m, t, p).unwrap();
    assert_eq!(
        roots.len(),
        1,
        "the interior multiparameter loop is excluded"
    );
    assert!(roots[0].state.molar_density > 13000.);
    let r = gerg_tp_flash(&m, t, p, FlashOptions::default()).unwrap();
    let permuted = model(&["propane", "nitrogen", "methane"], &[0.5, 0., 0.5]);
    let other = gerg_tp_flash(&permuted, t, p, FlashOptions::default()).unwrap();
    assert!((r.phases[0].phase_fraction - other.phases[0].phase_fraction).abs() < 2e-8);
    for (a, b) in r.phases.iter().zip(&other.phases) {
        assert_eq!(b.composition[1], 0.);
        assert!((a.composition[0] - b.composition[2]).abs() < 2e-8);
        assert!((a.state.enthalpy - b.state.enthalpy).abs() < 2e-5);
    }
    for options in [
        FlashOptions {
            max_iterations: 0,
            ..Default::default()
        },
        FlashOptions {
            tolerance: f64::NAN,
            ..Default::default()
        },
        FlashOptions {
            tolerance: 1e-6,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            gerg_tp_flash(&m, t, p, options),
            Err(Error::UnsupportedCapability(_))
        ));
    }
    assert!(matches!(
        gerg_tp_flash(
            &m,
            t,
            p,
            FlashOptions {
                max_iterations: 1,
                ..Default::default()
            }
        ),
        Err(Error::NonConvergence(_))
    ));
    assert!(
        gerg_assess_stability(
            &m,
            Temperature::new(330.).unwrap(),
            Pressure::new(2e6).unwrap(),
            StabilityOptions {
                max_iterations: 1,
                ..Default::default()
            }
        )
        .is_err()
    );
    for (t, p) in [(229., 3e6), (351., 3e6), (250., 99999.), (250., 1.01e7)] {
        assert!(matches!(
            gerg_tp_flash(
                &m,
                Temperature::new(t).unwrap(),
                Pressure::new(p).unwrap(),
                FlashOptions::default()
            ),
            Err(Error::OutsideValidityDomain(_))
        ));
    }
    assert!(matches!(
        gerg_tp_flash(
            &model(&["methane", "water"], &[0.99, 0.01]),
            t,
            p,
            FlashOptions::default()
        ),
        Err(Error::UnsupportedCapability(_))
    ));
    let unqualified = model(
        &["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
        &[0.22228531, 0.10768145, 0.02121619, 0.1428685, 0.50594855],
    );
    assert!(matches!(
        gerg_outer_roots(
            &unqualified,
            Temperature::new(240.).unwrap(),
            Pressure::new(5e6).unwrap()
        ),
        Err(Error::OutsideValidityDomain(_))
    ));
}

#[test]
fn ph_absolute_energy_floor_and_explicit_failures() {
    let methane = model(&["methane"], &[1.]);
    let p = Pressure::new(1e5).unwrap();
    let zero = gerg_ph_flash(&methane, p, 0., GergPhFlashOptions::default()).unwrap();
    assert!(zero.enthalpy.abs() <= 1e-10);
    assert!(matches!(
        gerg_ph_flash(&methane, p, f64::NAN, Default::default()),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        gerg_ph_flash(&methane, p, 1e10, Default::default()),
        Err(Error::OutsideValidityDomain(_))
    ));
    for options in [
        GergPhFlashOptions {
            minimum_temperature_k: 200.,
            ..Default::default()
        },
        GergPhFlashOptions {
            max_iterations: 0,
            ..Default::default()
        },
        GergPhFlashOptions {
            energy_tolerance: 0.,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            gerg_ph_flash(&methane, p, 0., options),
            Err(Error::UnsupportedCapability(_))
        ));
    }
    assert!(matches!(
        gerg_ph_flash(
            &methane,
            p,
            0.,
            GergPhFlashOptions {
                max_iterations: 1,
                ..Default::default()
            }
        ),
        Err(Error::NonConvergence(_))
    ));
    let f = fixture();
    let c = &f["ph"][6];
    let propane = model(&["nitrogen", "propane"], &[0., 1.]);
    assert!(matches!(
        gerg_tp_flash(
            &propane,
            Temperature::new(c["T"].as_f64().unwrap()).unwrap(),
            Pressure::new(c["p"].as_f64().unwrap()).unwrap(),
            Default::default()
        ),
        Err(Error::UnsupportedCapability(_))
    ));
    let pure = gerg_ph_flash(
        &propane,
        Pressure::new(c["p"].as_f64().unwrap()).unwrap(),
        c["target_h"].as_f64().unwrap(),
        Default::default(),
    )
    .unwrap();
    assert!(pure.pure_saturation.is_some());
    assert_eq!(pure.phases[0].composition[0], 0.);
}
