//! A7 gamma–phi solver behavior: equilibrium identities, TP→PH recovery, component
//! permutation, inactive components, liquid–liquid detection and explicit failures.
use thermo_data::{AlphaTerms, ComponentDatabase, LiquidDatabase, TauTerms};
use thermo_models::{ActivityModel, NrtlParameters, VaporModel};
use thermo_solvers::{
    ActivityOptions, activity_bubble_pressure, activity_bubble_temperature, activity_dew_pressure,
    activity_dew_temperature, activity_ph_flash, activity_tp_flash, liquid_stability,
};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, InputError, NormalizationPolicy,
    Phase, Pressure, Temperature,
};

fn composition(names: [&str; 2], z: [f64; 2]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        names
            .iter()
            .zip(z)
            .map(|(n, x)| (ComponentId::new(*n).unwrap(), x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

fn model_ordered(names: [&str; 2], z: [f64; 2]) -> ActivityModel {
    ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(names, z),
        VaporModel::IdealGas,
    )
    .unwrap()
}

fn model(z: [f64; 2]) -> ActivityModel {
    model_ordered(["methanol", "water"], z)
}

fn t(v: f64) -> Temperature {
    Temperature::new(v).unwrap()
}
fn p(v: f64) -> Pressure {
    Pressure::new(v).unwrap()
}

#[test]
fn two_phase_split_satisfies_equilibrium_balances_and_gibbs_decrease() {
    let m = model([0.4, 0.6]);
    let flash = activity_tp_flash(&m, t(355.0), p(101_325.0), ActivityOptions::default()).unwrap();
    assert_eq!(flash.phases.len(), 2);
    assert!(flash.max_ln_fugacity_residual <= 1e-11);
    assert!(flash.max_component_balance_residual <= 1e-12);
    let (v, l) = (&flash.phases[0], &flash.phases[1]);
    assert_eq!((v.phase, l.phase), (Phase::Vapor, Phase::Liquid));
    // ln f equality through the phase states: ln(y φ^V p) = ln(x φ^L p).
    for i in 0..2 {
        let lv = v.composition[i].ln() + v.state.ln_fugacity_coefficients[i];
        let ll = l.composition[i].ln() + l.state.ln_fugacity_coefficients[i];
        assert!((lv - ll).abs() < 1e-10);
    }
    // The split lowers G relative to the homogeneous liquid and vapor feeds.
    let g_split = v.phase_fraction * v.state.gibbs_energy + l.phase_fraction * l.state.gibbs_energy;
    let gl = m.liquid(t(355.0), p(101_325.0)).unwrap().gibbs_energy;
    let gv = m.vapor(t(355.0), p(101_325.0)).unwrap().gibbs_energy;
    assert!(g_split < gl && g_split < gv);
    assert!(flash.liquid_stability.as_ref().unwrap().stable_in_trials);
    assert!(
        flash
            .diagnostics
            .iter()
            .any(|d| d.code == "A7-LIQUID-TPD-FINITE-TRIALS")
    );
}

#[test]
fn saturation_points_bound_the_two_phase_region() {
    let m = model([0.4, 0.6]);
    let o = ActivityOptions::default();
    let bubble = activity_bubble_pressure(&m, t(350.0)).unwrap();
    let dew = activity_dew_pressure(&m, t(350.0), o).unwrap();
    assert!(dew.pressure_pa < bubble.pressure_pa);
    let at = |pressure| activity_tp_flash(&m, t(350.0), p(pressure), o).unwrap();
    assert_eq!(at(bubble.pressure_pa).phases[0].phase, Phase::Liquid);
    assert_eq!(at(dew.pressure_pa).phases[0].phase, Phase::Vapor);
    let near_bubble = at(bubble.pressure_pa * (1.0 - 1e-7));
    assert!(near_bubble.vapor_fraction() < 1e-5);
    // Incipient vapor at the bubble point agrees with the split limit.
    for i in 0..2 {
        assert!((near_bubble.phases[0].composition[i] - bubble.vapor[i]).abs() < 1e-5);
    }
    // Bubble/dew temperature invert bubble/dew pressure.
    let bt = activity_bubble_temperature(&m, p(bubble.pressure_pa), o).unwrap();
    let dt = activity_dew_temperature(&m, p(dew.pressure_pa), o).unwrap();
    assert!((bt.temperature_k - 350.0).abs() < 1e-8);
    assert!((dt.temperature_k - 350.0).abs() < 1e-8);
}

#[test]
fn tp_to_ph_round_trips_recover_temperature_and_phase_split() {
    let o = ActivityOptions::default();
    for (z, temperature, pressure) in [
        ([0.4, 0.6], 355.0, 101_325.0),
        ([0.4, 0.6], 320.0, 101_325.0),
        ([0.4, 0.6], 390.0, 101_325.0),
        ([0.9, 0.1], 340.0, 101_325.0),
        ([0.05, 0.95], 395.0, 250_000.0),
        ([0.5, 0.5], 290.0, 10_000.0),
    ] {
        let m = model(z);
        let tp = activity_tp_flash(&m, t(temperature), p(pressure), o).unwrap();
        let ph = activity_ph_flash(&m, p(pressure), tp.enthalpy(), o).unwrap();
        assert!(
            (ph.temperature_k - temperature).abs() < 1e-7,
            "{z:?} {temperature}"
        );
        assert_eq!(ph.phases.len(), tp.phases.len());
        assert!((ph.phases[0].phase_fraction - tp.phases[0].phase_fraction).abs() < 1e-8);
        assert!(ph.scaled_energy_residual.abs() <= 1e-10);
        assert!(ph.bubble_temperature_k.unwrap() <= ph.dew_temperature_k.unwrap());
    }
}

/// Review regression: two-phase targets near h = 0 J/mol (about 95 % vaporized at
/// 1 atm) must converge; a 1 J/mol energy floor is unattainable for the split.
#[test]
fn ph_two_phase_targets_near_zero_enthalpy_converge() {
    let o = ActivityOptions::default();
    // z = 0.95 methanol hits the dew-temperature endpoint, where the split converges to
    // β = 1 + O(1e-15); that roundoff must be reported as saturated vapor, not a failure.
    for z in [[0.1, 0.9], [0.25, 0.75], [0.9, 0.1], [0.95, 0.05]] {
        let m = model(z);
        for target in [-0.01, 0.0, 0.1, 0.3, 1.1, 2.2] {
            let ph = activity_ph_flash(&m, p(101_325.0), target, o)
                .unwrap_or_else(|e| panic!("{z:?} {target}: {e:?}"));
            assert_eq!(ph.phases.len(), 2, "{z:?} {target}");
            assert!(ph.energy_scale > 1e4, "{z:?} {target}");
            assert!(ph.scaled_energy_residual.abs() <= o.energy_tolerance);
            let tp = activity_tp_flash(&m, t(ph.temperature_k), p(101_325.0), o).unwrap();
            assert!((tp.enthalpy() - target).abs() <= o.energy_tolerance * ph.energy_scale);
        }
    }
}

/// Review regression: at low pressure the feed bubble temperature lies below the
/// bracket; single-phase vapor and two-phase targets must still be solved.
#[test]
fn ph_with_saturation_outside_the_bracket_solves_the_reachable_branches() {
    let o = ActivityOptions::default();
    let m = model([0.4, 0.6]);
    for temperature in [350.0, 300.0, 280.0, 274.0] {
        let tp = activity_tp_flash(&m, t(temperature), p(1000.0), o).unwrap();
        let ph = activity_ph_flash(&m, p(1000.0), tp.enthalpy(), o).unwrap();
        assert!(
            (ph.temperature_k - temperature).abs() < 1e-7,
            "{temperature}"
        );
        assert_eq!(ph.phases.len(), tp.phases.len(), "{temperature}");
        assert_eq!(ph.bubble_temperature_k, None);
        assert!(
            ph.diagnostics
                .iter()
                .any(|d| d.code == "A7-PH-SATURATION-OUTSIDE-BRACKET")
        );
    }
    // Below the two-phase enthalpy at the bracket minimum there is no liquid branch.
    let low = activity_tp_flash(&m, t(273.16), p(1000.0), o)
        .unwrap()
        .enthalpy();
    assert!(matches!(
        activity_ph_flash(&m, p(1000.0), low - 1000.0, o),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    // High pressure: the dew temperature lies above the bracket, liquid still solvable.
    let hp = model([0.05, 0.95]);
    let tp = activity_tp_flash(&hp, t(400.0), p(500_000.0), o).unwrap();
    let ph = activity_ph_flash(&hp, p(500_000.0), tp.enthalpy(), o).unwrap();
    assert!((ph.temperature_k - 400.0).abs() < 1e-7);
}

#[test]
fn component_order_does_not_change_physical_results() {
    let o = ActivityOptions::default();
    let a = activity_tp_flash(&model([0.3, 0.7]), t(360.0), p(101_325.0), o).unwrap();
    let b = activity_tp_flash(
        &model_ordered(["water", "methanol"], [0.7, 0.3]),
        t(360.0),
        p(101_325.0),
        o,
    )
    .unwrap();
    assert!((a.vapor_fraction() - b.vapor_fraction()).abs() < 1e-12);
    for (pa, pb) in a.phases.iter().zip(&b.phases) {
        assert!((pa.composition[0] - pb.composition[1]).abs() < 1e-12);
        assert!((pa.state.enthalpy - pb.state.enthalpy).abs() < 1e-8);
    }
    let ha = activity_ph_flash(&model([0.3, 0.7]), p(101_325.0), -20_000.0, o).unwrap();
    let hb = activity_ph_flash(
        &model_ordered(["water", "methanol"], [0.7, 0.3]),
        p(101_325.0),
        -20_000.0,
        o,
    )
    .unwrap();
    assert!((ha.temperature_k - hb.temperature_k).abs() < 1e-8);
}

#[test]
fn inactive_components_are_retained_and_pure_ph_uses_the_lever_rule() {
    let o = ActivityOptions::default();
    let m = model([0.0, 1.0]);
    let bubble = activity_bubble_pressure(&m, t(373.0)).unwrap();
    assert_eq!(bubble.vapor[0], 0.0);
    let liquid = activity_tp_flash(&m, t(373.0), p(bubble.pressure_pa * 1.01), o).unwrap();
    let vapor = activity_tp_flash(&m, t(373.0), p(bubble.pressure_pa * 0.99), o).unwrap();
    assert_eq!(liquid.phases[0].phase, Phase::Liquid);
    assert_eq!(vapor.phases[0].phase, Phase::Vapor);
    let sat = activity_bubble_temperature(&m, p(101_325.0), o).unwrap();
    let l = m
        .liquid(t(sat.temperature_k), p(101_325.0))
        .unwrap()
        .enthalpy;
    let v = m
        .vapor(t(sat.temperature_k), p(101_325.0))
        .unwrap()
        .enthalpy;
    for beta in [0.0, 0.4, 1.0] {
        let ph = activity_ph_flash(&m, p(101_325.0), l + beta * (v - l), o).unwrap();
        assert!((ph.temperature_k - sat.temperature_k).abs() < 1e-12);
        let vf: f64 = ph
            .phases
            .iter()
            .filter(|x| x.phase == Phase::Vapor)
            .map(|x| x.phase_fraction)
            .sum();
        assert!((vf - beta).abs() < 1e-9, "beta {beta}: {vf}");
    }
}

#[test]
fn narrow_boiling_mixture_resolves_ph_inside_the_interval() {
    let o = ActivityOptions::default();
    let m = model([0.002, 0.998]);
    let bt = activity_bubble_temperature(&m, p(101_325.0), o)
        .unwrap()
        .temperature_k;
    let dt = activity_dew_temperature(&m, p(101_325.0), o)
        .unwrap()
        .temperature_k;
    assert!(dt > bt && dt - bt < 1.0);
    let mid = activity_tp_flash(&m, t(0.5 * (bt + dt)), p(101_325.0), o).unwrap();
    assert_eq!(mid.phases.len(), 2);
    let ph = activity_ph_flash(&m, p(101_325.0), mid.enthalpy(), o).unwrap();
    assert!((ph.temperature_k - 0.5 * (bt + dt)).abs() < 1e-8);
}

fn lle_model() -> ActivityModel {
    // Strongly positive, explicitly supplied interactions: a liquid–liquid split.
    let z = TauTerms {
        t0: 0.0,
        t1: 0.0,
        t2: 0.0,
        t3: 0.0,
    };
    let s = TauTerms { t0: 3.0, ..z };
    let a = AlphaTerms { a0: 0.2, a1: 0.0 };
    let za = AlphaTerms { a0: 0.0, a1: 0.0 };
    let parameters = NrtlParameters::new(
        vec![
            ComponentId::new("methanol").unwrap(),
            ComponentId::new("water").unwrap(),
        ],
        vec![vec![z, s], vec![s, z]],
        vec![vec![za, a], vec![a, za]],
        "hypothetical immiscible set for LLE detection tests only",
    )
    .unwrap();
    ActivityModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(["methanol", "water"], [0.5, 0.5]),
        VaporModel::IdealGas,
        parameters,
    )
    .unwrap()
}

#[test]
fn liquid_liquid_split_is_detected_and_reported_unsupported() {
    let o = ActivityOptions::default();
    let m = lle_model();
    let s = liquid_stability(&m, t(300.0), &[0.5, 0.5], o).unwrap();
    assert!(!s.stable_in_trials);
    assert!(matches!(
        activity_tp_flash(&m, t(300.0), p(101_325.0), o),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    // The real dataset is miscible in the same trials.
    let real = liquid_stability(&model([0.5, 0.5]), t(300.0), &[0.5, 0.5], o).unwrap();
    assert!(real.stable_in_trials && real.trials.len() == 2);
}

#[test]
fn iteration_limits_and_invalid_options_are_explicit() {
    let m = model([0.4, 0.6]);
    let limited = ActivityOptions {
        max_iterations: 1,
        ..ActivityOptions::default()
    };
    assert!(matches!(
        activity_dew_pressure(&m, t(350.0), limited),
        Err(CalculationError::NonConvergence(_))
    ));
    assert!(matches!(
        activity_tp_flash(&m, t(355.0), p(101_325.0), limited),
        Err(CalculationError::NonConvergence(_))
    ));
    let outer = ActivityOptions {
        max_temperature_iterations: 1,
        ..ActivityOptions::default()
    };
    assert!(matches!(
        activity_bubble_temperature(&m, p(101_325.0), outer),
        Err(CalculationError::NonConvergence(_))
    ));
    let zero = ActivityOptions {
        max_iterations: 0,
        ..ActivityOptions::default()
    };
    assert!(matches!(
        activity_tp_flash(&m, t(355.0), p(101_325.0), zero),
        Err(CalculationError::NonConvergence(_))
    ));
    for bad in [
        ActivityOptions {
            tolerance: 0.0,
            ..ActivityOptions::default()
        },
        ActivityOptions {
            energy_tolerance: 1e-3,
            ..ActivityOptions::default()
        },
        ActivityOptions {
            minimum_temperature_k: 200.0,
            ..ActivityOptions::default()
        },
        ActivityOptions {
            tm_threshold: 0.0,
            ..ActivityOptions::default()
        },
    ] {
        assert!(matches!(
            activity_tp_flash(&m, t(355.0), p(101_325.0), bad),
            Err(CalculationError::UnsupportedCapability(_))
        ));
    }
}

#[test]
fn unsupported_states_and_targets_fail_explicitly() {
    let o = ActivityOptions::default();
    let m = model([0.4, 0.6]);
    assert!(matches!(
        activity_tp_flash(&m, t(500.0), p(101_325.0), o),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        activity_tp_flash(&m, t(350.0), p(2e6), o),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        activity_ph_flash(&m, p(101_325.0), f64::NAN, o),
        Err(CalculationError::InvalidInput(InputError::NonFinite(_)))
    ));
    // Far below the liquid enthalpy at the lower bracket and far above the vapor.
    for h in [-1e6, 1e6] {
        assert!(matches!(
            activity_ph_flash(&m, p(101_325.0), h, o),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
    }
    // Bubble temperature outside the bracket (water at 100 Pa boils below 273.16 K).
    assert!(matches!(
        activity_bubble_temperature(&model([0.0, 1.0]), p(100.0), o),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    // Bubble pressure outside the package pressure envelope.
    assert!(matches!(
        activity_bubble_pressure(&model([1.0, 0.0]), t(450.0)),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}
