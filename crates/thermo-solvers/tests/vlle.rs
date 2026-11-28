//! Multiphase (VLLE/LLE) stability and flash (#224): structure, consistency,
//! robustness and failure checks. The independent Clapeyron.jl comparisons are in
//! `vlle_clapeyron.rs`.

use thermo_data::{AlphaTerms, ComponentDatabase, LiquidDatabase, TauTerms};
use thermo_models::{ActivityModel, CubicModel, CubicVariant, NrtlParameters, VaporModel};
use thermo_solvers::*;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase,
    Pressure, Temperature,
};

/// Explicit water–hydrocarbon interaction for these checks (an approximation, not
/// a qualified dataset); hydrocarbon pairs are zero by explicit choice.
const K_WATER_HC: f64 = 0.5;

fn composition(ids: &[&str], x: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(x)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::Normalize,
    )
    .unwrap()
}

fn cubic(ids: &[&str], x: &[f64], variant: CubicVariant) -> CubicModel {
    let n = ids.len();
    let kij = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i != j && (ids[i] == "water") != (ids[j] == "water") {
                        K_WATER_HC
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition(ids, x),
        variant,
        kij,
        "explicit kij 0.5 for water-hydrocarbon pairs (VLLE verification approximation)",
    )
    .unwrap()
}

/// Synthetic partially miscible NRTL pair on the methanol/water pure-liquid data:
/// tau_12 = 350 K/T, tau_21 = 520 K/T, alpha = 0.2. Not a physical mixture.
fn synthetic_lle(x: &[f64]) -> ActivityModel {
    let ids: Vec<ComponentId> = ["methanol", "water"]
        .iter()
        .map(|i| ComponentId::new(*i).unwrap())
        .collect();
    let tau = |t1: f64| TauTerms {
        t0: 0.0,
        t1,
        t2: 0.0,
        t3: 0.0,
    };
    let alpha = AlphaTerms { a0: 0.2, a1: 0.0 };
    let none = AlphaTerms { a0: 0.0, a1: 0.0 };
    let parameters = NrtlParameters::new(
        ids,
        vec![vec![tau(0.0), tau(350.0)], vec![tau(520.0), tau(0.0)]],
        vec![vec![none, alpha], vec![alpha, none]],
        "synthetic partially miscible NRTL pair for VLLE/LLE implementation checks; not a physical mixture",
    )
    .unwrap();
    ActivityModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(&["methanol", "water"], x),
        VaporModel::IdealGas,
        parameters,
    )
    .unwrap()
}

fn tp(t: f64, p: f64) -> (Temperature, Pressure) {
    (Temperature::new(t).unwrap(), Pressure::new(p).unwrap())
}

const C1C3W: [&str; 3] = ["methane", "propane", "water"];
const FEED: [f64; 3] = [0.1, 0.4, 0.5];

fn roles<S>(r: &MultiphaseTpFlash<S>) -> Vec<PhaseRole> {
    r.phases.iter().map(|p| p.role).collect()
}

fn roles_ph<S>(r: &MultiphasePhFlash<S>) -> Vec<PhaseRole> {
    r.phases.iter().map(|p| p.role).collect()
}

/// Material balance, fractions and fugacity equality of a converged split.
fn check_split<S>(r: &MultiphaseTpFlash<S>, z: &[f64]) {
    assert!(
        r.max_ln_fugacity_residual <= 1e-10,
        "{}",
        r.max_ln_fugacity_residual
    );
    assert!(r.max_component_balance_residual <= 1e-12);
    let total: f64 = r.phases.iter().map(|p| p.phase_fraction).sum();
    assert!((total - 1.0).abs() < 1e-12);
    for (i, zi) in z.iter().enumerate() {
        let sum: f64 = r
            .phases
            .iter()
            .map(|p| p.phase_fraction * p.composition[i])
            .sum();
        assert!((sum - zi).abs() < 1e-12);
    }
    if r.phases.len() > 1 {
        assert!(r.reduced_gibbs_change < 0.0);
        assert!(r.final_stability.as_ref().unwrap().stable_in_trials);
    }
}

#[test]
fn water_hydrocarbon_feed_splits_into_vapor_and_two_liquids() {
    let m = cubic(&C1C3W, &FEED, CubicVariant::PengRobinson);
    let (t, p) = tp(300.0, 20e5);
    let r = vlle_tp_flash(&m, t, p, MultiphaseOptions::default()).unwrap();
    check_split(&r, &FEED);
    assert_eq!(
        roles(&r),
        [
            PhaseRole::Vapor,
            PhaseRole::LightLiquid,
            PhaseRole::HeavyLiquid
        ]
    );
    assert_eq!(
        r.phases.iter().map(|p| p.phase).collect::<Vec<_>>(),
        [Phase::Vapor, Phase::Liquid, Phase::Liquid]
    );
    assert_eq!(r.liquid_count(), 2);
    // The hydrocarbon liquid is lighter than the aqueous liquid.
    assert!(r.phases[2].composition[2] > 0.999);
    assert!(r.phases[1].composition[1] > 0.9);
    assert!(r.phases[1].state.mass_density < r.phases[2].state.mass_density);
    // The aqueous and hydrocarbon trial phases were run; the aqueous one proves
    // instability of the feed.
    let trials = &r.feed_stability.trials;
    for kind in [MultiphaseTrial::Aqueous, MultiphaseTrial::Hydrocarbon] {
        assert!(trials.iter().any(|t| t.trial == kind), "{kind:?}");
    }
    let aqueous = trials
        .iter()
        .find(|t| t.trial == MultiphaseTrial::Aqueous)
        .unwrap();
    assert!(aqueous.tangent_plane_distance < -1e-3 && aqueous.composition[2] > 0.99);
    assert!(
        r.diagnostics
            .iter()
            .any(|d| d.code == "VLLE-STABILITY-FINITE-TRIALS")
    );
    // The standalone stability test agrees.
    let s = vlle_stability(&m, t, p, MultiphaseOptions::default()).unwrap();
    assert!(!s.stable_in_trials);
}

#[test]
fn hydrocarbon_water_without_vapor_is_a_liquid_liquid_split() {
    let m = cubic(
        &["propane", "water"],
        &[0.5, 0.5],
        CubicVariant::PengRobinson,
    );
    let (t, p) = tp(300.0, 50e5);
    let r = vlle_tp_flash(&m, t, p, MultiphaseOptions::default()).unwrap();
    check_split(&r, &[0.5, 0.5]);
    assert_eq!(roles(&r), [PhaseRole::LightLiquid, PhaseRole::HeavyLiquid]);
    assert!(r.phases.iter().all(|p| p.phase == Phase::Liquid));
}

#[test]
fn dry_vapor_is_a_stable_single_phase() {
    let z = [0.6, 0.39, 0.01];
    let m = cubic(&C1C3W, &z, CubicVariant::PengRobinson);
    let (t, p) = tp(400.0, 5e5);
    let r = vlle_tp_flash(&m, t, p, MultiphaseOptions::default()).unwrap();
    assert!(r.feed_stability.stable_in_trials && r.final_stability.is_none());
    assert_eq!(roles(&r), [PhaseRole::Vapor]);
    assert_eq!(r.phases[0].composition, m.fractions());
}

#[test]
fn component_order_does_not_change_the_split() {
    let a = cubic(&C1C3W, &FEED, CubicVariant::PengRobinsonBostonMathias);
    let b = cubic(
        &["water", "propane", "methane"],
        &[0.5, 0.4, 0.1],
        CubicVariant::PengRobinsonBostonMathias,
    );
    let (t, p) = tp(300.0, 25e5);
    let ra = vlle_tp_flash(&a, t, p, MultiphaseOptions::default()).unwrap();
    let rb = vlle_tp_flash(&b, t, p, MultiphaseOptions::default()).unwrap();
    assert_eq!(ra.phases.len(), 3);
    assert_eq!(roles(&ra), roles(&rb));
    for (pa, pb) in ra.phases.iter().zip(&rb.phases) {
        assert!((pa.phase_fraction - pb.phase_fraction).abs() < 1e-9);
        for (i, j) in [(0, 2), (1, 1), (2, 0)] {
            let (x, y) = (pa.composition[i], pb.composition[j]);
            assert!((x - y).abs() <= 1e-9 + 1e-7 * x.abs(), "{x} vs {y}");
        }
    }
}

#[test]
fn split_converges_across_the_cubic_temperature_envelope() {
    // The PH search visits temperatures across its bracket; every TP flash of
    // this feed must converge there.
    for i in 0..=80 {
        let t = 200.0 + 5.0 * i as f64;
        let m = cubic(&C1C3W, &FEED, CubicVariant::PengRobinson);
        let (t_, p_) = tp(t, 20e5);
        let r = vlle_tp_flash(&m, t_, p_, MultiphaseOptions::default())
            .unwrap_or_else(|e| panic!("{t} K: {e}"));
        check_split(&r, &FEED);
    }
}

#[test]
fn ph_reproduces_tp_states() {
    for (ids, z, t, p) in [
        (&C1C3W[..], &FEED[..], 300.0, 20e5),
        (&C1C3W[..], &FEED[..], 320.0, 5e5),
        (&["propane", "water"][..], &[0.5, 0.5][..], 300.0, 50e5),
    ] {
        let m = cubic(ids, z, CubicVariant::PengRobinson);
        let (t_, p_) = tp(t, p);
        let r = vlle_tp_flash(&m, t_, p_, MultiphaseOptions::default()).unwrap();
        let h: f64 = r
            .phases
            .iter()
            .map(|p| p.phase_fraction * p.state.enthalpy)
            .sum();
        let ph = vlle_ph_flash(&m, p_, h, MultiphaseOptions::default()).unwrap();
        assert!(
            (ph.temperature_k - t).abs() < 1e-6,
            "{} vs {t}",
            ph.temperature_k
        );
        assert!(ph.scaled_energy_residual.abs() <= 1e-10);
        assert!(!ph.invariant_temperature);
        assert_eq!(ph.phases.len(), r.phases.len());
    }
}

#[test]
fn propane_rich_ph_inside_the_three_phase_enthalpy_step() {
    let model = cubic(
        &["propane", "water"],
        &[0.9957, 0.0043],
        CubicVariant::PengRobinson,
    );
    let pressure = Pressure::new(20e5).unwrap();
    // Clapeyron.jl reference: validation/julia/check_vlle_propane_rich.jl.
    let options = MultiphaseOptions::default();
    let result = vlle_ph_flash(&model, pressure, -6795.525837714403, options).unwrap();
    assert!(result.invariant_temperature);
    assert_eq!(
        roles_ph(&result),
        [
            PhaseRole::Vapor,
            PhaseRole::LightLiquid,
            PhaseRole::HeavyLiquid
        ]
    );
    assert!((result.temperature_k - 329.8018791069918).abs() < 2e-6);
    assert!(result.scaled_energy_residual.abs() <= options.energy_tolerance);
    assert!(result.max_component_balance_residual <= 1e-9);
    assert!(result.max_ln_fugacity_residual <= 1e-7);
    for (phase, (fraction, propane)) in result.phases.iter().zip([
        (0.43, 0.9928975992411447),
        (0.5697877700651349, 0.998185749513824),
        (0.00021222993486513878, 1.4484158740005578e-10),
    ]) {
        assert!((phase.phase_fraction - fraction).abs() < 2e-8);
        assert!((phase.composition[0] - propane).abs() < 1e-8 + 1e-7 * propane);
    }
}

#[test]
fn pure_component_ph_inside_the_enthalpy_step_uses_the_lever_rule() {
    let m = cubic(&["propane"], &[1.0], CubicVariant::PengRobinson);
    let p = Pressure::new(10e5).unwrap();
    // Reference: the existing two-phase PH flash on its pure saturation path.
    let reference = ph_flash(&m, p, -12_000.0, PhFlashOptions::default()).unwrap();
    assert!(reference.pure_saturation.is_some());
    let r = vlle_ph_flash(&m, p, -12_000.0, MultiphaseOptions::default()).unwrap();
    assert!(r.invariant_temperature);
    assert!((r.temperature_k - reference.temperature_k).abs() < 1e-6);
    assert!(r.scaled_energy_residual.abs() <= 1e-10);
    assert_eq!(roles_ph(&r), [PhaseRole::Vapor, PhaseRole::Liquid]);
    let vapor = reference
        .phases
        .iter()
        .find(|p| p.phase == Phase::Vapor)
        .unwrap();
    assert!((r.phases[0].phase_fraction - vapor.phase_fraction).abs() < 1e-7);
    assert!(r.max_ln_fugacity_residual <= 1e-7);
}

#[test]
fn nrtl_partially_miscible_pair_splits_into_two_liquids() {
    let z = [0.5, 0.5];
    let m = synthetic_lle(&z);
    let (t, p) = tp(300.0, 1e5);
    let r = activity_vlle_tp_flash(&m, t, p, MultiphaseOptions::default()).unwrap();
    check_split(&r, &z);
    assert_eq!(roles(&r), [PhaseRole::LightLiquid, PhaseRole::HeavyLiquid]);
    assert!(r.phases.iter().all(|p| p.kind == FluidKind::Liquid));
    // Equal activities x_i gamma_i in both liquids.
    let activity = |k: usize| -> Vec<f64> {
        let ex = r.phases[k].state.excess.as_ref().unwrap();
        (0..2)
            .map(|i| r.phases[k].composition[i] * ex.ln_activity_coefficients[i].exp())
            .collect()
    };
    let (a0, a1) = (activity(0), activity(1));
    for i in 0..2 {
        assert!((a0[i] / a1[i] - 1.0).abs() < 1e-9);
    }
    // The two-phase gamma-phi flash refuses this liquid explicitly.
    let error = activity_tp_flash(&m, t, p, ActivityOptions::default()).unwrap_err();
    assert!(matches!(error, CalculationError::UnsupportedCapability(_)));
    let s = activity_vlle_stability(&m, t, p, MultiphaseOptions::default()).unwrap();
    assert!(!s.stable_in_trials);
}

#[test]
fn nrtl_liquid_above_the_critical_solution_temperature_is_stable() {
    let m = synthetic_lle(&[0.5, 0.5]);
    let (t, p) = tp(420.0, 5e5);
    let r = activity_vlle_tp_flash(&m, t, p, MultiphaseOptions::default()).unwrap();
    assert_eq!(r.phases.len(), 1);
}

#[test]
fn failures_near_the_plait_point_are_explicit() {
    // Just below the critical solution temperature the TPD trials converge too
    // slowly to decide stability: the flash fails explicitly instead of
    // reporting one or two phases.
    let m = synthetic_lle(&[0.5, 0.5]);
    let (t, p) = tp(378.0, 5e5);
    match activity_vlle_tp_flash(&m, t, p, MultiphaseOptions::default()) {
        Err(CalculationError::NonConvergence(message)) => {
            assert!(message.contains("plait"), "{message}")
        }
        other => panic!("unexpected {other:?}"),
    }
}

#[test]
fn invalid_options_and_domains_are_rejected() {
    let m = cubic(&C1C3W, &FEED, CubicVariant::PengRobinson);
    let (t, p) = tp(300.0, 20e5);
    for options in [
        MultiphaseOptions {
            tolerance: 0.0,
            ..Default::default()
        },
        MultiphaseOptions {
            tm_threshold: 0.0,
            ..Default::default()
        },
        MultiphaseOptions {
            minimum_temperature_k: 400.0,
            maximum_temperature_k: 300.0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            vlle_tp_flash(&m, t, p, options),
            Err(CalculationError::UnsupportedCapability(_))
        ));
    }
    let zero = MultiphaseOptions {
        max_iterations: 0,
        ..Default::default()
    };
    assert!(matches!(
        vlle_tp_flash(&m, t, p, zero),
        Err(CalculationError::NonConvergence(_))
    ));
    // Gamma-phi: pressure outside the ideal-gas vapor envelope, and a PH bracket
    // outside the activity temperature range.
    let a = synthetic_lle(&[0.5, 0.5]);
    let (t, p) = tp(300.0, 20e5);
    assert!(matches!(
        activity_vlle_tp_flash(&a, t, p, MultiphaseOptions::default()),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        activity_vlle_ph_flash(
            &a,
            Pressure::new(1e5).unwrap(),
            0.0,
            MultiphaseOptions::default()
        ),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    // A PH target beyond the bracket's enthalpy range is outside the domain.
    assert!(matches!(
        vlle_ph_flash(
            &m,
            Pressure::new(20e5).unwrap(),
            1e7,
            MultiphaseOptions::default()
        ),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}
