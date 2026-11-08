//! A7 NRTL and gamma–phi package: identities, derivatives, permutation, directed
//! parameters, explicit missing-data and domain failures. Engine-only checks; the
//! independent references live in `thermo-solvers/tests/nrtl_clapeyron.rs`.
use thermo_data::{AlphaTerms, ComponentDatabase, LiquidDatabase, TauTerms};
use thermo_models::{ActivityModel, NrtlModel, NrtlParameters, VaporModel};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, InputError,
    MOLAR_GAS_CONSTANT as R, NormalizationPolicy, Pressure, Temperature, ValidityStatus,
};

const ZT: TauTerms = TauTerms {
    t0: 0.0,
    t1: 0.0,
    t2: 0.0,
    t3: 0.0,
};
const ZA: AlphaTerms = AlphaTerms { a0: 0.0, a1: 0.0 };

fn ids(names: &[&str]) -> Vec<ComponentId> {
    names
        .iter()
        .map(|n| ComponentId::new(*n).unwrap())
        .collect()
}

fn tau(t0: f64, t1: f64, t2: f64, t3: f64) -> TauTerms {
    TauTerms { t0, t1, t2, t3 }
}

/// Directed ternary with every temperature term active.
fn ternary() -> NrtlParameters {
    let t = vec![
        vec![
            ZT,
            tau(0.8, -150.0, 0.05, -0.001),
            tau(-0.4, 220.0, -0.02, 0.0015),
        ],
        vec![
            tau(1.3, 310.0, -0.03, 0.002),
            ZT,
            tau(0.25, -75.0, 0.04, -0.0005),
        ],
        vec![
            tau(0.6, -40.0, 0.01, 0.0007),
            tau(-0.9, 500.0, -0.06, 0.0012),
            ZT,
        ],
    ];
    let a = |a0, a1| AlphaTerms { a0, a1 };
    let alpha = vec![
        vec![ZA, a(0.2, 0.0004), a(0.35, -0.0002)],
        vec![a(0.2, 0.0004), ZA, a(0.3, 0.0001)],
        vec![a(0.35, -0.0002), a(0.3, 0.0001), ZA],
    ];
    NrtlParameters::new(ids(&["a", "b", "c"]), t, alpha, "synthetic test set").unwrap()
}

fn mixture(z: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        vec![
            (ComponentId::new("methanol").unwrap(), z[0]),
            (ComponentId::new("water").unwrap(), z[1]),
        ],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

fn package(z: &[f64]) -> ActivityModel {
    ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &mixture(z),
        VaporModel::IdealGas,
    )
    .unwrap()
}

#[test]
fn zero_interactions_give_ideal_solution() {
    let a = AlphaTerms { a0: 0.3, a1: 0.0 };
    let p = NrtlParameters::new(
        ids(&["a", "b"]),
        vec![vec![ZT; 2]; 2],
        vec![vec![ZA, a], vec![a, ZA]],
        "ideal",
    )
    .unwrap();
    let ex = NrtlModel::new(p).excess(350.0, &[0.3, 0.7]).unwrap();
    assert!(ex.ln_activity_coefficients.iter().all(|g| *g == 0.0));
    assert_eq!(ex.excess_gibbs, 0.0);
    assert_eq!(ex.excess_enthalpy, 0.0);
    assert_eq!(ex.excess_heat_capacity, 0.0);
}

#[test]
fn excess_gibbs_is_the_mole_weighted_sum_and_obeys_gibbs_duhem() {
    let m = NrtlModel::new(ternary());
    for t in [280.0, 330.0, 400.0] {
        for x in [[0.2, 0.3, 0.5], [0.6, 0.4, 0.0], [0.05, 0.9, 0.05]] {
            let ex = m.excess(t, &x).unwrap();
            let sum: f64 = x
                .iter()
                .zip(&ex.ln_activity_coefficients)
                .map(|(x, g)| x * g)
                .sum();
            assert!((sum - ex.reduced_excess_gibbs).abs() < 1e-13);
            // Gibbs–Duhem at fixed T: Σ x_i d ln γ_i = 0 along any direction keeping Σx = 1.
            let h = 1e-6;
            let dir = [1.0, -0.4, -0.6];
            let shift =
                |s: f64| -> Vec<f64> { x.iter().zip(dir).map(|(x, d)| x + s * d).collect() };
            if x.iter().any(|v| *v < 2.0 * h) {
                continue;
            }
            let up = m.excess(t, &shift(h)).unwrap().ln_activity_coefficients;
            let dn = m.excess(t, &shift(-h)).unwrap().ln_activity_coefficients;
            let gd: f64 = (0..3).map(|i| x[i] * (up[i] - dn[i]) / (2.0 * h)).sum();
            assert!(gd.abs() < 1e-8, "Gibbs–Duhem {gd:e} at T={t} x={x:?}");
        }
    }
}

#[test]
fn ln_gamma_is_the_partial_molar_derivative_of_n_ge() {
    let m = NrtlModel::new(ternary());
    let t = 330.0;
    let n = [0.2, 0.3, 0.5];
    let nge = |n: &[f64]| -> f64 {
        let total: f64 = n.iter().sum();
        let x: Vec<f64> = n.iter().map(|v| v / total).collect();
        total * m.excess(t, &x).unwrap().reduced_excess_gibbs
    };
    let ex = m.excess(t, &n).unwrap();
    for i in 0..3 {
        let h = 1e-6;
        let mut up = n;
        let mut dn = n;
        up[i] += h;
        dn[i] -= h;
        let fd = (nge(&up) - nge(&dn)) / (2.0 * h);
        assert!(
            (fd - ex.ln_activity_coefficients[i]).abs() < 1e-8,
            "component {i}"
        );
    }
}

#[test]
fn temperature_derivatives_match_finite_differences() {
    let m = NrtlModel::new(ternary());
    let x = [0.2, 0.3, 0.5];
    let t = 330.0;
    let h = 1e-3;
    let (up, mid, dn) = (
        m.excess(t + h, &x).unwrap(),
        m.excess(t, &x).unwrap(),
        m.excess(t - h, &x).unwrap(),
    );
    for i in 0..3 {
        let fd = (up.ln_activity_coefficients[i] - dn.ln_activity_coefficients[i]) / (2.0 * h);
        assert!((fd - mid.d_ln_activity_coefficients_dt[i]).abs() < 1e-9);
    }
    let fd = (up.reduced_excess_gibbs - dn.reduced_excess_gibbs) / (2.0 * h);
    assert!((-R * t * t * fd - mid.excess_enthalpy).abs() < 1e-5);
    let fd_cp = (up.excess_enthalpy - dn.excess_enthalpy) / (2.0 * h);
    assert!((fd_cp - mid.excess_heat_capacity).abs() < 1e-5);
    assert!((mid.excess_gibbs - (mid.excess_enthalpy - t * mid.excess_entropy)).abs() < 1e-9);
}

#[test]
fn component_permutation_permutes_results() {
    let base = NrtlModel::new(ternary());
    let order = [2, 0, 1];
    let permuted = NrtlModel::new(ternary().permuted(&order).unwrap());
    let x = [0.2, 0.3, 0.5];
    let px: Vec<f64> = order.iter().map(|&i| x[i]).collect();
    let a = base.excess(345.0, &x).unwrap();
    let b = permuted.excess(345.0, &px).unwrap();
    for (k, &i) in order.iter().enumerate() {
        assert!((a.ln_activity_coefficients[i] - b.ln_activity_coefficients[k]).abs() < 1e-14);
    }
    assert!((a.excess_enthalpy - b.excess_enthalpy).abs() < 1e-10);
    assert!(ternary().permuted(&[0, 0, 1]).is_err());
}

#[test]
fn parameters_are_directed() {
    // Swapping τ12 and τ21 is not the same model unless components are swapped too.
    let make = |t12: f64, t21: f64| {
        let t = vec![
            vec![ZT, tau(0.0, t12, 0.0, 0.0)],
            vec![tau(0.0, t21, 0.0, 0.0), ZT],
        ];
        let a = AlphaTerms { a0: 0.3, a1: 0.0 };
        NrtlModel::new(
            NrtlParameters::new(
                ids(&["a", "b"]),
                t,
                vec![vec![ZA, a], vec![a, ZA]],
                "directed",
            )
            .unwrap(),
        )
    };
    let forward = make(-95.0, 399.0).excess(350.0, &[0.3, 0.7]).unwrap();
    let swapped = make(399.0, -95.0).excess(350.0, &[0.3, 0.7]).unwrap();
    let mirrored = make(399.0, -95.0).excess(350.0, &[0.7, 0.3]).unwrap();
    assert!(
        (forward.ln_activity_coefficients[0] - swapped.ln_activity_coefficients[0]).abs() > 1e-3
    );
    assert!(
        (forward.ln_activity_coefficients[0] - mirrored.ln_activity_coefficients[1]).abs() < 1e-14
    );
    assert!(
        (forward.ln_activity_coefficients[1] - mirrored.ln_activity_coefficients[0]).abs() < 1e-14
    );
}

#[test]
fn database_orientation_and_missing_pairs_are_explicit() {
    let db = LiquidDatabase::embedded().unwrap();
    let mw = NrtlParameters::from_database(&db, &ids(&["methanol", "water"])).unwrap();
    let wm = NrtlParameters::from_database(&db, &ids(&["water", "methanol"])).unwrap();
    assert_eq!(mw.tau_terms(0, 1), wm.tau_terms(1, 0));
    assert_eq!(mw.tau_terms(1, 0), wm.tau_terms(0, 1));
    // ChemSep/DECHEMA: A12 = -189.0469 cal/mol (methanol→water), tau = A·4.184/(R T).
    assert!((mw.tau_terms(0, 1).t1 - (-189.0469 * 4.184 / R)).abs() < 1e-12);
    assert_eq!(mw.alpha_terms(0, 1).a0, 0.2999);
    for pair in [["methanol", "ethane"], ["water", "methane"]] {
        assert!(matches!(
            NrtlParameters::from_database(&db, &ids(&pair)),
            Err(CalculationError::MissingParameter(_))
        ));
    }
    assert!(matches!(
        db.pure("propane"),
        Err(CalculationError::MissingParameter(_))
    ));
    // A pure-liquid correlation is required even when interactions exist.
    let c = Composition::new(
        CompositionBasis::Mole,
        vec![(ComponentId::new("ethane").unwrap(), 1.0)],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let result = ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &db,
        &c,
        VaporModel::IdealGas,
    );
    assert!(matches!(result, Err(CalculationError::MissingParameter(_))));
}

#[test]
fn invalid_parameter_sets_are_rejected() {
    let two = ids(&["a", "b"]);
    let a = AlphaTerms { a0: 0.3, a1: 0.0 };
    let alpha = vec![vec![ZA, a], vec![a, ZA]];
    let good = vec![
        vec![ZT, tau(0.0, 100.0, 0.0, 0.0)],
        vec![tau(0.0, 50.0, 0.0, 0.0), ZT],
    ];
    let mut diagonal = good.clone();
    diagonal[0][0] = tau(0.1, 0.0, 0.0, 0.0);
    assert!(matches!(
        NrtlParameters::new(two.clone(), diagonal, alpha.clone(), "x"),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    let asym = vec![vec![ZA, a], vec![AlphaTerms { a0: 0.2, a1: 0.0 }, ZA]];
    assert!(matches!(
        NrtlParameters::new(two.clone(), good.clone(), asym, "x"),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    let mut nan = good.clone();
    nan[0][1].t1 = f64::NAN;
    assert!(matches!(
        NrtlParameters::new(two.clone(), nan, alpha.clone(), "x"),
        Err(CalculationError::InvalidInput(InputError::NonFinite(_)))
    ));
    assert!(matches!(
        NrtlParameters::new(two.clone(), good.clone(), alpha.clone(), "  "),
        Err(CalculationError::MissingParameter(_))
    ));
    assert!(matches!(
        NrtlParameters::new(two.clone(), vec![good[0].clone()], alpha.clone(), "x"),
        Err(CalculationError::MissingParameter(_))
    ));
    assert!(matches!(
        NrtlParameters::new(ids(&["a", "a"]), good.clone(), alpha.clone(), "x"),
        Err(CalculationError::InvalidInput(
            InputError::DuplicateComponent
        ))
    ));
    // alpha(T) = 0.3 - 0.001 T is nonpositive at 350 K.
    let falling = AlphaTerms {
        a0: 0.3,
        a1: -0.001,
    };
    let m = NrtlModel::new(
        NrtlParameters::new(two, good, vec![vec![ZA, falling], vec![falling, ZA]], "x").unwrap(),
    );
    assert!(matches!(
        m.excess(350.0, &[0.5, 0.5]),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(m.excess(250.0, &[0.5, 0.5]).is_ok());
    for bad in [[0.5, 0.6], [f64::NAN, 1.0], [-0.1, 1.1]] {
        assert!(matches!(
            m.excess(250.0, &bad),
            Err(CalculationError::InvalidInput(_))
        ));
    }
    assert!(m.excess(250.0, &[1.0]).is_err());
}

#[test]
fn validity_reports_regression_interval() {
    let m = package(&[0.5, 0.5]);
    let within = m.nrtl().excess(350.0, &[0.5, 0.5]).unwrap();
    let below = m.nrtl().excess(320.0, &[0.5, 0.5]).unwrap();
    let pure = m.nrtl().excess(320.0, &[1.0, 0.0]).unwrap();
    assert_eq!(within.validity, ValidityStatus::WithinDeclaredDomain);
    assert_eq!(below.validity, ValidityStatus::Extrapolated);
    assert_eq!(pure.validity, ValidityStatus::WithinDeclaredDomain);
    let explicit = NrtlModel::new(ternary())
        .excess(330.0, &[0.2, 0.3, 0.5])
        .unwrap();
    assert_eq!(explicit.validity, ValidityStatus::Unassessed);
}

#[test]
fn package_caloric_properties_are_thermodynamically_consistent() {
    let m = package(&[0.3, 0.7]);
    let p = Pressure::new(101_325.0).unwrap();
    let t = 340.0;
    let h = 1e-3;
    let at = |t: f64| Temperature::new(t).unwrap();
    for liquid in [true, false] {
        let s = |t: f64| {
            if liquid {
                m.liquid(at(t), p)
            } else {
                m.vapor(at(t), p)
            }
            .unwrap()
        };
        let (up, mid, dn) = (s(t + h), s(t), s(t - h));
        let cp = (up.enthalpy - dn.enthalpy) / (2.0 * h);
        assert!((cp - mid.cp).abs() < 1e-5 * mid.cp, "cp liquid={liquid}");
        let ds = (up.entropy - dn.entropy) / (2.0 * h);
        assert!(
            (ds - mid.cp / t).abs() < 1e-7,
            "dS/dT = cp/T liquid={liquid}"
        );
        // Gibbs–Helmholtz: d(g/T)/dT = −h/T².
        let gh = (up.gibbs_energy / (t + h) - dn.gibbs_energy / (t - h)) / (2.0 * h);
        assert!(
            (gh + mid.enthalpy / (t * t)).abs() < 1e-8,
            "Gibbs–Helmholtz liquid={liquid}"
        );
    }
    // Liquid h, s are pressure independent; vapor entropy carries −R ln p.
    let p2 = Pressure::new(200_000.0).unwrap();
    let (l1, l2) = (m.liquid(at(t), p).unwrap(), m.liquid(at(t), p2).unwrap());
    assert_eq!(l1.enthalpy, l2.enthalpy);
    assert_eq!(l1.entropy, l2.entropy);
    let (v1, v2) = (m.vapor(at(t), p).unwrap(), m.vapor(at(t), p2).unwrap());
    assert!((v1.entropy - v2.entropy - R * (2e5f64 / 101_325.0).ln()).abs() < 1e-10);
}

#[test]
fn pure_liquid_standard_state_matches_saturated_vapor_gibbs_energy() {
    // At p = psat the pure liquid and ideal-gas vapor Gibbs energies coincide.
    let m = package(&[0.0, 1.0]);
    for t in [300.0, 373.0, 420.0] {
        let temperature = Temperature::new(t).unwrap();
        let ps = m.pure_liquid(1, temperature).unwrap().vapor_pressure;
        let p = Pressure::new(ps).unwrap();
        let (l, v) = (
            m.liquid(temperature, p).unwrap(),
            m.vapor(temperature, p).unwrap(),
        );
        assert!((l.gibbs_energy - v.gibbs_energy).abs() < 1e-8 * v.gibbs_energy.abs().max(1.0));
        assert!(l.ln_fugacity_coefficients[1].abs() < 1e-14);
        let latent = v.enthalpy - l.enthalpy;
        assert!(
            (latent
                - m.pure_liquid(1, temperature)
                    .unwrap()
                    .enthalpy_of_vaporization)
                .abs()
                < 1e-9
        );
    }
}

#[test]
fn package_rejects_unsupported_inputs_and_domains() {
    let m = package(&[0.5, 0.5]);
    let p = Pressure::new(101_325.0).unwrap();
    for t in [273.0, 450.1] {
        assert!(matches!(
            m.liquid(Temperature::new(t).unwrap(), p),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
    }
    for p in [99.0, 600_000.0] {
        assert!(matches!(
            m.vapor(Temperature::new(350.0).unwrap(), Pressure::new(p).unwrap()),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
    }
    let mass = Composition::new(
        CompositionBasis::Mass,
        mixture(&[0.5, 0.5]).entries().to_vec(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    assert!(matches!(
        ActivityModel::from_databases(
            &ComponentDatabase::embedded().unwrap(),
            &LiquidDatabase::embedded().unwrap(),
            &mass,
            VaporModel::IdealGas
        ),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    // Parameters in a different order than the composition are a mismatch, not reordered.
    let reversed = NrtlParameters::from_database(
        &LiquidDatabase::embedded().unwrap(),
        &ids(&["water", "methanol"]),
    )
    .unwrap();
    assert!(matches!(
        ActivityModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &LiquidDatabase::embedded().unwrap(),
            &mixture(&[0.5, 0.5]),
            VaporModel::IdealGas,
            reversed
        ),
        Err(CalculationError::InvalidInput(
            InputError::ComponentMismatch
        ))
    ));
    assert!(m.with_fractions(&[0.5, 0.6]).is_err());
    assert!(m.with_fractions(&[0.5]).is_err());
    let state = m.liquid(Temperature::new(350.0).unwrap(), p).unwrap();
    assert!(state.model.formulation.contains("no Poynting"));
    assert!(state.data.identifier.contains("chemsep"));
    assert!(state.excess.is_some());
}
