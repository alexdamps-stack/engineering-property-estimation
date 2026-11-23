use proptest::prelude::*;
use thermo_data::ComponentDatabase;
use thermo_models::*;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, MOLAR_GAS_CONSTANT as R,
    NormalizationPolicy, Pressure, Temperature,
};

fn mixture(ids: &[&str], x: &[f64], variant: CubicVariant) -> CubicModel {
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
        variant,
        vec![vec![0.0; ids.len()]; ids.len()],
        "explicit synthetic zero kij",
    )
    .unwrap()
}
fn eval(m: &CubicModel, t: f64, p: f64) -> HomogeneousState {
    m.evaluate(
        Temperature::new(t).unwrap(),
        Pressure::new(p).unwrap(),
        RootSelection::LargestVolume,
    )
    .unwrap()
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= atol + rtol * b.abs(),
        "{a:e} != {b:e}, delta={:e}",
        (a - b).abs()
    );
}

#[test]
fn derivatives_and_thermodynamic_identities_on_reference_grid() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        for (ids, x) in [
            (vec!["methane", "propane"], vec![0.5, 0.5]),
            (
                vec!["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
                vec![0.88, 0.06, 0.03, 0.02, 0.01],
            ),
        ] {
            let m = mixture(&ids, &x, variant);
            for t in [320., 360., 400.] {
                for p in [1e5, 1e6, 5e6] {
                    let s = eval(&m, t, p);
                    let v = s.molar_volume;
                    // Independent of the closed forms in `evaluate`: finite differences
                    // at constant V use the volumetric EOS to find p(T±dT, V).
                    let pressure = |temp: f64, vol: f64| {
                        m.pressure_derivatives(Temperature::new(temp).unwrap(), vol)
                            .unwrap()
                            .0
                    };
                    let (dt, dv) = (1e-3, v * 1e-6);
                    let pt_fd = (pressure(t + dt, v) - pressure(t - dt, v)) / (2. * dt);
                    let pv_fd = (pressure(t, v + dv) - pressure(t, v - dv)) / (2. * dv);
                    close(s.dp_dt, pt_fd, 1e-6, 1e-7);
                    close(s.dp_dv, pv_fd, 1e-3, 1e-7);
                    let at_volume = |temp: f64| {
                        let state = eval(&m, temp, pressure(temp, v));
                        close(state.molar_volume, v, 0., 1e-9);
                        state
                    };
                    let (up, down) = (at_volume(t + dt), at_volume(t - dt));
                    let cv_fd = (up.internal_energy - down.internal_energy) / (2. * dt);
                    close(s.cv, cv_fd, 1e-5, 1e-7);
                    close(s.cp - s.cv, -t * pt_fd * pt_fd / pv_fd, 1e-5, 1e-6);
                    let g: f64 = s
                        .ln_fugacity_coefficients
                        .iter()
                        .zip(&x)
                        .map(|(l, x)| l * x)
                        .sum();
                    close(
                        g,
                        (s.residual_enthalpy - t * s.residual_entropy) / (R * t),
                        1e-11,
                        1e-10,
                    );
                    for dt in [0.01, 0.003] {
                        let plus = eval(&m, t + dt, p);
                        let minus = eval(&m, t - dt, p);
                        close(
                            (plus.enthalpy - minus.enthalpy) / (2. * dt),
                            s.cp,
                            2e-6,
                            2e-7,
                        );
                        let dp = p * 1e-5;
                        let sp = eval(&m, t, p + dp);
                        let sm = eval(&m, t, p - dp);
                        close(
                            (sp.entropy - sm.entropy) / (2. * dp),
                            -(plus.molar_volume - minus.molar_volume) / (2. * dt),
                            1e-10,
                            1e-6,
                        );
                    }
                    // Numerical amount derivatives at fixed T,p: symmetry and Gibbs-Duhem.
                    for step in [1e-4, 3e-5] {
                        let mut jac = vec![vec![0.; x.len()]; x.len()];
                        for j in 0..x.len() {
                            let mut xp = x.clone();
                            let mut xm = x.clone();
                            xp[j] += step;
                            xm[j] -= step;
                            let plus = eval(&mixture(&ids, &xp, variant), t, p);
                            let minus = eval(&mixture(&ids, &xm, variant), t, p);
                            for (i, row) in jac.iter_mut().enumerate() {
                                row[j] = (plus.ln_fugacity_coefficients[i]
                                    - minus.ln_fugacity_coefficients[i])
                                    / (2. * step);
                            }
                        }
                        for j in 0..x.len() {
                            close(
                                jac.iter().zip(&x).map(|(row, x)| row[j] * x).sum(),
                                0.,
                                2e-7,
                                0.,
                            );
                            for (i, row) in jac.iter().enumerate() {
                                close(row[j], jac[j][i], 2e-7, 1e-6);
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn exact_factors_satisfy_all_pure_critical_constraints() {
    let db = ComponentDatabase::embedded().unwrap();
    for c in db.components() {
        for variant in [
            CubicVariant::PengRobinson,
            CubicVariant::PengRobinsonBostonMathias,
        ] {
            let model = mixture(&[c.id()], &[1.0], variant);
            let tc = c.critical_temperature();
            let pc = c.critical_pressure();
            let vc = 0.307_401_308_698_703_8 * R * tc / pc;
            let (p, first, second) = model
                .pressure_derivatives(Temperature::new(tc).unwrap(), vc)
                .unwrap();
            close(p, pc, 1e-6, 1e-12);
            close(first * vc / pc, 0.0, 1e-12, 0.0);
            close(second * vc * vc / pc, 0.0, 1e-11, 0.0);
        }
    }
}

#[test]
fn alpha_derivatives_and_bm_transition_are_explicit() {
    let db = ComponentDatabase::embedded().unwrap();
    for c in db.components() {
        let w = c.acentric_factor();
        let low = alpha_derivatives(CubicVariant::PengRobinsonBostonMathias, 1., w).unwrap();
        let high =
            alpha_derivatives(CubicVariant::PengRobinsonBostonMathias, 1. + 1e-10, w).unwrap();
        close(low.0, high.0, 1e-9, 1e-9);
        close(low.1, high.1, 1e-9, 1e-9);
        for variant in [
            CubicVariant::PengRobinson,
            CubicVariant::PengRobinsonBostonMathias,
        ] {
            for tr in [0.7, 0.99, 1.01, 1.5, 2.5] {
                let (a, d, d2) = alpha_derivatives(variant, tr, w).unwrap();
                for step in [1e-3, 3e-4] {
                    let plus = alpha_derivatives(variant, tr + step, w).unwrap();
                    let minus = alpha_derivatives(variant, tr - step, w).unwrap();
                    close((plus.0 - minus.0) / (2. * step), d, 2e-6, 2e-6);
                    close((plus.0 - 2. * a + minus.0) / (step * step), d2, 3e-6, 3e-6);
                }
            }
        }
    }
    let m = mixture(&["propane"], &[1.], CubicVariant::PengRobinsonBostonMathias);
    assert!(matches!(
        m.evaluate(
            Temperature::new(369.89).unwrap(),
            Pressure::new(1e5).unwrap(),
            RootSelection::LargestVolume
        ),
        Err(CalculationError::UnsupportedCapability(_))
    ));
}

#[test]
fn zero_trace_pure_limit_and_low_density() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        let pure = eval(&mixture(&["methane"], &[1.], variant), 400., 100.);
        let zero = eval(
            &mixture(&["methane", "propane"], &[1., 0.], variant),
            400.,
            100.,
        );
        let trace = eval(
            &mixture(&["methane", "propane"], &[1., 1e-14], variant),
            400.,
            100.,
        );
        close(zero.enthalpy, pure.enthalpy, 1e-10, 1e-12);
        close(trace.enthalpy, pure.enthalpy, 1e-8, 1e-12);
        close(pure.compressibility, 1., 2e-5, 0.);
        close(pure.ln_fugacity_coefficients[0], 0., 2e-5, 0.);
        assert_eq!(zero.ln_fugacity_coefficients.len(), 2);
    }
}

#[test]
fn invalid_missing_and_singular_requests_fail_explicitly() {
    let db = ComponentDatabase::embedded().unwrap();
    assert!(matches!(
        db.get("unknown"),
        Err(CalculationError::MissingParameter(_))
    ));
    let c = Composition::new(
        CompositionBasis::Mole,
        vec![(ComponentId::new("methane").unwrap(), 1.)],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    assert!(CubicModel::new(&db, &c, CubicVariant::PengRobinson, vec![], "missing").is_err());
    assert!(
        CubicModel::new(
            &db,
            &c,
            CubicVariant::PengRobinson,
            vec![vec![f64::NAN]],
            "invalid"
        )
        .is_err()
    );
    let m = mixture(&["propane"], &[1.], CubicVariant::PengRobinson);
    assert!(matches!(
        m.evaluate(
            Temperature::new(199.).unwrap(),
            Pressure::new(1e5).unwrap(),
            RootSelection::LargestVolume
        ),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        m.compressibility_roots(
            Temperature::new(369.89).unwrap(),
            Pressure::new(4251200.).unwrap()
        ),
        Err(CalculationError::NumericalFailure(_))
    ));
    assert!(
        m.residual_helmholtz(Temperature::new(300.).unwrap(), 1., f64::NAN)
            .is_err()
    );
    assert!(matches!(
        m.residual_helmholtz(Temperature::new(300.).unwrap(), f64::MAX, f64::MIN_POSITIVE),
        Err(CalculationError::NumericalFailure(_))
    ));
    for c in db.components() {
        assert!(c.ideal_gas(Temperature::new(199.).unwrap()).is_err());
        assert!(c.ideal_gas(Temperature::new(6001.).unwrap()).is_err());
        let lo = c
            .ideal_gas(Temperature::new(1000. - 1e-6).unwrap())
            .unwrap();
        let hi = c
            .ideal_gas(Temperature::new(1000. + 1e-6).unwrap())
            .unwrap();
        close(lo.enthalpy, hi.enthalpy, 1e-3, 0.);
        close(lo.entropy, hi.entropy, 1e-6, 0.);
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: std::env::var("PROPTEST_CASES").ok().and_then(|v| v.parse().ok()).unwrap_or(64), rng_seed: proptest::test_runner::RngSeed::Fixed(0xA1A2A3), .. ProptestConfig::default() })]
    #[test]
    fn permutation_scaling_normalization_and_euler(t in 420f64..550.,p in 1e5f64..3e6,a in 0.05f64..0.95,scale in 0.1f64..100.,bm in any::<bool>()) {
        let variant=if bm { CubicVariant::PengRobinsonBostonMathias } else { CubicVariant::PengRobinson };
        let m=mixture(&["methane","propane"],&[a,1.-a],variant);
        let permuted=mixture(&["propane","methane"],&[(1.-a)*scale,a*scale],variant);
        let s=eval(&m,t,p); let other=eval(&permuted,t,p);
        close(s.enthalpy,other.enthalpy,1e-8,1e-11); close(s.cp,other.cp,1e-8,1e-11);
        close(s.ln_fugacity_coefficients[0],other.ln_fugacity_coefficients[1],1e-10,1e-10);
        let temp=Temperature::new(t).unwrap();
        let energy=m.residual_helmholtz(temp,s.molar_volume,1.).unwrap();
        close(m.residual_helmholtz(temp,s.molar_volume*scale,scale).unwrap(),energy*scale,1e-7,1e-11);
        // p = RT/v − (∂A_res/∂V)_T,n from the independent Helmholtz entry point.
        let dv=s.molar_volume*1e-6;
        let slope=(m.residual_helmholtz(temp,s.molar_volume+dv,1.).unwrap()-m.residual_helmholtz(temp,s.molar_volume-dv,1.).unwrap())/(2.*dv);
        close(R*t/s.molar_volume-slope,p,1e-6*p,1e-6);
    }
}
