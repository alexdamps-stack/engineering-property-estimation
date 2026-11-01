use thermo_models::alpha;
use thermo_types::{CubicModel, Temperature};

fn value(model: CubicModel, t: f64) -> (f64, f64) {
    alpha(
        model,
        Temperature::new(t).unwrap(),
        Temperature::new(190.564).unwrap(),
        0.01142,
    )
    .unwrap()
}

#[test]
fn alpha_derivatives_match_central_differences_in_smooth_regions() {
    for model in [
        CubicModel::PengRobinson,
        CubicModel::PengRobinsonBostonMathias,
    ] {
        for t in [140.0, 300.0, 700.0] {
            for step in [1e-2, 1e-3] {
                let derivative =
                    (value(model, t + step).0 - value(model, t - step).0) / (2.0 * step);
                assert!((derivative - value(model, t).1).abs() < 1e-9);
            }
        }
        assert_eq!(value(model, 190.564).0, 1.0);
    }
}

#[test]
fn bm_first_derivative_is_continuous_at_critical_temperature() {
    let model = CubicModel::PengRobinsonBostonMathias;
    let t = 190.564;
    // BM is C1, not C2, at the join. Central differences crossing the join
    // have first-order truncation error; use second-order one-sided checks.
    for step in [0.01, 0.001] {
        let left = (3.0 * value(model, t).0 - 4.0 * value(model, t - step).0
            + value(model, t - 2.0 * step).0)
            / (2.0 * step);
        let right = (-3.0 * value(model, t).0 + 4.0 * value(model, t + step).0
            - value(model, t + 2.0 * step).0)
            / (2.0 * step);
        assert!((left - value(model, t).1).abs() < 1e-9);
        assert!((right - value(model, t).1).abs() < 1e-9);
    }
}

#[test]
fn models_share_subcritical_alpha_but_diverge_supercritically() {
    assert_eq!(
        value(CubicModel::PengRobinson, 140.0),
        value(CubicModel::PengRobinsonBostonMathias, 140.0)
    );
    assert!(
        (value(CubicModel::PengRobinson, 700.0).0
            - value(CubicModel::PengRobinsonBostonMathias, 700.0).0)
            .abs()
            > 0.01
    );
}
