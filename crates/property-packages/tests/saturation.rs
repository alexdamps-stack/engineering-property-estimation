use property_packages::CubicPackage;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::{
    CubicSaturationKind as Kind, EnvelopeStop, PhFlashOptions, SaturationOptions,
};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};
#[test]
fn saturation_uses_package_model_composition_and_temperature_limits() {
    let composition = Composition::new(
        CompositionBasis::Mole,
        vec![(ComponentId::new("propane").unwrap(), 1.)],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition,
        CubicVariant::PengRobinsonBostonMathias,
        vec![vec![0.]],
        "explicit pure interaction",
    )
    .unwrap();
    let package = CubicPackage::new(
        model,
        PhFlashOptions {
            minimum_temperature_k: 290.,
            maximum_temperature_k: 360.,
            ..Default::default()
        },
    );
    let bounds = SaturationOptions {
        minimum_temperature_k: 300.,
        maximum_temperature_k: 360.,
        ..Default::default()
    };
    let bubble = package
        .saturation_pressure(&[1.], 300., Kind::Bubble, bounds)
        .unwrap();
    let dew = package
        .saturation_temperature(&[1.], bubble.pressure_pa, Kind::Dew, bounds)
        .unwrap();
    assert!((dew.temperature_k - 300.).abs() < 2e-6);
    let trace = package
        .phase_envelope(&[1.], Kind::Bubble, 5., 20, bounds)
        .unwrap();
    assert!(matches!(trace.stop, EnvelopeStop::NearCritical { .. }));
    assert!(
        package
            .saturation_pressure(&[1.], 300., Kind::Bubble, SaturationOptions::default())
            .is_err()
    );
    assert!(
        package
            .saturation_pressure(&[0.], 300., Kind::Bubble, bounds)
            .is_err()
    );
}

#[test]
fn package_accepts_near_unit_pure_fractions_without_panicking() {
    let composition = Composition::new(
        CompositionBasis::Mole,
        vec![
            (ComponentId::new("methane").unwrap(), 0.),
            (ComponentId::new("propane").unwrap(), 1.),
        ],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition,
            variant,
            vec![vec![0.; 2]; 2],
            "explicit pure interactions",
        )
        .unwrap();
        let package = CubicPackage::new(model, PhFlashOptions::default());
        let bounds = SaturationOptions::default();
        for kind in [Kind::Bubble, Kind::Dew] {
            let exact_p = package
                .saturation_pressure(&[0., 1.], 300., kind, bounds)
                .unwrap();
            let exact_t = package
                .saturation_temperature(&[0., 1.], 1e6, kind, bounds)
                .unwrap();
            for amount in [1. - 5e-11, 1. + 5e-11] {
                let p = package
                    .saturation_pressure(&[0., amount], 300., kind, bounds)
                    .unwrap();
                let t = package
                    .saturation_temperature(&[0., amount], 1e6, kind, bounds)
                    .unwrap();
                assert_eq!(p.pressure_pa, exact_p.pressure_pa);
                assert_eq!(t.temperature_k, exact_t.temperature_k);
                assert_eq!(p.liquid_composition, vec![0., 1.]);
                assert_eq!(t.vapor_composition, vec![0., 1.]);
                let exhausted = SaturationOptions {
                    max_iterations: 1,
                    ..bounds
                };
                assert!(matches!(
                    package.saturation_temperature(&[0., amount], 1e6, kind, exhausted),
                    Err(thermo_types::CalculationError::NonConvergence(_))
                ));
            }
            for invalid in [0.99, 1.01, f64::NAN] {
                assert!(matches!(
                    package.saturation_pressure(&[0., invalid], 300., kind, bounds),
                    Err(thermo_types::CalculationError::InvalidInput(_))
                ));
                assert!(matches!(
                    package.saturation_temperature(&[0., invalid], 1e6, kind, bounds),
                    Err(thermo_types::CalculationError::InvalidInput(_))
                ));
            }
        }
    }
}
