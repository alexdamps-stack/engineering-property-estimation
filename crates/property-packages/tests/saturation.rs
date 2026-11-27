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
    // A TP state outside the package interval is a state outside the validity
    // domain, the same kind the homogeneous model reports below 200 K.
    for t in [280., 370.] {
        assert!(matches!(
            package.tp(&[1.], t, 1e5),
            Err(thermo_types::CalculationError::OutsideValidityDomain(_))
        ));
    }
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

/// #227: the package's bubble/dew bounds are its own temperature interval and
/// 1 kPa..10 MPa; water-bearing and three-phase boundaries are package results.
#[test]
fn package_bounds_cover_the_interval_wet_components_and_three_phase_bubbles() {
    use property_packages::{
        CUBIC_PACKAGE_SATURATION_PRESSURE_RANGE_PA, FlowsheetPackage, SaturationKind,
    };
    let composition = Composition::new(
        CompositionBasis::Mole,
        vec![
            (ComponentId::new("propane").unwrap(), 0.5),
            (ComponentId::new("water").unwrap(), 0.5),
        ],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition,
        CubicVariant::PengRobinson,
        vec![vec![0., 0.5], vec![0.5, 0.]],
        "explicit water/hydrocarbon kij approximation",
    )
    .unwrap();
    let package = CubicPackage::new(
        model,
        PhFlashOptions {
            minimum_temperature_k: 200.,
            maximum_temperature_k: 600.,
            ..Default::default()
        },
    );
    assert_eq!(
        FlowsheetPackage::saturation_temperature_interval(&package).unwrap(),
        [200., 600.]
    );
    assert_eq!(CUBIC_PACKAGE_SATURATION_PRESSURE_RANGE_PA, [1e3, 1e7]);
    // Pure water near ambient temperature (vapor pressure about 3 kPa).
    let water =
        FlowsheetPackage::saturation_pressure(&package, &[0., 1.], 300., SaturationKind::Bubble)
            .unwrap();
    assert!(water.pressure_pa > 2e3 && water.pressure_pa < 5e3);
    assert!(water.second_liquid.is_none());
    // Propane and water are immiscible liquids: the bubble point has three phases.
    let bubble =
        FlowsheetPackage::saturation_pressure(&package, &[0.5, 0.5], 320., SaturationKind::Bubble)
            .unwrap();
    let second = bubble.second_liquid.as_ref().unwrap();
    assert!(second.composition[1] > 0.99 && bubble.liquid_composition[0] > 0.9);
    assert!(second.material_fraction > 0.4 && second.material_fraction < 0.6);
    assert!(
        bubble
            .diagnostics
            .iter()
            .any(|d| d.code == "A13-THREE-PHASE-BUBBLE")
    );
    // Its dew point is an ordinary two-phase boundary with an aqueous liquid.
    let dew =
        FlowsheetPackage::saturation_pressure(&package, &[0.5, 0.5], 320., SaturationKind::Dew)
            .unwrap();
    assert!(dew.second_liquid.is_none());
    assert!(dew.liquid_composition[1] > 0.99);
    assert!(dew.pressure_pa < bubble.pressure_pa);
}
