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
