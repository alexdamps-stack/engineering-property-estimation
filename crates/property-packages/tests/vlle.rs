//! Three-phase flash through the `FlowsheetPackage` interface (#224).

use property_packages::{FlowsheetPackage, GergPackage, NrtlPackage, PhaseRole};
use thermo_data::{AlphaTerms, ComponentDatabase, LiquidDatabase, TauTerms};
use thermo_models::{ActivityModel, GergModel, NrtlParameters, VaporModel};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy,
};

fn composition(ids: &[&str], x: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(x)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

/// Synthetic partially miscible NRTL pair on the methanol/water pure data.
fn synthetic_nrtl_package() -> NrtlPackage {
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
    let model = ActivityModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(&["methanol", "water"], &[0.5, 0.5]),
        VaporModel::IdealGas,
        parameters,
    )
    .unwrap();
    NrtlPackage::new(model, NrtlPackage::TEMPERATURE_RANGE_K).unwrap()
}

#[test]
fn nrtl_package_reports_two_liquids_with_roles_and_metadata() {
    let package = synthetic_nrtl_package();
    let r = package.vlle_tp(&[0.5, 0.5], 300.0, 1e5).unwrap();
    assert_eq!(r.roles, [PhaseRole::LightLiquid, PhaseRole::HeavyLiquid]);
    assert_eq!(r.state.phases.len(), 2);
    let light = r.phase_with_role(PhaseRole::LightLiquid).unwrap();
    let heavy = r.phase_with_role(PhaseRole::HeavyLiquid).unwrap();
    assert!(light.state.mass_density < heavy.state.mass_density);
    assert!(r.state.convergence.max_ln_fugacity_residual <= 1e-10);
    assert!(!light.state.model.identifier.is_empty());
    // PH at the TP enthalpy reproduces the state.
    let ph = package
        .vlle_ph(&[0.5, 0.5], 1e5, r.state.enthalpy())
        .unwrap();
    assert!((ph.state.temperature_k - 300.0).abs() < 1e-6);
    assert_eq!(ph.roles, r.roles);
    // Outside the package interval the request is outside the domain.
    assert!(matches!(
        package.vlle_tp(&[0.5, 0.5], 500.0, 1e5),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}

#[test]
fn packages_without_the_capability_return_unsupported() {
    let model = GergModel::new(&composition(&["methane", "ethane"], &[0.9, 0.1])).unwrap();
    let package = GergPackage::new(model, [250.0, 300.0]).unwrap();
    for result in [
        package.vlle_tp(&[0.9, 0.1], 280.0, 50e5),
        package.vlle_ph(&[0.9, 0.1], 50e5, 0.0),
    ] {
        assert!(matches!(
            result,
            Err(CalculationError::UnsupportedCapability(m)) if m.contains("VLLE")
        ));
    }
}
