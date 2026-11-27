//! #228 identity tests: attaching the COSTALD engineering liquid density to a
//! cubic package leaves every equilibrium, fugacity and caloric result
//! bit-for-bit unchanged, and the engineering density is reported separately
//! from the EOS density.
use property_packages::{CostaldLiquidDensity, CubicPackage, FlowsheetPackage, GergPackage};
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, GergModel};
use thermo_solvers::PhFlashOptions;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, InputError, NormalizationPolicy,
    Phase, ValidityStatus,
};

const IDS: [&str; 2] = ["ethane", "propane"];

fn composition(ids: &[&str]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .map(|id| (ComponentId::new(*id).unwrap(), 1.0))
            .collect(),
        NormalizationPolicy::Normalize,
    )
    .unwrap()
}

fn package(variant: CubicVariant) -> CubicPackage {
    let n = IDS.len();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition(&IDS),
        variant,
        vec![vec![0.0; n]; n],
        "explicit zero interaction for the identity test",
    )
    .unwrap();
    CubicPackage::new(
        model,
        PhFlashOptions {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 400.0,
            ..Default::default()
        },
    )
}

const VARIANTS: [CubicVariant; 3] = [
    CubicVariant::PengRobinson,
    CubicVariant::PengRobinsonBostonMathias,
    CubicVariant::SoaveRedlichKwong,
];

/// Liquid, two-phase, vapor and supercritical-pressure states.
const STATES: [(f64, f64); 5] = [
    (250.0, 5e6),
    (260.0, 1e6),
    (300.0, 1.2e6),
    (350.0, 1e5),
    (380.0, 8e6),
];

#[test]
fn engineering_density_leaves_equilibrium_and_caloric_results_unchanged() {
    let x = [0.4, 0.6];
    for variant in VARIANTS {
        let plain = package(variant);
        let with = package(variant)
            .with_engineering_liquid_density(CostaldLiquidDensity::embedded(&IDS).unwrap())
            .unwrap();
        for (t, p) in STATES {
            let a = plain.tp(&x, t, p).unwrap();
            let b = with.tp(&x, t, p).unwrap();
            // Debug output prints every f64 in round-trip form: equality is
            // bit-for-bit equality of all phase properties, including ln φ,
            // h, s, cp, cv, w, Z, molar volume and EOS density.
            assert_eq!(format!("{a:?}"), format!("{b:?}"), "{variant:?} TP {t} {p}");
            let h = a.enthalpy();
            let a = plain.ph(&x, p, h).unwrap();
            let b = with.ph(&x, p, h).unwrap();
            assert_eq!(format!("{a:?}"), format!("{b:?}"), "{variant:?} PH {t} {p}");
            let s = a.entropy();
            let a = plain.ps(&x, p, s).unwrap();
            let b = with.ps(&x, p, s).unwrap();
            assert_eq!(format!("{a:?}"), format!("{b:?}"), "{variant:?} PS {t} {p}");
            let a = plain.ln_fugacity_coefficients(&x, t, p).unwrap();
            let b = with.ln_fugacity_coefficients(&x, t, p).unwrap();
            assert_eq!(a, b);
        }
    }
}

#[test]
fn engineering_density_is_reported_beside_the_eos_density() {
    let correlation = CostaldLiquidDensity::embedded(&IDS).unwrap();
    for variant in VARIANTS {
        let pkg = package(variant)
            .with_engineering_liquid_density(correlation.clone())
            .unwrap();
        let state = FlowsheetPackage::tp(&pkg, &[0.4, 0.6], 250.0, 5e6).unwrap();
        let liquid = state
            .phases
            .iter()
            .find(|p| p.phase == Phase::Liquid)
            .expect("compressed liquid");
        let engineering = FlowsheetPackage::engineering_liquid_density(
            &pkg,
            &liquid.composition,
            state.temperature_k,
            state.pressure_pa,
        )
        .unwrap();
        let direct = correlation
            .evaluate(&liquid.composition, state.temperature_k, state.pressure_pa)
            .unwrap();
        assert_eq!(engineering, direct);
        assert_eq!(engineering.validity, ValidityStatus::WithinDeclaredDomain);
        assert!(engineering.pressure_corrected);
        // Untranslated cubic liquid density is not the engineering density;
        // both remain available and distinct.
        let relative = engineering.mass_density / liquid.state.mass_density - 1.0;
        assert!(relative.abs() > 1e-3, "{variant:?}: {relative}");
        // The phase's own EOS properties are those of the plain package.
        let plain = package(variant).tp(&[0.4, 0.6], 250.0, 5e6).unwrap();
        assert_eq!(
            plain.phases[0].state.mass_density.to_bits(),
            liquid.state.mass_density.to_bits()
        );
    }
}

#[test]
fn unsupported_and_mismatched_selections_fail_explicitly() {
    let plain = package(CubicVariant::PengRobinson);
    assert!(matches!(
        FlowsheetPackage::engineering_liquid_density(&plain, &[0.4, 0.6], 250.0, 5e6),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    let reversed = CostaldLiquidDensity::embedded(&["propane", "ethane"]).unwrap();
    assert!(matches!(
        package(CubicVariant::PengRobinson).with_engineering_liquid_density(reversed),
        Err(CalculationError::InvalidInput(
            InputError::ComponentMismatch
        ))
    ));
    let gerg = GergPackage::new(
        GergModel::new(&composition(&["methane", "ethane"])).unwrap(),
        [230.0, 350.0],
    )
    .unwrap();
    assert!(matches!(
        gerg.engineering_liquid_density(&[0.5, 0.5], 250.0, 5e6),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    // Outside the correlation domain: no extrapolated value.
    let with = package(CubicVariant::PengRobinson)
        .with_engineering_liquid_density(CostaldLiquidDensity::embedded(&IDS).unwrap())
        .unwrap();
    assert!(matches!(
        with.engineering_liquid_density(&[1.0, 0.0], 310.0, 6e6),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}
