//! Cubic packages with petroleum pseudo-components (#229).
//!
//! The kerosene surrogate case compares two representations of the same
//! three-species cut with PR: the *surrogate* (defined compounds with
//! CoolProp HEOS critical constants, acentric factors and molar masses) and
//! the *pseudo-component* characterisation from each species' Tb and SG
//! alone. Species data, an independent evaluation of the correlations and the
//! budgets declared before the comparison are in
//! `validation/fixtures/pseudo-components-v1.json`. This is a model
//! comparison, not experimental validation.
use property_packages::{
    CubicPackage, EquilibriumState, FlowsheetPackage, PSEUDO_CUBIC_TEMPERATURE_RANGE_K,
    cubic_package_with_pseudo_components,
};
use serde_json::Value;
use thermo_data::{PseudoComponentInput, SpecifiedConstants};
use thermo_models::CubicVariant;
use thermo_solvers::{CubicSaturationKind, PhFlashOptions, SaturationOptions};
use thermo_types::{CalculationError, MOLAR_GAS_CONSTANT as R, Phase};

const FIXTURE: &str = include_str!("../../../validation/fixtures/pseudo-components-v1.json");
const KIJ_SOURCE: &str = "explicit approximation: kij = 0 for every hydrocarbon/hydrocarbon pseudo-component pair (Fahim et al. 2010, sec. 3.6.1); no fitted values";

fn fixture() -> Value {
    serde_json::from_str(FIXTURE).unwrap()
}

fn budget(name: &str) -> f64 {
    fixture()["tolerances"][name].as_f64().unwrap()
}

fn f(v: &Value, key: &str) -> f64 {
    v[key]
        .as_f64()
        .unwrap_or_else(|| panic!("{key} missing in {v}"))
}

/// (surrogate inputs, pseudo-component inputs, mole fractions)
fn representations() -> (
    Vec<PseudoComponentInput>,
    Vec<PseudoComponentInput>,
    Vec<f64>,
) {
    let fx = fixture();
    let species = fx["species"].as_array().unwrap();
    let surrogate = species
        .iter()
        .map(|s| PseudoComponentInput {
            id: s["id"].as_str().unwrap().into(),
            normal_boiling_point_k: f(s, "normal_boiling_point_K"),
            specific_gravity: f(s, "specific_gravity_60F"),
            molar_mass_kg_per_mol: Some(f(s, "molar_mass_kg_per_mol")),
            specified_constants: Some(SpecifiedConstants {
                critical_temperature_k: f(s, "critical_temperature_K"),
                critical_pressure_pa: f(s, "critical_pressure_Pa"),
                acentric_factor: f(s, "acentric_factor"),
                source: format!(
                    "CoolProp 7.2.0 HEOS {}",
                    s["coolprop_fluid"].as_str().unwrap()
                ),
            }),
        })
        .collect();
    let pseudo = species
        .iter()
        .map(|s| {
            PseudoComponentInput::new(
                s["id"].as_str().unwrap(),
                f(s, "normal_boiling_point_K"),
                f(s, "specific_gravity_60F"),
            )
        })
        .collect();
    let x = species
        .iter()
        .map(|s| f(s, "surrogate_mole_fraction"))
        .collect();
    (surrogate, pseudo, x)
}

fn options() -> PhFlashOptions {
    PhFlashOptions {
        minimum_temperature_k: PSEUDO_CUBIC_TEMPERATURE_RANGE_K[0],
        maximum_temperature_k: PSEUDO_CUBIC_TEMPERATURE_RANGE_K[1],
        ..Default::default()
    }
}

fn package(inputs: &[PseudoComponentInput]) -> CubicPackage {
    let ids: Vec<String> = inputs.iter().map(|p| p.id.clone()).collect();
    let n = ids.len();
    cubic_package_with_pseudo_components(
        CubicVariant::PengRobinson,
        &ids,
        inputs,
        vec![vec![0.0; n]; n],
        KIJ_SOURCE,
        options(),
    )
    .unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// Checks a deviation against its declared budget, or, for an exceedance
/// recorded in the fixture, requires that it still exceeds the budget and
/// reproduces the recorded value. `index` selects a temperature entry.
#[track_caller]
fn check_budget(species: &str, quantity: &str, index: Option<usize>, deviation: f64) {
    let fx = fixture();
    let limit = budget(quantity);
    let recorded = fx["recorded_exceedances"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["species"] == species && r["quantity"] == quantity);
    match recorded {
        None => assert!(
            deviation <= limit,
            "{species} {quantity}[{index:?}]: {deviation}"
        ),
        Some(r) => {
            let observed = match index {
                Some(i) => r["observed"][i].as_f64().unwrap(),
                None => r["observed"].as_f64().unwrap(),
            };
            assert!(
                deviation > limit,
                "{species} {quantity}: no longer exceeds; remove the record"
            );
            assert!(
                (deviation - observed).abs() <= budget("recorded_exceedance_reproduction_absolute"),
                "{species} {quantity}[{index:?}]: {deviation} vs recorded {observed}"
            );
        }
    }
}

#[test]
fn correlations_agree_with_the_independent_python_evaluation() {
    let tol = budget("correlation_relative");
    let (_, pseudo, _) = representations();
    let fx = fixture();
    for (input, s) in pseudo.iter().zip(fx["species"].as_array().unwrap()) {
        let c = input.characterise().unwrap();
        let r = &s["correlations"];
        for (actual, key) in [
            (c.watson_k, "watson_k"),
            (c.critical_temperature_k, "critical_temperature_K"),
            (c.critical_pressure_pa, "critical_pressure_Pa"),
            (c.acentric_factor, "acentric_factor"),
            (c.molar_mass_kg_per_mol, "molar_mass_kg_per_mol"),
        ] {
            assert!(rel(actual, f(r, key)) <= tol, "{} {key}", input.id);
        }
        for (actual, expected) in c
            .ideal_gas_cp_over_r
            .iter()
            .zip(r["ideal_gas_cp_over_R"].as_array().unwrap())
        {
            assert!(
                rel(*actual, expected.as_f64().unwrap()) <= tol,
                "{}",
                input.id
            );
        }
    }
}

/// Model comparison of the characterised constants and Kesler–Lee cp with
/// the defined compounds (CoolProp HEOS).
#[test]
fn characterised_constants_and_ideal_gas_cp_against_the_defined_compounds() {
    let (surrogate, pseudo, _) = representations();
    let fx = fixture();
    let temperatures: Vec<f64> = fx["ideal_gas_cp_temperatures_K"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t.as_f64().unwrap())
        .collect();
    for ((s_in, p_in), s) in surrogate
        .iter()
        .zip(&pseudo)
        .zip(fx["species"].as_array().unwrap())
    {
        let p = p_in.characterise().unwrap();
        let id = &p_in.id;
        let deviations = [
            (
                rel(p.critical_temperature_k, f(s, "critical_temperature_K")),
                "critical_temperature_vs_surrogate_relative",
            ),
            (
                rel(p.critical_pressure_pa, f(s, "critical_pressure_Pa")),
                "critical_pressure_vs_surrogate_relative",
            ),
            (
                rel(p.acentric_factor, f(s, "acentric_factor")),
                "acentric_factor_vs_surrogate_relative",
            ),
            (
                rel(p.molar_mass_kg_per_mol, f(s, "molar_mass_kg_per_mol")),
                "molar_mass_vs_surrogate_relative",
            ),
        ];
        for (deviation, name) in deviations {
            println!("{id} {name}: {deviation:.4}");
            check_budget(id, name, None, deviation);
        }
        let cp_budget = budget("kesler_lee_cp0_vs_coolprop_relative");
        // The pseudo representation is subject to recorded exceedances; the
        // surrogate representation (HEOS omega and M) must meet the budget.
        for (pseudo_representation, c) in [(true, p), (false, s_in.characterise().unwrap())] {
            for (i, (t, reference)) in temperatures
                .iter()
                .zip(s["ideal_gas_cp_J_per_mol_K"].as_array().unwrap())
                .enumerate()
            {
                let [c0, c1, c2] = c.ideal_gas_cp_over_r;
                let cp = R * (c0 + t * (c1 + t * c2));
                let deviation = rel(cp, reference.as_f64().unwrap());
                println!("{id} cp0 at {t} K (pseudo {pseudo_representation}): {deviation:.4}");
                if pseudo_representation {
                    check_budget(
                        id,
                        "kesler_lee_cp0_vs_coolprop_relative",
                        Some(i),
                        deviation,
                    );
                } else {
                    assert!(deviation <= cp_budget, "{id} cp0 at {t} K: {deviation}");
                }
            }
        }
    }
}

fn vapor_fraction(state: &EquilibriumState) -> f64 {
    state
        .phases
        .iter()
        .filter(|p| p.phase == Phase::Vapor)
        .map(|p| p.phase_fraction)
        .sum()
}

/// Temperature at which the 1 atm TP flash changes between `below` and not
/// `below`, by bisection to 1e-4 K.
fn transition(
    package: &CubicPackage,
    x: &[f64],
    mut lo: f64,
    mut hi: f64,
    below: impl Fn(f64) -> bool,
) -> f64 {
    let p = 101_325.0;
    let at = |t: f64| below(vapor_fraction(&package.tp(x, t, p).unwrap()));
    assert!(at(lo) && !at(hi));
    while hi - lo > 1e-4 {
        let mid = 0.5 * (lo + hi);
        if at(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

struct CaseResult {
    bubble_k: f64,
    dew_k: f64,
    liquid_density: f64,
    specific_enthalpy_change: f64,
}

fn evaluate(package: &CubicPackage, inputs: &[PseudoComponentInput], x: &[f64]) -> CaseResult {
    let molar_mass: f64 = inputs
        .iter()
        .zip(x)
        .map(|(c, x)| x * c.characterise().unwrap().molar_mass_kg_per_mol)
        .sum();
    let bubble_k = transition(package, x, 350.0, 560.0, |v| v == 0.0);
    let dew_k = transition(package, x, 350.0, 560.0, |v| v < 1.0);
    let liquid = package.tp(x, 300.0, 101_325.0).unwrap();
    assert_eq!(liquid.phases.len(), 1);
    assert_eq!(liquid.phases[0].phase, Phase::Liquid);
    let vapor = package.tp(x, 560.0, 101_325.0).unwrap();
    assert_eq!(vapor_fraction(&vapor), 1.0);
    CaseResult {
        bubble_k,
        dew_k,
        liquid_density: liquid.phases[0].state.mass_density,
        specific_enthalpy_change: (vapor.enthalpy() - liquid.enthalpy()) / molar_mass,
    }
}

#[test]
fn kerosene_surrogate_against_its_pseudo_component_characterisation() {
    let (surrogate_inputs, pseudo_inputs, x) = representations();
    let surrogate = evaluate(&package(&surrogate_inputs), &surrogate_inputs, &x);
    let pseudo = evaluate(&package(&pseudo_inputs), &pseudo_inputs, &x);
    assert!(surrogate.bubble_k < surrogate.dew_k && pseudo.bubble_k < pseudo.dew_k);
    let dt = budget("bubble_dew_temperature_vs_surrogate_K");
    println!(
        "bubble {:.3} vs {:.3} K, dew {:.3} vs {:.3} K, liquid density {:.2} vs {:.2} kg/m3, dh {:.1} vs {:.1} J/kg",
        pseudo.bubble_k,
        surrogate.bubble_k,
        pseudo.dew_k,
        surrogate.dew_k,
        pseudo.liquid_density,
        surrogate.liquid_density,
        pseudo.specific_enthalpy_change,
        surrogate.specific_enthalpy_change
    );
    assert!((pseudo.bubble_k - surrogate.bubble_k).abs() <= dt);
    assert!((pseudo.dew_k - surrogate.dew_k).abs() <= dt);
    assert!(
        rel(pseudo.liquid_density, surrogate.liquid_density)
            <= budget("liquid_mass_density_vs_surrogate_relative")
    );
    assert!(
        rel(
            pseudo.specific_enthalpy_change,
            surrogate.specific_enthalpy_change
        ) <= budget("specific_enthalpy_change_vs_surrogate_relative")
    );
}

/// PH and PS reproduce a two-phase TP state of the pseudo-component package,
/// and the homogeneous cp equals the temperature derivative of h.
#[test]
fn ph_ps_round_trips_and_cp_identity_with_pseudo_components() {
    let (_, inputs, x) = representations();
    let package = package(&inputs);
    let tp = package.tp(&x, 455.0, 101_325.0).unwrap();
    assert_eq!(tp.phases.len(), 2);
    let ph = package.ph(&x, 101_325.0, tp.enthalpy()).unwrap();
    assert!(
        (ph.temperature_k - 455.0).abs() < 1e-6,
        "{}",
        ph.temperature_k
    );
    let ps = package.ps(&x, 101_325.0, tp.entropy()).unwrap();
    assert!(
        (ps.temperature_k - 455.0).abs() < 1e-6,
        "{}",
        ps.temperature_k
    );
    for (t, p) in [(300.0, 101_325.0), (560.0, 101_325.0), (500.0, 20_000.0)] {
        let h = |t: f64| package.tp(&x, t, p).unwrap().enthalpy();
        let step = 1e-3;
        let numeric = (h(t + step) - h(t - step)) / (2.0 * step);
        let state = package.tp(&x, t, p).unwrap();
        assert_eq!(state.phases.len(), 1);
        let cp = state.phases[0].state.cp;
        assert!(rel(cp, numeric) < 1e-6, "{t} K: {cp} vs {numeric}");
    }
    // Pseudo-components carry the characterisation version in data provenance.
    let data = &tp.phases[0].state.data;
    assert!(data.version.contains("+pseudo-"), "{}", data.version);
    assert!(data.source.contains("kesler-lee"), "{}", data.source);
}

#[test]
fn pseudo_component_packages_fail_explicitly() {
    let (_, inputs, x) = representations();
    let ids: Vec<String> = inputs.iter().map(|p| p.id.clone()).collect();
    let zero = |n: usize| vec![vec![0.0; n]; n];
    let build = |ids: &[String], inputs: &[PseudoComponentInput], kij, source: &str, options| {
        cubic_package_with_pseudo_components(
            CubicVariant::PengRobinson,
            ids,
            inputs,
            kij,
            source,
            options,
        )
    };
    // A selected pseudo-component without a definition.
    assert!(matches!(
        build(&ids, &inputs[..2], zero(3), KIJ_SOURCE, options()),
        Err(CalculationError::MissingParameter(_))
    ));
    // A definition that is not selected.
    assert!(matches!(
        build(&ids[..2], &inputs, zero(2), KIJ_SOURCE, options()),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    // A pseudo-component colliding with another definition.
    let mut duplicate = inputs.clone();
    duplicate[1].id = duplicate[0].id.clone();
    assert!(build(&ids, &duplicate, zero(3), KIJ_SOURCE, options()).is_err());
    // Interactions are never filled in: the matrix must cover every pair.
    assert!(matches!(
        build(&ids, &inputs, zero(2), KIJ_SOURCE, options()),
        Err(CalculationError::MissingParameter(_))
    ));
    assert!(matches!(
        build(&ids, &inputs, zero(3), " ", options()),
        Err(CalculationError::MissingParameter(_))
    ));
    // The search interval must stay inside the Kesler–Lee ideal-gas domain.
    let mut wide = options();
    wide.minimum_temperature_k = 200.0;
    assert!(matches!(
        build(&ids, &inputs, zero(3), KIJ_SOURCE, wide),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    // Bubble/dew solvers are qualified for the A13 dry components only.
    let package = package(&inputs);
    assert!(matches!(
        package.saturation_temperature(
            &x,
            101_325.0,
            CubicSaturationKind::Bubble,
            SaturationOptions {
                minimum_temperature_k: 300.0,
                maximum_temperature_k: 350.0,
                ..Default::default()
            }
        ),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    // A pseudo-component next to database components in one package.
    let mut mixed = vec!["methane".to_string()];
    mixed.extend(ids.iter().cloned());
    let mixed = build(&mixed, &inputs, zero(4), KIJ_SOURCE, options()).unwrap();
    assert_eq!(mixed.component_ids()[0], "methane");
    let state = FlowsheetPackage::tp(&mixed, &[0.1, 0.18, 0.36, 0.36], 400.0, 5e5).unwrap();
    assert!(state.phases.iter().all(|p| p.state.enthalpy.is_finite()));
}
