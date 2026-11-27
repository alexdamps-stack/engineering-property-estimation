//! Package operations beyond TP/PH (#226) through `FlowsheetPackage`:
//! PS, homogeneous phases, fugacity coefficients, K-values and bubble/dew
//! points for the cubic, NRTL-IG, GERG-2008 and EOS-CG-2021 packages.
//!
//! Expected values come from committed independent fixtures, never from Rust:
//! `nrtl-v1.json` (Clapeyron.jl, BigFloat modified Raoult) and
//! `gerg-equilibrium-v1.json` (NIST GERG C++ calorics, teqp, SciPy). PS targets
//! are the reference entropies of independently solved TP states and of pure
//! saturation endpoints (lever rule in entropy); the expected PS result is the
//! reference temperature and phase split. Budgets are the fixtures' own
//! declared tolerances, except where stated below.

use property_packages::{
    CubicPackage, EosCgPackage, FlowsheetPackage, GergPackage, NrtlPackage, RootSelection,
    SaturationKind,
};
use serde_json::Value;
use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_models::{ActivityModel, CubicModel, CubicVariant, GergModel, VaporModel};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase,
};

const NRTL: &str = include_str!("../../../validation/fixtures/nrtl-v1.json");
const GERG: &str = include_str!("../../../validation/fixtures/gerg-equilibrium-v1.json");

/// The A14 PS budgets of `validation/rotating/ps-v1.json`, fixed before any
/// comparison. The PS entropy tolerance (1e-10 scaled) limits the attainable
/// temperature precision to a few 1e-8 K, which the 1e-8 K NRTL PH temperature
/// budget does not cover, so NRTL PS uses the A14 temperature and fraction budgets.
const PS_TEMPERATURE: [f64; 2] = [2e-6, 0.0];
const PS_FRACTION: [f64; 2] = [2e-8, 0.0];

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(number).collect()
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap().to_owned())
        .collect()
}

fn budget(fixture: &Value, key: &str) -> [f64; 2] {
    let t = numbers(&fixture["tolerances"][key]);
    [t[0], t[1]]
}

fn close(what: &str, actual: f64, expected: f64, [atol, rtol]: [f64; 2]) {
    let allowed = atol + rtol * expected.abs();
    assert!(
        (actual - expected).abs() <= allowed,
        "{what}: {actual} vs reference {expected} (allowed {allowed})"
    );
}

fn composition(ids: &[String], z: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(z)
            .map(|(id, x)| (ComponentId::new(id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

fn nrtl() -> NrtlPackage {
    let ids = ["methanol".to_owned(), "water".to_owned()];
    let model = ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(&ids, &[0.5, 0.5]),
        VaporModel::IdealGas,
    )
    .unwrap();
    NrtlPackage::new(model, NrtlPackage::TEMPERATURE_RANGE_K).unwrap()
}

fn gerg(ids: &[String]) -> GergPackage {
    let mut z = vec![0.0; ids.len()];
    z[0] = 1.0;
    let model = GergModel::new(&composition(ids, &z)).unwrap();
    GergPackage::new(model, GergPackage::TEMPERATURE_RANGE_K).unwrap()
}

fn vapor_fraction(phases: &[property_packages::EquilibriumPhase]) -> f64 {
    phases
        .iter()
        .filter(|p| p.phase == Phase::Vapor)
        .map(|p| p.phase_fraction)
        .sum()
}

#[test]
fn nrtl_ps_reproduces_independent_tp_states_and_pure_saturation() {
    let fixture: Value = serde_json::from_str(NRTL).unwrap();
    let package = nrtl();
    for case in fixture["tp"].as_array().unwrap() {
        let (t, p, z) = (number(&case["T"]), number(&case["p"]), numbers(&case["z"]));
        let phases = case["phases"].as_array().unwrap();
        let entropy: f64 = phases
            .iter()
            .map(|ph| number(&ph["fraction"]) * number(&ph["state"]["s"]))
            .sum();
        let beta: f64 = phases
            .iter()
            .filter(|ph| ph["phase"] == "vapor")
            .map(|ph| number(&ph["fraction"]))
            .sum();
        let r = package.ps(&z, p, entropy).unwrap();
        let what = format!("NRTL PS T={t} p={p} z={z:?}");
        close(&format!("{what} T"), r.temperature_k, t, PS_TEMPERATURE);
        close(
            &format!("{what} beta"),
            vapor_fraction(&r.phases),
            beta,
            PS_FRACTION,
        );
        close(
            &format!("{what} s"),
            r.entropy(),
            entropy,
            budget(&fixture, "s"),
        );
        assert!(r.convergence.scaled_entropy_residual.unwrap().abs() <= 1e-10);
        // PH at the same state's enthalpy agrees (section 3.1 consistency).
        close(
            &format!("{what} h"),
            r.enthalpy(),
            number(&case["h"]),
            budget(&fixture, "h"),
        );
    }
    // Saturated pure liquid (lever-rule endpoint): s_liquid at p = psat(T).
    for pure in fixture["pure"].as_array().unwrap() {
        let t = number(&pure["T"]);
        let p = number(&pure["psat"]);
        if !(1e3..=1e6).contains(&p) || !(274.0..=449.0).contains(&t) {
            continue;
        }
        let z = if pure["id"] == "methanol" {
            [1.0, 0.0]
        } else {
            [0.0, 1.0]
        };
        let r = package.ps(&z, p, number(&pure["s_liquid"])).unwrap();
        close(
            &format!("NRTL pure {} PS T", pure["id"]),
            r.temperature_k,
            t,
            PS_TEMPERATURE,
        );
        assert!(vapor_fraction(&r.phases) <= PS_FRACTION[0]);
    }
}

#[test]
fn nrtl_saturation_matches_the_independent_modified_raoult_points() {
    let fixture: Value = serde_json::from_str(NRTL).unwrap();
    let package = nrtl();
    for s in fixture["saturation"].as_array().unwrap() {
        let kind = s["kind"].as_str().unwrap();
        let (t, p) = (number(&s["T"]), number(&s["p"]));
        let (x, y) = (numbers(&s["liquid"]), numbers(&s["vapor"]));
        let r = match kind {
            "bubble-pressure" => package.saturation_pressure(&x, t, SaturationKind::Bubble),
            "dew-pressure" => package.saturation_pressure(&y, t, SaturationKind::Dew),
            "bubble-temperature" => package.saturation_temperature(&x, p, SaturationKind::Bubble),
            "dew-temperature" => package.saturation_temperature(&y, p, SaturationKind::Dew),
            other => panic!("unknown kind {other}"),
        }
        .unwrap();
        let what = format!("NRTL {kind} T={t} p={p}");
        close(
            &format!("{what} T"),
            r.temperature_k,
            t,
            budget(&fixture, "temperature"),
        );
        close(
            &format!("{what} p"),
            r.pressure_pa,
            p,
            budget(&fixture, "pressure"),
        );
        for (a, e) in r.liquid_composition.iter().zip(&x) {
            close(
                &format!("{what} x"),
                *a,
                *e,
                budget(&fixture, "composition"),
            );
        }
        for (a, e) in r.vapor_composition.iter().zip(&y) {
            close(
                &format!("{what} y"),
                *a,
                *e,
                budget(&fixture, "composition"),
            );
        }
        // K-values from the package's own phases reproduce y/x.
        let k = package
            .k_values(&x, &y, r.temperature_k, r.pressure_pa)
            .unwrap();
        for i in 0..2 {
            if x[i] > 0.0 && y[i] > 0.0 {
                close(&format!("{what} K"), k[i], y[i] / x[i], [0.0, 1e-8]);
            }
        }
        assert!(r.max_ln_fugacity_residual <= 1e-10, "{what}");
    }
}

#[test]
fn nrtl_homogeneous_phases_are_the_gamma_phi_liquid_and_ideal_gas() {
    let package = nrtl();
    let z = [0.4, 0.6];
    let liquid = package
        .homogeneous(&z, 330.0, 101_325.0, RootSelection::SmallestVolume)
        .unwrap();
    let vapor = package
        .homogeneous(&z, 330.0, 101_325.0, RootSelection::LargestVolume)
        .unwrap();
    assert!(vapor.ln_fugacity_coefficients.iter().all(|v| *v == 0.0));
    assert!((vapor.compressibility - 1.0).abs() < 1e-12);
    assert!(liquid.molar_volume < vapor.molar_volume / 100.0);
    // Subcooled at 330 K and 1 atm: minimum Gibbs selects the liquid.
    let stable = package
        .ln_fugacity_coefficients(&z, 330.0, 101_325.0, RootSelection::MinimumGibbs)
        .unwrap();
    assert_eq!(stable, liquid.ln_fugacity_coefficients);
    // ln φ^L = ln γ + ln p^sat − ln p, so K = γ p^sat / p.
    let tp = package.tp(&z, 330.0, 101_325.0).unwrap();
    assert_eq!(tp.phases.len(), 1);
    assert_eq!(
        tp.phases[0].state.ln_fugacity_coefficients,
        liquid.ln_fugacity_coefficients
    );
    assert!(matches!(
        package.homogeneous(&z, 500.0, 101_325.0, RootSelection::SmallestVolume),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    assert!(matches!(
        package.saturation_pressure(&z, 460.0, SaturationKind::Bubble),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}

/// Phase-weighted reference entropy of a GERG TP case.
fn gerg_tp_entropy(case: &Value) -> f64 {
    if case["phases"] == 1 {
        number(&case["s"])
    } else {
        let beta = number(&case["vapor_fraction"]);
        beta * number(&case["vapor_properties"]["s"])
            + (1.0 - beta) * number(&case["liquid_properties"]["s"])
    }
}

/// Representative mixture states (vapor, liquid, binary and natural-gas splits,
/// narrow boiling); each GERG PS solve runs many sampled-root TP flashes.
const GERG_PS_CASES: [&str; 6] = [
    "binary-split",
    "binary-vapor",
    "binary-liquid",
    "ng-vapor",
    "ng-rich-split",
    "narrow-ethane-propane",
];

#[test]
fn gerg_ps_reproduces_independent_tp_states() {
    let fixture: Value = serde_json::from_str(GERG).unwrap();
    for case in fixture["tp"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        if !GERG_PS_CASES.contains(&id) {
            continue;
        }
        let ids = strings(&case["ids"]);
        let (t, p, z) = (number(&case["T"]), number(&case["p"]), numbers(&case["z"]));
        let package = gerg(&ids);
        let entropy = gerg_tp_entropy(case);
        let r = package
            .ps(&z, p, entropy)
            .unwrap_or_else(|e| panic!("{id}: {e}"));
        close(
            &format!("GERG PS {id} T"),
            r.temperature_k,
            t,
            budget(&fixture, "temperature"),
        );
        if case["phases"] == 2 {
            assert_eq!(r.phases.len(), 2, "{id}");
            close(
                &format!("GERG PS {id} beta"),
                r.phases[0].phase_fraction,
                number(&case["vapor_fraction"]),
                budget(&fixture, "vapor_fraction"),
            );
        } else {
            assert_eq!(r.phases.len(), 1, "{id}");
        }
        close(
            &format!("GERG PS {id} h"),
            r.enthalpy(),
            number(&case["enthalpy"]),
            budget(&fixture, "h"),
        );
    }
}

#[test]
fn gerg_pure_ps_and_saturation_reproduce_independent_coexistence() {
    let fixture: Value = serde_json::from_str(GERG).unwrap();
    let propane = gerg(&["propane".to_owned()]);
    for case in fixture["ph"].as_array().unwrap() {
        let Some(sat) = case.get("saturation") else {
            continue;
        };
        let id = case["id"].as_str().unwrap();
        let beta = number(&case["vapor_fraction"]);
        let entropy = beta * number(&sat["vapor_properties"]["s"])
            + (1.0 - beta) * number(&sat["liquid_properties"]["s"]);
        let r = propane.ps(&[1.0], number(&sat["p"]), entropy).unwrap();
        close(
            &format!("GERG PS {id} T"),
            r.temperature_k,
            number(&sat["T"]),
            budget(&fixture, "temperature"),
        );
        close(
            &format!("GERG PS {id} beta"),
            vapor_fraction(&r.phases),
            beta,
            budget(&fixture, "vapor_fraction"),
        );
        if beta != 0.25 {
            continue;
        }
        // Pure bubble pressure and dew temperature, once per reference point.
        let bubble = propane
            .saturation_pressure(&[1.0], number(&sat["T"]), SaturationKind::Bubble)
            .unwrap();
        close(
            &format!("GERG {id} psat"),
            bubble.pressure_pa,
            number(&sat["p"]),
            [0.0, 2e-9],
        );
        close(
            &format!("GERG {id} liquid s"),
            bubble.liquid.entropy,
            number(&sat["liquid_properties"]["s"]),
            budget(&fixture, "s"),
        );
        let dew = propane
            .saturation_temperature(&[1.0], number(&sat["p"]), SaturationKind::Dew)
            .unwrap();
        close(
            &format!("GERG {id} Tsat"),
            dew.temperature_k,
            number(&sat["T"]),
            budget(&fixture, "temperature"),
        );
    }
}

#[test]
fn gerg_homogeneous_and_fugacity_match_independent_tp_phases() {
    let fixture: Value = serde_json::from_str(GERG).unwrap();
    let case = fixture["tp"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "binary-split")
        .unwrap();
    let ids = strings(&case["ids"]);
    let package = gerg(&ids);
    let (t, p) = (number(&case["T"]), number(&case["p"]));
    let x = numbers(&case["liquid"]);
    let y = numbers(&case["vapor"]);
    let liquid = package
        .homogeneous(&x, t, p, RootSelection::SmallestVolume)
        .unwrap();
    let vapor = package
        .homogeneous(&y, t, p, RootSelection::LargestVolume)
        .unwrap();
    for (state, expected) in [
        (&liquid, &case["liquid_properties"]),
        (&vapor, &case["vapor_properties"]),
    ] {
        for (key, actual) in [
            ("h", state.enthalpy),
            ("s", state.entropy),
            ("cp", state.cp),
            ("Z", state.compressibility),
        ] {
            close(
                &format!("GERG {key}"),
                actual,
                number(&expected[key]),
                budget(&fixture, key),
            );
        }
        close(
            "GERG rho",
            1.0 / state.molar_volume,
            number(&expected["rho"]),
            budget(&fixture, "rho"),
        );
    }
    // At the reference split the package K-values reproduce y/x.
    let k = package.k_values(&x, &y, t, p).unwrap();
    for ((k, x), y) in k.iter().zip(&x).zip(&y) {
        close("GERG K", *k, y / x, [0.0, 1e-7]);
    }
    assert!(matches!(
        package.formation_components(),
        Err(CalculationError::UnsupportedCapability(_))
    ));
    assert!(matches!(
        package.homogeneous(&x, 400.0, p, RootSelection::SmallestVolume),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}

#[test]
fn gerg_mixture_saturation_bounds_the_independent_two_phase_states() {
    // The independent split at 250 K and 3 MPa lies inside the bubble/dew
    // envelope of its feed: bubble pressure above, dew pressure below 3 MPa.
    let fixture: Value = serde_json::from_str(GERG).unwrap();
    let case = fixture["tp"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "binary-split")
        .unwrap();
    let ids = strings(&case["ids"]);
    let package = gerg(&ids);
    let (t, p, z) = (number(&case["T"]), number(&case["p"]), numbers(&case["z"]));
    let bubble = package
        .saturation_pressure(&z, t, SaturationKind::Bubble)
        .unwrap();
    let dew = package
        .saturation_pressure(&z, t, SaturationKind::Dew)
        .unwrap();
    assert!(dew.pressure_pa < p && p < bubble.pressure_pa);
    // The reference liquid is at its own bubble point at (T, p), with the
    // reference vapor as incipient phase.
    let liquid = numbers(&case["liquid"]);
    let vapor = numbers(&case["vapor"]);
    let b = package
        .saturation_pressure(&liquid, t, SaturationKind::Bubble)
        .unwrap();
    close(
        "GERG bubble p of reference liquid",
        b.pressure_pa,
        p,
        [0.0, 1e-7],
    );
    for (a, e) in b.vapor_composition.iter().zip(&vapor) {
        close(
            "GERG incipient vapor",
            *a,
            *e,
            budget(&fixture, "composition"),
        );
    }
    let d = package
        .saturation_temperature(&vapor, p, SaturationKind::Dew)
        .unwrap();
    close(
        "GERG dew T of reference vapor",
        d.temperature_k,
        t,
        [2e-6, 0.0],
    );
    for (a, e) in d.liquid_composition.iter().zip(&liquid) {
        close(
            "GERG incipient liquid",
            *a,
            *e,
            budget(&fixture, "composition"),
        );
    }
}

fn eos_cg() -> EosCgPackage {
    EosCgPackage::new(
        vec!["carbon-dioxide".into(), "hydrogen".into()],
        EosCgPackage::TEMPERATURE_RANGE_K,
    )
    .unwrap()
}

#[test]
fn eos_cg_ps_round_trips_on_its_gas_branch_and_rejects_other_operations() {
    let package = eos_cg();
    let z = [0.3, 0.7];
    for (t, p) in [(335.0, 500_000.0), (348.0, 1_800_000.0)] {
        let tp = package.tp(&z, t, p).unwrap();
        let ps = package.ps(&z, p, tp.entropy()).unwrap();
        assert!((ps.temperature_k - t).abs() < 1e-6, "{t} {p}");
        assert!(ps.convergence.scaled_entropy_residual.unwrap().abs() <= 1e-10);
        let gas = package
            .homogeneous(&z, t, p, RootSelection::LargestVolume)
            .unwrap();
        assert_eq!(gas.enthalpy, tp.enthalpy());
        assert_eq!(
            package
                .ln_fugacity_coefficients(&z, t, p, RootSelection::LargestVolume)
                .unwrap(),
            gas.ln_fugacity_coefficients
        );
    }
    let unsupported = |r: Result<(), CalculationError>| {
        assert!(
            matches!(r, Err(CalculationError::UnsupportedCapability(_))),
            "{r:?}"
        );
    };
    unsupported(
        package
            .homogeneous(&z, 340.0, 1e6, RootSelection::SmallestVolume)
            .map(|_| ()),
    );
    unsupported(
        package
            .homogeneous(&z, 340.0, 1e6, RootSelection::MinimumGibbs)
            .map(|_| ()),
    );
    unsupported(package.k_values(&z, &z, 340.0, 1e6).map(|_| ()));
    unsupported(
        package
            .saturation_pressure(&z, 340.0, SaturationKind::Bubble)
            .map(|_| ()),
    );
    unsupported(
        package
            .saturation_temperature(&z, 1e6, SaturationKind::Dew)
            .map(|_| ()),
    );
    unsupported(package.saturation_temperature_interval().map(|_| ()));
    unsupported(package.formation_components().map(|_| ()));
    // Outside the gas-branch domain.
    assert!(matches!(
        package.ps(&z, 5e6, 0.0),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
}

#[test]
fn cubic_package_operations_delegate_to_the_cubic_model() {
    let ids = ["methane".to_owned(), "propane".to_owned()];
    let z = [0.5, 0.5];
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition(&ids, &z),
        CubicVariant::PengRobinsonBostonMathias,
        vec![vec![0.0; 2]; 2],
        "explicit zero interactions",
    )
    .unwrap();
    let package = CubicPackage::new(model, Default::default());
    let dynamic: &dyn FlowsheetPackage = &package;
    let tp = package.tp(&z, 250.0, 3e6).unwrap();
    let ps = dynamic.ps(&z, 3e6, tp.entropy()).unwrap();
    assert!((ps.temperature_k - 250.0).abs() < 1e-6);
    let bubble = dynamic
        .saturation_temperature(&z, 3e6, SaturationKind::Bubble)
        .unwrap();
    // The package operation searches its own interval and the package pressure
    // bounds (#227).
    let [t_lo, t_hi] = package.temperature_interval();
    let [p_lo, p_hi] = property_packages::CUBIC_PACKAGE_SATURATION_PRESSURE_RANGE_PA;
    let direct = package
        .saturation_temperature(
            &z,
            3e6,
            thermo_solvers::CubicSaturationKind::Bubble,
            thermo_solvers::SaturationOptions {
                minimum_temperature_k: t_lo,
                maximum_temperature_k: t_hi,
                minimum_pressure_pa: p_lo,
                maximum_pressure_pa: p_hi,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(bubble.temperature_k, direct.temperature_k);
    assert_eq!(bubble.vapor_composition, direct.vapor_composition);
    let k = dynamic
        .k_values(
            &bubble.liquid_composition,
            &bubble.vapor_composition,
            bubble.temperature_k,
            3e6,
        )
        .unwrap();
    for ((k, x), y) in k
        .iter()
        .zip(&bubble.liquid_composition)
        .zip(&bubble.vapor_composition)
    {
        close("cubic K", *k, y / x, [0.0, 1e-8]);
    }
    let state = dynamic
        .homogeneous(&z, 300.0, 1e6, RootSelection::MinimumGibbs)
        .unwrap();
    assert_eq!(
        state.ln_fugacity_coefficients,
        package.ln_fugacity_coefficients(&z, 300.0, 1e6).unwrap()
    );
    assert_eq!(dynamic.formation_components().unwrap().len(), 2);
}
