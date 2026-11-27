//! Liquid–liquid extractor (#189) through the flowsheet engine: independent
//! NRTL references, single stage against the LLE flash, adiabatic energy
//! balance, a cubic package and explicit failure paths.

use flowsheet_core::units::multiphase::{
    ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorSpec,
};
use flowsheet_core::units::{MultiphaseSpec, columns::packing};
use flowsheet_core::*;
use property_packages::{CubicPackage, FlowsheetPackage, GergPackage, NrtlPackage, PhaseRole};
use serde_json::Value;
use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_models::{ActivityModel, CubicModel, CubicVariant, GergModel, VaporModel};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

const FIXTURE: &str = include_str!("../../../validation/multiphase/extractor-v1.json");
const IDS: [&str; 3] = ["water", "methanol", "1-butanol"];

fn composition(ids: &[&str], z: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(z)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::Normalize,
    )
    .unwrap()
}

fn nrtl() -> NrtlPackage {
    let model = ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(&IDS, &[1.0, 1.0, 1.0]),
        VaporModel::IdealGas,
    )
    .unwrap();
    NrtlPackage::new(model, NrtlPackage::TEMPERATURE_RANGE_K).unwrap()
}

/// Feed and solvent sources at (T, p, flow, z).
type Inlet<'a> = (f64, f64, f64, &'a [f64]);

struct Rig {
    flowsheet: Flowsheet,
    feeds: FeedCompositions,
    unit: UnitId,
    extract: StreamId,
    raffinate: StreamId,
}

fn rig(ids: &[&str], feed: Inlet, solvent: Inlet, spec: LiquidExtractorSpec) -> Rig {
    let mut flowsheet = Flowsheet::new();
    let unit = flowsheet.add_unit(UnitKind::LiquidExtractor);
    flowsheet
        .try_set_spec(
            unit,
            UnitSpec::Multiphase(MultiphaseSpec::LiquidExtractor(spec)),
        )
        .unwrap();
    let mut feeds = FeedCompositions::new();
    for (port, (t, p, flow, z)) in [feed, solvent].into_iter().enumerate() {
        let source = flowsheet.add_unit(UnitKind::Feed);
        assert!(flowsheet.set_spec(
            source,
            UnitSpec::Feed {
                temperature_k: Some(t),
                pressure_pa: Some(p),
                molar_flow_mol_s: Some(flow),
            }
        ));
        flowsheet
            .connect(PortRef::new(source, 0), PortRef::new(unit, port))
            .unwrap();
        feeds.insert(source, composition(ids, z));
    }
    let mut outlet = |port| {
        let product = flowsheet.add_unit(UnitKind::Product);
        flowsheet
            .connect(PortRef::new(unit, port), PortRef::new(product, 0))
            .unwrap()
    };
    let extract = outlet(2);
    let raffinate = outlet(3);
    Rig {
        flowsheet,
        feeds,
        unit,
        extract,
        raffinate,
    }
}

fn isothermal(stages: u16, t: f64, p: f64) -> LiquidExtractorSpec {
    LiquidExtractorSpec {
        internals: ExtractorInternals::Stages(stages),
        pressure_pa: Some(p),
        thermal: ExtractorThermal::Isothermal {
            temperature_k: Some(t),
        },
        feed_phase: FeedPhase::Detect,
    }
}

fn closed(calc: &Calculation, unit: UnitId) {
    let balance = calc.units[&unit].balance.as_ref().unwrap();
    assert!(
        balance.max_scaled_component_residual <= 1e-8,
        "{}",
        balance.max_scaled_component_residual
    );
    assert!(
        balance.scaled_energy_residual <= 2e-9,
        "{}",
        balance.scaled_energy_residual
    );
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn within(tolerance: &Value, actual: f64, expected: f64, what: &str) {
    let [atol, rtol] = [
        tolerance[0].as_f64().unwrap(),
        tolerance[1].as_f64().unwrap(),
    ];
    assert!(
        (actual - expected).abs() <= atol + rtol * expected.abs(),
        "{what}: {actual} vs {expected}"
    );
}

#[test]
fn nrtl_cascades_match_the_independent_references() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let package = nrtl();
    for case in fixture["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let (t, p) = (
            case["temperature_K"].as_f64().unwrap(),
            case["pressure_Pa"].as_f64().unwrap(),
        );
        let stages = case["stages"].as_u64().unwrap() as u16;
        let feed_z = numbers(&case["feed"]["z"]);
        let solvent_z = numbers(&case["solvent"]["z"]);
        let rig = rig(
            &IDS,
            (
                t,
                p,
                case["feed"]["molar_flow_mol_s"].as_f64().unwrap(),
                &feed_z,
            ),
            (
                t,
                p,
                case["solvent"]["molar_flow_mol_s"].as_f64().unwrap(),
                &solvent_z,
            ),
            isothermal(stages, t, p),
        );
        let calc = rig
            .flowsheet
            .solve(&package, &rig.feeds)
            .unwrap_or_else(|e| panic!("{id}: {e:?}"));
        closed(&calc, rig.unit);
        let profile = calc.units[&rig.unit].liquid_extractor.as_ref().unwrap();
        let tolerances = &fixture["tolerances"];
        let reference = case["profile"].as_array().expect("stage profile");
        assert_eq!(reference.len(), usize::from(stages), "{id}");
        for (j, stage) in reference.iter().enumerate() {
            let what = format!("{id} stage {}", j + 1);
            within(
                &tolerances["molar_flow_mol_s"],
                profile.extract_flows_mol_s[j],
                stage["extract_mol_s"].as_f64().unwrap(),
                &format!("{what} extract flow"),
            );
            within(
                &tolerances["molar_flow_mol_s"],
                profile.raffinate_flows_mol_s[j],
                stage["raffinate_mol_s"].as_f64().unwrap(),
                &format!("{what} raffinate flow"),
            );
            for (x, expected) in profile.extract_compositions[j]
                .iter()
                .zip(numbers(&stage["x_extract"]))
                .chain(
                    profile.raffinate_compositions[j]
                        .iter()
                        .zip(numbers(&stage["x_raffinate"])),
                )
            {
                within(&tolerances["composition"], *x, expected, &what);
            }
            assert_eq!(profile.stage_temperatures_k[j], t);
        }
        // Water-rich raffinate is the heavier liquid: the feed enters at the top.
        assert!(profile.feed_is_heavy, "{id}");
        let extract = &calc.streams[&rig.extract];
        let raffinate = &calc.streams[&rig.raffinate];
        assert!((extract.molar_flow_mol_s - profile.extract_flows_mol_s[0]).abs() < 1e-12);
        assert!(
            (raffinate.molar_flow_mol_s - profile.raffinate_flows_mol_s[usize::from(stages) - 1])
                .abs()
                < 1e-12
        );
        // The isothermal duty is the outlet-minus-inlet enthalpy flow.
        assert_eq!(calc.units[&rig.unit].duty_w, profile.duty_w);
    }
}

#[test]
fn single_stage_equals_the_lle_flash_of_the_mixed_inlets() {
    let package = nrtl();
    let (t, p) = (298.15, 101_325.0);
    let feed = [0.9, 0.1, 0.0];
    let solvent = [0.0, 0.0, 1.0];
    let rig = rig(
        &IDS,
        (t, p, 1.0, &feed),
        (t, p, 0.15, &solvent),
        isothermal(1, t, p),
    );
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    closed(&calc, rig.unit);
    let z: Vec<f64> = feed
        .iter()
        .zip(solvent)
        .map(|(a, b)| (a * 1.0 + b * 0.15) / 1.15)
        .collect();
    let flash = package.vlle_tp(&z, t, p).unwrap();
    // Butanol-rich light liquid = extract, water-rich heavy liquid = raffinate.
    let light = flash.phase_with_role(PhaseRole::LightLiquid).unwrap();
    let heavy = flash.phase_with_role(PhaseRole::HeavyLiquid).unwrap();
    for (stream, phase) in [(rig.extract, light), (rig.raffinate, heavy)] {
        let s = &calc.streams[&stream];
        let state = s.state.as_ref().unwrap();
        assert!((s.molar_flow_mol_s - 1.15 * phase.phase_fraction).abs() < 1e-8);
        for (a, b) in state.mole_fractions.iter().zip(&phase.composition) {
            assert!((a - b).abs() < 1e-8, "{a} vs {b}");
        }
        assert!((state.enthalpy() - phase.state.enthalpy).abs() < 1e-4);
    }
}

#[test]
fn adiabatic_cascade_closes_the_energy_balance_and_matches_a_ph_flash_on_one_stage() {
    let package = nrtl();
    let p = 150_000.0;
    let feed = [0.92, 0.08, 0.0];
    let solvent = [0.05, 0.0, 0.95];
    let adiabatic = |stages| LiquidExtractorSpec {
        internals: ExtractorInternals::Stages(stages),
        pressure_pa: Some(p),
        thermal: ExtractorThermal::Adiabatic,
        feed_phase: FeedPhase::Heavy,
    };
    let rig4 = rig(
        &IDS,
        (300.0, p, 1.0, &feed),
        (320.0, 2e5, 0.15, &solvent),
        adiabatic(4),
    );
    let calc = rig4.flowsheet.solve(&package, &rig4.feeds).unwrap();
    closed(&calc, rig4.unit);
    let unit = &calc.units[&rig4.unit];
    assert_eq!(unit.duty_w, 0.0);
    let profile = unit.liquid_extractor.as_ref().unwrap();
    for t in &profile.stage_temperatures_k {
        assert!((295.0..=325.0).contains(t), "{t}");
    }
    assert!(profile.stage_temperatures_k[0] != profile.stage_temperatures_k[3]);

    let rig1 = rig(
        &IDS,
        (300.0, p, 1.0, &feed),
        (320.0, 2e5, 0.15, &solvent),
        adiabatic(1),
    );
    let calc = rig1.flowsheet.solve(&package, &rig1.feeds).unwrap();
    closed(&calc, rig1.unit);
    let inlet_h: f64 = calc
        .streams
        .iter()
        .filter(|(id, _)| ![rig1.extract, rig1.raffinate].contains(id))
        .map(|(_, s)| s.enthalpy_flow_w())
        .sum();
    let z: Vec<f64> = feed
        .iter()
        .zip(solvent)
        .map(|(a, b)| (a * 1.0 + b * 0.15) / 1.15)
        .collect();
    let flash = package.vlle_ph(&z, p, inlet_h / 1.15).unwrap();
    let extract = calc.streams[&rig1.extract].state.as_ref().unwrap();
    assert!((extract.temperature_k - flash.state.temperature_k).abs() < 1e-6);
    let light = flash.phase_with_role(PhaseRole::LightLiquid).unwrap();
    for (a, b) in extract.mole_fractions.iter().zip(&light.composition) {
        assert!((a - b).abs() < 1e-7, "{a} vs {b}");
    }
}

#[test]
fn packed_section_resolves_its_stages_by_hetp() {
    let package = nrtl();
    let (t, p) = (298.15, 101_325.0);
    let mut spec = isothermal(1, t, p);
    spec.internals = ExtractorInternals::Packed(packing::PackedSectionSpec {
        hetp_m: Some(0.6),
        basis: packing::PackedSectionBasis::PackedHeight {
            height_m: Some(1.8),
            rounding: packing::StageRounding::Nearest,
        },
    });
    let rig = rig(
        &IDS,
        (t, p, 1.0, &[0.92, 0.08, 0.0]),
        (t, p, 0.15, &[0.0, 0.0, 1.0]),
        spec,
    );
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    closed(&calc, rig.unit);
    let profile = calc.units[&rig.unit].liquid_extractor.as_ref().unwrap();
    assert_eq!(profile.stage_temperatures_k.len(), 3);
    assert_eq!(profile.packed_height_m, Some(1.8));
}

#[test]
fn cubic_package_extracts_through_liquid_roots() {
    // Explicit water-propane kij 0.5 (approximation, as in the VLLE tests);
    // methanol pairs zero. Liquid propane at 300 K and 50 bar, water solvent.
    let ids = ["propane", "methanol", "water"];
    let kij = vec![
        vec![0.0, 0.0, 0.5],
        vec![0.0, 0.0, 0.0],
        vec![0.5, 0.0, 0.0],
    ];
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition(&ids, &[1.0, 1.0, 1.0]),
        CubicVariant::PengRobinson,
        kij,
        "explicit approximation for the extractor test",
    )
    .unwrap();
    let package = CubicPackage::new(model, Default::default());
    let (t, p) = (300.0, 50e5);
    for stages in [1, 4] {
        let rig = rig(
            &ids,
            (t, p, 1.0, &[0.9, 0.1, 0.0]),
            (t, p, 0.8, &[0.0, 0.0, 1.0]),
            isothermal(stages, t, p),
        );
        let calc = rig
            .flowsheet
            .solve(&package, &rig.feeds)
            .unwrap_or_else(|e| panic!("{stages}: {e:?}"));
        closed(&calc, rig.unit);
        if stages == 1 {
            let z = [0.9 / 1.8, 0.1 / 1.8, 0.8 / 1.8];
            let flash = package.vlle_tp(&z, t, p).unwrap();
            let raffinate = calc.streams[&rig.raffinate].state.as_ref().unwrap();
            let propane_rich = flash
                .state
                .phases
                .iter()
                .max_by(|a, b| a.composition[0].total_cmp(&b.composition[0]))
                .unwrap();
            for (a, b) in raffinate
                .mole_fractions
                .iter()
                .zip(&propane_rich.composition)
            {
                assert!((a - b).abs() < 1e-7, "{a} vs {b}");
            }
        }
    }
}

fn failure(package: &dyn FlowsheetPackage, rig: Rig) -> String {
    let error = rig.flowsheet.solve(package, &rig.feeds).unwrap_err();
    assert_eq!(error.unit, Some(rig.unit), "{error:?}");
    error.message
}

#[test]
fn failures_are_explicit() {
    let package = nrtl();
    let (t, p) = (298.15, 101_325.0);
    let water_methanol = [0.92, 0.08, 0.0];
    let butanol = [0.0, 0.0, 1.0];

    // Miscible: methanol as solvent gives one liquid.
    let message = failure(
        &package,
        rig(
            &IDS,
            (t, p, 1.0, &water_methanol),
            (t, p, 1.0, &[0.0, 1.0, 0.0]),
            isothermal(3, t, p),
        ),
    );
    assert!(message.contains("no liquid-liquid split"), "{message}");

    // Vapor forms at the column conditions.
    let message = failure(
        &package,
        rig(
            &IDS,
            (t, p, 1.0, &water_methanol),
            (t, p, 0.15, &butanol),
            isothermal(3, 365.0, 30_000.0),
        ),
    );
    assert!(message.contains("vapor forms"), "{message}");

    // A specified feed phase that contradicts the densities.
    let mut spec = isothermal(3, t, p);
    spec.feed_phase = FeedPhase::Light;
    let message = failure(
        &package,
        rig(
            &IDS,
            (t, p, 1.0, &water_methanol),
            (t, p, 0.15, &butanol),
            spec,
        ),
    );
    assert!(message.contains("specified light"), "{message}");

    // Column pressure above an inlet pressure.
    let message = failure(
        &package,
        rig(
            &IDS,
            (t, p, 1.0, &water_methanol),
            (t, p, 0.15, &butanol),
            isothermal(3, t, 2e5),
        ),
    );
    assert!(message.contains("pressure exceeds"), "{message}");

    // Zero solvent flow.
    let message = failure(
        &package,
        rig(
            &IDS,
            (t, p, 1.0, &water_methanol),
            (t, p, 0.0, &butanol),
            isothermal(3, t, p),
        ),
    );
    assert!(message.contains("solvent"), "{message}");

    // Temperature outside the package interval.
    let message = failure(
        &package,
        rig(
            &IDS,
            (t, p, 1.0, &water_methanol),
            (t, p, 0.15, &butanol),
            isothermal(3, 500.0, p),
        ),
    );
    assert!(
        message.contains("outside the package interval"),
        "{message}"
    );

    // Missing specification values are spec issues before any calculation.
    let rig_missing = rig(
        &IDS,
        (t, p, 1.0, &water_methanol),
        (t, p, 0.15, &butanol),
        LiquidExtractorSpec::default(),
    );
    assert!(
        rig_missing
            .flowsheet
            .solve(&package, &rig_missing.feeds)
            .is_err()
    );
}

#[test]
fn unsupported_package_fails_explicitly() {
    let ids = ["methane", "ethane"];
    let package = GergPackage::new(
        GergModel::new(&composition(&ids, &[0.9, 0.1])).unwrap(),
        GergPackage::TEMPERATURE_RANGE_K,
    )
    .unwrap();
    let message = failure(
        &package,
        rig(
            &ids,
            (250.0, 5e6, 1.0, &[0.9, 0.1]),
            (250.0, 5e6, 1.0, &[0.1, 0.9]),
            isothermal(2, 250.0, 5e6),
        ),
    );
    assert!(
        message.contains("VLLE") && message.contains("not supported"),
        "{message}"
    );
}

#[test]
fn liquid_extractor_is_calculated_and_other_multiphase_kinds_stay_drawn_only() {
    assert_eq!(
        UnitKind::LiquidExtractor.calculation_status(),
        CalculationStatus::Calculated
    );
    for kind in [
        UnitKind::ThreePhaseSeparator,
        UnitKind::Decanter,
        UnitKind::ThreePhaseColumn,
    ] {
        assert!(matches!(
            kind.calculation_status(),
            CalculationStatus::DrawnOnly(_)
        ));
    }
}
