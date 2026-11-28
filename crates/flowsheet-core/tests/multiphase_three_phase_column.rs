//! Three-phase distillation column (#188) through the flowsheet engine:
//! independent NRTL references, the decanter against the LLE flash of the
//! condensate, energy ports, packed internals and explicit failure paths;
//! columns with two liquids on trays and reboiler against independent
//! NRTL-IG and PR (water/propane) references.

use flowsheet_core::units::MultiphaseSpec;
use flowsheet_core::units::columns::ColumnInternals;
use flowsheet_core::units::columns::packing::{PackedSectionBasis, PackedSectionSpec};
use flowsheet_core::units::multiphase::{DecanterReflux, PhaseReflux, ThreePhaseColumnSpec};
use flowsheet_core::*;
use property_packages::{CubicPackage, FlowsheetPackage, GergPackage, NrtlPackage, PhaseRole};
use serde_json::Value;
use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_models::{ActivityModel, CubicModel, CubicVariant, GergModel, VaporModel};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

const FIXTURE: &str = include_str!("../../../validation/multiphase/three-phase-column-v1.json");
const IDS: [&str; 2] = ["water", "1-butanol"];

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

fn nrtl(ids: &[&str]) -> NrtlPackage {
    let model = ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(ids, &vec![1.0; ids.len()]),
        VaporModel::IdealGas,
    )
    .unwrap();
    NrtlPackage::new(model, NrtlPackage::TEMPERATURE_RANGE_K).unwrap()
}

/// Feed (T, p, flow, z).
type Feed<'a> = (f64, f64, f64, &'a [f64]);

struct Rig {
    flowsheet: Flowsheet,
    feeds: FeedCompositions,
    unit: UnitId,
    /// Distillate, decant, bottoms.
    outlets: Vec<StreamId>,
    condenser: StreamId,
    reboiler: StreamId,
}

fn rig(ids: &[&str], feed: Feed, spec: ThreePhaseColumnSpec) -> Rig {
    let mut flowsheet = Flowsheet::new();
    let unit = flowsheet.add_unit(UnitKind::ThreePhaseColumn);
    flowsheet
        .try_set_spec(
            unit,
            UnitSpec::Multiphase(MultiphaseSpec::ThreePhaseColumn(spec)),
        )
        .unwrap();
    let source = flowsheet.add_unit(UnitKind::Feed);
    let (t, p, flow, z) = feed;
    assert!(flowsheet.set_spec(
        source,
        UnitSpec::Feed {
            temperature_k: Some(t),
            pressure_pa: Some(p),
            molar_flow_mol_s: Some(flow),
        }
    ));
    flowsheet
        .connect(PortRef::new(source, 0), PortRef::new(unit, 0))
        .unwrap();
    let mut feeds = FeedCompositions::new();
    feeds.insert(source, composition(ids, z));
    let outlets = (1..=3)
        .map(|port| {
            let product = flowsheet.add_unit(UnitKind::Product);
            flowsheet
                .connect(PortRef::new(unit, port), PortRef::new(product, 0))
                .unwrap()
        })
        .collect();
    let sink = flowsheet.add_unit(UnitKind::EnergySink);
    let condenser = flowsheet
        .connect(PortRef::new(unit, 4), PortRef::new(sink, 0))
        .unwrap();
    let duty = flowsheet.add_unit(UnitKind::EnergySource);
    assert!(flowsheet.set_spec(
        duty,
        UnitSpec::EnergySource {
            duty: EnergySourceDuty::Calculated,
        }
    ));
    let reboiler = flowsheet
        .connect(PortRef::new(duty, 0), PortRef::new(unit, 5))
        .unwrap();
    Rig {
        flowsheet,
        feeds,
        unit,
        outlets,
        condenser,
        reboiler,
    }
}

fn base_spec() -> ThreePhaseColumnSpec {
    ThreePhaseColumnSpec {
        internals: ColumnInternals::Trays(6),
        feed_stage: 2,
        pressure_pa: Some(101_325.0),
        decanter_temperature_k: Some(320.0),
        reflux: DecanterReflux::Light(PhaseReflux::Total),
        bottoms_fraction: Some(0.71),
    }
}

const BUTANOL_FEED: Feed = (330.0, 2e5, 1.0, &[0.3, 0.7]);

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

fn reflux_of(r: &Value) -> DecanterReflux {
    let ratio = || PhaseReflux::Ratio(r["ratio"].as_f64());
    match r["kind"].as_str().unwrap() {
        "light_total" => DecanterReflux::Light(PhaseReflux::Total),
        "heavy_total" => DecanterReflux::Heavy(PhaseReflux::Total),
        "light_ratio" => DecanterReflux::Light(ratio()),
        "heavy_ratio" => DecanterReflux::Heavy(ratio()),
        "fractions" => DecanterReflux::Fractions {
            light: r["light"].as_f64(),
            heavy: r["heavy"].as_f64(),
        },
        other => panic!("unknown reflux kind {other}"),
    }
}

#[test]
fn nrtl_columns_match_the_independent_references() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let tolerances = &fixture["tolerances"];
    let package = nrtl(&IDS);
    for case in fixture["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let feed = &case["feed"];
        let z = numbers(&feed["z"]);
        let spec = ThreePhaseColumnSpec {
            internals: ColumnInternals::Trays(case["internal_stages"].as_u64().unwrap() as u16),
            feed_stage: case["feed_stage"].as_u64().unwrap() as u16,
            pressure_pa: case["pressure_Pa"].as_f64(),
            decanter_temperature_k: case["decanter_temperature_K"].as_f64(),
            reflux: reflux_of(&case["reflux"]),
            bottoms_fraction: case["bottoms_fraction"].as_f64(),
        };
        let rig = rig(
            &IDS,
            (
                feed["temperature_K"].as_f64().unwrap(),
                feed["pressure_Pa"].as_f64().unwrap(),
                feed["molar_flow_mol_s"].as_f64().unwrap(),
                &z,
            ),
            spec,
        );
        let calc = rig
            .flowsheet
            .solve(&package, &rig.feeds)
            .unwrap_or_else(|e| panic!("{id}: {e:?}"));
        closed(&calc, rig.unit);
        let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
        let t = &tolerances["temperature_K"];
        let x = &tolerances["composition"];
        let f = &tolerances["molar_flow_mol_s"];
        let q = &tolerances["duty_W"];
        let profile = case["profile"].as_array().unwrap();
        assert_eq!(r.stage_temperatures_k.len(), profile.len() + 1, "{id}");
        for stage in profile {
            let j = stage["stage"].as_u64().unwrap() as usize - 1;
            let what = format!("{id} stage {}", j + 1);
            within(
                t,
                r.stage_temperatures_k[j],
                stage["temperature_K"].as_f64().unwrap(),
                &what,
            );
            within(
                f,
                r.liquid_flows_mol_s[j],
                stage["liquid_mol_s"].as_f64().unwrap(),
                &format!("{what} L"),
            );
            within(
                f,
                r.vapor_flows_mol_s[j],
                stage["vapor_mol_s"].as_f64().unwrap(),
                &format!("{what} V"),
            );
            for (a, b) in r.liquid_compositions[j]
                .iter()
                .zip(numbers(&stage["x"]))
                .chain(r.vapor_compositions[j].iter().zip(numbers(&stage["y"])))
            {
                within(x, *a, b, &format!("{what} composition"));
            }
        }
        let decanter = &case["decanter"];
        within(
            f,
            r.light_liquid_mol_s,
            decanter["light_mol_s"].as_f64().unwrap(),
            &format!("{id} light"),
        );
        within(
            f,
            r.heavy_liquid_mol_s,
            decanter["heavy_mol_s"].as_f64().unwrap(),
            &format!("{id} heavy"),
        );
        for (a, b) in r
            .light_liquid_composition
            .iter()
            .zip(numbers(&decanter["x_light"]))
            .chain(
                r.heavy_liquid_composition
                    .iter()
                    .zip(numbers(&decanter["x_heavy"])),
            )
        {
            within(x, *a, b, &format!("{id} decanter composition"));
        }
        let qc = case["condenser_duty_W"].as_f64().unwrap();
        let qr = case["reboiler_duty_W"].as_f64().unwrap();
        within(q, r.condenser_duty_w, qc, &format!("{id} condenser"));
        within(q, r.reboiler_duty_w, qr, &format!("{id} reboiler"));
        // Energy ports carry the duties, positive in the port direction.
        within(
            q,
            calc.energy_streams[&rig.condenser].duty_w,
            qc,
            &format!("{id} condenser port"),
        );
        within(
            q,
            calc.energy_streams[&rig.reboiler].duty_w,
            qr,
            &format!("{id} reboiler port"),
        );
        for (stream, name) in rig.outlets.iter().zip(["distillate", "decant", "bottoms"]) {
            let product = &case["products"][name];
            let s = &calc.streams[stream];
            within(
                f,
                s.molar_flow_mol_s,
                product["molar_flow_mol_s"].as_f64().unwrap(),
                &format!("{id} {name} flow"),
            );
            for (a, b) in s
                .state
                .as_ref()
                .unwrap()
                .mole_fractions
                .iter()
                .zip(numbers(&product["x"]))
            {
                within(x, *a, b, &format!("{id} {name} composition"));
            }
        }
    }
}

#[test]
fn decanter_equals_the_liquid_liquid_flash_of_the_condensate() {
    let package = nrtl(&IDS);
    let spec = ThreePhaseColumnSpec {
        reflux: DecanterReflux::Fractions {
            light: Some(0.9),
            heavy: Some(0.2),
        },
        bottoms_fraction: Some(0.68),
        ..base_spec()
    };
    let rig = rig(&IDS, BUTANOL_FEED, spec);
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    closed(&calc, rig.unit);
    let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
    // The condensate is the vapor leaving stage 2, condensed completely.
    let condensate = &r.vapor_compositions[1];
    assert_eq!(&r.vapor_compositions[0], condensate);
    assert!((r.vapor_flows_mol_s[1] - r.light_liquid_mol_s - r.heavy_liquid_mol_s).abs() < 1e-12);
    let flash = package.vlle_tp(condensate, 320.0, 101_325.0).unwrap();
    assert_eq!(
        flash.roles,
        [PhaseRole::LightLiquid, PhaseRole::HeavyLiquid]
    );
    let light = flash.phase_with_role(PhaseRole::LightLiquid).unwrap();
    let heavy = flash.phase_with_role(PhaseRole::HeavyLiquid).unwrap();
    let total = r.light_liquid_mol_s + r.heavy_liquid_mol_s;
    assert!((light.phase_fraction - r.light_liquid_mol_s / total).abs() < 1e-8);
    for (a, b) in light
        .composition
        .iter()
        .zip(&r.light_liquid_composition)
        .chain(heavy.composition.iter().zip(&r.heavy_liquid_composition))
    {
        assert!((a - b).abs() < 1e-8, "{a} vs {b}");
    }
    // Reflux = 0.9 light + 0.2 heavy; the rest leaves as distillate and decant.
    let reflux = r.liquid_flows_mol_s[0];
    assert!((reflux - 0.9 * r.light_liquid_mol_s - 0.2 * r.heavy_liquid_mol_s).abs() < 1e-12);
    let distillate = &calc.streams[&rig.outlets[0]];
    let decant = &calc.streams[&rig.outlets[1]];
    assert!((distillate.molar_flow_mol_s - 0.1 * r.light_liquid_mol_s).abs() < 1e-12);
    assert!((decant.molar_flow_mol_s - 0.8 * r.heavy_liquid_mol_s).abs() < 1e-12);
    // Products leave at the decanter temperature as liquids.
    for s in [distillate, decant] {
        let state = s.state.as_ref().unwrap();
        assert_eq!(state.temperature_k, 320.0);
        assert_eq!(state.phases.len(), 1);
    }
    assert!(r.light_refluxed);
    assert_eq!(r.stage_temperatures_k[0], 320.0);
}

#[test]
fn heavy_phase_reflux_returns_the_aqueous_liquid() {
    let package = nrtl(&IDS);
    let spec = ThreePhaseColumnSpec {
        reflux: DecanterReflux::Heavy(PhaseReflux::Ratio(Some(9.0))),
        bottoms_fraction: Some(0.992),
        ..base_spec()
    };
    let rig = rig(&IDS, (330.0, 2e5, 1.0, &[0.996, 0.004]), spec);
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    closed(&calc, rig.unit);
    let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
    assert!(!r.light_refluxed);
    // Distillate: a tenth of the heavy (aqueous) liquid; decant: the light one.
    let distillate = &calc.streams[&rig.outlets[0]];
    let decant = &calc.streams[&rig.outlets[1]];
    assert!((distillate.molar_flow_mol_s - 0.1 * r.heavy_liquid_mol_s).abs() < 1e-12);
    assert!((decant.molar_flow_mol_s - r.light_liquid_mol_s).abs() < 1e-12);
    assert!(distillate.state.as_ref().unwrap().mole_fractions[0] > 0.98);
    assert!(decant.state.as_ref().unwrap().mole_fractions[1] > 0.3);
    // Bottoms are nearly pure water.
    let bottoms = calc.streams[&rig.outlets[2]].state.as_ref().unwrap();
    assert!(bottoms.mole_fractions[0] > 0.999);
}

#[test]
fn packed_internals_reproduce_the_equivalent_tray_column() {
    let package = nrtl(&IDS);
    let trays = rig(&IDS, BUTANOL_FEED, base_spec());
    let packed = rig(
        &IDS,
        BUTANOL_FEED,
        ThreePhaseColumnSpec {
            internals: ColumnInternals::Packed(PackedSectionSpec {
                hetp_m: Some(0.5),
                basis: PackedSectionBasis::TheoreticalStages(6),
            }),
            ..base_spec()
        },
    );
    let a = trays.flowsheet.solve(&package, &trays.feeds).unwrap();
    let b = packed.flowsheet.solve(&package, &packed.feeds).unwrap();
    let ra = a.units[&trays.unit].three_phase_column.as_ref().unwrap();
    let rb = b.units[&packed.unit].three_phase_column.as_ref().unwrap();
    assert_eq!(ra.packed_height_m, None);
    assert_eq!(rb.packed_height_m, Some(3.0));
    for (x, y) in ra.stage_temperatures_k.iter().zip(&rb.stage_temperatures_k) {
        assert!((x - y).abs() < 1e-9);
    }
}

fn failure(
    ids: &[&str],
    package: &dyn FlowsheetPackage,
    feed: Feed,
    spec: ThreePhaseColumnSpec,
) -> String {
    let rig = rig(ids, feed, spec);
    match rig.flowsheet.solve(package, &rig.feeds) {
        Ok(_) => panic!("the column should fail"),
        Err(e) => e.message,
    }
}

#[test]
fn failures_are_explicit() {
    let package = nrtl(&IDS);
    // Decanter above the condensate bubble point: vapor forms.
    let message = failure(
        &IDS,
        &package,
        BUTANOL_FEED,
        ThreePhaseColumnSpec {
            decanter_temperature_k: Some(375.0),
            ..base_spec()
        },
    );
    assert!(message.contains("vapor forms in the decanter"), "{message}");
    // A superheated vapor feed at the reboiler that must largely condense:
    // the reboiler would have to remove heat.
    let message = failure(
        &IDS,
        &package,
        (449.0, 1.2e5, 1.0, &[0.3, 0.7]),
        ThreePhaseColumnSpec {
            feed_stage: 8,
            bottoms_fraction: Some(0.8),
            ..base_spec()
        },
    );
    assert!(
        message.contains("heat input at the reboiler") && message.contains("infeasible"),
        "{message}"
    );
    // More aqueous decant than the feed water allows: no solution.
    let message = failure(
        &IDS,
        &package,
        BUTANOL_FEED,
        ThreePhaseColumnSpec {
            bottoms_fraction: Some(0.69),
            ..base_spec()
        },
    );
    assert!(
        message.contains("did not converge") || message.contains("stalled"),
        "{message}"
    );
    // Column pressure above the feed pressure.
    let message = failure(
        &IDS,
        &package,
        (330.0, 0.9e5, 1.0, &[0.3, 0.7]),
        base_spec(),
    );
    assert!(message.contains("exceeds the feed pressure"), "{message}");
    // Unset specification.
    let message = failure(
        &IDS,
        &package,
        BUTANOL_FEED,
        ThreePhaseColumnSpec {
            decanter_temperature_k: None,
            ..base_spec()
        },
    );
    assert!(message.contains("decanter temperature"), "{message}");
}

#[test]
fn a_miscible_condensate_is_an_explicit_failure() {
    // Methanol/water condenses to one liquid: use the ordinary column.
    let ids = ["water", "methanol", "1-butanol"];
    let package = nrtl(&ids);
    let message = failure(
        &ids,
        &package,
        (330.0, 2e5, 1.0, &[0.5, 0.5, 0.0]),
        ThreePhaseColumnSpec {
            bottoms_fraction: Some(0.5),
            ..base_spec()
        },
    );
    assert!(
        message.contains("no liquid-liquid split in the decanter")
            && message.contains("distillation column"),
        "{message}"
    );
}

#[test]
fn a_package_without_vlle_is_rejected() {
    let ids = ["ethane", "propane"];
    let package = GergPackage::new(
        GergModel::new(&composition(&ids, &[0.5, 0.5])).unwrap(),
        [250.0, 300.0],
    )
    .unwrap();
    let message = failure(
        &ids,
        &package,
        (280.0, 30e5, 1.0, &[0.5, 0.5]),
        ThreePhaseColumnSpec {
            pressure_pa: Some(20e5),
            decanter_temperature_k: Some(260.0),
            bottoms_fraction: Some(0.5),
            ..base_spec()
        },
    );
    assert!(message.contains("needs the VLLE/LLE flash"), "{message}");
}

#[test]
fn three_phase_column_is_calculated() {
    assert!(UnitKind::ThreePhaseColumn.calculation_planned());
    assert!(!UnitKind::ThreePhaseSeparator.calculation_planned());
    assert!(!UnitKind::Decanter.calculation_planned());
}

const VLLE_FIXTURE: &str =
    include_str!("../../../validation/multiphase/three-phase-column-vlle-v1.json");

/// PR water/propane with the explicit water-hydrocarbon kij = 0.5
/// (approximation of the #224 VLLE references).
fn water_propane_pr() -> CubicPackage {
    let ids = ["water", "propane"];
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition(&ids, &[1.0, 1.0]),
        CubicVariant::PengRobinson,
        vec![vec![0.0, 0.5], vec![0.5, 0.0]],
        "explicit water-hydrocarbon kij 0.5 (VLLE verification approximation)",
    )
    .unwrap();
    CubicPackage::new(model, Default::default())
}

fn fixture_spec(case: &Value) -> ThreePhaseColumnSpec {
    ThreePhaseColumnSpec {
        internals: ColumnInternals::Trays(case["internal_stages"].as_u64().unwrap() as u16),
        feed_stage: case["feed_stage"].as_u64().unwrap() as u16,
        pressure_pa: case["pressure_Pa"].as_f64(),
        decanter_temperature_k: case["decanter_temperature_K"].as_f64(),
        reflux: reflux_of(&case["reflux"]),
        bottoms_fraction: case["bottoms_fraction"].as_f64(),
    }
}

#[test]
fn vlle_stage_columns_match_the_independent_references() {
    let fixture: Value = serde_json::from_str(VLLE_FIXTURE).unwrap();
    let tolerances = &fixture["tolerances"];
    let (t, x, f, q) = (
        &tolerances["temperature_K"],
        &tolerances["composition"],
        &tolerances["molar_flow_mol_s"],
        &tolerances["duty_W"],
    );
    let nrtl = nrtl(&IDS);
    let pr = water_propane_pr();
    let mut packages_seen = [false; 2];
    for case in fixture["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let ids: Vec<&str> = case["components"]
            .as_array()
            .unwrap()
            .iter()
            .map(|c| c.as_str().unwrap())
            .collect();
        let package: &dyn FlowsheetPackage = match case["package"].as_str().unwrap() {
            "NRTL-IG" => {
                packages_seen[0] = true;
                &nrtl
            }
            "PR" => {
                packages_seen[1] = true;
                assert_eq!(numbers(&case["kij"][0]), [0.0, 0.5]);
                &pr
            }
            other => panic!("unknown package {other}"),
        };
        let feed = &case["feed"];
        let z = numbers(&feed["z"]);
        let rig = rig(
            &ids,
            (
                feed["temperature_K"].as_f64().unwrap(),
                feed["pressure_Pa"].as_f64().unwrap(),
                feed["molar_flow_mol_s"].as_f64().unwrap(),
                &z,
            ),
            fixture_spec(case),
        );
        let calc = rig
            .flowsheet
            .solve(package, &rig.feeds)
            .unwrap_or_else(|e| panic!("{id}: {e:?}"));
        closed(&calc, rig.unit);
        let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
        let profile = case["profile"].as_array().unwrap();
        assert_eq!(r.stage_phase_counts.len(), profile.len() + 1, "{id}");
        assert_eq!(r.stage_phase_counts[0], 2, "{id}");
        let mut three_phase = 0;
        for stage in profile {
            let j = stage["stage"].as_u64().unwrap() as usize - 1;
            let what = format!("{id} stage {}", j + 1);
            let count = stage["phase_count"].as_u64().unwrap() as u8;
            assert_eq!(r.stage_phase_counts[j], count, "{what} phase count");
            within(
                t,
                r.stage_temperatures_k[j],
                stage["temperature_K"].as_f64().unwrap(),
                &what,
            );
            within(
                f,
                r.liquid_flows_mol_s[j],
                stage["liquid_mol_s"].as_f64().unwrap(),
                &format!("{what} L"),
            );
            within(
                f,
                r.vapor_flows_mol_s[j],
                stage["vapor_mol_s"].as_f64().unwrap(),
                &format!("{what} V"),
            );
            for (a, b) in r.liquid_compositions[j]
                .iter()
                .zip(numbers(&stage["x"]))
                .chain(r.vapor_compositions[j].iter().zip(numbers(&stage["y"])))
            {
                within(x, *a, b, &format!("{what} composition"));
            }
            match (&r.stage_liquid_splits[j], count) {
                (None, 2) => {}
                (Some(split), 3) => {
                    three_phase += 1;
                    within(
                        f,
                        split.light_mol_s,
                        stage["light_mol_s"].as_f64().unwrap(),
                        &format!("{what} light"),
                    );
                    within(
                        f,
                        split.heavy_mol_s,
                        stage["heavy_mol_s"].as_f64().unwrap(),
                        &format!("{what} heavy"),
                    );
                    for (a, b) in split
                        .light_composition
                        .iter()
                        .zip(numbers(&stage["x_light"]))
                        .chain(
                            split
                                .heavy_composition
                                .iter()
                                .zip(numbers(&stage["x_heavy"])),
                        )
                    {
                        within(x, *a, b, &format!("{what} liquid compositions"));
                    }
                }
                (split, count) => panic!("{what}: {count} phases, split {split:?}"),
            }
        }
        // Every reference column carries two liquids on some stages and one
        // on others.
        assert!(three_phase > 0 && three_phase < profile.len(), "{id}");
        let decanter = &case["decanter"];
        within(
            f,
            r.light_liquid_mol_s,
            decanter["light_mol_s"].as_f64().unwrap(),
            &format!("{id} light"),
        );
        within(
            f,
            r.heavy_liquid_mol_s,
            decanter["heavy_mol_s"].as_f64().unwrap(),
            &format!("{id} heavy"),
        );
        for (a, b) in r
            .light_liquid_composition
            .iter()
            .zip(numbers(&decanter["x_light"]))
            .chain(
                r.heavy_liquid_composition
                    .iter()
                    .zip(numbers(&decanter["x_heavy"])),
            )
        {
            within(x, *a, b, &format!("{id} decanter composition"));
        }
        let qc = case["condenser_duty_W"].as_f64().unwrap();
        let qr = case["reboiler_duty_W"].as_f64().unwrap();
        within(q, r.condenser_duty_w, qc, &format!("{id} condenser"));
        within(q, r.reboiler_duty_w, qr, &format!("{id} reboiler"));
        within(
            q,
            calc.energy_streams[&rig.condenser].duty_w,
            qc,
            &format!("{id} condenser port"),
        );
        within(
            q,
            calc.energy_streams[&rig.reboiler].duty_w,
            qr,
            &format!("{id} reboiler port"),
        );
        for (stream, name) in rig.outlets.iter().zip(["distillate", "decant", "bottoms"]) {
            let product = &case["products"][name];
            let s = &calc.streams[stream];
            within(
                f,
                s.molar_flow_mol_s,
                product["molar_flow_mol_s"].as_f64().unwrap(),
                &format!("{id} {name} flow"),
            );
            for (a, b) in s
                .state
                .as_ref()
                .unwrap()
                .mole_fractions
                .iter()
                .zip(numbers(&product["x"]))
            {
                within(x, *a, b, &format!("{id} {name} composition"));
            }
        }
    }
    assert_eq!(packages_seen, [true, true]);
}

#[test]
fn a_column_without_internal_splits_keeps_single_liquid_stages() {
    // The decanter-only reference columns report one liquid on every stage
    // below the decanter, with the decanter split as stage 1.
    let package = nrtl(&IDS);
    let rig = rig(&IDS, BUTANOL_FEED, base_spec());
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
    assert_eq!(r.stage_phase_counts, [2; 8]);
    let decanter = r.stage_liquid_splits[0].as_ref().unwrap();
    assert_eq!(decanter.light_mol_s, r.light_liquid_mol_s);
    assert_eq!(decanter.heavy_composition, r.heavy_liquid_composition);
    assert!(r.stage_liquid_splits[1..].iter().all(Option::is_none));
}

#[test]
fn three_phase_stages_satisfy_both_equilibria_and_the_package_flash() {
    // Most of both liquids refluxed: two liquids on the upper trays at the
    // water/1-butanol heteroazeotropic temperature.
    let package = nrtl(&IDS);
    let spec = ThreePhaseColumnSpec {
        reflux: DecanterReflux::Fractions {
            light: Some(0.95),
            heavy: Some(0.9),
        },
        ..base_spec()
    };
    let rig = rig(&IDS, BUTANOL_FEED, spec);
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    closed(&calc, rig.unit);
    let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
    assert_eq!(r.stage_phase_counts[1], 3);
    let (p, t) = (101_325.0, r.stage_temperatures_k[1]);
    let split = r.stage_liquid_splits[1].as_ref().unwrap();
    let ln_f = |x: &[f64], root| {
        let phi = package.ln_fugacity_coefficients(x, t, p, root).unwrap();
        x.iter()
            .zip(phi)
            .map(|(x, phi)| x.ln() + phi)
            .collect::<Vec<_>>()
    };
    let light = ln_f(
        &split.light_composition,
        thermo_models::RootSelection::SmallestVolume,
    );
    let heavy = ln_f(
        &split.heavy_composition,
        thermo_models::RootSelection::SmallestVolume,
    );
    let vapor = ln_f(
        &r.vapor_compositions[1],
        thermo_models::RootSelection::LargestVolume,
    );
    for i in 0..2 {
        assert!((light[i] - heavy[i]).abs() < 1e-9);
        assert!((light[i] - vapor[i]).abs() < 1e-9);
    }
    // The total liquid is the sum of both liquids.
    let total = split.light_mol_s + split.heavy_mol_s;
    assert!((total - r.liquid_flows_mol_s[1]).abs() < 1e-12);
    for i in 0..2 {
        let x = (split.light_mol_s * split.light_composition[i]
            + split.heavy_mol_s * split.heavy_composition[i])
            / total;
        assert!((x - r.liquid_compositions[1][i]).abs() < 1e-12);
    }
    // The package's VLLE PH flash of the stage contents finds the same
    // three phases.
    let vapor_flow = r.vapor_flows_mol_s[1];
    let stage = total + vapor_flow;
    let z: Vec<f64> = (0..2)
        .map(|i| {
            (total * r.liquid_compositions[1][i] + vapor_flow * r.vapor_compositions[1][i]) / stage
        })
        .collect();
    let props = |x: &[f64], root| package.homogeneous(x, t, p, root).unwrap().enthalpy;
    let h = (split.light_mol_s
        * props(
            &split.light_composition,
            thermo_models::RootSelection::SmallestVolume,
        )
        + split.heavy_mol_s
            * props(
                &split.heavy_composition,
                thermo_models::RootSelection::SmallestVolume,
            )
        + vapor_flow
            * props(
                &r.vapor_compositions[1],
                thermo_models::RootSelection::LargestVolume,
            ))
        / stage;
    let flash = package.vlle_ph(&z, p, h).unwrap();
    assert_eq!(
        flash.roles,
        [
            PhaseRole::Vapor,
            PhaseRole::LightLiquid,
            PhaseRole::HeavyLiquid
        ]
    );
    assert!((flash.state.temperature_k - t).abs() < 1e-6);
    let heavy_phase = flash.phase_with_role(PhaseRole::HeavyLiquid).unwrap();
    assert!((heavy_phase.phase_fraction - split.heavy_mol_s / stage).abs() < 1e-5);
}

#[test]
fn nonconvergence_with_two_liquid_stages_names_them() {
    // Light-liquid reflux ratio 3 on the wet-propane column: the stage
    // equations with free water on the trays find no solution.
    let ids = ["water", "propane"];
    let message = failure(
        &ids,
        &water_propane_pr(),
        (400.0, 25e5, 1.0, &[0.05, 0.95]),
        ThreePhaseColumnSpec {
            feed_stage: 3,
            pressure_pa: Some(20e5),
            decanter_temperature_k: Some(300.0),
            reflux: DecanterReflux::Light(PhaseReflux::Ratio(Some(3.0))),
            bottoms_fraction: Some(0.9),
            ..base_spec()
        },
    );
    assert!(
        message.contains("three-phase column") && message.contains("two liquids on stages"),
        "{message}"
    );
}

#[test]
fn a_three_phase_reboiler_delivers_both_liquids_as_bottoms() {
    // Wet propane with more water left in the bottoms than propane
    // dissolves: free water at the reboiler.
    let ids = ["water", "propane"];
    let package = water_propane_pr();
    let rig = rig(
        &ids,
        (400.0, 25e5, 1.0, &[0.05, 0.95]),
        ThreePhaseColumnSpec {
            feed_stage: 3,
            pressure_pa: Some(20e5),
            decanter_temperature_k: Some(300.0),
            bottoms_fraction: Some(0.952),
            ..base_spec()
        },
    );
    let calc = rig.flowsheet.solve(&package, &rig.feeds).unwrap();
    closed(&calc, rig.unit);
    let r = calc.units[&rig.unit].three_phase_column.as_ref().unwrap();
    assert_eq!(*r.stage_phase_counts.last().unwrap(), 3);
    let split = r.stage_liquid_splits.last().unwrap().as_ref().unwrap();
    let bottoms = calc.streams[&rig.outlets[2]].state.as_ref().unwrap();
    assert_eq!(bottoms.phases.len(), 2);
    assert!(
        bottoms
            .phases
            .iter()
            .all(|p| p.phase == thermo_types::Phase::Liquid)
    );
    assert!((calc.streams[&rig.outlets[2]].molar_flow_mol_s - 0.952).abs() < 1e-12);
    let free_water = &bottoms.phases[1];
    assert!((free_water.phase_fraction - split.heavy_mol_s / 0.952).abs() < 1e-12);
    assert!(free_water.composition[0] > 0.999_999);
}
