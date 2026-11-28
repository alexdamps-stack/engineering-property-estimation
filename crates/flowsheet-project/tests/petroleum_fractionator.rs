//! Petroleum fractionator (#190): crude-like pseudo-component columns with
//! two or three side products (molar and standard-liquid-volume flows), a
//! reboiler or stripping steam, against independent references; explicit
//! failures; project persistence and the drawn-only legacy record.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{
    BottomSpecification, ColumnInternals, ColumnPressureDrop, DistillationCondenser, DrawRate,
    FractionatorBottom, PetroleumFractionatorSpec, SideProduct, TopSpecification,
    WATER_DENSITY_60F_KG_M3,
};
use flowsheet_core::{EnergySourceDuty, PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    PackageChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};
use thermo_data::{PseudoComponentInput, cut_mole_fractions, cut_tbp_curve};

const FIXTURE: &str = include_str!("../../../validation/columns/petroleum-v1.json");

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn count(value: &Value) -> u16 {
    u16::try_from(value.as_u64().unwrap()).unwrap()
}

fn numbers(value: &Value) -> Vec<f64> {
    value.as_array().unwrap().iter().map(number).collect()
}

fn fixture() -> Value {
    serde_json::from_str(FIXTURE).unwrap()
}

fn case(id: &str) -> Value {
    fixture()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == id)
        .unwrap_or_else(|| panic!("fixture case {id}"))
        .clone()
}

fn pseudo_inputs(case: &Value) -> Vec<PseudoComponentInput> {
    case["components"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|c| c["id"].as_str().unwrap().starts_with("pseudo:"))
        .map(|c| {
            PseudoComponentInput::new(
                c["id"].as_str().unwrap(),
                number(&c["normal_boiling_point_K"]),
                number(&c["specific_gravity"]),
            )
        })
        .collect()
}

fn selection(case: &Value) -> PackageSelection {
    assert_eq!(case["variant"], "PR");
    let database = if case["water"].as_bool().unwrap() {
        vec!["water".to_owned()]
    } else {
        Vec::new()
    };
    PackageSelection::with_pseudo_components(
        PackageChoice::PengRobinson,
        database,
        pseudo_inputs(case),
        case["kij"]
            .as_array()
            .unwrap()
            .iter()
            .map(numbers)
            .collect(),
        case["kij_note"].as_str().unwrap(),
    )
}

/// The fractionator specification of a fixture case.
fn spec_for(case: &Value) -> PetroleumFractionatorSpec {
    let draws = case["draws"].as_array().unwrap();
    PetroleumFractionatorSpec {
        side_draws: u16::try_from(draws.len()).unwrap(),
        feeds: u16::try_from(case["feeds"].as_array().unwrap().len()).unwrap(),
        feed_stages: case["feeds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| count(&f["stage"]))
            .collect(),
        products: draws
            .iter()
            .map(|d| SideProduct {
                stage: count(&d["stage"]),
                rate: match d.get("standard_liquid_volume_flow_m3_s") {
                    Some(v) => DrawRate::StandardLiquidVolume(Some(number(v))),
                    None => DrawRate::Molar(Some(number(&d["molar_flow_mol_s"]))),
                },
            })
            .collect(),
        internals: ColumnInternals::Trays(count(&case["internal_stages"])),
        pressure_pa: Some(number(&case["pressure_Pa"])),
        condenser_pressure_drop_pa: number(&case["condenser_pressure_drop_Pa"]),
        pressure_drop: ColumnPressureDrop::PerStage(number(&case["stage_pressure_drop_Pa"])),
        condenser: match case["condenser"].as_str().unwrap() {
            "total" => DistillationCondenser::Total,
            "partial_vapor" => DistillationCondenser::PartialVapor,
            other => panic!("unexpected condenser {other}"),
        },
        top: TopSpecification::RefluxRatio(Some(number(&case["top"]["reflux_ratio"]))),
        bottom: if case["steam"].is_null() {
            FractionatorBottom::Reboiler(BottomSpecification::DistillateFraction(Some(number(
                &case["bottom"]["distillate_fraction"],
            ))))
        } else {
            FractionatorBottom::StrippingSteam
        },
        ..Default::default()
    }
}

struct Rig {
    project: Project,
    column: UnitId,
    /// Overhead, bottoms and side products, in output order.
    products: Vec<StreamId>,
    condenser: StreamId,
    reboiler: Option<StreamId>,
}

fn add_feed(project: &mut Project, inlet: &Value) -> UnitId {
    let feed = project.flowsheet.add_unit(UnitKind::Feed);
    project.flowsheet.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(number(&inlet["temperature_K"])),
            pressure_pa: Some(number(&inlet["pressure_Pa"])),
            molar_flow_mol_s: Some(number(&inlet["molar_flow_mol_s"])),
        },
    );
    project.feeds.insert(feed, numbers(&inlet["z"]));
    feed
}

impl Rig {
    fn new(selection: PackageSelection, spec: PetroleumFractionatorSpec, case: &Value) -> Self {
        let mut project = Project::new(selection);
        let column = project.flowsheet.add_unit(UnitKind::PetroleumFractionator);
        project.flowsheet.set_spec(
            column,
            UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec.clone())),
        );
        let ports = project.flowsheet.unit(column).unwrap().ports();
        let index = |label: &str| {
            ports
                .iter()
                .position(|p| p.label() == label)
                .unwrap_or_else(|| panic!("port {label}"))
        };
        for (i, inlet) in case["feeds"].as_array().unwrap().iter().enumerate() {
            let feed = add_feed(&mut project, inlet);
            let port = if i == 0 {
                0
            } else {
                index(&format!("feed {}", i + 1))
            };
            project
                .flowsheet
                .connect(PortRef::new(feed, 0), PortRef::new(column, port))
                .unwrap();
        }
        if !case["steam"].is_null() {
            let steam = add_feed(&mut project, &case["steam"]);
            project
                .flowsheet
                .connect(
                    PortRef::new(steam, 0),
                    PortRef::new(column, index("stripping steam")),
                )
                .unwrap();
        }
        let outlets: Vec<usize> = [1, 2]
            .into_iter()
            .chain((1..=spec.side_draws).map(|k| index(&format!("side draw {k}"))))
            .collect();
        let fs = &mut project.flowsheet;
        let products = outlets
            .into_iter()
            .map(|i| {
                let product = fs.add_unit(UnitKind::Product);
                fs.connect(PortRef::new(column, i), PortRef::new(product, 0))
                    .unwrap()
            })
            .collect();
        let sink = fs.add_unit(UnitKind::EnergySink);
        let condenser = fs
            .connect(PortRef::new(column, 3), PortRef::new(sink, 0))
            .unwrap();
        let reboiler = (!spec.steam()).then(|| {
            let source = fs.add_unit(UnitKind::EnergySource);
            fs.set_spec(
                source,
                UnitSpec::EnergySource {
                    duty: EnergySourceDuty::Calculated,
                },
            );
            fs.connect(
                PortRef::new(source, 0),
                PortRef::new(column, index("reboiler duty")),
            )
            .unwrap()
        });
        Self {
            project,
            column,
            products,
            condenser,
            reboiler,
        }
    }

    fn from_case(case: &Value) -> Self {
        Self::new(selection(case), spec_for(case), case)
    }

    fn set(&mut self, spec: PetroleumFractionatorSpec) {
        self.project.flowsheet.set_spec(
            self.column,
            UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec)),
        );
    }

    fn solve(&self) -> Result<ResultSet, SolveFailure> {
        solve_project(
            &self.project,
            &self.project.package.prepare().unwrap(),
            |_| ControlFlow::Continue(()),
            0,
        )
    }

    fn error(&self) -> String {
        let error = self.solve().unwrap_err();
        assert_eq!(error.unit, Some(self.column), "{}", error.message);
        error.message
    }
}

/// Standard liquid molar volume of each fixture component, m³/mol, from the
/// engine's own characterisation and water density.
fn standard_volumes(case: &Value) -> Vec<f64> {
    pseudo_inputs(case)
        .iter()
        .map(|p| {
            let c = p.characterise().unwrap();
            c.molar_mass_kg_per_mol / (p.specific_gravity * WATER_DENSITY_60F_KG_M3)
        })
        .collect()
}

#[test]
fn assay_cuts_and_water_density_match_the_reference_definitions() {
    let fixture = fixture();
    assert_eq!(
        number(&fixture["water_density_60F_kg_m3"]),
        WATER_DENSITY_60F_KG_M3
    );
    for case in fixture["cases"].as_array().unwrap() {
        let tbp: Vec<(f64, f64)> = case["tbp_curve_K"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| (number(&p[0]), number(&p[1])))
            .collect();
        let cuts = cut_tbp_curve(
            &tbp,
            usize::from(count(&case["cuts"])),
            number(&case["watson_k"]),
            "pseudo:crude-",
        )
        .unwrap();
        let fractions = cut_mole_fractions(&cuts).unwrap();
        let expected = numbers(&case["cut_mole_fractions"]);
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * b.abs().max(1e-300);
        for ((cut, x), (reference, x_ref)) in cuts
            .iter()
            .zip(&fractions)
            .zip(pseudo_inputs(case).iter().zip(&expected))
        {
            assert_eq!(cut.input.id, reference.id);
            assert!(close(
                cut.input.normal_boiling_point_k,
                reference.normal_boiling_point_k
            ));
            assert!(close(
                cut.input.specific_gravity,
                reference.specific_gravity
            ));
            assert!(close(*x, *x_ref), "{} {x} {x_ref}", cut.input.id);
        }
        let volumes = standard_volumes(case);
        let pseudo = case["components"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|c| c["id"].as_str().unwrap().starts_with("pseudo:"));
        for (v, c) in volumes.iter().zip(pseudo) {
            assert!(close(*v, number(&c["standard_liquid_volume_m3_per_mol"])));
        }
    }
}

#[test]
fn crude_fractionators_match_independent_references() {
    let fixture = fixture();
    let compare = |id: &Value, field: &str, actual: f64, expected: f64| {
        let range = &fixture["tolerances"][field];
        let tolerance = number(&range[0]) + number(&range[1]) * expected.abs();
        assert!(
            (actual - expected).abs() <= tolerance,
            "{id} {field}: {actual} differs from {expected} by more than {tolerance}"
        );
    };
    for case in fixture["cases"].as_array().unwrap() {
        let id = &case["id"];
        let mut rig = Rig::from_case(case);
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{id}: {}", e.message));
        let unit = &result.units[&rig.column];
        assert!(unit.warnings.is_empty(), "{id}: {:?}", unit.warnings);
        let column = unit.distillation_column.as_ref().unwrap();
        let temperatures = numbers(&case["stage_temperatures_K"]);
        assert_eq!(column.stage_temperatures_k.len(), temperatures.len());
        for (actual, expected) in column.stage_temperatures_k.iter().zip(&temperatures) {
            compare(id, "temperature_K", *actual, *expected);
        }
        for (actual, expected) in column
            .stage_pressures_pa
            .iter()
            .zip(numbers(&case["stage_pressures_Pa"]))
        {
            assert!(
                (actual - expected).abs() <= 1e-9 * expected,
                "{id} pressure"
            );
        }
        compare(
            id,
            "duty_W",
            column.condenser_duty_w,
            number(&case["condenser_duty_W"]),
        );
        assert_eq!(
            result.energy_streams[&rig.condenser].duty_w,
            column.condenser_duty_w
        );
        match rig.reboiler {
            Some(port) => {
                compare(
                    id,
                    "duty_W",
                    column.reboiler_duty_w,
                    number(&case["reboiler_duty_W"]),
                );
                assert_eq!(result.energy_streams[&port].duty_w, column.reboiler_duty_w);
            }
            None => assert_eq!(column.reboiler_duty_w, 0.0),
        }
        for (actual, field) in [
            (column.distillate_mol_s, "distillate_mol_s"),
            (column.reflux_mol_s, "reflux_mol_s"),
            (column.bottoms_mol_s, "bottoms_mol_s"),
            (column.boilup_mol_s, "boilup_mol_s"),
        ] {
            compare(
                id,
                "molar_flow_mol_s",
                actual.unwrap(),
                number(&case[field]),
            );
        }
        // Component and energy balances close, including the steam feed.
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9, "{id}");
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9, "{id}");
        let volumes = standard_volumes(case);
        let outlets = case["outlets"].as_array().unwrap();
        assert_eq!(rig.products.len(), outlets.len(), "{id}");
        for (k, (stream, expected)) in rig.products.iter().zip(outlets).enumerate() {
            let actual = &result.streams[stream];
            compare(
                id,
                "molar_flow_mol_s",
                actual.molar_flow_mol_s,
                number(&expected["molar_flow_mol_s"]),
            );
            let state = actual.state.as_ref().unwrap();
            let vapor = if expected["phase"] == "vapor" {
                1.0
            } else {
                0.0
            };
            assert_eq!(state.vapor_fraction, Some(vapor), "{id}");
            for (actual, expected) in state.mole_fractions.iter().zip(numbers(&expected["z"])) {
                compare(id, "composition", *actual, expected);
            }
            // Side products specified by standard liquid volume meet it.
            if let Some(DrawRate::StandardLiquidVolume(Some(target))) = k
                .checked_sub(2)
                .and_then(|k| spec_for(case).products.get(k).map(|p| p.rate))
            {
                let volume: f64 = state
                    .mole_fractions
                    .iter()
                    .zip(&volumes)
                    .map(|(x, v)| x * v)
                    .sum::<f64>()
                    * actual.molar_flow_mol_s;
                assert!((volume - target).abs() <= 1e-10 * target, "{id}");
                compare(
                    id,
                    "standard_liquid_volume_flow_m3_s",
                    volume,
                    number(&expected["standard_liquid_volume_flow_m3_s"]),
                );
            }
        }
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        let reopened = Project::from_json(&saved).unwrap();
        assert_eq!(reopened.to_json().unwrap(), saved, "{id}");
        assert_eq!(
            reopened.flowsheet.unit(rig.column).unwrap().spec,
            rig.project.flowsheet.unit(rig.column).unwrap().spec
        );
    }
}

#[test]
fn infeasible_or_unsupported_fractionators_fail_explicitly() {
    let case = case("crude-two-draws-reboiled-PR");
    let base = spec_for(&case);
    let mut rig = Rig::from_case(&case);

    // Side products that take more than the feed less the distillate.
    let mut spec = base.clone();
    spec.products[1].rate = DrawRate::Molar(Some(60.0));
    rig.set(spec);
    assert!(
        rig.error().contains("must leave positive bottoms flow"),
        "{}",
        rig.error()
    );

    // A volume flow the column cannot deliver above the bottoms.
    let mut spec = base.clone();
    spec.products[0].rate = DrawRate::StandardLiquidVolume(Some(0.05));
    rig.set(spec);
    let message = rig.error();
    assert!(message.contains("petroleum fractionator"), "{message}");

    // A side product on the condenser or reboiler stage.
    let mut spec = base.clone();
    spec.products[0].stage = 12;
    rig.set(spec);
    let message = rig.solve().unwrap_err().message;
    assert!(message.contains("side-draw stage"), "{message}");

    // A missing side-product flow is a specification issue.
    let mut spec = base.clone();
    spec.products[0].rate = DrawRate::StandardLiquidVolume(None);
    rig.set(spec);
    assert!(
        rig.project
            .flowsheet
            .unit(rig.column)
            .unwrap()
            .issues()
            .iter()
            .any(|issue| issue.field == "side-product standard liquid volume flow")
    );

    // Volume-specified products need a standard liquid volume for every
    // component: water has no specific gravity.
    let steam_case = case_with_volume_draw_and_steam();
    let rig = Rig::new(selection(&steam_case), spec_for(&steam_case), &steam_case);
    let message = rig.error();
    assert!(
        message.contains("water has no specific gravity"),
        "{message}"
    );

    // Missing characterisation: a selected pseudo-component without input.
    let mut missing = selection(&case);
    missing.pseudo_components.pop();
    assert!(missing.prepare().is_err());

    // A package without pseudo-components (PR water/methanol), and a package
    // that cannot supply the component characterisation (GERG-2008).
    let mut methanol = case.clone();
    methanol["feeds"][0]["z"] = json!([0.5, 0.5]);
    methanol["feeds"][0]["temperature_K"] = json!(300.0);
    let package = PackageSelection::default_for(
        PackageChoice::PengRobinson,
        vec!["water".into(), "methanol".into()],
        "zero kij",
    );
    let none = Rig::new(package, base.clone(), &methanol);
    let message = none.error();
    assert!(
        message.contains("the selected package has none"),
        "{message}"
    );
    let gerg = PackageSelection::default_for(
        PackageChoice::Gerg2008,
        vec!["ethane".into(), "propane".into()],
        "no interaction matrix",
    );
    let mut gas = methanol;
    gas["feeds"][0]["pressure_Pa"] = json!(2.0e6);
    gas["feeds"][0]["temperature_K"] = json!(295.0);
    let mut spec = base;
    spec.pressure_pa = Some(1.5e6);
    let message = Rig::new(gerg, spec, &gas).error();
    assert!(
        message.contains("needs a cubic package with characterised pseudo-components"),
        "{message}"
    );

    // A pseudo-component package whose feed is only water.
    let mut steam_only = self::case("crude-steam-stripped-PR");
    let mut z = vec![0.0; 9];
    z[0] = 1.0;
    steam_only["feeds"][0]["z"] = json!(z);
    steam_only["feeds"][0]["temperature_K"] = json!(300.0);
    let rig = Rig::new(selection(&steam_only), spec_for(&steam_only), &steam_only);
    let message = rig.error();
    assert!(
        message.contains("feed carries no characterised pseudo-component"),
        "{message}"
    );
}

/// The steam case with its first side product specified by volume.
fn case_with_volume_draw_and_steam() -> Value {
    let mut steam = case("crude-steam-stripped-PR");
    steam["draws"][0] = json!({"stage": 4, "standard_liquid_volume_flow_m3_s": 2.5e-3});
    steam
}

#[test]
fn fractionator_ports_follow_the_specification() {
    let case = case("crude-steam-stripped-PR");
    let rig = Rig::from_case(&case);
    let labels: Vec<String> = rig
        .project
        .flowsheet
        .unit(rig.column)
        .unwrap()
        .ports()
        .iter()
        .map(|p| p.label())
        .collect();
    assert_eq!(
        labels,
        [
            "feed",
            "overhead",
            "bottoms",
            "condenser duty",
            "side draw 1",
            "side draw 2",
            "stripping steam"
        ]
    );
    let mut spec = spec_for(&case);
    spec.feeds = 2;
    spec.condenser = DistillationCondenser::PartialVaporLiquid {
        liquid_distillate_mol_s: Some(1.0),
    };
    spec.bottom = FractionatorBottom::default();
    let ports = flowsheet_core::ports_for(
        UnitKind::PetroleumFractionator,
        &UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec)),
    );
    let labels: Vec<String> = ports.iter().map(|p| p.label()).collect();
    assert_eq!(
        &labels[4..],
        [
            "side draw 1",
            "side draw 2",
            "feed 2",
            "liquid distillate",
            "reboiler duty"
        ]
    );
}

#[test]
fn records_round_trip_and_drawn_only_records_load() {
    let case = case("crude-volume-draws-partial-pressure-drop-PR");
    let mut rig = Rig::from_case(&case);
    let mut spec = spec_for(&case);
    spec.top = TopSpecification::RefluxRate(Some(70.0));
    spec.bottom = FractionatorBottom::Reboiler(BottomSpecification::ReboilerDuty(Some(1.0e6)));
    spec.internals = ColumnInternals::Trays(12);
    rig.set(spec.clone());
    let saved: Value = serde_json::from_str(&rig.project.to_json().unwrap()).unwrap();
    let record = saved["units"]
        .as_array()
        .unwrap()
        .iter()
        .find(|u| u["kind"] == "petroleum_fractionator")
        .unwrap()["spec"]["unit"]
        .clone();
    assert_eq!(
        record["products"][0],
        json!({"stage": 4, "rate": {"standard_liquid_volume_flow_m3_s": 2.5e-3}})
    );
    assert_eq!(
        record["top_specification"],
        json!({"reflux_rate_mol_s": 70.0})
    );
    let reopened = Project::from_json(&rig.project.to_json().unwrap()).unwrap();
    assert_eq!(
        reopened.flowsheet.unit(rig.column).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec.clone()))
    );

    // Stripping steam round trips and refuses a reboiler specification.
    let mut steam = spec.clone();
    steam.bottom = FractionatorBottom::StrippingSteam;
    // A connected reboiler-duty port cannot be removed; use a fresh column.
    assert!(!rig.project.flowsheet.set_spec(
        rig.column,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(steam.clone()))
    ));
    let mut fresh = Project::new(selection(&case));
    let id = fresh.flowsheet.add_unit(UnitKind::PetroleumFractionator);
    assert!(fresh.flowsheet.set_spec(
        id,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(steam.clone()))
    ));
    let text = fresh.to_json().unwrap();
    assert!(text.contains("\"stripping_steam\": true"));
    let reopened = Project::from_json(&text).unwrap();
    assert_eq!(
        reopened.flowsheet.unit(id).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(steam))
    );
    let conflicting = text.replace(
        "\"stripping_steam\": true",
        "\"stripping_steam\": true, \"distillate_fraction\": 0.3",
    );
    assert!(Project::from_json(&conflicting).is_err());

    // A drawn-only fractionator record (side draws only) with connected side
    // draws loads with the defaults and keeps its connections.
    let mut legacy = Project::new(selection(&case));
    let fs = &mut legacy.flowsheet;
    let column = fs.add_unit(UnitKind::PetroleumFractionator);
    fs.set_spec(
        column,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(
            PetroleumFractionatorSpec {
                side_draws: 2,
                ..Default::default()
            },
        )),
    );
    let product = fs.add_unit(UnitKind::Product);
    let draw = fs
        .connect(PortRef::new(column, 5), PortRef::new(product, 0))
        .unwrap();
    let mut value: Value = serde_json::from_str(&legacy.to_json().unwrap()).unwrap();
    let units = value["units"].as_array_mut().unwrap();
    let unit = units
        .iter_mut()
        .find(|u| u["kind"] == "petroleum_fractionator")
        .unwrap();
    assert_eq!(
        unit["spec"]["unit"],
        json!({"kind": "petroleum_fractionator", "side_draws": 2})
    );
    unit["spec"]["unit"] = json!({"kind": "petroleum_fractionator", "side_draws": 2});
    let text = serde_json::to_string(&value).unwrap();
    let loaded = Project::from_json(&text).unwrap();
    assert_eq!(
        loaded.flowsheet.unit(column).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(
            PetroleumFractionatorSpec {
                side_draws: 2,
                ..Default::default()
            }
        ))
    );
    let stream = loaded.flowsheet.stream(draw).unwrap();
    assert_eq!(stream.from, PortRef::new(column, 5));
    assert_eq!(
        loaded.flowsheet.unit(column).unwrap().ports()[5].label(),
        "side draw 2"
    );
}
