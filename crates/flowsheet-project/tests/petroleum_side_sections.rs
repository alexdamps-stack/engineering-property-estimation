//! Petroleum fractionator pump-arounds and side strippers (#190): crude-like
//! pseudo-component columns against the independent references of
//! `validation/columns/petroleum-v2.json`, balance closure, explicit
//! failures, ports and project persistence.

use std::ops::ControlFlow;

use flowsheet_core::units::columns::{
    BottomSpecification, ColumnInternals, ColumnPressureDrop, DistillationCondenser, DrawRate,
    FractionatorBottom, PetroleumFractionatorSpec, PumpAround, PumpAroundCooling, SideProduct,
    SideStripper, SideStripperBottom, TopSpecification, WATER_DENSITY_60F_KG_M3,
};
use flowsheet_core::units::{ColumnSpec, GroupSpec};
use flowsheet_core::{
    EnergySourceDuty, PortDirection, PortMedium, PortRef, StreamId, UnitId, UnitKind, UnitSpec,
};
use flowsheet_project::{
    PackageChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};
use thermo_data::PseudoComponentInput;

const FIXTURE: &str = include_str!("../../../validation/columns/petroleum-v2.json");

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
        pump_arounds: case["pump_arounds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| PumpAround {
                draw_stage: count(&p["draw_stage"]),
                return_stage: count(&p["return_stage"]),
                molar_flow_mol_s: Some(number(&p["molar_flow_mol_s"])),
                cooling: match p.get("return_temperature_K") {
                    Some(t) => PumpAroundCooling::ReturnTemperature(Some(number(t))),
                    None => PumpAroundCooling::Duty(Some(number(&p["duty_W"]))),
                },
            })
            .collect(),
        side_strippers: draws
            .iter()
            .enumerate()
            .filter_map(|(k, d)| {
                let s = d.get("stripper")?;
                Some(SideStripper {
                    product: u16::try_from(k + 1).unwrap(),
                    stages: count(&s["stages"]),
                    return_stage: count(&s["return_stage"]),
                    bottom: if let Some(q) = s.get("reboiler_duty_W") {
                        SideStripperBottom::ReboilerDuty(Some(number(q)))
                    } else if let Some(r) = s.get("boilup_ratio") {
                        SideStripperBottom::BoilupRatio(Some(number(r)))
                    } else {
                        SideStripperBottom::Steam
                    },
                })
            })
            .collect(),
    }
}

/// Steam inlets of the case's steam-stripped side strippers, in order.
fn stripper_steam(case: &Value) -> Vec<Value> {
    case["draws"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|d| d.get("stripper")?.get("steam").cloned())
        .collect()
}

struct Rig {
    project: Project,
    column: UnitId,
    /// Stream at each port of the column, by port index.
    streams: Vec<StreamId>,
    labels: Vec<String>,
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
    /// Connects every port: feeds from the case's `feeds`, `steam` and
    /// stripper steam, products, energy sinks and calculated energy sources.
    fn new(selection: PackageSelection, spec: PetroleumFractionatorSpec, case: &Value) -> Self {
        let mut project = Project::new(selection);
        let column = project.flowsheet.add_unit(UnitKind::PetroleumFractionator);
        assert!(project.flowsheet.set_spec(
            column,
            UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec)),
        ));
        let ports = project.flowsheet.unit(column).unwrap().ports();
        let labels: Vec<String> = ports.iter().map(|p| p.label()).collect();
        let steam = stripper_steam(case);
        let mut streams = Vec::new();
        for (i, port) in ports.iter().enumerate() {
            let label = port.label();
            let at = PortRef::new(column, i);
            let fs = &mut project.flowsheet;
            let stream = match (port.medium, port.direction) {
                (PortMedium::Material, PortDirection::Inlet) => {
                    let inlet = if label == "feed" {
                        case["feeds"][0].clone()
                    } else if let Some(k) = label.strip_prefix("feed ") {
                        case["feeds"][k.parse::<usize>().unwrap() - 1].clone()
                    } else if label == "stripping steam" {
                        case["steam"].clone()
                    } else if let Some(k) = label.strip_prefix("stripper steam ") {
                        steam[k.parse::<usize>().unwrap() - 1].clone()
                    } else {
                        panic!("unexpected inlet {label}")
                    };
                    let feed = add_feed(&mut project, &inlet);
                    project
                        .flowsheet
                        .connect(PortRef::new(feed, 0), at)
                        .unwrap()
                }
                (PortMedium::Material, PortDirection::Outlet) => {
                    let product = fs.add_unit(UnitKind::Product);
                    fs.connect(at, PortRef::new(product, 0)).unwrap()
                }
                (PortMedium::Energy, PortDirection::Outlet) => {
                    let sink = fs.add_unit(UnitKind::EnergySink);
                    fs.connect(at, PortRef::new(sink, 0)).unwrap()
                }
                (PortMedium::Energy, PortDirection::Inlet) => {
                    let source = fs.add_unit(UnitKind::EnergySource);
                    fs.set_spec(
                        source,
                        UnitSpec::EnergySource {
                            duty: EnergySourceDuty::Calculated,
                        },
                    );
                    fs.connect(PortRef::new(source, 0), at).unwrap()
                }
            };
            streams.push(stream);
        }
        Self {
            project,
            column,
            streams,
            labels,
        }
    }

    fn from_case(case: &Value) -> Self {
        Self::new(selection(case), spec_for(case), case)
    }

    fn find(&self, label: &str) -> Option<StreamId> {
        let i = self.labels.iter().position(|l| l == label)?;
        Some(self.streams[i])
    }

    fn port(&self, label: &str) -> StreamId {
        self.find(label).unwrap_or_else(|| panic!("port {label}"))
    }

    fn set(&mut self, spec: PetroleumFractionatorSpec) {
        assert!(self.project.flowsheet.set_spec(
            self.column,
            UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec)),
        ));
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

    fn issues(&self) -> Vec<&'static str> {
        self.project
            .flowsheet
            .unit(self.column)
            .unwrap()
            .issues()
            .iter()
            .map(|i| i.field)
            .collect()
    }
}

/// Specification issue fields of a fractionator.
fn issues(spec: &PetroleumFractionatorSpec) -> Vec<&'static str> {
    ColumnSpec::PetroleumFractionator(spec.clone())
        .issues()
        .iter()
        .map(|i| i.field)
        .collect()
}

/// Standard liquid molar volume of each fixture component, m³/mol (water,
/// when present, has none).
fn standard_volumes(case: &Value) -> Vec<Option<f64>> {
    let water = case["water"].as_bool().unwrap();
    (water.then_some(None).into_iter())
        .chain(pseudo_inputs(case).iter().map(|p| {
            let c = p.characterise().unwrap();
            Some(c.molar_mass_kg_per_mol / (p.specific_gravity * WATER_DENSITY_60F_KG_M3))
        }))
        .collect()
}

#[test]
fn pump_arounds_and_side_strippers_match_independent_references() {
    let fixture = fixture();
    let compare = |id: &Value, field: &str, actual: f64, expected: f64| {
        let range = &fixture["tolerances"][field];
        let tolerance = number(&range[0]) + number(&range[1]) * expected.abs();
        assert!(
            (actual - expected).abs() <= tolerance,
            "{id} {field}: {actual} differs from {expected} by more than {tolerance}"
        );
    };
    assert_eq!(
        number(&fixture["water_density_60F_kg_m3"]),
        WATER_DENSITY_60F_KG_M3
    );
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
            assert!((actual - expected).abs() <= 1e-9 * expected, "{id}");
        }
        compare(
            id,
            "duty_W",
            column.condenser_duty_w,
            number(&case["condenser_duty_W"]),
        );
        assert_eq!(
            result.energy_streams[&rig.port("condenser duty")].duty_w,
            column.condenser_duty_w
        );
        if case["steam"].is_null() {
            compare(
                id,
                "duty_W",
                column.reboiler_duty_w,
                number(&case["reboiler_duty_W"]),
            );
            assert_eq!(
                result.energy_streams[&rig.port("reboiler duty")].duty_w,
                column.reboiler_duty_w
            );
        } else {
            assert_eq!(column.reboiler_duty_w, 0.0);
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
        // Pump-around duties: result, energy stream and reference.
        let duties = numbers(&case["pump_around_duties_W"]);
        assert_eq!(column.pump_around_duties_w.len(), duties.len(), "{id}");
        for (k, (actual, expected)) in column.pump_around_duties_w.iter().zip(&duties).enumerate() {
            compare(id, "duty_W", *actual, *expected);
            let port = rig.port(&format!("pump-around duty {}", k + 1));
            assert_eq!(result.energy_streams[&port].duty_w, *actual, "{id}");
            if let Some(q) = case["pump_arounds"][k].get("duty_W") {
                assert!((actual - number(q)).abs() <= 1e-9 * number(q), "{id}");
            }
        }
        // Side strippers: stage temperatures, draws, vapor returns, duties.
        let strippers = case["side_strippers"].as_array().unwrap();
        assert_eq!(column.side_strippers.len(), strippers.len(), "{id}");
        for (k, (actual, expected)) in column.side_strippers.iter().zip(strippers).enumerate() {
            let temperatures = numbers(&expected["stage_temperatures_K"]);
            assert_eq!(actual.stage_temperatures_k.len(), temperatures.len());
            for (t, t_ref) in actual.stage_temperatures_k.iter().zip(&temperatures) {
                compare(id, "temperature_K", *t, *t_ref);
            }
            compare(
                id,
                "molar_flow_mol_s",
                actual.draw_mol_s,
                number(&expected["draw_mol_s"]),
            );
            compare(
                id,
                "molar_flow_mol_s",
                actual.vapor_return_mol_s,
                number(&expected["vapor_return_mol_s"]),
            );
            compare(
                id,
                "duty_W",
                actual.reboiler_duty_w,
                number(&expected["reboiler_duty_W"]),
            );
            if let Some(port) = rig.find(&format!("stripper reboiler duty {}", k + 1)) {
                assert_eq!(
                    result.energy_streams[&port].duty_w, actual.reboiler_duty_w,
                    "{id}"
                );
            } else {
                assert_eq!(actual.reboiler_duty_w, 0.0, "{id}");
            }
        }
        // Component and energy balances close around the whole system,
        // including steam to the column and to the strippers.
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9, "{id}");
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9, "{id}");
        // Products in port order: overhead, bottoms, side draw 1..k.
        let volumes = standard_volumes(case);
        let spec = spec_for(case);
        let outlets = case["outlets"].as_array().unwrap();
        let products: Vec<StreamId> = ["overhead".to_owned(), "bottoms".to_owned()]
            .into_iter()
            .chain((1..=spec.side_draws).map(|k| format!("side draw {k}")))
            .map(|label| rig.port(&label))
            .collect();
        assert_eq!(products.len(), outlets.len(), "{id}");
        for (k, (stream, expected)) in products.iter().zip(outlets).enumerate() {
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
            if let Some(DrawRate::StandardLiquidVolume(Some(target))) =
                k.checked_sub(2).map(|k| spec.products[k].rate)
            {
                let volume: f64 = state
                    .mole_fractions
                    .iter()
                    .zip(&volumes)
                    .map(|(x, v)| x * v.unwrap())
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
        // A stripped product is the stripper's bottoms at its bottom-stage
        // temperature.
        for (stripper, expected) in spec.side_strippers.iter().zip(strippers) {
            let stream = &result.streams[&products[1 + usize::from(stripper.product)]];
            let t = stream.state.as_ref().unwrap().temperature_k;
            let bottom = *numbers(&expected["stage_temperatures_K"]).last().unwrap();
            compare(id, "temperature_K", t, bottom);
            compare(
                id,
                "molar_flow_mol_s",
                stream.molar_flow_mol_s,
                number(&expected["bottoms_mol_s"]),
            );
        }
        // The project with results round trips.
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        let reopened = Project::from_json(&saved).unwrap();
        assert_eq!(reopened.to_json().unwrap(), saved, "{id}");
        assert_eq!(
            reopened.flowsheet.unit(rig.column).unwrap().spec,
            rig.project.flowsheet.unit(rig.column).unwrap().spec
        );
        assert_eq!(reopened.results, rig.project.results, "{id}");
    }
}

#[test]
fn infeasible_pump_arounds_and_side_strippers_fail_explicitly() {
    let case = case("crude-pump-arounds-reboiled-PR");
    let base = spec_for(&case);
    let mut rig = Rig::from_case(&case);

    // A return temperature above the draw-stage temperature heats the
    // column instead of cooling it.
    let mut spec = base.clone();
    spec.pump_arounds[0].cooling = PumpAroundCooling::ReturnTemperature(Some(455.0));
    rig.set(spec);
    let message = rig.error();
    assert!(
        message.contains("pump-around 1 return temperature 455 K is not below"),
        "{message}"
    );

    // A pump-around returning below its draw, missing values, and a
    // negative duty are specification issues.
    let mut spec = base.clone();
    spec.pump_arounds[0].return_stage = 8;
    spec.pump_arounds[1].molar_flow_mol_s = None;
    spec.pump_arounds.push(PumpAround {
        draw_stage: 5,
        return_stage: 3,
        molar_flow_mol_s: Some(10.0),
        cooling: PumpAroundCooling::Duty(Some(-1.0)),
    });
    rig.set(spec);
    let found = rig.issues();
    for field in [
        "pump-around stage",
        "pump-around molar flow",
        "pump-around duty",
    ] {
        assert!(found.contains(&field), "{field}: {found:?}");
    }
    // Specification issues stop the solve before the calculation.
    let message = rig.solve().unwrap_err().message;
    assert!(message.contains("pump-around stage"), "{message}");

    // A side stripper on a product that does not exist, two strippers on
    // one product, and a vapor return below the draw.
    let mut spec = base.clone();
    spec.pump_arounds.clear();
    let stripper = SideStripper {
        product: 3,
        stages: 2,
        return_stage: 3,
        bottom: SideStripperBottom::ReboilerDuty(Some(1.0e5)),
    };
    spec.side_strippers = vec![stripper];
    assert_eq!(issues(&spec), ["side-stripper product"]);
    let message = Rig::new(selection(&case), spec.clone(), &case)
        .solve()
        .unwrap_err()
        .message;
    assert!(message.contains("side-stripper product"), "{message}");
    spec.side_strippers = vec![
        SideStripper {
            product: 1,
            ..stripper
        },
        SideStripper {
            product: 1,
            return_stage: 5,
            ..stripper
        },
    ];
    assert_eq!(
        issues(&spec),
        ["side-stripper product", "side-stripper return stage"]
    );
    spec.side_strippers = vec![SideStripper {
        product: 1,
        stages: 0,
        bottom: SideStripperBottom::BoilupRatio(None),
        ..stripper
    }];
    assert_eq!(
        issues(&spec),
        ["side-stripper stages", "side-stripper boilup ratio"]
    );

    // A boilup too small to carry vapor back up the stripper drives the
    // main column dry below the draw: the coupled Newton solve fails.
    spec.side_strippers = vec![SideStripper {
        product: 1,
        stages: 2,
        bottom: SideStripperBottom::BoilupRatio(Some(0.3)),
        ..stripper
    }];
    let rig = Rig::new(selection(&case), spec, &case);
    let message = rig.error();
    assert!(message.contains("column MESH"), "{message}");

    // A steam-stripped side stripper without steam flow.
    let steam_case = self::case("crude-steam-side-stripper-pump-around-PR");
    let mut no_steam = steam_case.clone();
    no_steam["draws"][0]["stripper"]["steam"]["molar_flow_mol_s"] = json!(0.0);
    let rig = Rig::from_case(&no_steam);
    let message = rig.error();
    assert!(
        message.contains("side stripper 1 needs a steam flow"),
        "{message}"
    );
    // Steam below the stripper pressure.
    let mut low = steam_case;
    low["draws"][0]["stripper"]["steam"]["pressure_Pa"] = json!(1.0e5);
    let message = Rig::from_case(&low).error();
    assert!(message.contains("steam pressure is below"), "{message}");
}

#[test]
fn side_section_ports_follow_the_specification() {
    let case = case("crude-steam-side-stripper-pump-around-PR");
    let rig = Rig::from_case(&case);
    assert_eq!(
        rig.labels,
        [
            "feed",
            "overhead",
            "bottoms",
            "condenser duty",
            "side draw 1",
            "side draw 2",
            "stripping steam",
            "pump-around duty 1",
            "stripper steam 1"
        ]
    );
    let mut spec = spec_for(&self::case("crude-reboiled-side-strippers-volume-PR"));
    spec.pump_arounds.push(PumpAround::default());
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
            "reboiler duty",
            "pump-around duty 1",
            "stripper reboiler duty 1",
            "stripper reboiler duty 2"
        ]
    );
}

#[test]
fn side_section_records_round_trip_and_older_records_load() {
    let case = case("crude-steam-side-stripper-pump-around-PR");
    let rig = Rig::from_case(&case);
    let spec = spec_for(&case);
    let saved: Value = serde_json::from_str(&rig.project.to_json().unwrap()).unwrap();
    let record = saved["units"]
        .as_array()
        .unwrap()
        .iter()
        .find(|u| u["kind"] == "petroleum_fractionator")
        .unwrap()["spec"]["unit"]
        .clone();
    assert_eq!(
        record["pump_arounds"],
        json!([{"draw_stage": 7, "return_stage": 6, "molar_flow_mol_s": 25.0,
                "cooling": {"duty_w": 2.0e5}}])
    );
    assert_eq!(
        record["side_strippers"],
        json!([{"product": 1, "stages": 3, "return_stage": 3, "bottom": "steam"}])
    );
    let reopened = Project::from_json(&rig.project.to_json().unwrap()).unwrap();
    assert_eq!(
        reopened.flowsheet.unit(rig.column).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec.clone()))
    );
    // Other cooling and stripping choices round trip.
    let mut other = spec_for(&self::case("crude-reboiled-side-strippers-volume-PR"));
    other.pump_arounds = vec![PumpAround {
        cooling: PumpAroundCooling::ReturnTemperature(Some(400.0)),
        ..PumpAround::default()
    }];
    let mut fresh = Project::new(selection(&case));
    let id = fresh.flowsheet.add_unit(UnitKind::PetroleumFractionator);
    assert!(fresh.flowsheet.set_spec(
        id,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(other.clone()))
    ));
    let text = fresh.to_json().unwrap();
    assert!(text.contains("\"return_temperature_k\": 400.0"), "{text}");
    assert!(text.contains("\"boilup_ratio\": 0.5"), "{text}");
    assert!(text.contains("\"reboiler_duty_w\": 150000.0"), "{text}");
    let reopened = Project::from_json(&text).unwrap();
    assert_eq!(
        reopened.flowsheet.unit(id).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(other))
    );
    // More side strippers than side products is refused on load.
    let mut value: Value = serde_json::from_str(&text).unwrap();
    let unit = value["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "petroleum_fractionator")
        .unwrap();
    let strippers = unit["spec"]["unit"]["side_strippers"].clone();
    let mut three = strippers.as_array().unwrap().clone();
    three.push(three[0].clone());
    unit["spec"]["unit"]["side_strippers"] = json!(three);
    assert!(Project::from_json(&serde_json::to_string(&value).unwrap()).is_err());

    // A fractionator record written before pump-arounds and side strippers
    // (no such keys) loads without them and is written back unchanged.
    let mut plain = spec;
    plain.pump_arounds.clear();
    plain.side_strippers.clear();
    plain.bottom = FractionatorBottom::default();
    let mut older = Project::new(selection(&case));
    let id = older.flowsheet.add_unit(UnitKind::PetroleumFractionator);
    assert!(older.flowsheet.set_spec(
        id,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(plain.clone()))
    ));
    let text = older.to_json().unwrap();
    assert!(!text.contains("pump_arounds") && !text.contains("side_strippers"));
    let loaded = Project::from_json(&text).unwrap();
    assert_eq!(
        loaded.flowsheet.unit(id).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::PetroleumFractionator(plain))
    );
    assert_eq!(loaded.to_json().unwrap(), text);
}
