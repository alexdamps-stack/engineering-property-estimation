//! Petroleum fractionator product curves and cut-point side products (#190)
//! against the independent references of `petroleum-curves-v1.json` (Python
//! standard library on the teqp columns of `petroleum-v1.json`); explicit
//! failures; persistence and results saved before the curves.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{
    BottomSpecification, CUT_POINT_TOLERANCE_K, ColumnInternals, ColumnPressureDrop, CurveBasis,
    CutPoint, DistillationCondenser, DrawRate, FractionatorBottom, PetroleumFractionatorSpec,
    SideProduct, TopSpecification,
};
use flowsheet_core::{EnergySourceDuty, PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    PackageChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};
use thermo_data::PseudoComponentInput;

const COLUMNS: &str = include_str!("../../../validation/columns/petroleum-v1.json");
const CURVES: &str = include_str!("../../../validation/columns/petroleum-curves-v1.json");

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn count(value: &Value) -> u16 {
    u16::try_from(value.as_u64().unwrap()).unwrap()
}

fn numbers(value: &Value) -> Vec<f64> {
    value.as_array().unwrap().iter().map(number).collect()
}

fn columns() -> Value {
    serde_json::from_str(COLUMNS).unwrap()
}

fn curves() -> Value {
    serde_json::from_str(CURVES).unwrap()
}

fn case(id: &str) -> Value {
    columns()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == id)
        .unwrap()
        .clone()
}

fn selection(case: &Value) -> PackageSelection {
    let database = if case["water"].as_bool().unwrap() {
        vec!["water".to_owned()]
    } else {
        Vec::new()
    };
    let pseudo = case["components"]
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
        .collect();
    PackageSelection::with_pseudo_components(
        PackageChoice::PengRobinson,
        database,
        pseudo,
        case["kij"]
            .as_array()
            .unwrap()
            .iter()
            .map(numbers)
            .collect(),
        case["kij_note"].as_str().unwrap(),
    )
}

fn spec_for(case: &Value) -> PetroleumFractionatorSpec {
    let draws = case["draws"].as_array().unwrap();
    PetroleumFractionatorSpec {
        side_draws: u16::try_from(draws.len()).unwrap(),
        feeds: 1,
        feed_stages: vec![count(&case["feeds"][0]["stage"])],
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
        pump_arounds: Vec::new(),
        side_strippers: Vec::new(),
        internals: ColumnInternals::Trays(count(&case["internal_stages"])),
        pressure_pa: Some(number(&case["pressure_Pa"])),
        condenser_pressure_drop_pa: number(&case["condenser_pressure_drop_Pa"]),
        pressure_drop: ColumnPressureDrop::PerStage(number(&case["stage_pressure_drop_Pa"])),
        condenser: match case["condenser"].as_str().unwrap() {
            "total" => DistillationCondenser::Total,
            _ => DistillationCondenser::PartialVapor,
        },
        top: TopSpecification::RefluxRatio(Some(number(&case["top"]["reflux_ratio"]))),
        bottom: if case["steam"].is_null() {
            FractionatorBottom::Reboiler(BottomSpecification::DistillateFraction(Some(number(
                &case["bottom"]["distillate_fraction"],
            ))))
        } else {
            FractionatorBottom::StrippingSteam
        },
    }
}

/// The reference case with its side products specified by the fixture's
/// cut points, starting from `estimates` (mol/s).
fn cut_point_spec(estimates: [Option<f64>; 2]) -> PetroleumFractionatorSpec {
    let fixture = curves();
    let base = &fixture["cut_point_case"];
    let mut spec = spec_for(&case(base["base_case"].as_str().unwrap()));
    for (cut, estimate) in base["cut_points"].as_array().unwrap().iter().zip(estimates) {
        let k = usize::from(count(&cut["side_draw"])) - 1;
        spec.products[k].rate = DrawRate::CutPoint(CutPoint {
            basis: match cut["basis"].as_str().unwrap() {
                "tbp" => CurveBasis::Tbp,
                _ => CurveBasis::AstmD86,
            },
            volume_percent: number(&cut["volume_percent"]),
            temperature_k: Some(number(&cut["temperature_K"])),
            estimate_mol_s: estimate,
        });
    }
    spec
}

struct Rig {
    project: Project,
    column: UnitId,
    products: Vec<StreamId>,
}

impl Rig {
    fn new(case: &Value, spec: PetroleumFractionatorSpec) -> Self {
        let mut project = Project::new(selection(case));
        let fs = &mut project.flowsheet;
        let column = fs.add_unit(UnitKind::PetroleumFractionator);
        fs.set_spec(
            column,
            UnitSpec::Columns(ColumnSpec::PetroleumFractionator(spec.clone())),
        );
        let ports = fs.unit(column).unwrap().ports();
        let index = |label: &str| ports.iter().position(|p| p.label() == label).unwrap();
        let add_feed = |project: &mut Project, inlet: &Value, port: usize| {
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
            project
                .flowsheet
                .connect(PortRef::new(feed, 0), PortRef::new(column, port))
                .unwrap();
        };
        add_feed(&mut project, &case["feeds"][0], 0);
        if !case["steam"].is_null() {
            add_feed(&mut project, &case["steam"], index("stripping steam"));
        }
        let fs = &mut project.flowsheet;
        let outlets: Vec<usize> = [1, 2]
            .into_iter()
            .chain((1..=spec.side_draws).map(|k| index(&format!("side draw {k}"))))
            .collect();
        let products = outlets
            .into_iter()
            .map(|i| {
                let product = fs.add_unit(UnitKind::Product);
                fs.connect(PortRef::new(column, i), PortRef::new(product, 0))
                    .unwrap()
            })
            .collect();
        let sink = fs.add_unit(UnitKind::EnergySink);
        fs.connect(PortRef::new(column, 3), PortRef::new(sink, 0))
            .unwrap();
        if !spec.steam() {
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
            .unwrap();
        }
        Self {
            project,
            column,
            products,
        }
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

#[test]
fn reported_product_curves_match_the_references() {
    let fixture = curves();
    for id in ["crude-two-draws-reboiled-PR", "crude-steam-stripped-PR"] {
        let case = case(id);
        let rig = Rig::new(&case, spec_for(&case));
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{id}: {}", e.message));
        let unit = &result.units[&rig.column];
        assert!(unit.warnings.is_empty(), "{id}: {:?}", unit.warnings);
        let curves = &unit.distillation_column.as_ref().unwrap().product_curves;
        let expected: Vec<&Value> = fixture["products"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|p| p["case"] == id)
            .collect();
        assert_eq!(curves.len(), expected.len(), "{id}");
        for (actual, expected) in curves.iter().zip(expected) {
            assert_eq!(actual.product, expected["product"].as_str().unwrap());
            assert!(
                (actual.excluded_mole_fraction - number(&expected["excluded_mole_fraction"])).abs()
                    <= number(&expected["excluded_budget"]),
                "{id} {}",
                actual.product
            );
            for (points, reference, budget) in [
                (&actual.tbp_k, "tbp_K", "tbp_budget_K"),
                (&actual.d86_k, "d86_K", "d86_budget_K"),
            ] {
                for ((p, r), b) in points
                    .iter()
                    .zip(expected[reference].as_array().unwrap())
                    .zip(numbers(&expected[budget]))
                {
                    assert_eq!(p[0], number(&r[0]));
                    assert!(
                        (p[1] - number(&r[1])).abs() <= b,
                        "{id} {} {reference} {}: {} vs {} (budget {b})",
                        actual.product,
                        p[0],
                        p[1],
                        r[1]
                    );
                }
            }
        }
    }
}

#[test]
fn cut_point_products_reproduce_the_reference_column() {
    let fixture = curves();
    let cut_case = &fixture["cut_point_case"];
    let base = case(cut_case["base_case"].as_str().unwrap());
    let range = |field: &str| number(&columns()["tolerances"][field][0]);
    let mut rig = Rig::new(&base, cut_point_spec([None, None]));
    let result = rig.solve().unwrap_or_else(|e| panic!("{}", e.message));
    let unit = &result.units[&rig.column];
    assert!(unit.warnings.is_empty(), "{:?}", unit.warnings);
    assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
    assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
    let column = unit.distillation_column.as_ref().unwrap();
    for (actual, expected) in column
        .stage_temperatures_k
        .iter()
        .zip(numbers(&base["stage_temperatures_K"]))
    {
        assert!((actual - expected).abs() <= range("temperature_K"));
    }
    for (stream, expected) in rig.products.iter().zip(base["outlets"].as_array().unwrap()) {
        let actual = &result.streams[stream];
        let flow = number(&expected["molar_flow_mol_s"]);
        assert!(
            (actual.molar_flow_mol_s - flow).abs() <= range("molar_flow_mol_s"),
            "{} vs {flow}",
            actual.molar_flow_mol_s
        );
        let state = actual.state.as_ref().unwrap();
        for (x, z) in state.mole_fractions.iter().zip(numbers(&expected["z"])) {
            assert!((x - z).abs() <= range("composition"));
        }
    }
    // The reported curves meet the cut points.
    for cut in cut_case["cut_points"].as_array().unwrap() {
        let label = format!("side draw {}", cut["side_draw"]);
        let curve = column
            .product_curves
            .iter()
            .find(|c| c.product == label)
            .unwrap();
        let points = if cut["basis"] == "tbp" {
            &curve.tbp_k
        } else {
            &curve.d86_k
        };
        let point = points
            .iter()
            .find(|p| p[0] == number(&cut["volume_percent"]))
            .unwrap();
        assert!((point[1] - number(&cut["temperature_K"])).abs() <= CUT_POINT_TOLERANCE_K);
    }

    // Saved with its results, the project reloads unchanged.
    rig.project.results = Some(result);
    let saved = rig.project.to_json().unwrap();
    let reopened = Project::from_json(&saved).unwrap();
    assert_eq!(reopened.to_json().unwrap(), saved);
    let value: Value = serde_json::from_str(&saved).unwrap();
    let record = value["units"]
        .as_array()
        .unwrap()
        .iter()
        .find(|u| u["kind"] == "petroleum_fractionator")
        .unwrap()["spec"]["unit"]
        .clone();
    assert_eq!(
        record["products"][0]["rate"],
        json!({"cut_point": {
            "basis": "tbp",
            "volume_percent": 95.0,
            "temperature_k": number(&cut_case["cut_points"][0]["temperature_K"]),
        }})
    );
    assert_eq!(
        record["products"][1]["rate"]["cut_point"]["basis"],
        "astm_d86"
    );

    // Results saved before the curves existed load without them.
    let mut legacy = value.clone();
    let units = legacy["results"]["units"].as_object_mut().unwrap();
    let mut removed = false;
    for (_, u) in units.iter_mut() {
        if let Some(c) = u
            .get_mut("distillation_column")
            .and_then(Value::as_object_mut)
        {
            removed |= c.remove("product_curves").is_some();
        }
    }
    assert!(removed);
    let loaded = Project::from_json(&serde_json::to_string(&legacy).unwrap()).unwrap();
    let estimate = loaded.results.as_ref().unwrap().units[&rig.column]
        .distillation_column
        .as_ref()
        .unwrap();
    assert!(estimate.product_curves.is_empty());
}

#[test]
fn unattainable_or_invalid_cut_points_fail_explicitly() {
    let base = case("crude-two-draws-reboiled-PR");
    let mut rig = Rig::new(&base, spec_for(&base));
    let molar = spec_for(&base);
    let cut = |basis, volume_percent, temperature_k: Option<f64>| {
        DrawRate::CutPoint(CutPoint {
            basis,
            volume_percent,
            temperature_k,
            estimate_mol_s: Some(15.0),
        })
    };

    // Outside the pseudo-components' TBP range.
    let mut spec = molar.clone();
    spec.products[0].rate = cut(CurveBasis::Tbp, 50.0, Some(600.0));
    rig.set(spec);
    let message = rig.error();
    assert!(
        message.contains("is not attainable") && message.contains("TBP curve spans"),
        "{message}"
    );

    // Inside the range but heavier than a side product above the feed can
    // be: the draw grows until the column leaves no bottoms.
    let mut spec = molar.clone();
    spec.products[1].rate = cut(CurveBasis::Tbp, 50.0, Some(530.0));
    rig.set(spec);
    let message = rig.error();
    assert!(
        message.contains("cut-point specification could not be met"),
        "{message}"
    );

    // Specification issues: D86 above 95 vol% and a missing temperature.
    let mut spec = molar.clone();
    spec.products[0].rate = cut(CurveBasis::AstmD86, 97.0, None);
    rig.set(spec);
    let fields: Vec<_> = rig
        .project
        .flowsheet
        .unit(rig.column)
        .unwrap()
        .issues()
        .iter()
        .map(|i| i.field)
        .collect();
    assert_eq!(
        fields,
        [
            "side-product cut-point volume percent",
            "side-product cut-point temperature"
        ]
    );
}
