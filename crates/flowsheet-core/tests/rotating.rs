//! Pump, compressor, expander and multistage compressor (A14, #182) through
//! the flowsheet solver, against the independent isentropic-efficiency
//! references in `validation/rotating/ps-v1.json` (teqp/Cantera/SciPy brentq
//! outlet solves, no Rust calls; see `validation/rotating/generate-ps-references.py`).

use flowsheet_core::units::RotatingSpec;
use flowsheet_core::units::rotating::MultistageCompressorSpec;
use flowsheet_core::{
    Calculation, EnergySourceDuty, FeedCompositions, Flowsheet, PortRef, PressureRise,
    PressureSpec, SolveError, StreamId, UnitId, UnitKind, UnitSpec,
};
use property_packages::CubicPackage;
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase};

const FIXTURE: &str = include_str!("../../../validation/rotating/ps-v1.json");
/// Energy ports follow the two material ports.
const POWER: usize = 2;
const INTERCOOLER: usize = 3;

fn fixture() -> Value {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(fixture["schema"], "a14-ps-reference-v1");
    fixture
}

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn ids(case: &Value) -> Vec<&str> {
    case["ids"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect()
}

fn composition(ids: &[&str], z: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(z)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

struct Machine {
    flowsheet: Flowsheet,
    package: CubicPackage,
    feeds: FeedCompositions,
    unit: UnitId,
    outlet: StreamId,
    power: StreamId,
    intercooler: Option<StreamId>,
}

impl Machine {
    /// Feed → machine → product, with the power port (and the intercooler
    /// port of a multistage compressor) on calculated energy boundaries.
    fn new(
        ids: &[&str],
        z: &[f64],
        variant: CubicVariant,
        inlet: (f64, f64, f64),
        kind: UnitKind,
        spec: UnitSpec,
    ) -> Self {
        let mut flowsheet = Flowsheet::new();
        let source = flowsheet.add_unit(UnitKind::Feed);
        let unit = flowsheet.add_unit(kind);
        let product = flowsheet.add_unit(UnitKind::Product);
        assert!(flowsheet.set_spec(
            source,
            UnitSpec::Feed {
                temperature_k: Some(inlet.0),
                pressure_pa: Some(inlet.1),
                molar_flow_mol_s: Some(inlet.2),
            }
        ));
        flowsheet.try_set_spec(unit, spec).unwrap();
        flowsheet
            .connect(PortRef::new(source, 0), PortRef::new(unit, 0))
            .unwrap();
        let outlet = flowsheet
            .connect(PortRef::new(unit, 1), PortRef::new(product, 0))
            .unwrap();
        let boundary = |flowsheet: &mut Flowsheet, port: usize, into_unit: bool| {
            if into_unit {
                let e = flowsheet.add_unit(UnitKind::EnergySource);
                assert!(flowsheet.set_spec(
                    e,
                    UnitSpec::EnergySource {
                        duty: EnergySourceDuty::Calculated
                    }
                ));
                flowsheet
                    .connect(PortRef::new(e, 0), PortRef::new(unit, port))
                    .unwrap()
            } else {
                let e = flowsheet.add_unit(UnitKind::EnergySink);
                flowsheet
                    .connect(PortRef::new(unit, port), PortRef::new(e, 0))
                    .unwrap()
            }
        };
        let power = boundary(&mut flowsheet, POWER, kind != UnitKind::Expander);
        let intercooler = (kind == UnitKind::MultistageCompressor)
            .then(|| boundary(&mut flowsheet, INTERCOOLER, false));
        let n = ids.len();
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition(ids, z),
            variant,
            vec![vec![0.0; n]; n],
            "explicit zero interactions, as in the rotating reference",
        )
        .unwrap();
        let mut feeds = FeedCompositions::new();
        feeds.insert(source, composition(ids, z));
        Self {
            flowsheet,
            package: CubicPackage::new(model, Default::default()),
            feeds,
            unit,
            outlet,
            power,
            intercooler,
        }
    }

    fn solve(&self) -> Result<Calculation, SolveError> {
        self.flowsheet.solve(&self.package, &self.feeds)
    }
}

fn check(fixture: &Value, case: &str, quantity: &str, tolerance: &str, actual: f64, expected: f64) {
    let t: Vec<f64> = fixture["tolerances"][tolerance]
        .as_array()
        .unwrap()
        .iter()
        .map(number)
        .collect();
    let allowed = t[0] + t[1] * expected.abs();
    assert!(
        (actual - expected).abs() <= allowed,
        "{case} {quantity}: {actual} vs reference {expected} (allowed {allowed})"
    );
}

fn variant(case: &Value) -> CubicVariant {
    match case["variant"].as_str().unwrap() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        other => panic!("unknown variant {other}"),
    }
}

fn machine_for(case: &Value) -> Machine {
    let ids = ids(case);
    let z: Vec<f64> = case["z"].as_array().unwrap().iter().map(number).collect();
    let inlet = &case["inlet"];
    let flow = number(&case["molar_flow_mol_s"]);
    let efficiency = number(&case["efficiency"]);
    let p_out = number(&case["outlet"]["p"]);
    let (kind, spec) = match case["kind"].as_str().unwrap() {
        "pump" => (
            UnitKind::Pump,
            UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(Some(p_out)),
                efficiency,
            },
        ),
        "compressor" => (
            UnitKind::Compressor,
            UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(Some(p_out)),
                efficiency,
            },
        ),
        "expander" => (
            UnitKind::Expander,
            UnitSpec::Expand {
                pressure: PressureSpec::OutletPressure(Some(p_out)),
                efficiency,
            },
        ),
        "multistage_compressor" => (
            UnitKind::MultistageCompressor,
            UnitSpec::Rotating(RotatingSpec::MultistageCompressor(
                MultistageCompressorSpec {
                    stages: case["stages"].as_u64().unwrap() as u16,
                    outlet_pressure_pa: Some(p_out),
                    efficiency,
                    intercooler_temperature_k: Some(number(&case["intercooler_T"])),
                },
            )),
        ),
        other => panic!("unknown kind {other}"),
    };
    Machine::new(
        &ids,
        &z,
        variant(case),
        (number(&inlet["T"]), number(&inlet["p"]), flow),
        kind,
        spec,
    )
}

#[test]
fn rotating_units_match_independent_isentropic_efficiency_references() {
    let fixture = fixture();
    let cases = fixture["units"].as_array().unwrap();
    assert_eq!(cases.len(), 10);
    for case in cases {
        let id = format!(
            "{} {}",
            case["variant"].as_str().unwrap(),
            case["id"].as_str().unwrap()
        );
        let m = machine_for(case);
        let calc = m.solve().unwrap_or_else(|e| panic!("{id}: {e:?}"));
        let unit = &calc.units[&m.unit];
        let out = calc.streams[&m.outlet].state.as_ref().unwrap();
        let reference = &case["outlet"];
        check(
            &fixture,
            &id,
            "outlet T",
            "temperature_K",
            out.temperature_k,
            number(&reference["T"]),
        );
        check(
            &fixture,
            &id,
            "outlet h",
            "enthalpy_J_per_mol",
            out.enthalpy(),
            number(&reference["target_h"]),
        );
        check(
            &fixture,
            &id,
            "outlet s",
            "entropy_J_per_mol_K",
            out.entropy(),
            number(&reference["target_s"]),
        );
        assert_eq!(out.pressure_pa, number(&reference["p"]), "{id}");
        let phases = reference["phases"].as_u64().unwrap() as usize;
        assert_eq!(out.phases.len(), phases, "{id} outlet phases");
        if phases == 2 {
            check(
                &fixture,
                &id,
                "outlet vapor fraction",
                "phase_fraction",
                out.vapor_fraction().unwrap(),
                number(&reference["vapor_fraction"]),
            );
        }

        let power = number(&case["power_W"]);
        let result = unit.rotating.as_ref().expect("rotating result");
        check(
            &fixture,
            &id,
            "shaft power",
            "power_W",
            result.power_w(),
            power,
        );
        // Power leaves an expander through its energy outlet.
        check(
            &fixture,
            &id,
            "power stream",
            "power_W",
            calc.energy_streams[&m.power].duty_w,
            power.abs(),
        );
        let balance = unit.balance.as_ref().unwrap();
        assert!(balance.max_scaled_component_residual <= 1e-8, "{id}");
        assert!(balance.scaled_energy_residual <= 2e-9, "{id}");

        if let Some(stages) = case["stage_outlets"].as_array() {
            assert_eq!(result.stages.len(), stages.len(), "{id}");
            for (k, (stage, reference)) in result.stages.iter().zip(stages).enumerate() {
                let label = format!("stage {} discharge T", k + 1);
                check(
                    &fixture,
                    &id,
                    &label,
                    "temperature_K",
                    stage.outlet_temperature_k,
                    number(&reference["T"]),
                );
                check(
                    &fixture,
                    &id,
                    "stage pressure",
                    "power_W",
                    stage.outlet_pressure_pa,
                    number(&reference["p"]),
                );
            }
            for (stage, reference) in result
                .stages
                .iter()
                .zip(case["stage_power_W"].as_array().unwrap())
            {
                check(
                    &fixture,
                    &id,
                    "stage power",
                    "power_W",
                    stage.power_w,
                    number(reference),
                );
            }
            let removed = -number(&case["intercooler_duty_W"]);
            check(
                &fixture,
                &id,
                "intercooler duty",
                "power_W",
                result.intercooler_duty_w,
                removed,
            );
            check(
                &fixture,
                &id,
                "intercooler stream",
                "power_W",
                calc.energy_streams[&m.intercooler.unwrap()].duty_w,
                removed,
            );
            check(
                &fixture,
                &id,
                "net duty",
                "power_W",
                unit.duty_w,
                power - removed,
            );
        } else {
            let ideal = &case["isentropic_outlet"];
            assert_eq!(result.stages.len(), 1);
            let stage = &result.stages[0];
            check(
                &fixture,
                &id,
                "isentropic outlet T",
                "temperature_K",
                stage.isentropic_outlet_temperature_k,
                number(&ideal["T"]),
            );
            let isentropic_power = number(&case["molar_flow_mol_s"])
                * (number(&ideal["target_h"]) - number(&case["inlet"]["target_h"]));
            check(
                &fixture,
                &id,
                "isentropic power",
                "power_W",
                stage.isentropic_power_w,
                isentropic_power,
            );
            check(&fixture, &id, "net duty", "power_W", unit.duty_w, power);
        }
    }
}

fn ng() -> (Vec<&'static str>, Vec<f64>) {
    (
        vec!["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
        vec![0.7, 0.15, 0.13, 0.01, 0.01],
    )
}

fn compressor(p_out: Option<f64>, efficiency: f64, inlet: (f64, f64, f64)) -> Machine {
    let (ids, z) = ng();
    Machine::new(
        &ids,
        &z,
        CubicVariant::PengRobinsonBostonMathias,
        inlet,
        UnitKind::Compressor,
        UnitSpec::Pressurize {
            pressure: PressureRise::OutletPressure(p_out),
            efficiency,
        },
    )
}

#[test]
fn unit_efficiency_reproduces_the_isentropic_outlet() {
    let m = compressor(Some(6e6), 1.0, (300.0, 2e6, 10.0));
    let calc = m.solve().unwrap();
    let stage = &calc.units[&m.unit].rotating.as_ref().unwrap().stages[0];
    let out = calc.streams[&m.outlet].state.as_ref().unwrap();
    assert_eq!(
        stage.outlet_temperature_k,
        stage.isentropic_outlet_temperature_k
    );
    assert_eq!(stage.power_w, stage.isentropic_power_w);
    assert_eq!(out.temperature_k, stage.outlet_temperature_k);
}

#[test]
fn pressure_increase_and_unchanged_pressure() {
    let (ids, z) = ng();
    let spec = |dp| UnitSpec::Pressurize {
        pressure: PressureRise::PressureIncrease(Some(dp)),
        efficiency: 0.75,
    };
    let m = Machine::new(
        &ids,
        &z,
        CubicVariant::PengRobinson,
        (300.0, 2e6, 10.0),
        UnitKind::Compressor,
        spec(4e6),
    );
    let by_increase = m.solve().unwrap();
    let reference = fixture()["units"][0].clone();
    assert_eq!(reference["id"], "ng-compressor");
    assert_eq!(reference["variant"], "PR");
    let out = by_increase.streams[&m.outlet].state.as_ref().unwrap();
    assert_eq!(out.pressure_pa, 6e6);
    check(
        &fixture(),
        "increase",
        "outlet T",
        "temperature_K",
        out.temperature_k,
        number(&reference["outlet"]["T"]),
    );

    let m = Machine::new(
        &ids,
        &z,
        CubicVariant::PengRobinson,
        (300.0, 2e6, 10.0),
        UnitKind::Compressor,
        spec(0.0),
    );
    let calc = m.solve().unwrap();
    let inlet = calc
        .streams
        .values()
        .find(|s| s.state.as_ref().is_some_and(|s| s.pressure_pa == 2e6))
        .unwrap();
    let out = calc.streams[&m.outlet].state.as_ref().unwrap();
    assert_eq!(
        out.temperature_k,
        inlet.state.as_ref().unwrap().temperature_k
    );
    assert_eq!(calc.units[&m.unit].duty_w, 0.0);
}

#[test]
fn inlet_phase_rules_and_pressure_direction_fail_explicitly() {
    // Natural gas at 230 K and 2 MPa splits into vapor and liquid.
    let wet = compressor(Some(6e6), 0.75, (230.0, 2e6, 10.0));
    let e = wet.solve().unwrap_err();
    assert_eq!(e.unit, Some(wet.unit));
    assert!(
        e.message.contains("compressor inlet contains liquid"),
        "{}",
        e.message
    );

    let (ids, z) = ng();
    let pump = Machine::new(
        &ids,
        &z,
        CubicVariant::PengRobinson,
        (300.0, 2e6, 10.0),
        UnitKind::Pump,
        UnitSpec::Pressurize {
            pressure: PressureRise::OutletPressure(Some(6e6)),
            efficiency: 0.7,
        },
    );
    assert!(
        pump.solve()
            .unwrap_err()
            .message
            .contains("pump inlet contains vapor")
    );

    let backwards = compressor(Some(1e6), 0.75, (300.0, 2e6, 10.0));
    assert!(
        backwards
            .solve()
            .unwrap_err()
            .message
            .contains("below the inlet pressure")
    );

    let expander = Machine::new(
        &ids,
        &z,
        CubicVariant::PengRobinson,
        (300.0, 2e6, 10.0),
        UnitKind::Expander,
        UnitSpec::Expand {
            pressure: PressureSpec::OutletPressure(Some(3e6)),
            efficiency: 0.8,
        },
    );
    assert!(
        expander
            .solve()
            .unwrap_err()
            .message
            .contains("exceeds the inlet pressure")
    );

    let missing = compressor(None, 0.75, (300.0, 2e6, 10.0));
    assert!(missing.solve().is_err());
    let invalid = compressor(Some(6e6), 1.2, (300.0, 2e6, 10.0));
    assert!(invalid.solve().unwrap_err().message.contains("efficiency"));
}

#[test]
fn zero_flow_gives_an_empty_outlet_and_zero_power() {
    let m = compressor(Some(6e6), 0.75, (300.0, 2e6, 0.0));
    let calc = m.solve().unwrap();
    let unit = &calc.units[&m.unit];
    assert_eq!(unit.duty_w, 0.0);
    assert_eq!(calc.energy_streams[&m.power].duty_w, 0.0);
    assert!(calc.streams[&m.outlet].state.is_none());
    assert!(unit.warnings.iter().any(|w| w.contains("Zero inlet flow")));
}

#[test]
fn intercooler_condensation_and_single_stage_operation() {
    let (ids, z) = ng();
    let multistage = |stages, t| {
        Machine::new(
            &ids,
            &z,
            CubicVariant::PengRobinsonBostonMathias,
            (300.0, 2e6, 10.0),
            UnitKind::MultistageCompressor,
            UnitSpec::Rotating(RotatingSpec::MultistageCompressor(
                MultistageCompressorSpec {
                    stages,
                    outlet_pressure_pa: Some(6e6),
                    efficiency: 0.75,
                    intercooler_temperature_k: t,
                },
            )),
        )
    };
    // Cooling a 3.5 MPa rich gas to 230 K condenses liquid ahead of stage 2.
    let wet = multistage(2, Some(230.0));
    let e = wet.solve().unwrap_err();
    assert!(e.message.contains("condenses liquid"), "{}", e.message);
    let hot = multistage(2, Some(500.0));
    assert!(
        hot.solve()
            .unwrap_err()
            .message
            .contains("exceeds the stage discharge")
    );

    // One stage equals the single compressor.
    let single = multistage(1, None).solve().unwrap();
    let reference = compressor(Some(6e6), 0.75, (300.0, 2e6, 10.0))
        .solve()
        .unwrap();
    let power = |c: &Calculation| {
        c.units
            .values()
            .find_map(|u| u.rotating.as_ref())
            .unwrap()
            .power_w()
    };
    assert!((power(&single) - power(&reference)).abs() <= 1e-9 * power(&reference));
    let m = multistage(1, None);
    let calc = m.solve().unwrap();
    assert_eq!(calc.energy_streams[&m.intercooler.unwrap()].duty_w, 0.0);
}

#[test]
fn expander_reports_a_two_phase_outlet() {
    let fixture = fixture();
    let case = fixture["units"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "propane-two-phase-expander")
        .unwrap();
    let m = machine_for(case);
    let calc = m.solve().unwrap();
    let out = calc.streams[&m.outlet].state.as_ref().unwrap();
    assert!(out.phases.iter().any(|p| p.phase == Phase::Liquid));
    assert!(
        calc.units[&m.unit]
            .warnings
            .iter()
            .any(|w| w.contains("two-phase"))
    );
}
