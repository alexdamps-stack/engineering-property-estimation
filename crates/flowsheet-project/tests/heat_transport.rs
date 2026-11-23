//! #179 and #180 through the normal solver and version 4 records.
use std::ops::ControlFlow;

use flowsheet_core::units::HeatTransportSpec;
use flowsheet_core::units::heat_transport::{AirCoolerSpec, AirCoolerThermal, FiredHeaterSpec};
use flowsheet_core::{
    EnergySourceDuty, PortRef, StreamId, ThermalSpec, UnitId, UnitKind, UnitSpec,
};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};

struct ProcessLine {
    project: Project,
    unit: UnitId,
    source: UnitId,
    inlet: StreamId,
    outlet: StreamId,
    energy: StreamId,
}

impl ProcessLine {
    fn new(kind: UnitKind, package: PackageSelection, z: Vec<f64>, t: f64, p: f64, n: f64) -> Self {
        let mut project = Project::new(package);
        let fs = &mut project.flowsheet;
        let feed = fs.add_unit(UnitKind::Feed);
        fs.set_spec(
            feed,
            UnitSpec::Feed {
                temperature_k: Some(t),
                pressure_pa: Some(p),
                molar_flow_mol_s: Some(n),
            },
        );
        project.feeds.insert(feed, z);
        let unit = fs.add_unit(kind);
        let product = fs.add_unit(UnitKind::Product);
        let source = fs.add_unit(UnitKind::EnergySource);
        let inlet = fs
            .connect(PortRef::new(feed, 0), PortRef::new(unit, 0))
            .unwrap();
        let outlet = fs
            .connect(PortRef::new(unit, 1), PortRef::new(product, 0))
            .unwrap();
        let energy = fs
            .connect(PortRef::new(source, 0), PortRef::new(unit, 2))
            .unwrap();
        Self {
            project,
            unit,
            source,
            inlet,
            outlet,
            energy,
        }
    }

    fn air(n: f64) -> Self {
        Self::new(
            UnitKind::AirCooler,
            PackageSelection::zero_interactions(
                CubicVariantChoice::PengRobinson,
                vec!["methane".into()],
                "single-component air-cooler test",
            ),
            vec![1.0],
            360.0,
            2e6,
            n,
        )
    }

    fn set(&mut self, spec: HeatTransportSpec) {
        self.project
            .flowsheet
            .set_spec(self.unit, UnitSpec::HeatTransport(spec));
    }

    fn source(&mut self, duty: EnergySourceDuty) {
        self.project
            .flowsheet
            .set_spec(self.source, UnitSpec::EnergySource { duty });
    }

    fn solve(&self) -> Result<ResultSet, SolveFailure> {
        solve_project(
            &self.project,
            &self.project.package.prepare().unwrap(),
            |_| ControlFlow::Continue(()),
            0,
        )
    }

    fn check_balances(&self, result: &ResultSet) {
        let inlet = &result.streams[&self.inlet];
        let outlet = &result.streams[&self.outlet];
        assert_eq!(outlet.molar_flow_mol_s, inlet.molar_flow_mol_s);
        if inlet.molar_flow_mol_s > 0.0 {
            assert_eq!(
                outlet.state.as_ref().unwrap().mole_fractions,
                inlet.state.as_ref().unwrap().mole_fractions,
            );
        }
        let unit = &result.units[&self.unit];
        assert!(unit.max_scaled_component_residual.unwrap() <= 1e-8);
        assert!(unit.scaled_energy_residual.unwrap().abs() <= 2e-9);
    }

    fn reopen(&mut self) {
        let saved = self.project.to_json().unwrap();
        self.project = Project::from_json(&saved).unwrap();
        assert_eq!(self.project.to_json().unwrap(), saved);
    }
}

#[test]
fn air_cooler_keeps_process_duty_and_fan_work_separate() {
    let mut line = ProcessLine::air(25.0);
    let mut spec = AirCoolerSpec {
        thermal: AirCoolerThermal::OutletTemperature(Some(320.0)),
        pressure_drop_pa: 20_000.0,
        ambient_temperature_k: Some(300.0),
        air_mass_flow_kg_s: Some(5.0),
        fan_power_w: 250.0,
    };
    line.set(HeatTransportSpec::AirCooler(spec.clone()));
    line.reopen();
    let result = line.solve().unwrap();
    line.check_balances(&result);
    let outlet = result.streams[&line.outlet].state.as_ref().unwrap();
    assert_eq!(outlet.temperature_k, 320.0);
    assert_eq!(outlet.pressure_pa, 1_980_000.0);
    let duty = result.units[&line.unit].duty_w;
    assert!(duty < 0.0);
    assert_eq!(result.energy_streams[&line.energy].duty_w, 250.0);
    let air_t = spec.air_outlet_temperature_k(duty).unwrap();
    let air_gain = 5_025.0 * (air_t - 300.0);
    assert!((air_gain + duty - 250.0).abs() < 1e-8);
    line.project.results = Some(result);
    line.reopen();
    let saved = line.project.results.take().unwrap();
    assert_eq!(saved.units[&line.unit].duty_w, duty);
    assert_eq!(saved.energy_streams[&line.energy].duty_w, 250.0);

    spec.thermal = AirCoolerThermal::Duty(Some(duty));
    line.set(HeatTransportSpec::AirCooler(spec.clone()));
    line.reopen();
    let result = line.solve().unwrap();
    line.check_balances(&result);
    let t = result.streams[&line.outlet]
        .state
        .as_ref()
        .unwrap()
        .temperature_k;
    assert!((t - 320.0).abs() <= 2e-6, "{t}");
    line.source(EnergySourceDuty::Specified(Some(250.0)));
    assert!(line.solve().is_err());

    spec.thermal = AirCoolerThermal::OutletTemperature(Some(290.0));
    line.source(EnergySourceDuty::Calculated);
    line.set(HeatTransportSpec::AirCooler(spec.clone()));
    assert_eq!(line.solve().unwrap_err().unit, Some(line.unit));
    spec.thermal = AirCoolerThermal::OutletTemperature(Some(320.0));
    spec.air_mass_flow_kg_s = Some(0.001);
    line.set(HeatTransportSpec::AirCooler(spec.clone()));
    assert_eq!(line.solve().unwrap_err().unit, Some(line.unit));
    spec.fan_power_w = f64::INFINITY;
    line.set(HeatTransportSpec::AirCooler(spec));
    assert!(line.project.to_json().is_err());
}

#[test]
fn zero_process_flow_has_no_process_duty_but_fan_can_run() {
    let mut line = ProcessLine::air(0.0);
    line.set(HeatTransportSpec::AirCooler(AirCoolerSpec {
        thermal: AirCoolerThermal::OutletTemperature(Some(320.0)),
        ambient_temperature_k: Some(300.0),
        air_mass_flow_kg_s: Some(1.0),
        fan_power_w: 502.5,
        ..Default::default()
    }));
    let result = line.solve().unwrap();
    line.check_balances(&result);
    assert_eq!(result.units[&line.unit].duty_w, 0.0);
    assert_eq!(result.streams[&line.outlet].state, None);
    assert_eq!(result.energy_streams[&line.energy].duty_w, 502.5);
}

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

#[test]
fn fired_heater_reverses_independent_cooling_reference() {
    // Reverse the independent process fixture's constant-pressure cooler.
    // Composition and total molar flow stay fixed across both directions.
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/process-v1.json")).unwrap();
    for c in fixture["cases"].as_array().unwrap() {
        let streams = c["streams"].as_array().unwrap();
        let cold = streams.iter().find(|s| s["id"] == "cool").unwrap();
        let hot = streams.iter().find(|s| s["id"] == "mix").unwrap();
        let variant = if c["variant"] == "PR" {
            CubicVariantChoice::PengRobinson
        } else {
            CubicVariantChoice::PengRobinsonBostonMathias
        };
        let package = PackageSelection::zero_interactions(
            variant,
            cold["ids"]
                .as_array()
                .unwrap()
                .iter()
                .map(|id| id.as_str().unwrap().into())
                .collect(),
            "matched zero-interaction approximation from the independent process fixture",
        );
        let n = number(&cold["flow_mol_per_s"]);
        let mut line = ProcessLine::new(
            UnitKind::FiredHeater,
            package,
            cold["z"].as_array().unwrap().iter().map(number).collect(),
            number(&cold["T"]),
            number(&cold["p"]),
            n,
        );
        let expected_duty = n * (number(&hot["target_h"]) - number(&cold["target_h"]));
        let expected_t = number(&hot["T"]);
        for thermal in [
            ThermalSpec::OutletTemperature(Some(expected_t)),
            ThermalSpec::Duty(Some(expected_duty)),
            ThermalSpec::EnergyStream,
        ] {
            line.set(HeatTransportSpec::FiredHeater(FiredHeaterSpec {
                thermal,
                pressure_drop_pa: 0.0,
            }));
            if thermal == ThermalSpec::EnergyStream {
                // The connection is retained when switching its determining end.
                assert!(line.solve().is_err());
                line.source(EnergySourceDuty::Specified(Some(expected_duty)));
            }
            line.reopen();
            let result = line.solve().unwrap();
            line.check_balances(&result);
            let duty = result.units[&line.unit].duty_w;
            assert!((duty - expected_duty).abs() / n <= 4e-5);
            assert_eq!(duty, result.energy_streams[&line.energy].duty_w);
            let outlet = result.streams[&line.outlet].state.as_ref().unwrap();
            assert!((outlet.temperature_k - expected_t).abs() <= 2e-6);
            let expected_h = number(&hot["target_h"]);
            assert!((outlet.enthalpy_j_mol - expected_h).abs() <= 2e-5 + 2e-9 * expected_h.abs());
        }
    }
}

#[test]
fn legacy_v4_heat_transport_records_load_with_inputs_unset() {
    for (kind, name) in [
        (UnitKind::AirCooler, "air_cooler"),
        (UnitKind::FiredHeater, "fired_heater"),
    ] {
        let mut line = ProcessLine::air(25.0);
        if kind == UnitKind::FiredHeater {
            line = ProcessLine::new(kind, line.project.package, vec![1.0], 360.0, 2e6, 25.0);
        }
        let mut saved: Value = serde_json::from_str(&line.project.to_json().unwrap()).unwrap();
        assert_eq!(saved["format_version"], 4);
        let unit = saved["units"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|u| u["kind"] == name)
            .unwrap();
        unit["spec"]["unit"] = json!({"kind": name});
        let original_ports = line.project.flowsheet.unit(line.unit).unwrap().ports();
        line.project = Project::from_json(&saved.to_string()).unwrap();
        let loaded = line.project.flowsheet.unit(line.unit).unwrap();
        assert_eq!(loaded.spec, kind.default_spec());
        assert_eq!(loaded.ports(), original_ports);
        assert!(
            line.solve()
                .unwrap_err()
                .message
                .contains("outlet temperature")
        );
    }
}
