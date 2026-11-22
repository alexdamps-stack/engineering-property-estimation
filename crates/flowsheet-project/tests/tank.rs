//! #178 steady-state tank through project persistence and the normal solver.
use std::ops::ControlFlow;

use flowsheet_core::units::SimpleLogicalSpec;
use flowsheet_core::units::simple_logical::TankSpec;
use flowsheet_core::{PortRef, SeparatorThermal, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::Value;

struct TankLine {
    project: Project,
    tank: UnitId,
    inlet: StreamId,
    outlet: StreamId,
    energy: StreamId,
}

fn line(package: PackageSelection, z: Vec<f64>, t: f64, p: f64, n: f64) -> TankLine {
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
    let tank = fs.add_unit(UnitKind::Tank);
    let product = fs.add_unit(UnitKind::Product);
    let source = fs.add_unit(UnitKind::EnergySource);
    let inlet = fs
        .connect(PortRef::new(feed, 0), PortRef::new(tank, 0))
        .unwrap();
    let outlet = fs
        .connect(PortRef::new(tank, 1), PortRef::new(product, 0))
        .unwrap();
    let energy = fs
        .connect(PortRef::new(source, 0), PortRef::new(tank, 2))
        .unwrap();
    TankLine {
        project,
        tank,
        inlet,
        outlet,
        energy,
    }
}

fn binary_line(n: f64) -> TankLine {
    line(
        PackageSelection::zero_interactions(
            CubicVariantChoice::PengRobinson,
            vec!["methane".into(), "propane".into()],
            "explicit zero-interaction approximation for the binary tank test",
        ),
        vec![0.5, 0.5],
        250.0,
        3e6,
        n,
    )
}

impl TankLine {
    fn set(&mut self, thermal: SeparatorThermal, pressure_drop_pa: f64) {
        self.project.flowsheet.set_spec(
            self.tank,
            UnitSpec::SimpleLogical(SimpleLogicalSpec::Tank(TankSpec {
                thermal,
                pressure_drop_pa,
            })),
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

    fn check_balances(&self, r: &ResultSet) {
        let inlet = &r.streams[&self.inlet];
        let outlet = &r.streams[&self.outlet];
        assert_eq!(outlet.molar_flow_mol_s, inlet.molar_flow_mol_s);
        assert_eq!(
            outlet.state.as_ref().map(|s| &s.mole_fractions),
            inlet.state.as_ref().map(|s| &s.mole_fractions)
        );
        let u = &r.units[&self.tank];
        assert_eq!(r.energy_streams[&self.energy].duty_w, u.duty_w);
        assert!(u.max_scaled_component_residual.unwrap() <= 1e-8);
        assert!(u.scaled_energy_residual.unwrap().abs() <= 2e-9);
    }
}

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

#[test]
fn adiabatic_tank_with_pressure_drop_matches_independent_isenthalpic_outlet() {
    let fixture: Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/process-v1.json")).unwrap();
    for c in fixture["cases"].as_array().unwrap() {
        // The reference valve and adiabatic tank use the same balances for
        // a single outlet and assigned pressure loss.
        let streams = c["streams"].as_array().unwrap();
        let inlet = streams.iter().find(|s| s["id"] == "cool").unwrap();
        let outlet = streams.iter().find(|s| s["id"] == "expand").unwrap();
        let variant = if c["variant"] == "PR" {
            CubicVariantChoice::PengRobinson
        } else {
            CubicVariantChoice::PengRobinsonBostonMathias
        };
        let package = PackageSelection::zero_interactions(
            variant,
            inlet["ids"]
                .as_array()
                .unwrap()
                .iter()
                .map(|id| id.as_str().unwrap().into())
                .collect(),
            "matched zero-interaction approximation from the independent process fixture",
        );
        let mut tank = line(
            package,
            inlet["z"].as_array().unwrap().iter().map(number).collect(),
            number(&inlet["T"]),
            number(&inlet["p"]),
            number(&inlet["flow_mol_per_s"]),
        );
        tank.set(
            SeparatorThermal::Adiabatic,
            number(&inlet["p"]) - number(&outlet["p"]),
        );
        let r = tank.solve().unwrap();
        tank.check_balances(&r);
        assert_eq!(r.units[&tank.tank].duty_w, 0.0);
        let state = r.streams[&tank.outlet].state.as_ref().unwrap();
        let close = |key: &str, actual: f64, expected: f64| {
            let b = &fixture["tolerances"][key];
            let limit = number(&b[0]) + number(&b[1]) * expected.abs();
            assert!(
                (actual - expected).abs() <= limit,
                "{key}: {actual} vs {expected}"
            );
        };
        assert_eq!(state.pressure_pa, number(&outlet["p"]));
        close("temperature_K", state.temperature_k, number(&outlet["T"]));
        close(
            "enthalpy_J_per_mol",
            state.enthalpy_j_mol,
            number(&outlet["target_h"]),
        );
        close(
            "phase_fraction",
            state.vapor_fraction.unwrap(),
            number(&outlet["vapor_fraction"]),
        );
    }
}

#[test]
fn outlet_temperature_and_signed_duty_round_trip_in_v4() {
    let mut tank = binary_line(10.0);
    for setpoint in [225.0, 340.0] {
        tank.set(SeparatorThermal::OutletTemperature(Some(setpoint)), 0.0);
        let r = tank.solve().unwrap();
        tank.check_balances(&r);
        let q = r.units[&tank.tank].duty_w;
        assert_eq!(q.signum(), (setpoint - 250.0).signum());
        assert_eq!(
            r.streams[&tank.outlet]
                .state
                .as_ref()
                .unwrap()
                .temperature_k,
            setpoint
        );
        tank.set(SeparatorThermal::Duty(Some(q)), 0.0);
        let json = tank.project.to_json().unwrap();
        tank.project = Project::from_json(&json).unwrap();
        assert_eq!(tank.project.to_json().unwrap(), json);
        let r = tank.solve().unwrap();
        tank.check_balances(&r);
        assert_eq!(r.units[&tank.tank].duty_w, q);
        let t = r.streams[&tank.outlet]
            .state
            .as_ref()
            .unwrap()
            .temperature_k;
        assert!((t - setpoint).abs() < 2e-6, "{t} vs {setpoint}");
    }
}

#[test]
fn legacy_tank_keeps_inlet_state_and_zero_flow_gives_zero_duty() {
    let mut tank = binary_line(10.0);
    let mut json: Value = serde_json::from_str(&tank.project.to_json().unwrap()).unwrap();
    let legacy = json["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "tank")
        .unwrap();
    legacy["spec"]["unit"] = serde_json::json!({"kind": "tank"});
    tank.project = Project::from_json(&json.to_string()).unwrap();
    let r = tank.solve().unwrap();
    tank.check_balances(&r);
    assert_eq!(r.streams[&tank.outlet].state, r.streams[&tank.inlet].state);

    let mut empty = binary_line(0.0);
    empty.set(SeparatorThermal::Duty(Some(100.0)), 0.0);
    let r = empty.solve().unwrap();
    assert_eq!(r.streams[&empty.outlet].molar_flow_mol_s, 0.0);
    assert_eq!(r.streams[&empty.outlet].state, None);
    assert_eq!(r.energy_streams[&empty.energy].duty_w, 0.0);
    assert!(r.units[&empty.tank].warnings[0].contains("Zero inlet flow"));
}

#[test]
fn invalid_pressure_drop_and_unset_thermal_spec_fail() {
    let mut tank = binary_line(10.0);
    tank.set(SeparatorThermal::Adiabatic, 3e6);
    let error = tank.solve().unwrap_err();
    assert_eq!(error.unit, Some(tank.tank));
    assert!(
        error
            .message
            .contains("nonpositive absolute outlet pressure")
    );

    tank.set(SeparatorThermal::OutletTemperature(None), 0.0);
    assert!(
        tank.solve()
            .unwrap_err()
            .message
            .contains("outlet temperature")
    );

    tank.set(SeparatorThermal::Duty(Some(f64::NAN)), 0.0);
    assert!(tank.project.to_json().is_err());
}
