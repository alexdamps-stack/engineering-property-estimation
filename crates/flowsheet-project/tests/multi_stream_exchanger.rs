//! Fixed terminal approaches, v4 records and normal project solving for #181.

use std::ops::ControlFlow;

use flowsheet_core::units::HeatTransportSpec;
use flowsheet_core::units::heat_transport::{ExchangerPassThermal, MultiStreamExchangerSpec};
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};

struct ExchangerProject {
    project: Project,
    exchanger: UnitId,
    inlets: [StreamId; 3],
    outlets: [StreamId; 3],
}

impl ExchangerProject {
    fn new() -> Self {
        let mut project = Project::new(PackageSelection::zero_interactions(
            CubicVariantChoice::PengRobinson,
            vec!["methane".into(), "nitrogen".into()],
            "explicit zero interactions for the multi-stream project test",
        ));
        let fs = &mut project.flowsheet;
        let exchanger = fs.add_unit(UnitKind::MultiStreamExchanger);
        let connections: [(StreamId, StreamId); 3] = std::array::from_fn(|i| {
            let feed = fs.add_unit(UnitKind::Feed);
            let product = fs.add_unit(UnitKind::Product);
            let (t, n, z) = [
                (410.0, 2.0, [1.0, 0.0]),
                (300.0, 1.0, [0.0, 1.0]),
                (290.0, 1.5, [0.6, 0.4]),
            ][i];
            fs.set_spec(
                feed,
                UnitSpec::Feed {
                    temperature_k: Some(t),
                    pressure_pa: Some(2e6),
                    molar_flow_mol_s: Some(n),
                },
            );
            project.feeds.insert(feed, z.to_vec());
            let inlet = fs
                .connect(PortRef::new(feed, 0), PortRef::new(exchanger, i * 2))
                .unwrap();
            let outlet = fs
                .connect(PortRef::new(exchanger, i * 2 + 1), PortRef::new(product, 0))
                .unwrap();
            (inlet, outlet)
        });
        Self {
            project,
            exchanger,
            inlets: connections.map(|(inlet, _)| inlet),
            outlets: connections.map(|(_, outlet)| outlet),
        }
    }

    fn set(&mut self, spec: MultiStreamExchangerSpec) {
        self.project.flowsheet.set_spec(
            self.exchanger,
            UnitSpec::HeatTransport(HeatTransportSpec::MultiStreamExchanger(spec)),
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

    fn reopen(&mut self) {
        let saved = self.project.to_json().unwrap();
        self.project = Project::from_json(&saved).unwrap();
        assert_eq!(self.project.to_json().unwrap(), saved);
    }
}

fn fixed_spec() -> MultiStreamExchangerSpec {
    let mut spec = MultiStreamExchangerSpec::default();
    spec.passes[0].thermal = ExchangerPassThermal::TerminalApproach {
        reference_pass: 2,
        approach_k: Some(65.0),
    };
    spec.passes[1].thermal = ExchangerPassThermal::TerminalApproach {
        reference_pass: 1,
        approach_k: Some(80.0),
    };
    spec.passes[0].pressure_drop_pa = 20_000.0;
    spec.passes[2].pressure_drop_pa = 5_000.0;
    spec.minimum_approach_k = 5.0;
    spec
}

#[test]
fn fixed_approaches_reach_both_terminals_and_balance_after_v4_reload() {
    let mut line = ExchangerProject::new();
    line.set(fixed_spec());
    line.reopen();
    let results = line.solve().unwrap();
    let mut duties = [0.0; 3];
    for (k, duty) in duties.iter_mut().enumerate() {
        let inlet = &results.streams[&line.inlets[k]];
        let outlet = &results.streams[&line.outlets[k]];
        assert_eq!(inlet.molar_flow_mol_s, outlet.molar_flow_mol_s);
        let a = inlet.state.as_ref().unwrap();
        let b = outlet.state.as_ref().unwrap();
        assert_eq!(a.mole_fractions, b.mole_fractions);
        assert_eq!(b.pressure_pa, [1_980_000.0, 2_000_000.0, 1_995_000.0][k]);
        *duty = inlet.molar_flow_mol_s * (b.enthalpy_j_mol - a.enthalpy_j_mol);
    }
    let t = |stream: StreamId| {
        results.streams[&stream]
            .state
            .as_ref()
            .unwrap()
            .temperature_k
    };
    assert!((t(line.outlets[0]) - t(line.inlets[1]) - 65.0).abs() <= 1e-9);
    assert!((t(line.inlets[0]) - t(line.outlets[1]) - 80.0).abs() <= 1e-9);
    assert!(duties[0] < 0.0 && duties[1] > 0.0 && duties[2] > 0.0);
    assert!(duties.iter().sum::<f64>().abs() <= 1e-5, "{duties:?}");
    let unit = &results.units[&line.exchanger];
    assert_eq!(unit.duty_w, 0.0);
    assert!((unit.transferred_duty_w.unwrap() + duties[0]).abs() <= 1e-5);
    assert!(unit.max_scaled_component_residual.unwrap() <= 1e-8);
    assert!(unit.scaled_energy_residual.unwrap().abs() <= 2e-9);
    line.project.results = Some(results);
    line.reopen();
    assert!(
        line.project.results.unwrap().units[&line.exchanger]
            .transferred_duty_w
            .unwrap()
            > 0.0
    );
}

#[test]
fn rejects_infeasible_and_over_or_under_specified_approaches() {
    let mut line = ExchangerProject::new();
    let mut spec = fixed_spec();
    spec.minimum_approach_k = 66.0;
    line.set(spec);
    let error = line.solve().unwrap_err();
    assert_eq!(error.unit, Some(line.exchanger));
    assert!(error.message.contains("approach"), "{error:?}");

    let mut spec = fixed_spec();
    spec.passes[1].thermal = ExchangerPassThermal::TerminalApproach {
        reference_pass: 1,
        approach_k: Some(500.0),
    };
    line.set(spec);
    assert_eq!(line.solve().unwrap_err().unit, Some(line.exchanger));

    for (reference_pass, approach_k) in [(3, Some(80.0)), (8, Some(80.0)), (1, None)] {
        let mut spec = fixed_spec();
        spec.passes[1].thermal = ExchangerPassThermal::TerminalApproach {
            reference_pass,
            approach_k,
        };
        line.set(spec);
        assert!(line.solve().is_err());
    }

    let mut spec = fixed_spec();
    spec.passes[2].thermal = ExchangerPassThermal::Duty(Some(1_000.0));
    line.set(spec);
    assert!(line.solve().is_err()); // No remaining degree of freedom for the heat balance.
}

#[test]
fn a_fixed_approach_requires_a_flowing_reference_pass() {
    let mut line = ExchangerProject::new();
    let feed = line
        .project
        .flowsheet
        .stream(line.inlets[1])
        .unwrap()
        .from
        .unit;
    line.project.flowsheet.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(300.0),
            pressure_pa: Some(2e6),
            molar_flow_mol_s: Some(0.0),
        },
    );
    let mut spec = fixed_spec();
    spec.passes[1].thermal = ExchangerPassThermal::Duty(Some(0.0));
    line.set(spec);
    let error = line.solve().unwrap_err();
    assert_eq!(error.unit, Some(line.exchanger));
    assert!(error.message.contains("reference"), "{error:?}");
}

#[test]
fn legacy_v4_count_and_connections_survive_without_thermal_settings() {
    let mut line = ExchangerProject::new();
    line.set(fixed_spec());
    let mut value: Value = serde_json::from_str(&line.project.to_json().unwrap()).unwrap();
    assert_eq!(value["format_version"], 4);
    let unit = value["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "multi_stream_exchanger")
        .unwrap();
    let record = &mut unit["spec"]["unit"];
    // A payload cannot specify an approach and an outlet temperature together.
    record["pass_specs"][0]["thermal"]["terminal_approach"]["outlet_temperature_k"] = json!(365.0);
    assert!(Project::from_json(&value.to_string()).is_err());
    let record = &mut value["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "multi_stream_exchanger")
        .unwrap()["spec"]["unit"];
    *record = json!({"kind": "multi_stream_exchanger", "passes": 3});
    line.project = Project::from_json(&value.to_string()).unwrap();
    for (k, port) in line.outlets.iter().enumerate() {
        let stream = line.project.flowsheet.stream(*port).unwrap();
        assert_eq!(stream.from, PortRef::new(line.exchanger, k * 2 + 1));
    }
    assert_eq!(
        line.project.flowsheet.unit(line.exchanger).unwrap().spec,
        UnitKind::MultiStreamExchanger.default_spec(),
    );
    assert!(line.solve().is_err());
}
