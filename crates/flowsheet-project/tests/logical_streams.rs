//! A17 logical boundaries through the normal project solve and save paths.
use std::ops::ControlFlow;

use flowsheet_core::units::SimpleLogicalSpec;
use flowsheet_core::units::simple_logical::{MultiplierSpec, SelectorSpec};
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, Currency, PackageSelection, Project, ResultSet, SolveFailure, StreamResult,
    solve_project,
};

fn project() -> Project {
    Project::new(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        vec!["methane".into(), "propane".into()],
        "explicit zero-interaction approximation for logical routing tests",
    ))
}

fn feed(project: &mut Project, n: f64, t: f64, z: [f64; 2]) -> UnitId {
    let id = project.flowsheet.add_unit(UnitKind::Feed);
    project.flowsheet.set_spec(
        id,
        UnitSpec::Feed {
            temperature_k: Some(t),
            pressure_pa: Some(2e6),
            molar_flow_mol_s: Some(n),
        },
    );
    project.feeds.insert(id, z.to_vec());
    id
}

fn connect(project: &mut Project, from: UnitId, out: usize, to: UnitId, inlet: usize) -> StreamId {
    project
        .flowsheet
        .connect(PortRef::new(from, out), PortRef::new(to, inlet))
        .unwrap()
}

fn multiplier(factor: Option<f64>) -> UnitSpec {
    UnitSpec::SimpleLogical(SimpleLogicalSpec::Multiplier(MultiplierSpec { factor }))
}

fn selector(inlets: u16, selected_inlet: Option<u16>) -> UnitSpec {
    UnitSpec::SimpleLogical(SimpleLogicalSpec::Selector(SelectorSpec {
        inlets,
        selected_inlet,
    }))
}

fn solve(project: &Project) -> Result<ResultSet, SolveFailure> {
    solve_project(
        project,
        &project.package.prepare().unwrap(),
        |_| ControlFlow::Continue(()),
        0,
    )
}

struct LogicalLine {
    project: Project,
    multiplier: UnitId,
    selector: UnitId,
    fresh: StreamId,
    comparison: StreamId,
    scaled: StreamId,
    product: StreamId,
}

fn logical_line() -> LogicalLine {
    let mut project = project();
    let first = feed(&mut project, 2.0, 300.0, [0.75, 0.25]);
    let other = feed(&mut project, 3.0, 360.0, [0.0, 1.0]);
    let multiplier = project.flowsheet.add_unit(UnitKind::Multiplier);
    let selector = project.flowsheet.add_unit(UnitKind::Selector);
    let output = project.flowsheet.add_unit(UnitKind::Product);
    let fresh = connect(&mut project, first, 0, multiplier, 0);
    let comparison = connect(&mut project, other, 0, selector, 2);
    let scaled = connect(&mut project, multiplier, 1, selector, 1);
    let product = connect(&mut project, selector, 0, output, 0);
    LogicalLine {
        project,
        multiplier,
        selector,
        fresh,
        comparison,
        scaled,
        product,
    }
}

fn enthalpy_flow(stream: &StreamResult) -> f64 {
    stream
        .state
        .as_ref()
        .map_or(0.0, |s| stream.molar_flow_mol_s * s.enthalpy_j_mol)
}

fn component_flow(stream: &StreamResult, i: usize) -> f64 {
    stream
        .state
        .as_ref()
        .map_or(0.0, |s| stream.molar_flow_mol_s * s.mole_fractions[i])
}

fn close(actual: f64, expected: f64, tol: f64) {
    assert!((actual - expected).abs() <= tol, "{actual} ≠ {expected}");
}

#[test]
fn scaled_and_unselected_streams_close_the_whole_flowsheet() {
    let mut line = logical_line();
    for (factor, selected) in [(2.0, 1), (0.5, 1), (1.0, 2), (2.0, 2), (0.0, 1)] {
        let fs = &mut line.project.flowsheet;
        fs.set_spec(line.multiplier, multiplier(Some(factor)));
        fs.set_spec(line.selector, selector(2, Some(selected)));
        line.project.results = Some(solve(&line.project).unwrap());
        let loaded = Project::from_json(&line.project.to_json().unwrap()).unwrap();
        assert!(
            loaded
                .flowsheet
                .same_calculation_inputs(&line.project.flowsheet)
        );
        assert_eq!(loaded.results, line.project.results);
        let fingerprints = loaded.fingerprints();
        let result = loaded.results.unwrap();
        let fresh = &result.streams[&line.fresh];
        let other = &result.streams[&line.comparison];
        let scaled = &result.streams[&line.scaled];
        let product = &result.streams[&line.product];
        let (chosen, rejected) = if selected == 1 {
            (scaled, other)
        } else {
            (other, scaled)
        };
        close(scaled.molar_flow_mol_s, factor * 2.0, 1e-12);
        if factor == 0.0 {
            assert!(scaled.state.is_none());
        } else {
            assert_eq!(scaled.state, fresh.state);
        }
        assert_eq!(product.molar_flow_mol_s, chosen.molar_flow_mol_s);
        assert_eq!(product.state, chosen.state);

        let added = result.units[&line.multiplier]
            .material_exchange
            .as_ref()
            .unwrap();
        let removed = result.units[&line.selector]
            .material_exchange
            .as_ref()
            .unwrap();
        close(added.molar_flow_mol_s, 2.0 * (factor - 1.0), 1e-12);
        close(removed.molar_flow_mol_s, -rejected.molar_flow_mol_s, 1e-12);
        close(
            added.enthalpy_flow_w,
            (factor - 1.0) * enthalpy_flow(fresh),
            1e-9,
        );
        close(removed.enthalpy_flow_w, -enthalpy_flow(rejected), 1e-9);
        for (i, inlet_component) in [1.5, 0.5].into_iter().enumerate() {
            close(
                added.component_flows_mol_s[i],
                (factor - 1.0) * inlet_component,
                1e-12,
            );
            close(
                removed.component_flows_mol_s[i],
                -component_flow(rejected, i),
                1e-12,
            );
            close(
                component_flow(product, i),
                inlet_component
                    + component_flow(other, i)
                    + added.component_flows_mol_s[i]
                    + removed.component_flows_mol_s[i],
                1e-12,
            );
        }
        close(
            enthalpy_flow(product),
            enthalpy_flow(fresh)
                + enthalpy_flow(other)
                + added.enthalpy_flow_w
                + removed.enthalpy_flow_w,
            1e-9,
        );
        for id in [line.multiplier, line.selector] {
            let u = &result.units[&id];
            assert_eq!(result.unit_currency(id, &fingerprints), Currency::Current);
            assert_eq!(u.duty_w, 0.0);
            assert!(u.max_scaled_component_residual.unwrap() < 1e-12);
            assert!(u.scaled_energy_residual.unwrap().abs() < 1e-12);
        }
    }
}

#[test]
fn selector_accounts_for_multiple_unselected_inlets_at_their_own_states() {
    let mut line = logical_line();
    line.project
        .flowsheet
        .set_spec(line.multiplier, multiplier(Some(1.0)));
    line.project
        .flowsheet
        .try_set_spec(line.selector, selector(3, Some(2)))
        .unwrap();
    let third = feed(&mut line.project, 4.0, 270.0, [0.25, 0.75]);
    let third = connect(&mut line.project, third, 0, line.selector, 3);

    let result = solve(&line.project).unwrap();
    let rejected = [&result.streams[&line.scaled], &result.streams[&third]];
    let exchange = result.units[&line.selector]
        .material_exchange
        .as_ref()
        .unwrap();
    assert_eq!(exchange.molar_flow_mol_s, -6.0);
    assert_eq!(exchange.component_flows_mol_s, [-2.5, -3.5]);
    close(
        exchange.enthalpy_flow_w,
        -rejected.into_iter().map(enthalpy_flow).sum::<f64>(),
        1e-9,
    );
    assert_eq!(
        result.streams[&line.product].state,
        result.streams[&line.comparison].state
    );
}

#[test]
fn original_v4_specs_and_invalid_logical_values_have_explicit_diagnostics() {
    let mut line = logical_line();
    let json = line.project.to_json().unwrap();
    assert!(!json.contains("\"factor\""));
    assert!(!json.contains("\"selected_inlet\""));
    line.project = Project::from_json(&json).unwrap();
    let error = solve(&line.project).unwrap_err().message;
    assert!(
        error.contains("flow factor") && error.contains("selected inlet"),
        "{error}"
    );

    line.project
        .flowsheet
        .set_spec(line.selector, selector(2, Some(2)));
    for factor in [f64::NAN, -1.0, f64::MAX] {
        line.project
            .flowsheet
            .set_spec(line.multiplier, multiplier(Some(factor)));
        let error = solve(&line.project).unwrap_err().message;
        let expected = if factor == f64::MAX {
            "nonfinite extensive"
        } else {
            "flow factor"
        };
        assert!(error.contains(expected), "{error}");
    }
    line.project
        .flowsheet
        .set_spec(line.multiplier, multiplier(Some(1.0)));
    line.project
        .flowsheet
        .set_spec(line.selector, selector(2, Some(3)));
    assert!(
        solve(&line.project)
            .unwrap_err()
            .message
            .contains("selected inlet")
    );
}

#[test]
fn final_recycle_closure_retains_the_multiplier_exchange() {
    let mut project = project();
    let feed = feed(&mut project, 2.0, 330.0, [1.0, 0.0]);
    let mix = project.flowsheet.add_unit(UnitKind::Mixer);
    let scale = project.flowsheet.add_unit(UnitKind::Multiplier);
    let split = project.flowsheet.add_unit(UnitKind::Splitter);
    let product = project.flowsheet.add_unit(UnitKind::Product);
    project.flowsheet.set_spec(scale, multiplier(Some(2.0)));
    project.flowsheet.set_spec(
        split,
        UnitSpec::Splitter {
            fraction_to_first: Some(0.25),
        },
    );
    let fresh = connect(&mut project, feed, 0, mix, 0);
    connect(&mut project, mix, 2, scale, 0);
    connect(&mut project, scale, 1, split, 0);
    connect(&mut project, split, 1, mix, 1);
    let product = connect(&mut project, split, 2, product, 0);

    let r = solve(&project).unwrap();
    assert!(r.recycle.as_ref().unwrap().max_scaled_residual <= 1e-10);
    let unit = &r.units[&scale];
    let exchange = unit.material_exchange.as_ref().unwrap();
    close(exchange.molar_flow_mol_s, 4.0, 1e-8);
    close(r.streams[&product].molar_flow_mol_s, 6.0, 1e-8);
    close(
        enthalpy_flow(&r.streams[&product]),
        enthalpy_flow(&r.streams[&fresh]) + exchange.enthalpy_flow_w,
        1e-5,
    );
    assert!(unit.max_scaled_component_residual.unwrap() <= 1e-8);
    assert!(unit.scaled_energy_residual.unwrap().abs() <= 2e-9);
}
