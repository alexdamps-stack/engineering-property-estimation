//! Shortcut distillation column (#185) through the flowsheet solver, against
//! the independent reference `validation/columns/shortcut-v1.json`
//! (teqp/Cantera/SciPy, no Rust calls; see `validation/columns/README.md`).

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{ShortcutColumnSpec, ShortcutReflux};
use flowsheet_core::{
    Calculation, EnergySourceDuty, FeedCompositions, Flowsheet, PortRef, SolveError, StreamId,
    UnitId, UnitKind, UnitSpec,
};
use property_packages::CubicPackage;
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

const FIXTURE: &str = include_str!("../../../validation/columns/shortcut-v1.json");
// Port indices: feed, distillate, bottoms, condenser duty (out), reboiler duty (in).
const CONDENSER: usize = 3;
const REBOILER: usize = 4;

struct Column {
    flowsheet: Flowsheet,
    package: CubicPackage,
    feeds: FeedCompositions,
    column: UnitId,
    feed: StreamId,
    distillate: StreamId,
    bottoms: StreamId,
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

impl Column {
    fn new(
        ids: &[&str],
        z: &[f64],
        variant: CubicVariant,
        feed: (f64, f64, f64),
        spec: ShortcutColumnSpec,
    ) -> Self {
        let mut flowsheet = Flowsheet::new();
        let source = flowsheet.add_unit(UnitKind::Feed);
        let column = flowsheet.add_unit(UnitKind::ShortcutColumn);
        let top = flowsheet.add_unit(UnitKind::Product);
        let bottom = flowsheet.add_unit(UnitKind::Product);
        assert!(flowsheet.set_spec(
            source,
            UnitSpec::Feed {
                temperature_k: Some(feed.0),
                pressure_pa: Some(feed.1),
                molar_flow_mol_s: Some(feed.2),
            }
        ));
        flowsheet
            .try_set_spec(column, UnitSpec::Columns(ColumnSpec::Shortcut(spec)))
            .unwrap();
        let feed = flowsheet
            .connect(PortRef::new(source, 0), PortRef::new(column, 0))
            .unwrap();
        let distillate = flowsheet
            .connect(PortRef::new(column, 1), PortRef::new(top, 0))
            .unwrap();
        let bottoms = flowsheet
            .connect(PortRef::new(column, 2), PortRef::new(bottom, 0))
            .unwrap();
        let n = ids.len();
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition(ids, z),
            variant,
            vec![vec![0.0; n]; n],
            "explicit zero interactions, as in the shortcut reference",
        )
        .unwrap();
        let mut feeds = FeedCompositions::new();
        feeds.insert(source, composition(ids, z));
        Self {
            flowsheet,
            package: CubicPackage::new(model, Default::default()),
            feeds,
            column,
            feed,
            distillate,
            bottoms,
        }
    }

    fn solve(&self) -> Result<Calculation, SolveError> {
        self.flowsheet.solve(&self.package, &self.feeds)
    }
}

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(number).collect()
}

fn check(fixture: &Value, case: &str, quantity: &str, tolerance: &str, actual: f64, expected: f64) {
    let t = numbers(&fixture["tolerances"][tolerance]);
    let allowed = t[0] + t[1] * expected.abs();
    assert!(
        (actual - expected).abs() <= allowed,
        "{case} {quantity}: {actual} vs reference {expected} (allowed {allowed})"
    );
}

fn case_column(case: &Value) -> Column {
    let ids: Vec<&str> = case["ids"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    let variant = match case["variant"].as_str().unwrap() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        other => panic!("unknown variant {other}"),
    };
    let value = number(&case["reflux"]["value"]);
    let reflux = match case["reflux"]["kind"].as_str().unwrap() {
        "ratio" => ShortcutReflux::Ratio(Some(value)),
        "multiple_of_minimum" => ShortcutReflux::MinimumMultiple(Some(value)),
        other => panic!("unknown reflux {other}"),
    };
    let spec = ShortcutColumnSpec {
        light_key: Some(case["light_key"].as_str().unwrap().into()),
        heavy_key: Some(case["heavy_key"].as_str().unwrap().into()),
        light_key_recovery: Some(number(&case["r_LK"])),
        heavy_key_recovery: Some(number(&case["r_HK"])),
        pressure_pa: Some(number(&case["p"])),
        reflux,
    };
    Column::new(
        &ids,
        &numbers(&case["z"]),
        variant,
        (
            number(&case["feed_T"]),
            number(&case["feed_p"]),
            number(&case["F_mol_s"]),
        ),
        spec,
    )
}

#[test]
fn shortcut_column_matches_the_independent_reference() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let cases = fixture["cases"].as_array().unwrap();
    assert!(cases.len() >= 4);
    for case in cases {
        let name = case["id"].as_str().unwrap();
        let column = case_column(case);
        let calc = column.solve().unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let unit = &calc.units[&column.column];
        let r = unit.shortcut_column.as_ref().unwrap();
        let ids = case["ids"].as_array().unwrap();
        let lk = ids.iter().position(|v| *v == case["light_key"]).unwrap();
        let c = |quantity: &str, tolerance: &str, actual: f64, key: &str| {
            check(
                &fixture,
                name,
                quantity,
                tolerance,
                actual,
                number(&case[key]),
            )
        };
        c("q", "shortcut_scalar", r.feed_condition_q, "q");
        check(
            &fixture,
            name,
            "alpha_LK",
            "relative_volatility",
            r.relative_volatility,
            numbers(&case["alpha"])[lk],
        );
        c("N_min", "stages", r.minimum_stages, "N_min");
        c("theta", "shortcut_scalar", r.underwood_root, "theta");
        c("R_min", "shortcut_scalar", r.minimum_reflux_ratio, "R_min");
        c("R", "shortcut_scalar", r.reflux_ratio, "R");
        c("X", "shortcut_scalar", r.gilliland_x, "X");
        c("N", "stages", r.theoretical_stages, "N");
        c("N_R", "stages", r.rectifying_stages, "N_R");
        c("T_top", "temperature_K", r.top_vapor_temperature_k, "T_top");
        c("Q_C", "duty_W", r.condenser_duty_w, "Q_C");
        c("Q_R", "duty_W", r.reboiler_duty_w, "Q_R");
        let iterations = numbers(&fixture["tolerances"]["alpha_iterations"])[0];
        assert!(
            (r.volatility_iterations as f64 - number(&case["alpha_iterations"])).abs()
                <= iterations,
            "{name} iterations {}",
            r.volatility_iterations
        );

        let n = ids.len();
        for (stream, total, flows, x, t, h) in [
            (column.distillate, "D", "d", "x_D", "T_D", "h_D"),
            (column.bottoms, "B", "b", "x_B", "T_bot", "h_B"),
        ] {
            let s = &calc.streams[&stream];
            let state = s.state.as_ref().unwrap();
            c(total, "molar_flow_mol_s", s.molar_flow_mol_s, total);
            for (i, (actual, expected)) in s
                .component_flows_mol_s(n)
                .into_iter()
                .zip(numbers(&case[flows]))
                .enumerate()
            {
                check(
                    &fixture,
                    name,
                    &format!("{flows}[{i}]"),
                    "molar_flow_mol_s",
                    actual,
                    expected,
                );
            }
            for (i, (actual, expected)) in state
                .mole_fractions
                .iter()
                .zip(numbers(&case[x]))
                .enumerate()
            {
                check(
                    &fixture,
                    name,
                    &format!("{x}[{i}]"),
                    "composition",
                    *actual,
                    expected,
                );
            }
            c(t, "temperature_K", state.temperature_k, t);
            c(h, "saturated_enthalpy_J_per_mol", state.enthalpy(), h);
            assert_eq!(state.pressure_pa, number(&case["p"]));
        }
        let feed = calc.streams[&column.feed].state.as_ref().unwrap();
        c("h_F", "enthalpy_J_per_mol", feed.enthalpy(), "h_F");

        // The engine's own closure check passed; report it explicitly.
        let balance = unit.balance.as_ref().unwrap();
        assert!(balance.max_scaled_component_residual < 1e-10, "{name}");
        assert!(balance.scaled_energy_residual.abs() < 1e-9, "{name}");
        assert!((unit.duty_w - (r.reboiler_duty_w - r.condenser_duty_w)).abs() < 1e-6);
    }
}

#[test]
fn condenser_and_reboiler_ports_report_the_duties() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let case = &fixture["cases"][0];
    let mut column = case_column(case);
    let fs = &mut column.flowsheet;
    let sink = fs.add_unit(UnitKind::EnergySink);
    let source = fs.add_unit(UnitKind::EnergySource);
    assert!(fs.set_spec(
        source,
        UnitSpec::EnergySource {
            duty: EnergySourceDuty::Calculated
        }
    ));
    let condenser = fs
        .connect(
            PortRef::new(column.column, CONDENSER),
            PortRef::new(sink, 0),
        )
        .unwrap();
    let reboiler = fs
        .connect(
            PortRef::new(source, 0),
            PortRef::new(column.column, REBOILER),
        )
        .unwrap();
    assert!(fs.diagnostics().is_empty(), "{:?}", fs.diagnostics());
    let calc = column.solve().unwrap();
    let r = calc.units[&column.column].shortcut_column.unwrap();
    assert_eq!(calc.energy_streams[&condenser].duty_w, r.condenser_duty_w);
    assert_eq!(calc.energy_streams[&reboiler].duty_w, r.reboiler_duty_w);
    check(
        &fixture,
        "ports",
        "Q_C",
        "duty_W",
        r.condenser_duty_w,
        number(&case["Q_C"]),
    );

    // A specified reboiler duty would over-specify the column.
    assert!(column.flowsheet.set_spec(
        source,
        UnitSpec::EnergySource {
            duty: EnergySourceDuty::Specified(Some(1.0e6))
        }
    ));
    assert!(!column.flowsheet.diagnostics().is_empty());
}

#[test]
fn shortcut_failures_are_explicit() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let case = &fixture["cases"][0];
    let message = |mutate: &dyn Fn(&mut ShortcutColumnSpec)| {
        let column = case_column(case);
        let UnitSpec::Columns(ColumnSpec::Shortcut(mut spec)) =
            column.flowsheet.unit(column.column).unwrap().spec.clone()
        else {
            unreachable!()
        };
        mutate(&mut spec);
        let mut column = column;
        column
            .flowsheet
            .try_set_spec(column.column, UnitSpec::Columns(ColumnSpec::Shortcut(spec)))
            .unwrap();
        format!("{:?}", column.solve().unwrap_err())
    };
    for (mutate, expected) in [
        (
            &(|s: &mut ShortcutColumnSpec| s.reflux = ShortcutReflux::Ratio(Some(1.0)))
                as &dyn Fn(&mut ShortcutColumnSpec),
            "must exceed the minimum reflux",
        ),
        (
            &|s: &mut ShortcutColumnSpec| std::mem::swap(&mut s.light_key, &mut s.heavy_key),
            "not more volatile",
        ),
        (
            &|s: &mut ShortcutColumnSpec| s.pressure_pa = Some(3.0e6),
            "exceeds the feed pressure",
        ),
        (
            &|s: &mut ShortcutColumnSpec| s.light_key = Some("methane".into()),
            "not a package component",
        ),
        (
            &|s: &mut ShortcutColumnSpec| s.heavy_key = None,
            "heavy key",
        ),
    ] {
        let error = message(mutate);
        assert!(error.contains(expected), "{expected}: {error}");
    }
}
