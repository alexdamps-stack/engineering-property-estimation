//! Steady-state multi-stream exchanger (#181) via the flowsheet solver.

use flowsheet_core::units::heat_transport::{
    ExchangerPassDirection::{Cold, Hot},
    ExchangerPassSpec,
    ExchangerPassThermal::{Balance, Duty, OutletTemperature},
    HeatTransportSpec, MultiStreamExchangerSpec,
};
use flowsheet_core::{
    Calculation, DiagnosticKind, FeedCompositions, Flowsheet, PortDirection, PortMedium, PortRef,
    SolveError, SpecProblem, StreamId, UnitId, UnitKind, UnitSpec,
};
use property_packages::CubicPackage;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

#[derive(Clone, Copy)]
struct Feed {
    temperature_k: f64,
    pressure_pa: f64,
    flow_mol_s: f64,
    methane: f64,
}

fn composition(methane: f64) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        [("methane", methane), ("nitrogen", 1.0 - methane)]
            .into_iter()
            .map(|(id, x)| (ComponentId::new(id).unwrap(), x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

fn unit_spec(spec: MultiStreamExchangerSpec) -> UnitSpec {
    UnitSpec::HeatTransport(HeatTransportSpec::MultiStreamExchanger(spec))
}

struct ExchangerLine {
    flowsheet: Flowsheet,
    package: CubicPackage,
    compositions: FeedCompositions,
    exchanger: UnitId,
    streams: Vec<(StreamId, StreamId)>,
}

impl ExchangerLine {
    fn new(spec: MultiStreamExchangerSpec, feeds: &[Feed]) -> Self {
        assert_eq!(spec.passes.len(), feeds.len());
        let mut flowsheet = Flowsheet::new();
        let exchanger = flowsheet.add_unit(UnitKind::MultiStreamExchanger);
        flowsheet.try_set_spec(exchanger, unit_spec(spec)).unwrap();
        let ports = flowsheet.unit(exchanger).unwrap().ports();
        assert_eq!(ports.len(), 2 * feeds.len());
        assert!(ports.iter().all(|p| p.medium == PortMedium::Material));

        let mut compositions = FeedCompositions::new();
        let mut streams = Vec::new();
        for (index, feed) in feeds.iter().enumerate() {
            let port = |name, direction| {
                let label = format!("{name} {}", index + 1);
                let index = ports.iter().position(|p| p.label() == label).unwrap();
                assert_eq!(ports[index].direction, direction);
                PortRef::new(exchanger, index)
            };
            let inlet = flowsheet.add_unit(UnitKind::Feed);
            let outlet = flowsheet.add_unit(UnitKind::Product);
            assert!(flowsheet.set_spec(
                inlet,
                UnitSpec::Feed {
                    temperature_k: Some(feed.temperature_k),
                    pressure_pa: Some(feed.pressure_pa),
                    molar_flow_mol_s: Some(feed.flow_mol_s),
                }
            ));
            compositions.insert(inlet, composition(feed.methane));
            let input = flowsheet
                .connect(PortRef::new(inlet, 0), port("in", PortDirection::Inlet))
                .unwrap();
            let output = flowsheet
                .connect(port("out", PortDirection::Outlet), PortRef::new(outlet, 0))
                .unwrap();
            streams.push((input, output));
        }

        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition(0.5),
            CubicVariant::PengRobinson,
            vec![vec![0.0; 2]; 2],
            "synthetic zero binary interactions for the dilute-gas exchanger test",
        )
        .unwrap();
        Self {
            flowsheet,
            package: CubicPackage::new(model, Default::default()),
            compositions,
            exchanger,
            streams,
        }
    }

    fn solve(&self) -> Result<Calculation, SolveError> {
        self.flowsheet.solve(&self.package, &self.compositions)
    }

    fn failure(&self, expected: &str) -> String {
        let error = self.solve().unwrap_err();
        assert!(
            error.unit == Some(self.exchanger)
                || error.diagnostics.iter().any(|d| matches!(
                    d.kind,
                    DiagnosticKind::InvalidSpecification(unit, _)
                        | DiagnosticKind::MissingSpecification(unit, _)
                        if unit == self.exchanger
                )),
            "{error:?}"
        );
        let message = error.message.to_lowercase();
        assert!(message.contains(expected), "{error}");
        message
    }
}

const REFERENCE_FEEDS: [Feed; 4] = [
    Feed {
        temperature_k: 460.0,
        pressure_pa: 1.0,
        flow_mol_s: 6.0,
        methane: 0.8,
    },
    Feed {
        temperature_k: 285.0,
        pressure_pa: 1.5,
        flow_mol_s: 8.0,
        methane: 0.4,
    },
    Feed {
        temperature_k: 420.0,
        pressure_pa: 2.0,
        flow_mol_s: 5.0,
        methane: 0.2,
    },
    Feed {
        temperature_k: 320.0,
        pressure_pa: 2.5,
        flow_mol_s: 4.0,
        methane: 0.65,
    },
];

fn reference_spec() -> MultiStreamExchangerSpec {
    let passes = [
        (Hot, OutletTemperature(Some(390.0)), 0.125),
        (Cold, Balance, 0.25),
        (Hot, Duty(Some(-3_600.0)), 0.375),
        (Cold, OutletTemperature(Some(375.0)), 0.5),
    ]
    .into_iter()
    .map(|(direction, thermal, pressure_drop_pa)| ExchangerPassSpec {
        direction,
        thermal,
        pressure_drop_pa,
    })
    .collect();
    MultiStreamExchangerSpec {
        passes,
        // The limiting pair is pass 3 inlet (420 K) and pass 4 outlet (375 K).
        minimum_approach_k: 45.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64, quantity: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{quantity}: {actual} vs {expected} (absolute tolerance {tolerance})"
    );
}

#[test]
fn four_pass_exchanger_matches_independent_dilute_gas_reference() {
    // CH4/N2 mole mixtures, 0.875–2.5 Pa absolute, single-phase gas. The Rust
    // package uses Peng–Robinson with zero kij and ideal h=0 at 298.15 K.
    // Independently integrate the low-range NASA7 cp/R coefficients a1–a5
    // (pinned Cantera data in data/source-extracts/components.provisional.json,
    // version 2026-11-07.1). Each duty is flow * R * sum over species x_i times
    // sum over j=1..5 of a_j * (T_out^j - T_in^j) / j, with R=8.31446261815324.
    // Set pass 2 duty to the negative sum of the other duties and invert both
    // PH targets by bisection in [200, 600] K. Caloric offsets cancel. The
    // 0.015 W and 0.001 K budgets allow for finite-pressure PR departures.
    let expected_duties = [
        -16_566.470_697_948,
        12_510.418_852_392,
        -3_600.0,
        7_656.051_845_556,
    ];
    let expected_temperatures = [390.0, 333.958_448_182, 397.243_470_460, 375.0];
    let spec = reference_spec();
    let line = ExchangerLine::new(spec.clone(), &REFERENCE_FEEDS);
    let result = line.solve().unwrap();
    let mut duties = [0.0; 4];
    for (i, &(input, output)) in line.streams.iter().enumerate() {
        let input = &result.streams[&input];
        let output = &result.streams[&output];
        let input_state = input.state.as_ref().unwrap();
        let output_state = output.state.as_ref().unwrap();
        assert_eq!(output_state.phases.len(), 1);
        assert_eq!(output.molar_flow_mol_s, input.molar_flow_mol_s);
        assert_eq!(output_state.mole_fractions, input_state.mole_fractions);
        assert_eq!(
            output_state.pressure_pa,
            REFERENCE_FEEDS[i].pressure_pa - spec.passes[i].pressure_drop_pa
        );
        duties[i] = output.enthalpy_flow_w() - input.enthalpy_flow_w();
        let label = format!("pass {}", i + 1);
        close(duties[i], expected_duties[i], 0.015, &format!("{label} W"));
        close(
            output_state.temperature_k,
            expected_temperatures[i],
            0.001,
            &format!("{label} outlet K"),
        );
    }

    let unit = &result.units[&line.exchanger];
    let carried = unit.transferred_duty_w.unwrap();
    assert!(result.energy_streams.is_empty());
    assert_eq!(unit.duty_w, 0.0);
    close(
        carried,
        20_166.470_697_948,
        0.015,
        "reference transferred W",
    );
    close(carried, -(duties[0] + duties[2]), 1e-5, "hot-side W");
    close(carried, duties[1] + duties[3], 1e-5, "cold-side W");
    close(duties[2], -3_600.0, 1e-5, "specified pass 3 W");
    let balance = unit.balance.as_ref().unwrap();
    assert!(balance.max_scaled_component_residual <= 1e-12);
    close(balance.energy_residual_w, 0.0, 1e-5, "unit balance W");
}

#[test]
fn rejects_incomplete_and_overdetermined_thermal_specifications() {
    let mut spec = reference_spec();
    spec.passes[1].thermal = Duty(Some(12_500.0));
    ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("balanc");

    let mut spec = reference_spec();
    spec.passes[0].thermal = Balance;
    ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("balanc");

    for missing in [Duty(None), OutletTemperature(None)] {
        let mut spec = reference_spec();
        spec.passes[2].thermal = missing;
        let message = ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("pass 3");
        assert!(message.contains("missing"), "{message}");
    }
}

#[test]
fn rejects_duty_directions_for_each_thermal_mode() {
    for (index, thermal, failing_pass) in [
        (2, Duty(Some(3_600.0)), 3),
        (3, Duty(Some(-3_600.0)), 4),
        (0, OutletTemperature(Some(470.0)), 1),
        // Reduced hot-side removal leaves a negative duty for cold balance pass 2.
        (0, OutletTemperature(Some(450.0)), 2),
    ] {
        let mut spec = reference_spec();
        spec.passes[index].thermal = thermal;
        let message =
            ExchangerLine::new(spec, &REFERENCE_FEEDS).failure(&format!("pass {failing_pass}"));
        assert!(
            ["duty", "heat", "direction"]
                .into_iter()
                .any(|word| message.contains(word)),
            "{message}"
        );
    }
}

#[test]
fn checks_all_hot_and_cold_pairs_at_both_terminals() {
    let mut spec = reference_spec();
    spec.minimum_approach_k = 45.1;
    let message = ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("approach");
    assert!(
        message.contains("pass 3") && message.contains("pass 4"),
        "{message}"
    );

    // At zero specified approach, violate each terminal separately. The
    // other terminal remains feasible and the cold balancing duty stays positive.
    for (cold_in, cold_out) in [(320.0, 430.0), (395.0, 410.0)] {
        let mut spec = reference_spec();
        spec.minimum_approach_k = 0.0;
        spec.passes[3].thermal = OutletTemperature(Some(cold_out));
        let mut feeds = REFERENCE_FEEDS;
        feeds[3].temperature_k = cold_in;
        let message = ExchangerLine::new(spec, &feeds).failure("pass 4");
        assert!(
            message.contains("approach") || message.contains("cross"),
            "{message}"
        );
    }
}

#[test]
fn rejects_invalid_pressure_drops_and_nonfinite_specs() {
    for drop in [-0.125, f64::INFINITY, REFERENCE_FEEDS[2].pressure_pa] {
        let mut spec = reference_spec();
        spec.passes[2].pressure_drop_pa = drop;
        let message = ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("pass 3");
        assert!(message.contains("pressure"), "{message}");
    }
    for thermal in [Duty(Some(f64::NAN)), OutletTemperature(Some(f64::INFINITY))] {
        let mut spec = reference_spec();
        spec.passes[2].thermal = thermal;
        ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("pass 3");
    }
    for approach in [-1.0, f64::NAN] {
        let mut spec = reference_spec();
        spec.minimum_approach_k = approach;
        ExchangerLine::new(spec, &REFERENCE_FEEDS).failure("approach");
    }
}

fn zero_duty_spec() -> MultiStreamExchangerSpec {
    let mut spec = reference_spec();
    spec.passes.truncate(2);
    spec.passes[0].thermal = Duty(Some(0.0));
    spec.minimum_approach_k = 0.0;
    spec
}

#[test]
fn requires_two_to_twelve_passes_with_both_directions() {
    for count in [1, 2, 12, 13] {
        let mut spec = zero_duty_spec();
        if count == 1 {
            spec.passes.remove(0);
        } else {
            let extra = ExchangerPassSpec {
                direction: Hot,
                thermal: Duty(Some(0.0)),
                pressure_drop_pa: 0.0,
            };
            spec.passes.resize(count, extra);
        }
        let issues = unit_spec(spec).issues();
        if (2..=12).contains(&count) {
            assert!(issues.is_empty(), "{count} passes: {issues:?}");
        } else {
            assert!(
                issues
                    .iter()
                    .any(|i| i.problem == SpecProblem::CountOutOfRange),
                "{count} passes: {issues:?}"
            );
        }
    }
    for direction in [Hot, Cold] {
        let mut spec = zero_duty_spec();
        for pass in &mut spec.passes {
            pass.direction = direction;
        }
        ExchangerLine::new(spec, &REFERENCE_FEEDS[..2]).failure("pass");
    }
}

#[test]
fn explicit_zero_duty_bypasses_empty_pass_but_checks_flowing_zero_duty_pass() {
    let spec = zero_duty_spec();
    let mut feeds = [REFERENCE_FEEDS[0], REFERENCE_FEEDS[1]];
    feeds[0].temperature_k = 260.0;
    feeds[0].flow_mol_s = 0.0;
    let line = ExchangerLine::new(spec.clone(), &feeds);
    let result = line.solve().unwrap();
    let empty = &result.streams[&line.streams[0].1];
    assert_eq!(empty.molar_flow_mol_s, 0.0);
    assert!(empty.state.is_none());
    let unit = &result.units[&line.exchanger];
    assert_eq!(unit.transferred_duty_w, Some(0.0));
    assert!(
        unit.warnings
            .iter()
            .any(|w| w.to_lowercase().contains("pass 1")),
        "{:?}",
        unit.warnings
    );

    // Restore flow at 260 K, below the cold pass inlet of 285 K. Both duties stay zero.
    feeds[0].flow_mol_s = 6.0;
    ExchangerLine::new(spec, &feeds).failure("approach");
}

#[test]
fn rejects_empty_passes_that_require_a_thermal_calculation() {
    let mut feeds = REFERENCE_FEEDS;
    feeds[0].flow_mol_s = 0.0;
    for thermal in [Duty(Some(-1.0)), OutletTemperature(Some(390.0)), Balance] {
        let mut spec = reference_spec();
        spec.passes[0].thermal = thermal;
        if thermal == Balance {
            spec.passes[1].thermal = Duty(Some(1_000.0));
        }
        let message = ExchangerLine::new(spec, &feeds).failure("pass 1");
        assert!(
            message.contains("empty") || (message.contains("zero") && message.contains("flow")),
            "{message}"
        );
    }
}
