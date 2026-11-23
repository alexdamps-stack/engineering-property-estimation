//! Heat-transport group: air cooler (#179), fired heater (#180),
//! multi-stream exchanger (#181), pipe segment (#183), relief valve (#184).
//!
//! Owned by the heat-transport work package (`docs/a17-architecture.md`).

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info, numbered};
use crate::unit::{
    CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, SpecProblem, ThermalSpec,
    UnitKind,
};

mod air_cooler;
mod fired_heater;

pub use air_cooler::{AirCoolerSpec, AirCoolerThermal};
pub use fired_heater::FiredHeaterSpec;

const SINGLE: [PortSpec; 2] = [PortSpec::inlet("in"), PortSpec::outlet("out")];

/// Metadata of a heat-transport kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::{HeatTransfer as H, PressureChange as P};
    const FAN: [PortSpec; 1] = [PortSpec::energy_inlet("fan power")];
    const DUTY: [PortSpec; 1] = [PortSpec::energy_inlet("duty")];
    const HEAT: [PortSpec; 1] = [PortSpec::energy_inlet("heat")];
    match kind {
        UnitKind::AirCooler => kind_info("Air cooler", "AC", H, &SINGLE, &FAN, 179),
        UnitKind::FiredHeater => kind_info("Fired heater", "H", H, &SINGLE, &DUTY, 180),
        // Passes are variable ports (see `MultiStreamExchangerSpec`).
        UnitKind::MultiStreamExchanger => {
            kind_info("Multi-stream exchanger", "LNG", H, &[], &[], 181)
        }
        UnitKind::PipeSegment => kind_info("Pipe segment", "PIPE", P, &SINGLE, &HEAT, 183),
        UnitKind::ReliefValve => kind_info("Relief valve", "PSV", P, &SINGLE, &[], 184),
        _ => unreachable!("{kind:?} is not a heat-transport kind"),
    }
}

/// Calculation status of a heat-transport kind.
pub(crate) fn status(kind: UnitKind) -> CalculationStatus {
    match kind {
        UnitKind::AirCooler | UnitKind::FiredHeater => CalculationStatus::Calculated,
        _ => DRAWN_ONLY,
    }
}

/// Specifications of the heat-transport kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum HeatTransportSpec {
    /// Air cooler (#179).
    AirCooler(AirCoolerSpec),
    /// Fired heater (#180).
    FiredHeater(FiredHeaterSpec),
    /// Multi-stream exchanger (#181).
    MultiStreamExchanger(MultiStreamExchangerSpec),
    /// Pipe segment (#183).
    PipeSegment(PipeSegmentSpec),
    /// Relief valve (#184).
    ReliefValve(ReliefValveSpec),
}

/// Largest number of passes of a multi-stream exchanger.
pub const MAX_PASSES: u16 = 12;

/// Multi-stream exchanger specification (#181). Each pass `k` has the
/// variable ports `in k` and `out k`.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiStreamExchangerSpec {
    /// Number of passes (streams), 2 to [`MAX_PASSES`].
    pub passes: u16,
}

impl Default for MultiStreamExchangerSpec {
    fn default() -> Self {
        Self { passes: 3 }
    }
}

/// Pipe segment specification (#183).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PipeSegmentSpec {}

/// Relief valve specification (#184).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReliefValveSpec {}

/// Default specification of a heat-transport kind.
pub(crate) fn default_spec(kind: UnitKind) -> HeatTransportSpec {
    match kind {
        UnitKind::AirCooler => HeatTransportSpec::AirCooler(Default::default()),
        UnitKind::FiredHeater => HeatTransportSpec::FiredHeater(Default::default()),
        UnitKind::MultiStreamExchanger => {
            HeatTransportSpec::MultiStreamExchanger(Default::default())
        }
        UnitKind::PipeSegment => HeatTransportSpec::PipeSegment(Default::default()),
        UnitKind::ReliefValve => HeatTransportSpec::ReliefValve(Default::default()),
        _ => unreachable!("{kind:?} has no HeatTransportSpec"),
    }
}

impl GroupSpec for HeatTransportSpec {
    fn kind(&self) -> UnitKind {
        match self {
            HeatTransportSpec::AirCooler(_) => UnitKind::AirCooler,
            HeatTransportSpec::FiredHeater(_) => UnitKind::FiredHeater,
            HeatTransportSpec::MultiStreamExchanger(_) => UnitKind::MultiStreamExchanger,
            HeatTransportSpec::PipeSegment(_) => UnitKind::PipeSegment,
            HeatTransportSpec::ReliefValve(_) => UnitKind::ReliefValve,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        match self {
            HeatTransportSpec::AirCooler(s) => s.issues(),
            HeatTransportSpec::FiredHeater(s) => s.issues(),
            HeatTransportSpec::MultiStreamExchanger(s) if !(2..=MAX_PASSES).contains(&s.passes) => {
                vec![SpecIssue {
                    field: "passes",
                    problem: SpecProblem::CountOutOfRange,
                }]
            }
            _ => Vec::new(),
        }
    }

    fn variable_ports(&self) -> Vec<PortSpec> {
        match self {
            HeatTransportSpec::MultiStreamExchanger(s) => numbered(PortSpec::inlet("in"), s.passes)
                .zip(numbered(PortSpec::outlet("out"), s.passes))
                .flat_map(|(a, b)| [a, b])
                .collect(),
            _ => Vec::new(),
        }
    }

    fn consumes_energy(&self, port: &PortSpec) -> bool {
        matches!(self, Self::FiredHeater(s) if s.thermal == ThermalSpec::EnergyStream)
            && port.same_port(&PortSpec::energy_inlet("duty"))
    }
}

/// Calculates a heat-transport unit.
pub(crate) fn calculate(
    kind: UnitKind,
    spec: &HeatTransportSpec,
    inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    match spec {
        HeatTransportSpec::AirCooler(s) => air_cooler::calculate(s, &inputs.material[0], ctx),
        HeatTransportSpec::FiredHeater(s) => {
            fired_heater::calculate(s, &inputs.material[0], inputs.energy, ctx)
        }
        _ => Err(super::not_implemented(kind)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::UnitGroup;

    #[test]
    fn heat_transport_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::HeatTransport);
    }

    #[test]
    fn exchanger_passes_are_numbered_port_pairs() {
        let spec = HeatTransportSpec::MultiStreamExchanger(MultiStreamExchangerSpec { passes: 2 });
        let labels: Vec<_> = spec.variable_ports().iter().map(|p| p.label()).collect();
        assert_eq!(labels, ["in 1", "out 1", "in 2", "out 2"]);
        let one = HeatTransportSpec::MultiStreamExchanger(MultiStreamExchangerSpec { passes: 1 });
        assert_eq!(one.issues()[0].problem, SpecProblem::CountOutOfRange);
    }
}
