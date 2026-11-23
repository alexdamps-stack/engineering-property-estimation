//! Rotating group (A14, #182): pumps, compressors, expanders and the
//! multistage compressor with intercooling.
//!
//! Owned by the rotating work package (`docs/a17-architecture.md`). The
//! pump/compressor/expander specifications are the shared
//! [`UnitSpec::Pressurize`] and [`UnitSpec::Expand`] variants in `unit.rs`
//! (this group may extend them); their calculation is routed here through
//! `calculate_core`. A PS flash for the cubic package is a prerequisite
//! ([`CalcContext::ps`]).

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info};
use crate::UnitSpec;
use crate::unit::{CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, UnitKind};

const SINGLE: [PortSpec; 2] = [PortSpec::inlet("in"), PortSpec::outlet("out")];
const POWER_IN: [PortSpec; 1] = [PortSpec::energy_inlet("power")];
const POWER_OUT: [PortSpec; 1] = [PortSpec::energy_outlet("power")];
const MULTISTAGE_ENERGY: [PortSpec; 2] = [
    PortSpec::energy_inlet("power"),
    PortSpec::energy_outlet("intercooler duty"),
];

/// Metadata of a rotating kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::PressureChange as P;
    let core = |label, tag_series, energy_ports| KindInfo {
        issue: None,
        ..kind_info(label, tag_series, P, &SINGLE, energy_ports, 0)
    };
    match kind {
        UnitKind::Pump => core("Pump", "P", &POWER_IN),
        UnitKind::Compressor => core("Compressor", "K", &POWER_IN),
        UnitKind::Expander => core("Expander", "EX", &POWER_OUT),
        UnitKind::MultistageCompressor => kind_info(
            "Multistage compressor",
            "K",
            P,
            &SINGLE,
            &MULTISTAGE_ENERGY,
            182,
        ),
        _ => unreachable!("{kind:?} is not a rotating kind"),
    }
}

/// Calculation status of a rotating kind.
pub(crate) fn status(kind: UnitKind) -> CalculationStatus {
    match kind {
        UnitKind::Pump | UnitKind::Compressor | UnitKind::Expander => CalculationStatus::DrawnOnly(
            "requires a rotating-equipment efficiency model (deferred)",
        ),
        _ => DRAWN_ONLY,
    }
}

/// Specifications of the rotating kinds other than pump, compressor and
/// expander.
#[derive(Debug, Clone, PartialEq)]
pub enum RotatingSpec {
    /// Multistage compressor with intercooling (#182).
    MultistageCompressor(MultistageCompressorSpec),
}

/// Multistage compressor specification (#182). Fields follow with the
/// calculation.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MultistageCompressorSpec {}

/// Default specification of a rotating group kind.
pub(crate) fn default_spec(kind: UnitKind) -> RotatingSpec {
    match kind {
        UnitKind::MultistageCompressor => {
            RotatingSpec::MultistageCompressor(MultistageCompressorSpec::default())
        }
        _ => unreachable!("{kind:?} has no RotatingSpec"),
    }
}

impl GroupSpec for RotatingSpec {
    fn kind(&self) -> UnitKind {
        match self {
            RotatingSpec::MultistageCompressor(_) => UnitKind::MultistageCompressor,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        Vec::new()
    }
}

impl CalcContext<'_> {
    /// PS equilibrium at the specified absolute pressure and molar entropy.
    pub fn ps(
        &self,
        z: &[f64],
        p_pa: f64,
        s_j_mol_k: f64,
    ) -> Result<property_packages::EquilibriumState, String> {
        self.package()
            .ps(z, p_pa, s_j_mol_k)
            .map_err(|e| e.to_string())
    }
}

/// Pump, compressor and expander ([`UnitSpec::Pressurize`] /
/// [`UnitSpec::Expand`]).
pub(crate) fn calculate_core(
    kind: UnitKind,
    _spec: &UnitSpec,
    _inputs: &UnitInputs,
    _ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    Err(format!(
        "{kind} requires a rotating-equipment efficiency model (deferred)"
    ))
}

/// Calculates a [`RotatingSpec`] unit.
pub(crate) fn calculate(
    kind: UnitKind,
    _spec: &RotatingSpec,
    _inputs: &UnitInputs,
    _ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    Err(super::not_implemented(kind))
}

#[cfg(test)]
mod tests {
    use crate::UnitGroup;

    #[test]
    fn rotating_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::Rotating);
        assert!(matches!(
            crate::UnitKind::Pump.calculation_status(),
            crate::CalculationStatus::DrawnOnly(r) if r.contains("efficiency model")
        ));
    }
}
