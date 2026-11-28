//! Multiphase group: VLLE flash with three-phase separator (#176), decanter
//! (#177), three-phase distillation (#188) and liquid–liquid extractor
//! (#189).
//!
//! Owned by the multiphase work package (`docs/a17-architecture.md`). The
//! VLLE / LLE flash is a prerequisite this group adds (the package exposes
//! TP/PH two-phase flashes and fugacity coefficients through
//! [`CalcContext`]).

use property_packages::MultiphaseEquilibrium;

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info};
use crate::unit::{CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, UnitKind};

// Liquid–liquid extractor (#189).
pub mod extractor;
pub use extractor::{
    ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorResult, LiquidExtractorSpec,
    MAX_EXTRACTOR_STAGES,
};

// Three-phase distillation column (#188).
pub mod three_phase_column;
pub use three_phase_column::{
    DecanterReflux, PhaseReflux, StageLiquidSplit, ThreePhaseColumnResult, ThreePhaseColumnSpec,
};

/// Metadata of a multiphase kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::{Columns as C, Separators as S};
    const THREE_PHASE: [PortSpec; 4] = [
        PortSpec::inlet("in"),
        PortSpec::outlet("vapor"),
        PortSpec::outlet("light liquid"),
        PortSpec::outlet("heavy liquid"),
    ];
    const DECANTER: [PortSpec; 3] = [
        PortSpec::inlet("in"),
        PortSpec::outlet("light liquid"),
        PortSpec::outlet("heavy liquid"),
    ];
    const COLUMN: [PortSpec; 4] = [
        PortSpec::inlet("feed"),
        PortSpec::outlet("distillate"),
        PortSpec::outlet("decant"),
        PortSpec::outlet("bottoms"),
    ];
    const EXTRACTOR: [PortSpec; 4] = [
        PortSpec::inlet("feed"),
        PortSpec::inlet("solvent"),
        PortSpec::outlet("extract"),
        PortSpec::outlet("raffinate"),
    ];
    const DUTY: [PortSpec; 1] = [PortSpec::energy_inlet("duty")];
    const CONDENSER_REBOILER: [PortSpec; 2] = [
        PortSpec::energy_outlet("condenser duty"),
        PortSpec::energy_inlet("reboiler duty"),
    ];
    match kind {
        UnitKind::ThreePhaseSeparator => {
            kind_info("Three-phase separator", "D", S, &THREE_PHASE, &DUTY, 176)
        }
        UnitKind::Decanter => kind_info("Decanter", "D", S, &DECANTER, &DUTY, 177),
        UnitKind::ThreePhaseColumn => kind_info(
            "Three-phase column",
            "T",
            C,
            &COLUMN,
            &CONDENSER_REBOILER,
            188,
        ),
        UnitKind::LiquidExtractor => {
            kind_info("Liquid-liquid extractor", "T", C, &EXTRACTOR, &[], 189)
        }
        _ => unreachable!("{kind:?} is not a multiphase kind"),
    }
}

/// Calculation status of a multiphase kind.
pub(crate) fn status(kind: UnitKind) -> CalculationStatus {
    match kind {
        UnitKind::LiquidExtractor | UnitKind::ThreePhaseColumn => CalculationStatus::Calculated,
        _ => DRAWN_ONLY,
    }
}

/// Specifications of the multiphase kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum MultiphaseSpec {
    /// Three-phase separator (#176).
    ThreePhaseSeparator(ThreePhaseSeparatorSpec),
    /// Decanter (#177).
    Decanter(DecanterSpec),
    /// Three-phase distillation (#188).
    ThreePhaseColumn(ThreePhaseColumnSpec),
    /// Liquid–liquid extractor (#189).
    LiquidExtractor(LiquidExtractorSpec),
}

/// Three-phase separator specification (#176).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThreePhaseSeparatorSpec {}

/// Decanter specification (#177).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DecanterSpec {}

/// Default specification of a multiphase kind.
pub(crate) fn default_spec(kind: UnitKind) -> MultiphaseSpec {
    match kind {
        UnitKind::ThreePhaseSeparator => MultiphaseSpec::ThreePhaseSeparator(Default::default()),
        UnitKind::Decanter => MultiphaseSpec::Decanter(Default::default()),
        UnitKind::ThreePhaseColumn => MultiphaseSpec::ThreePhaseColumn(Default::default()),
        UnitKind::LiquidExtractor => MultiphaseSpec::LiquidExtractor(Default::default()),
        _ => unreachable!("{kind:?} has no MultiphaseSpec"),
    }
}

impl GroupSpec for MultiphaseSpec {
    fn kind(&self) -> UnitKind {
        match self {
            MultiphaseSpec::ThreePhaseSeparator(_) => UnitKind::ThreePhaseSeparator,
            MultiphaseSpec::Decanter(_) => UnitKind::Decanter,
            MultiphaseSpec::ThreePhaseColumn(_) => UnitKind::ThreePhaseColumn,
            MultiphaseSpec::LiquidExtractor(_) => UnitKind::LiquidExtractor,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        match self {
            MultiphaseSpec::LiquidExtractor(s) => s.issues(),
            MultiphaseSpec::ThreePhaseColumn(s) => s.issues(),
            _ => Vec::new(),
        }
    }
}

impl CalcContext<'_> {
    /// Three-phase (VLLE/LLE) TP equilibrium of overall composition `z` with
    /// multiphase stability: phases vapor first, then liquids by increasing mass
    /// density, each with its [`PhaseRole`](property_packages::PhaseRole).
    pub fn vlle(&self, z: &[f64], t_k: f64, p_pa: f64) -> Result<MultiphaseEquilibrium, String> {
        self.package()
            .vlle_tp(z, t_k, p_pa)
            .map_err(|e| e.to_string())
    }

    /// Three-phase (VLLE/LLE) PH equilibrium (molar enthalpy in J/mol on the
    /// package caloric reference).
    pub fn vlle_ph(
        &self,
        z: &[f64],
        p_pa: f64,
        h_j_mol: f64,
    ) -> Result<MultiphaseEquilibrium, String> {
        self.package()
            .vlle_ph(z, p_pa, h_j_mol)
            .map_err(|e| e.to_string())
    }
}

/// Calculates a multiphase unit.
pub(crate) fn calculate(
    kind: UnitKind,
    spec: &MultiphaseSpec,
    inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    match spec {
        MultiphaseSpec::LiquidExtractor(s) => {
            let (mut output, result) = extractor::calculate(s, inputs.material, ctx)?;
            output.liquid_extractor = Some(result);
            Ok(output)
        }
        MultiphaseSpec::ThreePhaseColumn(s) => {
            let (mut output, result) = three_phase_column::calculate(s, inputs.material, ctx)?;
            output.three_phase_column = Some(result);
            Ok(output)
        }
        _ => Err(super::not_implemented(kind)),
    }
}

#[cfg(test)]
mod tests {
    use crate::UnitGroup;

    #[test]
    fn multiphase_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::Multiphase);
    }

    fn water_hydrocarbon_package() -> property_packages::CubicPackage {
        use thermo_data::ComponentDatabase;
        use thermo_models::{CubicModel, CubicVariant};
        use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};
        let ids = ["methane", "propane", "water"];
        let composition = Composition::new(
            CompositionBasis::Mole,
            ids.iter()
                .zip([0.1, 0.4, 0.5])
                .map(|(id, x)| (ComponentId::new(*id).unwrap(), x))
                .collect(),
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let kij = vec![
            vec![0.0, 0.0, 0.5],
            vec![0.0, 0.0, 0.5],
            vec![0.5, 0.5, 0.0],
        ];
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition,
            CubicVariant::PengRobinson,
            kij,
            "explicit water-hydrocarbon kij 0.5 for the VLLE context test",
        )
        .unwrap();
        property_packages::CubicPackage::new(model, Default::default())
    }

    #[test]
    fn context_reaches_the_package_three_phase_flash() {
        use property_packages::PhaseRole;
        let package = water_hydrocarbon_package();
        let ctx = super::CalcContext::new(&package, None);
        let z = [0.1, 0.4, 0.5];
        let tp = ctx.vlle(&z, 300.0, 20e5).unwrap();
        assert_eq!(
            tp.roles,
            [
                PhaseRole::Vapor,
                PhaseRole::LightLiquid,
                PhaseRole::HeavyLiquid
            ]
        );
        assert!(tp.state.convergence.max_ln_fugacity_residual <= 1e-10);
        // The heavy liquid is the aqueous phase.
        assert!(
            tp.phase_with_role(PhaseRole::HeavyLiquid)
                .unwrap()
                .composition[2]
                > 0.999
        );
        let ph = ctx.vlle_ph(&z, 20e5, tp.state.enthalpy()).unwrap();
        assert!((ph.state.temperature_k - 300.0).abs() < 1e-6);
        assert_eq!(ph.roles, tp.roles);
        // Component and energy balances close.
        for (i, zi) in z.iter().enumerate() {
            let total: f64 = ph
                .state
                .phases
                .iter()
                .map(|p| p.phase_fraction * p.composition[i])
                .sum();
            assert!((total - zi).abs() < 1e-8);
        }
        assert!(
            (ph.state.enthalpy() - tp.state.enthalpy()).abs()
                <= 2e-9 * tp.state.enthalpy().abs().max(1.0)
        );
    }

    #[test]
    fn context_reports_an_unsupported_three_phase_flash_explicitly() {
        use thermo_models::GergModel;
        use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};
        let composition = Composition::new(
            CompositionBasis::Mole,
            vec![
                (ComponentId::new("methane").unwrap(), 0.9),
                (ComponentId::new("ethane").unwrap(), 0.1),
            ],
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let package = property_packages::GergPackage::new(
            GergModel::new(&composition).unwrap(),
            [250.0, 300.0],
        )
        .unwrap();
        let ctx = super::CalcContext::new(&package, None);
        let error = ctx.vlle(&[0.9, 0.1], 280.0, 50e5).unwrap_err();
        assert!(
            error.contains("unsupported") && error.contains("VLLE"),
            "{error}"
        );
        assert!(
            ctx.vlle_ph(&[0.9, 0.1], 50e5, 0.0)
                .unwrap_err()
                .contains("VLLE")
        );
    }
}
