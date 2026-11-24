//! Unit-operation groups of work package A17 and the calculation context
//! they share.
//!
//! Each group module (`rotating`, `heat_transport`, `simple_logical`,
//! `reactors`, `columns`, `multiphase`, `solids`) owns, for its kinds:
//!
//! - `info(kind) -> KindInfo`: label, tag series, palette category, fixed
//!   material and energy ports, sub-issue;
//! - `status(kind) -> CalculationStatus`: calculated or drawn only;
//! - its specification type (an enum with one variant per kind) with
//!   `default_spec(kind)` and a [`GroupSpec`] implementation (issues,
//!   variable ports, energy roles);
//! - `calculate(kind, spec, inputs, ctx)`, called by the engine.
//!
//! Groups extend [`CalcContext`] additively with their own `impl` blocks in
//! their own file (Rust allows inherent impls anywhere in the crate), so
//! property-access helpers do not collide. See `docs/a17-architecture.md`.

use property_packages::{
    CubicPackage, EquilibriumState, FlowsheetPackage, PhaseProperties, RootSelection,
    SaturationKind, SaturationState,
};
use thermo_types::Composition;

use crate::engine::{MaterialExchange, MaterialStream};
use crate::unit::{PortSpec, SpecIssue, UnitKind};

pub mod columns;
pub mod heat_transport;
pub mod multiphase;
pub mod reactors;
pub mod rotating;
pub mod simple_logical;
pub mod solids;

pub use columns::ColumnSpec;
pub use heat_transport::HeatTransportSpec;
pub use multiphase::MultiphaseSpec;
pub use rotating::RotatingSpec;
pub use simple_logical::SimpleLogicalSpec;
pub use solids::SolidsSpec;

/// Behaviour every group specification provides to the shared registry.
pub trait GroupSpec {
    /// The unit kind this specification belongs to.
    fn kind(&self) -> UnitKind;

    /// Missing and invalid fields, as [`crate::UnitSpec::issues`].
    fn issues(&self) -> Vec<SpecIssue>;

    /// Ports defined by the specification (numbered groups), appended after
    /// the kind's fixed ports. Each must be unique by
    /// [`PortSpec::same_port`]; connected ports keep their stream when the
    /// group grows or shrinks around them.
    fn variable_ports(&self) -> Vec<PortSpec> {
        Vec::new()
    }

    /// Whether the unit takes the duty of the energy stream at `port` from
    /// the stream (`true`) instead of determining it (`false`, default).
    fn consumes_energy(&self, _port: &PortSpec) -> bool {
        false
    }
}

/// Inlet values handed to a unit calculation.
#[derive(Debug, Clone, Copy)]
pub struct UnitInputs<'a> {
    /// One stream per material inlet port, in port order. An unconnected
    /// optional inlet is a zero-flow stream without a state.
    pub material: &'a [MaterialStream],
    /// One value per energy *inlet* port, in port order: the duty in W
    /// delivered into the unit, or `None` when the port is unconnected or
    /// the stream's source takes its value from this unit (an energy source
    /// set to "calculated").
    pub energy: &'a [Option<f64>],
}

/// Result of one unit calculation.
#[derive(Debug, Clone)]
pub struct UnitOutput {
    /// One stream per material outlet port, in port order.
    pub outlets: Vec<MaterialStream>,
    /// Net external duty in W, positive into the unit; it enters the energy
    /// balance. Energy streams at the unit's energy ports are part of it.
    pub duty_w: f64,
    /// Heat moved internally between sides of an exchanger, W.
    pub transferred_duty_w: Option<f64>,
    /// Design estimates for a shortcut distillation column.
    pub shortcut_column: Option<columns::ShortcutColumnResult>,
    /// Stage temperatures and duties of a rigorous distillation column.
    pub distillation_column: Option<columns::DistillationColumnResult>,
    /// Stage temperatures of an adiabatic absorber or stripper.
    pub absorber: Option<columns::AbsorberResult>,
    /// Stages, isentropic states and power of a rotating machine.
    pub rotating: Option<rotating::RotatingResult>,
    /// Material crossing the flowsheet boundary at this unit.
    pub material_exchange: Option<MaterialExchange>,
    /// Species generation and formation-enthalpy terms for reacting units.
    pub reaction: Option<crate::ReactionTerms>,
    /// Warnings shown with the result.
    pub warnings: Vec<String>,
    /// Duty in W at each energy port, in energy-port order, positive in the
    /// port's direction (into the unit for an inlet, out of it for an
    /// outlet). Leave empty for the default: a unit with exactly one energy
    /// port puts `duty_w` through it (negated for an outlet).
    pub energy: Vec<Option<f64>>,
}

impl UnitOutput {
    /// Output with outlets, duty and warnings; no transferred duty and the
    /// default energy-port values.
    pub fn new(outlets: Vec<MaterialStream>, duty_w: f64, warnings: Vec<String>) -> Self {
        Self {
            outlets,
            duty_w,
            transferred_duty_w: None,
            shortcut_column: None,
            distillation_column: None,
            absorber: None,
            rotating: None,
            material_exchange: None,
            reaction: None,
            warnings,
            energy: Vec::new(),
        }
    }
}

/// The explicit failure of a kind whose calculation is not implemented.
pub fn not_implemented(kind: UnitKind) -> String {
    match kind.issue() {
        Some(issue) => format!("{kind} calculation is not yet implemented (#{issue})"),
        None => format!("{kind} calculation is not yet implemented"),
    }
}

/// Builds [`KindInfo`](crate::KindInfo) in the group modules.
pub(crate) const fn kind_info(
    label: &'static str,
    tag_series: &'static str,
    category: crate::PaletteCategory,
    ports: &'static [PortSpec],
    energy_ports: &'static [PortSpec],
    issue: u16,
) -> crate::KindInfo {
    crate::KindInfo {
        label,
        tag_series,
        category,
        ports,
        energy_ports,
        issue: Some(issue),
    }
}

/// `count` numbered copies of `port` (1-based), for variable port groups.
pub fn numbered(port: PortSpec, count: u16) -> impl Iterator<Item = PortSpec> {
    (1..=count).map(move |n| port.numbered(n))
}

/// Status of a group kind that is drawn only.
pub(crate) const DRAWN_ONLY: crate::CalculationStatus =
    crate::CalculationStatus::DrawnOnly("is drawn only; its calculation is not implemented yet");

/// Property access for unit calculations: one immutable property package and,
/// for feeds, the feed composition. Methods return `Err(String)` with the
/// property error, so calculations can use `?` and fail explicitly.
///
/// Group modules may add methods in their own `impl CalcContext<'_>` block.
#[derive(Debug, Clone, Copy)]
pub struct CalcContext<'a> {
    package: &'a dyn FlowsheetPackage,
    feed: Option<&'a Composition>,
}

impl<'a> CalcContext<'a> {
    /// Context for one unit calculation.
    pub fn new(package: &'a dyn FlowsheetPackage, feed: Option<&'a Composition>) -> Self {
        Self { package, feed }
    }

    /// The property package.
    pub fn package(&self) -> &'a dyn FlowsheetPackage {
        self.package
    }

    /// Legacy access to the cubic package. No unit uses it: call the package
    /// operations below, which every package implements or rejects explicitly.
    pub fn cubic_package(&self) -> Result<&'a CubicPackage, String> {
        self.package.cubic().map_err(|e| e.to_string())
    }

    /// Number of components, in package order.
    pub fn component_count(&self) -> usize {
        self.package.component_ids().len()
    }

    /// Component identifiers in package order.
    pub fn component_ids(&self) -> Vec<&'a str> {
        self.package.component_ids()
    }

    /// Feed composition (feeds only).
    pub fn feed(&self) -> Option<&'a Composition> {
        self.feed
    }

    /// TP equilibrium of overall composition `z`.
    pub fn tp(&self, z: &[f64], t_k: f64, p_pa: f64) -> Result<EquilibriumState, String> {
        self.package.tp(z, t_k, p_pa).map_err(|e| e.to_string())
    }

    /// PH equilibrium of overall composition `z` (molar enthalpy in J/mol on
    /// the package caloric reference).
    pub fn ph(&self, z: &[f64], p_pa: f64, h_j_mol: f64) -> Result<EquilibriumState, String> {
        self.package.ph(z, p_pa, h_j_mol).map_err(|e| e.to_string())
    }

    /// PS equilibrium of overall composition `z` (molar entropy in J/(mol K)
    /// on the package caloric reference).
    pub fn ps(&self, z: &[f64], p_pa: f64, s_j_mol_k: f64) -> Result<EquilibriumState, String> {
        self.package
            .ps(z, p_pa, s_j_mol_k)
            .map_err(|e| e.to_string())
    }

    /// Homogeneous single-phase properties (enthalpy, entropy, cp, density,
    /// ln φ, …) of composition `x` at T and p on the selected root, from the
    /// package's own phase model. No stability claim.
    pub fn homogeneous(
        &self,
        x: &[f64],
        t_k: f64,
        p_pa: f64,
        root: RootSelection,
    ) -> Result<PhaseProperties, String> {
        self.package
            .homogeneous(x, t_k, p_pa, root)
            .map_err(|e| e.to_string())
    }

    /// Fugacity coefficients (as ln φ, package order) of composition `x` on
    /// the selected root, without phase-stability assessment.
    pub fn ln_fugacity_coefficients(
        &self,
        x: &[f64],
        t_k: f64,
        p_pa: f64,
        root: RootSelection,
    ) -> Result<Vec<f64>, String> {
        self.package
            .ln_fugacity_coefficients(x, t_k, p_pa, root)
            .map_err(|e| e.to_string())
    }

    /// Equilibrium ratios K_i = φ_i(liquid x) / φ_i(vapor y) at T and p,
    /// liquid on the smallest-volume root and vapor on the largest-volume
    /// root. A successive-substitution building block for column and
    /// multiphase solvers; it does not check phase stability.
    pub fn k_values(&self, x: &[f64], y: &[f64], t_k: f64, p_pa: f64) -> Result<Vec<f64>, String> {
        self.package
            .k_values(x, y, t_k, p_pa)
            .map_err(|e| e.to_string())
    }

    /// Bubble or dew temperature of `x` at `p_pa` from the package's
    /// saturation solver, inside its declared interval.
    pub fn saturation_temperature(
        &self,
        x: &[f64],
        p_pa: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, String> {
        self.package
            .saturation_temperature(x, p_pa, kind)
            .map_err(|e| e.to_string())
    }
}

/// Shared test helpers for the group modules: every group kind has
/// consistent metadata, a fitting default specification, and fails with an
/// explicit "not yet implemented" error while it is drawn only.
#[cfg(test)]
pub(crate) mod testing {
    use crate::{CalculationStatus, UnitGroup, UnitKind};

    /// Kinds of `group`.
    pub(crate) fn kinds(group: UnitGroup) -> Vec<UnitKind> {
        UnitKind::ALL
            .into_iter()
            .filter(|k| k.group() == group)
            .collect()
    }

    /// Checks the scaffolding contract for every drawn-only kind of `group`.
    pub(crate) fn check_drawn_only_contract(group: UnitGroup) {
        let kinds = kinds(group);
        assert!(!kinds.is_empty());
        for kind in kinds {
            let spec = kind.default_spec();
            assert!(spec.fits(kind), "{kind}");
            if let CalculationStatus::DrawnOnly(reason) = kind.calculation_status() {
                assert!(!reason.is_empty());
                assert!(
                    super::not_implemented(kind).contains("not yet implemented"),
                    "{kind}"
                );
            }
        }
    }
}
