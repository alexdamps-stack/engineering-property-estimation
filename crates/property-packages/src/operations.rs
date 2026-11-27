//! Result types of the package operations beyond TP/PH (#226): bubble and dew
//! points shared by every package that supports them.

use thermo_types::{Diagnostic, DiagnosticSeverity};

use crate::PhaseProperties;

/// Saturation boundary of a specified material phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaturationKind {
    /// The specified composition is a liquid at its bubble point; an incipient
    /// vapor of zero amount is in equilibrium with it.
    Bubble,
    /// The specified composition is a vapor at its dew point; an incipient
    /// liquid of zero amount is in equilibrium with it.
    Dew,
}

impl From<SaturationKind> for thermo_solvers::CubicSaturationKind {
    fn from(kind: SaturationKind) -> Self {
        match kind {
            SaturationKind::Bubble => Self::Bubble,
            SaturationKind::Dew => Self::Dew,
        }
    }
}

/// A bubble or dew point in package component order, with the properties of
/// both coexisting phases from the package's own phase models.
#[derive(Debug, Clone)]
pub struct SaturationState {
    pub kind: SaturationKind,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// Liquid mole fractions: the feed for a bubble point, incipient for a dew point.
    pub liquid_composition: Vec<f64>,
    /// Vapor mole fractions: incipient for a bubble point, the feed for a dew point.
    pub vapor_composition: Vec<f64>,
    pub liquid: PhaseProperties,
    pub vapor: PhaseProperties,
    /// max_i |ln(x_i φ_i^L) − ln(y_i φ_i^V)| over components present in both phases.
    pub max_ln_fugacity_residual: f64,
    pub iterations: usize,
    /// The heavier liquid of a three-phase (heterogeneous) bubble point (#227):
    /// the specified liquid has split, `liquid`/`liquid_composition` are then the
    /// lighter liquid and this is the heavier one. `None` for a two-phase boundary.
    /// Callers that assume one liquid must reject a state that carries it.
    pub second_liquid: Option<SaturationSecondLiquid>,
    pub diagnostics: Vec<Diagnostic>,
}

/// Heavier liquid of a three-phase bubble point.
#[derive(Debug, Clone)]
pub struct SaturationSecondLiquid {
    pub composition: Vec<f64>,
    pub properties: PhaseProperties,
    /// Moles of this liquid per mole of the specified (overall) liquid.
    pub material_fraction: f64,
}

impl SaturationState {
    /// A bounded A13 saturation point (cubic or GERG-2008) with its finite-trial
    /// stability evidence summarized as a diagnostic.
    pub(crate) fn from_point(point: thermo_solvers::SaturationPoint, model: &str) -> Self {
        let mut diagnostics = vec![Diagnostic {
            code: "A13-BOUNDED-SATURATION".into(),
            severity: DiagnosticSeverity::Information,
            message: format!(
                "{model} bubble/dew point from the bounded A13 solver: liquid and vapor roots, fugacity re-check 1e-9, phases stable in finite TPD trials; no critical-point or retrograde continuation."
            ),
        }];
        if !point.multiphase_stability.is_empty() {
            diagnostics.push(Diagnostic {
                code: "A13-MULTIPHASE-EVIDENCE".into(),
                severity: DiagnosticSeverity::Information,
                message: "Extended (non-dry) component set: the boundary phases also pass the multiphase (VLLE) aqueous, hydrocarbon and near-pure TPD trials.".into(),
            });
        }
        if point.second_liquid.is_some() {
            diagnostics.push(Diagnostic {
                code: "A13-THREE-PHASE-BUBBLE".into(),
                severity: DiagnosticSeverity::Warning,
                message: "Three-phase bubble point: the specified liquid splits into two liquids in equilibrium with the incipient vapor.".into(),
            });
        }
        Self {
            kind: match point.kind {
                thermo_solvers::CubicSaturationKind::Bubble => SaturationKind::Bubble,
                thermo_solvers::CubicSaturationKind::Dew => SaturationKind::Dew,
            },
            temperature_k: point.temperature_k,
            pressure_pa: point.pressure_pa,
            liquid_composition: point.liquid_composition,
            vapor_composition: point.vapor_composition,
            liquid: point.liquid.into(),
            vapor: point.vapor.into(),
            max_ln_fugacity_residual: point.max_ln_fugacity_residual,
            iterations: point.iterations,
            second_liquid: point.second_liquid.map(|l| SaturationSecondLiquid {
                composition: l.composition,
                properties: l.state.into(),
                material_fraction: l.material_fraction,
            }),
            diagnostics,
        }
    }
}
