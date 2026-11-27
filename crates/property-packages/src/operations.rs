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
    pub diagnostics: Vec<Diagnostic>,
}

impl SaturationState {
    /// A bounded A13 saturation point (cubic or GERG-2008) with its finite-trial
    /// stability evidence summarized as a diagnostic.
    pub(crate) fn from_point(point: thermo_solvers::SaturationPoint, model: &str) -> Self {
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
            diagnostics: vec![Diagnostic {
                code: "A13-BOUNDED-SATURATION".into(),
                severity: DiagnosticSeverity::Information,
                message: format!(
                    "{model} bubble/dew point from the bounded A13 solver: liquid and vapor roots, fugacity re-check 1e-9, both phases stable in finite TPD trials; no critical-point or retrograde continuation."
                ),
            }],
        }
    }
}
