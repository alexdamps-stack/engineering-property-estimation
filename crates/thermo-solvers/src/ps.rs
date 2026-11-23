//! Pressure/entropy flash (A14) with the structure of the PH flash: a
//! bracketed outer temperature iteration around TP flash, and a separate
//! phase-fraction (lever-rule) path at the saturation temperature for a pure
//! component. Narrow-boiling mixtures are resolved by the same bisection, since
//! their equilibrium entropy is continuous (but steep) in temperature.
//!
//! Molar entropy is the phase-weighted total EOS entropy on the model's
//! caloric reference (ideal gas at 298.15 K and 1 bar, including the ideal
//! mixing term), the same reference as [`crate::TpFlash`] phase states.
use crate::{FlashPhase, PhFlashOptions, PureSaturation, TpFlash, tp_flash};
use thermo_models::{CUBIC_TEMPERATURE_RANGE_K, CubicModel, RootSelection};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, Phase, Pressure, Temperature,
};

/// Options of the PS flash. The PH options are reused: the same explicit
/// temperature interval, iteration limit and inner TP options apply, and
/// `energy_tolerance` bounds the scaled entropy residual
/// |s − s_target| / max(|s_target|, 1 J/(mol K)).
pub type PsFlashOptions = PhFlashOptions;

#[derive(Debug, Clone)]
pub struct PsFlash {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// J/(mol K) of feed, on the model's caloric reference.
    pub target_entropy: f64,
    pub entropy: f64,
    /// Phase-weighted molar enthalpy of the result, J/mol.
    pub enthalpy: f64,
    /// (s − s_target) / max(|s_target|, 1 J/(mol K)).
    pub scaled_entropy_residual: f64,
    pub phases: Vec<FlashPhase>,
    /// Outer temperature iterations; pure saturation inversion counted separately.
    pub iterations: usize,
    /// TP evidence, absent only for a pure saturation lever-rule result.
    pub tp: Option<TpFlash>,
    pub pure_saturation: Option<PureSaturation>,
    pub diagnostics: Vec<Diagnostic>,
}

impl TpFlash {
    /// Phase-fraction-weighted molar entropy, J/(mol K) of feed.
    pub fn entropy(&self) -> f64 {
        mixture_entropy(&self.phases)
    }
}

fn mixture_entropy(phases: &[FlashPhase]) -> f64 {
    phases
        .iter()
        .map(|ph| ph.phase_fraction * ph.state.entropy)
        .sum()
}

fn scaled(entropy: f64, target: f64) -> f64 {
    (entropy - target) / target.abs().max(1.0)
}

fn result(
    target: f64,
    phases: Vec<FlashPhase>,
    iterations: usize,
    tp: Option<TpFlash>,
    sat: Option<PureSaturation>,
) -> Result<PsFlash, Error> {
    let entropy = mixture_entropy(&phases);
    let enthalpy = phases
        .iter()
        .map(|ph| ph.phase_fraction * ph.state.enthalpy)
        .sum();
    let scaled_entropy_residual = scaled(entropy, target);
    if !scaled_entropy_residual.is_finite() {
        return Err(Error::NumericalFailure(
            "nonfinite PS entropy residual".into(),
        ));
    }
    let state = &phases[0].state;
    let mut diagnostics = tp
        .as_ref()
        .map(|r| r.diagnostics.clone())
        .unwrap_or_default();
    diagnostics.push(Diagnostic {
        code: if sat.is_some() { "A14-PURE-SATURATION-LEVER-RULE" } else { "A14-BRACKETED-TP-ENTROPY" }.into(),
        severity: DiagnosticSeverity::Information,
        message: format!("PS entropy uses phase-weighted total EOS entropy on reference {}; scale max(|s_target|, 1 J/(mol K)). Computational envelope only; no global stability proof or critical-region qualification.", state.caloric_reference.identifier),
    });
    Ok(PsFlash {
        temperature_k: state.temperature_k,
        pressure_pa: state.pressure_pa,
        target_entropy: target,
        entropy,
        enthalpy,
        scaled_entropy_residual,
        phases,
        iterations,
        tp,
        pure_saturation: sat,
        diagnostics,
    })
}

/// PS equilibrium of a fixed feed composition. No temperature guess is required.
/// Errors from inner TP/stability solves propagate; they are never skipped.
/// The bracket assumes a regular equilibrium entropy branch (entropy strictly
/// increasing with temperature at fixed pressure, cp > 0). Pure coexistence is
/// solved explicitly at saturation; near-critical coexistence above 0.95 Tc is
/// outside this implementation's supported domain, as for PH.
pub fn ps_flash(
    model: &CubicModel,
    pressure: Pressure,
    target_entropy: f64,
    options: PsFlashOptions,
) -> Result<PsFlash, Error> {
    if !target_entropy.is_finite() {
        return Err(thermo_types::InputError::NonFinite("target entropy").into());
    }
    // A single active component stored as 1 − ε takes the pure lever-rule path.
    if let Some(canonical) = crate::saturation::canonical_single_component(model)? {
        return ps_flash(&canonical, pressure, target_entropy, options);
    }
    let mut lo = options.minimum_temperature_k;
    let mut hi = options.maximum_temperature_k;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo < CUBIC_TEMPERATURE_RANGE_K[0]
        || hi > CUBIC_TEMPERATURE_RANGE_K[1]
        || lo >= hi
        || !options.energy_tolerance.is_finite()
        || options.energy_tolerance <= 0.0
        || options.energy_tolerance > 1e-6
    {
        let [t_lo, t_hi] = CUBIC_TEMPERATURE_RANGE_K;
        return Err(Error::UnsupportedCapability(format!(
            "PS requires ordered bounds within {t_lo}..={t_hi} K and entropy tolerance in (0, 1e-6]"
        )));
    }
    if options.max_iterations == 0 {
        return Err(Error::NonConvergence(
            "PS temperature iteration limit is zero".into(),
        ));
    }
    // Validate pressure and the inner solver before the separate saturation path.
    let low = tp_flash(model, Temperature::new(lo)?, pressure, options.tp)?;
    let high = tp_flash(model, Temperature::new(hi)?, pressure, options.tp)?;
    if let Some(index) = crate::saturation::single_active_component(model.fractions()) {
        let top = hi.min(0.95 * model.components()[index].critical_temperature());
        if top >= lo
            && let Some(sat) = crate::ph::pure_coexistence_temperature(
                model,
                index,
                pressure,
                lo,
                top,
                options.max_iterations,
            )?
        {
            let t = Temperature::new(sat.temperature_k)?;
            let liquid = model.evaluate(t, pressure, RootSelection::SmallestVolume)?;
            let vapor = model.evaluate(t, pressure, RootSelection::LargestVolume)?;
            let latent_entropy = vapor.entropy - liquid.entropy;
            if latent_entropy <= 0.0 || !latent_entropy.is_finite() {
                return Err(Error::NumericalFailure(
                    "PS pure coexistence has nonpositive entropy of vaporization".into(),
                ));
            }
            let budget = options.energy_tolerance * target_entropy.abs().max(1.0);
            if target_entropy >= liquid.entropy - budget
                && target_entropy <= vapor.entropy + budget
            {
                let beta = ((target_entropy - liquid.entropy) / latent_entropy).clamp(0.0, 1.0);
                let phases = vec![
                    FlashPhase {
                        phase: Phase::Vapor,
                        phase_fraction: beta,
                        composition: model.fractions().to_vec(),
                        state: vapor,
                    },
                    FlashPhase {
                        phase: Phase::Liquid,
                        phase_fraction: 1.0 - beta,
                        composition: model.fractions().to_vec(),
                        state: liquid,
                    },
                ];
                let r = result(target_entropy, phases, 0, None, Some(sat))?;
                if r.scaled_entropy_residual.abs() > options.energy_tolerance {
                    return Err(Error::NonConvergence("PS pure lever rule cannot satisfy the requested entropy tolerance at floating-point precision".into()));
                }
                return Ok(r);
            }
        }
    }
    if target_entropy < low.entropy() || target_entropy > high.entropy() {
        return Err(Error::OutsideValidityDomain(format!(
            "PS target {target_entropy} J/(mol K) is outside endpoint entropies [{}, {}] J/(mol K) at [{lo}, {hi}] K",
            low.entropy(),
            high.entropy()
        )));
    }
    for endpoint in [low, high] {
        if scaled(endpoint.entropy(), target_entropy).abs() <= options.energy_tolerance {
            return result(
                target_entropy,
                endpoint.phases.clone(),
                0,
                Some(endpoint),
                None,
            );
        }
    }
    for iteration in 1..=options.max_iterations {
        let mid = lo + (hi - lo) / 2.0;
        let flash = tp_flash(model, Temperature::new(mid)?, pressure, options.tp)?;
        let residual = scaled(flash.entropy(), target_entropy);
        if residual.abs() <= options.energy_tolerance {
            return result(
                target_entropy,
                flash.phases.clone(),
                iteration,
                Some(flash),
                None,
            );
        }
        if mid == lo || mid == hi {
            return Err(Error::NonConvergence("PS temperature bracket collapsed without entropy convergence; unresolved coexistence or ill-conditioned state".into()));
        }
        if residual < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(format!(
        "PS entropy iteration exhausted {} iterations",
        options.max_iterations
    )))
}
