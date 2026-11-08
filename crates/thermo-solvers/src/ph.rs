//! Pressure/enthalpy flash with bracketed temperature and pure-fluid lever rule.
use crate::{
    FlashOptions, FlashPhase, PR_CRITICAL_VOLUME_OVER_COVOLUME, PureSaturation, TpFlash, tp_flash,
};
use thermo_models::{CubicModel, RootSelection};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, Phase, Pressure, Temperature,
};

#[derive(Debug, Clone, Copy)]
pub struct PhFlashOptions {
    /// Explicit search interval in K, within the homogeneous 200..=600 K envelope.
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub max_iterations: usize,
    /// |h - h_target| / max(|h_target|, 1 J/mol). Must be positive and <= 1e-6.
    pub energy_tolerance: f64,
    pub tp: FlashOptions,
}

impl Default for PhFlashOptions {
    fn default() -> Self {
        Self {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 600.0,
            max_iterations: 100,
            energy_tolerance: 1e-10,
            tp: FlashOptions::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct PhFlash {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// J/mol of feed, on the model's caloric reference.
    pub target_enthalpy: f64,
    pub enthalpy: f64,
    /// (h - h_target) / max(|h_target|, 1 J/mol).
    pub scaled_energy_residual: f64,
    pub phases: Vec<FlashPhase>,
    /// Outer temperature iterations; pure saturation inversion counted separately.
    pub iterations: usize,
    /// TP evidence, absent only for a pure saturation lever-rule result.
    pub tp: Option<TpFlash>,
    pub pure_saturation: Option<PureSaturation>,
    pub diagnostics: Vec<Diagnostic>,
}

impl TpFlash {
    /// Phase-fraction-weighted molar enthalpy, J/mol of feed.
    pub fn enthalpy(&self) -> f64 {
        mixture_enthalpy(&self.phases)
    }
}

fn mixture_enthalpy(phases: &[FlashPhase]) -> f64 {
    phases
        .iter()
        .map(|ph| ph.phase_fraction * ph.state.enthalpy)
        .sum()
}

pub(crate) fn result(
    target: f64,
    phases: Vec<FlashPhase>,
    iterations: usize,
    tp: Option<TpFlash>,
    sat: Option<PureSaturation>,
) -> Result<PhFlash, Error> {
    let enthalpy = mixture_enthalpy(&phases);
    let scaled_energy_residual = (enthalpy - target) / target.abs().max(1.0);
    if !scaled_energy_residual.is_finite() {
        return Err(Error::NumericalFailure(
            "nonfinite PH energy residual".into(),
        ));
    }
    let state = &phases[0].state;
    let mut diagnostics = tp
        .as_ref()
        .map(|r| r.diagnostics.clone())
        .unwrap_or_default();
    diagnostics.push(Diagnostic {
        code: if sat.is_some() { "A5-PURE-SATURATION-LEVER-RULE" } else { "A5-BRACKETED-TP-ENERGY" }.into(),
        severity: DiagnosticSeverity::Information,
        message: format!("PH energy uses phase-weighted total EOS enthalpy on reference {}; scale max(|h_target|, 1 J/mol). Computational envelope only; no global stability proof or critical-region qualification.",state.caloric_reference.identifier),
    });
    Ok(PhFlash {
        temperature_k: state.temperature_k,
        pressure_pa: state.pressure_pa,
        target_enthalpy: target,
        enthalpy,
        scaled_energy_residual,
        phases,
        iterations,
        tp,
        pure_saturation: sat,
        diagnostics,
    })
}

/// PH equilibrium of a fixed feed composition. No temperature guess is required.
/// Errors from inner TP/stability solves propagate; they are never skipped.
/// The temperature bracket assumes a regular equilibrium enthalpy branch. Pure
/// coexistence is solved explicitly at saturation; near-critical coexistence
/// above 0.95 Tc is outside this implementation's supported domain.
pub fn ph_flash(
    model: &CubicModel,
    pressure: Pressure,
    target_enthalpy: f64,
    options: PhFlashOptions,
) -> Result<PhFlash, Error> {
    if !target_enthalpy.is_finite() {
        return Err(thermo_types::InputError::NonFinite("target enthalpy").into());
    }
    let mut lo = options.minimum_temperature_k;
    let mut hi = options.maximum_temperature_k;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo < 200.0
        || hi > 600.0
        || lo >= hi
        || !options.energy_tolerance.is_finite()
        || options.energy_tolerance <= 0.0
        || options.energy_tolerance > 1e-6
    {
        return Err(Error::UnsupportedCapability(
            "PH requires ordered bounds within 200..=600 K and energy tolerance in (0, 1e-6]"
                .into(),
        ));
    }
    if options.max_iterations == 0 {
        return Err(Error::NonConvergence(
            "PH temperature iteration limit is zero".into(),
        ));
    }
    // Validate pressure and the inner solver before the separate saturation path.
    let low = tp_flash(model, Temperature::new(lo)?, pressure, options.tp)?;
    let high = tp_flash(model, Temperature::new(hi)?, pressure, options.tp)?;
    if let Some(index) = model.fractions().iter().position(|x| *x == 1.0) {
        let top = hi.min(0.95 * model.components()[index].critical_temperature());
        if top >= lo
            && let Some(sat) =
                saturation_temperature(model, index, pressure, lo, top, options.max_iterations)?
        {
            let t = Temperature::new(sat.temperature_k)?;
            let liquid = model.evaluate(t, pressure, RootSelection::SmallestVolume)?;
            let vapor = model.evaluate(t, pressure, RootSelection::LargestVolume)?;
            let latent_heat = vapor.enthalpy - liquid.enthalpy;
            if latent_heat <= 0.0 || !latent_heat.is_finite() {
                return Err(Error::NumericalFailure(
                    "PH pure coexistence has nonpositive latent heat".into(),
                ));
            }
            let beta = (target_enthalpy - liquid.enthalpy) / latent_heat;
            // Include endpoints only within the requested energy tolerance.
            let budget = options.energy_tolerance * target_enthalpy.abs().max(1.0);
            if target_enthalpy >= liquid.enthalpy - budget
                && target_enthalpy <= vapor.enthalpy + budget
            {
                let beta = beta.clamp(0.0, 1.0);
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
                let r = result(target_enthalpy, phases, 0, None, Some(sat))?;
                if r.scaled_energy_residual.abs() > options.energy_tolerance {
                    return Err(Error::NonConvergence("PH pure lever rule cannot satisfy the requested energy tolerance at floating-point precision".into()));
                }
                return Ok(r);
            }
        }
    }
    if target_enthalpy < low.enthalpy() || target_enthalpy > high.enthalpy() {
        return Err(Error::OutsideValidityDomain(format!(
            "PH target {target_enthalpy} J/mol is outside endpoint enthalpies [{}, {}] J/mol at [{lo}, {hi}] K",
            low.enthalpy(),
            high.enthalpy()
        )));
    }
    for endpoint in [low, high] {
        if ((endpoint.enthalpy() - target_enthalpy) / target_enthalpy.abs().max(1.0)).abs()
            <= options.energy_tolerance
        {
            return result(
                target_enthalpy,
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
        let residual = (flash.enthalpy() - target_enthalpy) / target_enthalpy.abs().max(1.0);
        if residual.abs() <= options.energy_tolerance {
            return result(
                target_enthalpy,
                flash.phases.clone(),
                iteration,
                Some(flash),
                None,
            );
        }
        if mid == lo || mid == hi {
            return Err(Error::NonConvergence("PH temperature bracket collapsed without energy convergence; unresolved coexistence or ill-conditioned state".into()));
        }
        if residual < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(format!(
        "PH energy iteration exhausted {} iterations",
        options.max_iterations
    )))
}

/// Below 0.95 Tc, a single root lies on the liquid or vapor side of the
/// critical volume. With three roots use the actual outer-root Gibbs difference.
/// This brackets coexistence without evaluating saturation pressures outside the
/// EOS pressure envelope (e.g. water at the 200 K lower temperature bound).
fn saturation_temperature(
    model: &CubicModel,
    index: usize,
    pressure: Pressure,
    mut lo: f64,
    mut hi: f64,
    max_iterations: usize,
) -> Result<Option<PureSaturation>, Error> {
    let sample = |t: f64| -> Result<(f64, Option<PureSaturation>), Error> {
        let temperature = Temperature::new(t)?;
        let l = model.ln_fugacity(temperature, pressure, RootSelection::SmallestVolume)?;
        let v = model.ln_fugacity(temperature, pressure, RootSelection::LargestVolume)?;
        if l.root_count < 3 {
            let sign = if l.molar_volume / model.covolume() < PR_CRITICAL_VOLUME_OVER_COVOLUME {
                -1.0
            } else {
                1.0
            };
            return Ok((sign, None));
        }
        let difference = l.ln_fugacity_coefficients[index] - v.ln_fugacity_coefficients[index];
        Ok((
            difference,
            Some(PureSaturation {
                temperature_k: t,
                pressure_pa: pressure.value(),
                liquid_compressibility: l.compressibility,
                vapor_compressibility: v.compressibility,
                liquid_molar_volume: l.molar_volume,
                vapor_molar_volume: v.molar_volume,
                ln_fugacity_coefficient: v.ln_fugacity_coefficients[index],
                fugacity_residual: difference.abs(),
                iterations: 0,
            }),
        ))
    };
    let (left, left_sat) = sample(lo)?;
    let (right, right_sat) = sample(hi)?;
    for sat in [left_sat, right_sat].into_iter().flatten() {
        if sat.fugacity_residual <= 1e-12 {
            return Ok(Some(sat));
        }
    }
    if left > 0.0 || right < 0.0 {
        return Ok(None);
    }
    for iteration in 1..=max_iterations {
        let mid = lo + (hi - lo) / 2.0;
        let (difference, sat) = sample(mid)?;
        if let Some(mut sat) = sat
            && difference.abs() <= 1e-12
        {
            sat.iterations = iteration;
            return Ok(Some(sat));
        }
        if difference < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(
        "PH pure saturation temperature iteration exhausted".into(),
    ))
}
