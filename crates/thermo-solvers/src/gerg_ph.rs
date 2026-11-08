//! GERG PH equilibrium with native calorics and pure-fluid saturation lever rule.
use crate::equilibrium_model::EquilibriumModel;
use crate::gerg_flash::{GergEquilibrium, validate, validate_flash_options};
use crate::{FlashOptions, FlashPhase, PhFlash, PureSaturation, gerg_outer_roots, gerg_tp_flash};
use thermo_models::{GergModel, RootSelection};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, InputError, Phase, Pressure,
    Temperature,
};

#[derive(Debug, Clone, Copy)]
pub struct GergPhFlashOptions {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub max_iterations: usize,
    /// |h-h_target|/max(|h_target|,1 J/mol), positive and <=1e-6.
    pub energy_tolerance: f64,
    pub tp: FlashOptions,
}

impl Default for GergPhFlashOptions {
    fn default() -> Self {
        Self {
            minimum_temperature_k: 230.,
            maximum_temperature_k: 350.,
            max_iterations: 100,
            energy_tolerance: 1e-10,
            tp: FlashOptions::default(),
        }
    }
}

fn finish(
    target: f64,
    phases: Vec<FlashPhase>,
    iterations: usize,
    tp: Option<crate::TpFlash>,
    sat: Option<PureSaturation>,
    tolerance: f64,
) -> Result<PhFlash, Error> {
    let mut r = crate::ph::result(target, phases, iterations, tp, sat)?;
    if r.scaled_energy_residual.abs() > tolerance {
        return Err(Error::NonConvergence(
            "GERG PH final energy residual exceeds tolerance".into(),
        ));
    }
    r.diagnostics.push(Diagnostic {code:"A6-GERG-PH-NATIVE-CALORICS".into(),severity:DiagnosticSeverity::Information,
        message:"Native GERG ideal and residual enthalpies on the GERG reference; bracketed equilibrium energy solve or pure saturation lever rule. Inner root/stability failures propagate.".into()});
    Ok(r)
}

/// PH flash in the same bounded dry-gas fluid domain as GERG TP. Pure
/// coexistence below 0.95 Tc uses a saturation-temperature solve and lever rule.
/// Every evaluated inner state must be supported; no failed TP state is skipped.
pub fn gerg_ph_flash(
    model: &GergModel,
    p: Pressure,
    target: f64,
    options: GergPhFlashOptions,
) -> Result<PhFlash, Error> {
    if !target.is_finite() {
        return Err(InputError::NonFinite("GERG target enthalpy").into());
    }
    let (mut lo, mut hi) = (options.minimum_temperature_k, options.maximum_temperature_k);
    if !lo.is_finite()
        || !hi.is_finite()
        || lo < 230.
        || hi > 350.
        || lo >= hi
        || options.max_iterations == 0
        || !options.energy_tolerance.is_finite()
        || options.energy_tolerance <= 0.
        || options.energy_tolerance > 1e-6
    {
        return Err(Error::UnsupportedCapability("GERG PH requires ordered 230..=350 K bounds, positive iterations and energy tolerance <=1e-6".into()));
    }
    validate(model, Temperature::new(lo)?, p)?;
    validate_flash_options(options.tp)?;
    let adapter = GergEquilibrium::new(model)?;
    if let Some(tc) = model.pure_critical_temperature() {
        let top = hi.min(0.95 * tc);
        if top >= lo
            && let Some(sat) = pure_saturation(model, p, lo, top, options.max_iterations)?
        {
            let t = Temperature::new(sat.temperature_k)?;
            let liquid = adapter.evaluate(t, p, RootSelection::SmallestVolume)?;
            let vapor = adapter.evaluate(t, p, RootSelection::LargestVolume)?;
            let latent = vapor.enthalpy - liquid.enthalpy;
            if latent <= 0. || !latent.is_finite() {
                return Err(Error::NumericalFailure(
                    "GERG pure saturation latent heat is nonpositive".into(),
                ));
            }
            let budget = options.energy_tolerance * target.abs().max(1.);
            if target >= liquid.enthalpy - budget && target <= vapor.enthalpy + budget {
                let beta = ((target - liquid.enthalpy) / latent).clamp(0., 1.);
                return finish(
                    target,
                    vec![
                        FlashPhase {
                            phase: Phase::Vapor,
                            phase_fraction: beta,
                            composition: model.fractions().to_vec(),
                            state: vapor,
                        },
                        FlashPhase {
                            phase: Phase::Liquid,
                            phase_fraction: 1. - beta,
                            composition: model.fractions().to_vec(),
                            state: liquid,
                        },
                    ],
                    0,
                    None,
                    Some(sat),
                    options.energy_tolerance,
                );
            }
        }
    }
    let low = gerg_tp_flash(model, Temperature::new(lo)?, p, options.tp)?;
    let high = gerg_tp_flash(model, Temperature::new(hi)?, p, options.tp)?;
    let (mut hlo, mut hhi) = (low.enthalpy(), high.enthalpy());
    let scale = target.abs().max(1.);
    if target < hlo - options.energy_tolerance * scale
        || target > hhi + options.energy_tolerance * scale
    {
        return Err(Error::OutsideValidityDomain(
            "GERG PH target outside the endpoint equilibrium enthalpy interval".into(),
        ));
    }
    for endpoint in [low, high] {
        if (endpoint.enthalpy() - target).abs() / scale <= options.energy_tolerance {
            return finish(
                target,
                endpoint.phases.clone(),
                0,
                Some(endpoint),
                None,
                options.energy_tolerance,
            );
        }
    }
    for iteration in 1..=options.max_iterations {
        let secant = lo + (target - hlo) / (hhi - hlo) * (hi - lo);
        let margin = 0.1 * (hi - lo);
        let mid = if secant > lo + margin && secant < hi - margin {
            secant
        } else {
            lo + 0.5 * (hi - lo)
        };
        if mid == lo || mid == hi {
            return Err(Error::NonConvergence(
                "GERG PH bracket collapsed without energy convergence".into(),
            ));
        }
        let flash = gerg_tp_flash(model, Temperature::new(mid)?, p, options.tp)?;
        let h = flash.enthalpy();
        if (h - target).abs() / scale <= options.energy_tolerance {
            return finish(
                target,
                flash.phases.clone(),
                iteration,
                Some(flash),
                None,
                options.energy_tolerance,
            );
        }
        if h < target {
            lo = mid;
            hlo = h;
        } else {
            hi = mid;
            hhi = h;
        }
    }
    Err(Error::NonConvergence(
        "GERG PH energy iteration limit exhausted".into(),
    ))
}

fn pure_saturation(
    model: &GergModel,
    p: Pressure,
    mut lo: f64,
    mut hi: f64,
    limit: usize,
) -> Result<Option<PureSaturation>, Error> {
    let index = model.fractions().iter().position(|x| *x > 0.).unwrap();
    let sample = |t: f64| -> Result<(f64, Option<PureSaturation>), Error> {
        let roots = gerg_outer_roots(model, Temperature::new(t)?, p)?;
        if roots.len() == 1 {
            let sign = if roots[0].state.molar_density > model.reducing_density() {
                -1.
            } else {
                1.
            };
            return Ok((sign, None));
        }
        let vapor = &roots[0];
        let liquid = roots.last().unwrap();
        let difference =
            liquid.ln_fugacity_coefficients[index] - vapor.ln_fugacity_coefficients[index];
        Ok((
            difference,
            Some(PureSaturation {
                temperature_k: t,
                pressure_pa: p.value(),
                liquid_compressibility: liquid.state.compressibility,
                vapor_compressibility: vapor.state.compressibility,
                liquid_molar_volume: 1. / liquid.state.molar_density,
                vapor_molar_volume: 1. / vapor.state.molar_density,
                ln_fugacity_coefficient: vapor.ln_fugacity_coefficients[index],
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
    if left > 0. || right < 0. {
        return Ok(None);
    }
    for iteration in 1..=limit {
        let mid = lo + 0.5 * (hi - lo);
        let (difference, sat) = sample(mid)?;
        if let Some(mut sat) = sat
            && difference.abs() <= 1e-12
        {
            sat.iterations = iteration;
            return Ok(Some(sat));
        }
        if mid == lo || mid == hi {
            break;
        }
        if difference < 0. {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(
        "GERG pure saturation temperature iteration exhausted".into(),
    ))
}
