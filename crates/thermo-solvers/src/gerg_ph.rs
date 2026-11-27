//! GERG PH and PS equilibrium with native calorics and pure-fluid saturation
//! lever rule. Both share one bracketed structure.
use crate::equilibrium_model::EquilibriumModel;
use crate::gerg_flash::{GergEquilibrium, validate, validate_flash_options};
use crate::{
    FlashOptions, FlashPhase, PhFlash, PsFlash, PureSaturation, TpFlash, gerg_outer_roots,
    gerg_tp_flash,
};
use thermo_models::{GergModel, HomogeneousState, RootSelection};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, InputError, Phase, Pressure,
    Temperature,
};

#[derive(Debug, Clone, Copy)]
pub struct GergPhFlashOptions {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub max_iterations: usize,
    /// |h-h_target|/max(|h_target|,1 J/mol), positive and <=1e-6. The PS flash
    /// applies it to |s-s_target|/max(|s_target|,1 J/(mol K)).
    pub energy_tolerance: f64,
    pub tp: FlashOptions,
}

/// Options of the GERG PS flash: the same bracket, limits and tolerances as PH.
pub type GergPsFlashOptions = GergPhFlashOptions;

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

/// The caloric quantity specified with pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Spec {
    Enthalpy,
    Entropy,
}

impl Spec {
    fn of(self, s: &HomogeneousState) -> f64 {
        match self {
            Self::Enthalpy => s.enthalpy,
            Self::Entropy => s.entropy,
        }
    }
    fn mixture(self, flash: &TpFlash) -> f64 {
        match self {
            Self::Enthalpy => flash.enthalpy(),
            Self::Entropy => flash.entropy(),
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::Enthalpy => "PH",
            Self::Entropy => "PS",
        }
    }
}

/// Phases and evidence of a converged GERG PH or PS solve.
struct Solved {
    phases: Vec<FlashPhase>,
    iterations: usize,
    tp: Option<TpFlash>,
    sat: Option<PureSaturation>,
}

const NATIVE_CALORICS: &str = "Native GERG ideal and residual calorics on the GERG reference; bracketed equilibrium solve or pure saturation lever rule. Inner root/stability failures propagate.";

/// PH flash in the same bounded dry-gas fluid domain as GERG TP. Pure
/// coexistence below 0.95 Tc uses a saturation-temperature solve and lever rule.
/// Every evaluated inner state must be supported; no failed TP state is skipped.
pub fn gerg_ph_flash(
    model: &GergModel,
    p: Pressure,
    target: f64,
    options: GergPhFlashOptions,
) -> Result<PhFlash, Error> {
    let s = solve(model, p, target, options, Spec::Enthalpy)?;
    let mut r = crate::ph::result(target, s.phases, s.iterations, s.tp, s.sat)?;
    if r.scaled_energy_residual.abs() > options.energy_tolerance {
        return Err(Error::NonConvergence(
            "GERG PH final energy residual exceeds tolerance".into(),
        ));
    }
    r.diagnostics.push(Diagnostic {
        code: "A6-GERG-PH-NATIVE-CALORICS".into(),
        severity: DiagnosticSeverity::Information,
        message: NATIVE_CALORICS.into(),
    });
    Ok(r)
}

/// PS flash with the structure of [`gerg_ph_flash`]: the same domain, pure
/// saturation lever rule (in entropy) and bracketed equilibrium solve, with
/// the total native GERG entropy (including ideal mixing) as the specification.
pub fn gerg_ps_flash(
    model: &GergModel,
    p: Pressure,
    target: f64,
    options: GergPsFlashOptions,
) -> Result<PsFlash, Error> {
    let s = solve(model, p, target, options, Spec::Entropy)?;
    let mut r = crate::ps::result(target, s.phases, s.iterations, s.tp, s.sat)?;
    if r.scaled_entropy_residual.abs() > options.energy_tolerance {
        return Err(Error::NonConvergence(
            "GERG PS final entropy residual exceeds tolerance".into(),
        ));
    }
    r.diagnostics.push(Diagnostic {
        code: "A6-GERG-PS-NATIVE-CALORICS".into(),
        severity: DiagnosticSeverity::Information,
        message: NATIVE_CALORICS.into(),
    });
    Ok(r)
}

fn solve(
    model: &GergModel,
    p: Pressure,
    target: f64,
    options: GergPhFlashOptions,
    spec: Spec,
) -> Result<Solved, Error> {
    let name = spec.name();
    if !target.is_finite() {
        return Err(InputError::NonFinite(match spec {
            Spec::Enthalpy => "GERG target enthalpy",
            Spec::Entropy => "GERG target entropy",
        })
        .into());
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
        return Err(Error::UnsupportedCapability(format!(
            "GERG {name} requires ordered 230..=350 K bounds, positive iterations and a tolerance <=1e-6"
        )));
    }
    validate(model, Temperature::new(lo)?, p)?;
    validate_flash_options(options.tp)?;
    let adapter = GergEquilibrium::new(model)?;
    if let Some(tc) = model.pure_critical_temperature() {
        let top = hi.min(0.95 * tc);
        if top >= lo
            && let Some(sat) =
                gerg_pure_saturation_temperature(model, p, lo, top, options.max_iterations)?
        {
            let t = Temperature::new(sat.temperature_k)?;
            let liquid = adapter.evaluate(t, p, RootSelection::SmallestVolume)?;
            let vapor = adapter.evaluate(t, p, RootSelection::LargestVolume)?;
            let latent = spec.of(&vapor) - spec.of(&liquid);
            if latent <= 0. || !latent.is_finite() {
                return Err(Error::NumericalFailure(format!(
                    "GERG pure saturation {} is nonpositive",
                    match spec {
                        Spec::Enthalpy => "latent heat",
                        Spec::Entropy => "entropy of vaporization",
                    }
                )));
            }
            let budget = options.energy_tolerance * target.abs().max(1.);
            if target >= spec.of(&liquid) - budget && target <= spec.of(&vapor) + budget {
                let beta = ((target - spec.of(&liquid)) / latent).clamp(0., 1.);
                return Ok(Solved {
                    phases: vec![
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
                    iterations: 0,
                    tp: None,
                    sat: Some(sat),
                });
            }
        }
    }
    let low = gerg_tp_flash(model, Temperature::new(lo)?, p, options.tp)?;
    let high = gerg_tp_flash(model, Temperature::new(hi)?, p, options.tp)?;
    let (mut vlo, mut vhi) = (spec.mixture(&low), spec.mixture(&high));
    let scale = target.abs().max(1.);
    if target < vlo - options.energy_tolerance * scale
        || target > vhi + options.energy_tolerance * scale
    {
        return Err(Error::OutsideValidityDomain(format!(
            "GERG {name} target outside the endpoint equilibrium {} interval",
            match spec {
                Spec::Enthalpy => "enthalpy",
                Spec::Entropy => "entropy",
            }
        )));
    }
    for endpoint in [low, high] {
        if (spec.mixture(&endpoint) - target).abs() / scale <= options.energy_tolerance {
            return Ok(Solved {
                phases: endpoint.phases.clone(),
                iterations: 0,
                tp: Some(endpoint),
                sat: None,
            });
        }
    }
    for iteration in 1..=options.max_iterations {
        let secant = lo + (target - vlo) / (vhi - vlo) * (hi - lo);
        let margin = 0.1 * (hi - lo);
        let mid = if secant > lo + margin && secant < hi - margin {
            secant
        } else {
            lo + 0.5 * (hi - lo)
        };
        if mid == lo || mid == hi {
            return Err(Error::NonConvergence(format!(
                "GERG {name} bracket collapsed without convergence"
            )));
        }
        let flash = gerg_tp_flash(model, Temperature::new(mid)?, p, options.tp)?;
        let value = spec.mixture(&flash);
        if (value - target).abs() / scale <= options.energy_tolerance {
            return Ok(Solved {
                phases: flash.phases.clone(),
                iterations: iteration,
                tp: Some(flash),
                sat: None,
            });
        }
        if value < target {
            lo = mid;
            vlo = value;
        } else {
            hi = mid;
            vhi = value;
        }
    }
    Err(Error::NonConvergence(format!(
        "GERG {name} iteration limit exhausted"
    )))
}

/// Saturation pressure and bisection iterations of a pure GERG fluid at `t`,
/// in the 0.1..=10 MPa equilibrium envelope. On the exterior branches
/// ln φ^L − ln φ^V decreases with p and changes sign at saturation; a single
/// root is classified by its density relative to the reducing density.
pub(crate) fn gerg_pure_saturation_pressure(
    model: &GergModel,
    t: f64,
    limit: usize,
) -> Result<(f64, usize), Error> {
    let index = model
        .fractions()
        .iter()
        .position(|x| *x > 0.)
        .ok_or_else(|| Error::NumericalFailure("GERG pure model has no active component".into()))?;
    let temperature = Temperature::new(t)?;
    // exp(ln p) can round just outside the envelope at its ends.
    let envelope = |ln_p: f64| ln_p.exp().clamp(1e5, 1e7);
    let sample = |ln_p: f64| -> Result<(f64, bool), Error> {
        let roots = gerg_outer_roots(model, temperature, Pressure::new(envelope(ln_p))?)?;
        if roots.len() == 1 {
            let liquid = roots[0].state.molar_density > model.reducing_density();
            return Ok((if liquid { -1. } else { 1. }, false));
        }
        let difference = roots.last().unwrap().ln_fugacity_coefficients[index]
            - roots[0].ln_fugacity_coefficients[index];
        Ok((difference, true))
    };
    let (mut lo, mut hi) = (1e5_f64.ln(), 1e7_f64.ln());
    let (left, left_two) = sample(lo)?;
    let (right, right_two) = sample(hi)?;
    for (value, two, ln_p) in [(left, left_two, lo), (right, right_two, hi)] {
        if two && value.abs() <= 1e-12 {
            return Ok((envelope(ln_p), 0));
        }
    }
    if left < 0. || right > 0. {
        return Err(Error::OutsideValidityDomain(format!(
            "GERG pure saturation pressure at {t} K lies outside 0.1..=10 MPa"
        )));
    }
    for iteration in 1..=limit {
        let mid = lo + 0.5 * (hi - lo);
        let (difference, two) = sample(mid)?;
        if two && difference.abs() <= 1e-12 {
            return Ok((envelope(mid), iteration));
        }
        if mid == lo || mid == hi {
            break;
        }
        if difference > 0. {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(
        "GERG pure saturation pressure iteration exhausted".into(),
    ))
}

/// Saturation state of a pure GERG fluid at `p` in `lo..=hi` K, or `None` when
/// the saturation temperature is not bracketed there.
pub(crate) fn gerg_pure_saturation_temperature(
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
