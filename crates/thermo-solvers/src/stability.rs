//! Michelsen tangent-plane-distance (TPD) stability assessment with a finite,
//! documented trial set. A negative TPD at any trial proves instability; the
//! absence of negative values among finite trials is evidence, not proof, of
//! global stability.
use crate::equilibrium_model::EquilibriumModel;
use thermo_models::{CubicModel, RootSelection};
use thermo_types::{CalculationError as Error, Pressure, Temperature};

/// Iteration limits and thresholds for each TPD trial.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StabilityOptions {
    /// Successive-substitution limit per trial.
    pub max_iterations: usize,
    /// Convergence: max |Δ ln W_i| between iterations.
    pub tolerance: f64,
    /// Modified TPD below this value classifies the feed as unstable.
    pub tm_threshold: f64,
    /// Add one near-pure trial per active component (n ≥ 2).
    pub pure_component_trials: bool,
}

impl Default for StabilityOptions {
    fn default() -> Self {
        Self {
            max_iterations: 2000,
            tolerance: 1e-10,
            tm_threshold: -1e-8,
            pure_component_trials: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrialInitialization {
    /// W_i = z_i K_i with Wilson K.
    WilsonVapor,
    /// W_i = z_i / K_i with Wilson K.
    WilsonLiquid,
    /// Near-pure trial of this component index (model order).
    PureComponent(usize),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrialOutcome {
    /// Nontrivial stationary point.
    Converged,
    /// Converged to the feed composition (no information).
    TrivialSolution,
    /// Limit reached after a TPD below the threshold was found: instability is shown,
    /// but the stationary point was not located.
    IterationLimitAfterNegativeTpd,
    /// Limit reached without a TPD below the threshold: inconclusive. The assessment
    /// returns `NonConvergence` unless another trial proved instability.
    IterationLimitInconclusive,
    /// Stopped once the iterate came within 1e-6 (mole fraction) of the nontrivial
    /// stationary point already located by the trial at this index of `trials`; it
    /// carries no new information and is not used to initialise a split.
    DuplicateOf(usize),
}

impl TrialOutcome {
    /// True for a trial whose composition is a usable nontrivial stationary point
    /// (or a negative-TPD iterate) for classifying the feed and initialising a split.
    pub fn locates_stationary_point(self) -> bool {
        matches!(self, Self::Converged | Self::IterationLimitAfterNegativeTpd)
    }
}

#[derive(Debug, Clone)]
pub struct TrialResult {
    pub initialization: TrialInitialization,
    pub outcome: TrialOutcome,
    /// Modified TPD tm = 1 + Σ W_i (ln W_i + ln φ_i(w) − ln z_i − ln φ_i(z) − 1).
    pub tangent_plane_distance: f64,
    /// Normalized trial composition w = W / ΣW in model order (zeros for inactive).
    pub composition: Vec<f64>,
    pub compressibility: f64,
    pub iterations: usize,
}

#[derive(Debug, Clone)]
pub struct StabilityAssessment {
    /// True when no trial produced tm below the threshold. Finite trials only.
    pub stable_in_trials: bool,
    /// Minimum-Gibbs root of the tested composition.
    pub feed_compressibility: f64,
    pub trials: Vec<TrialResult>,
    pub options: StabilityOptions,
}

impl StabilityAssessment {
    /// Trial with the lowest TPD, if any trial ran.
    pub fn most_negative(&self) -> Option<&TrialResult> {
        self.trials.iter().min_by(|a, b| {
            a.tangent_plane_distance
                .total_cmp(&b.tangent_plane_distance)
        })
    }
}

/// Wilson K_i = (Pc_i/p) exp(5.373 (1+ω_i)(1 − Tc_i/T)) in model order.
pub fn wilson_k(model: &CubicModel, temperature: Temperature, pressure: Pressure) -> Vec<f64> {
    model
        .components()
        .iter()
        .map(|c| {
            c.critical_pressure() / pressure.value()
                * (5.373
                    * (1.0 + c.acentric_factor())
                    * (1.0 - c.critical_temperature() / temperature.value()))
                .exp()
        })
        .collect()
}

/// Assess the stability of the model's composition at T and p. Each trial is a
/// successive-substitution search for a TPD stationary point using the minimum-Gibbs
/// root at the trial composition. Trial set: Wilson vapor-like, Wilson liquid-like
/// and (optionally) one near-pure trial per active component.
///
/// The first trial that proves instability (tm below the threshold) ends the
/// assessment; later trials are not run. A trial that exhausts its iteration limit
/// without a negative TPD is inconclusive: the assessment returns `NonConvergence`
/// only when no trial proved instability.
pub fn assess_stability(
    model: &CubicModel,
    temperature: Temperature,
    pressure: Pressure,
    options: StabilityOptions,
) -> Result<StabilityAssessment, Error> {
    assess_stability_impl(model, temperature, pressure, options)
}

pub(crate) fn assess_stability_impl<M: EquilibriumModel>(
    model: &M,
    temperature: Temperature,
    pressure: Pressure,
    options: StabilityOptions,
) -> Result<StabilityAssessment, Error> {
    let z = model.fractions();
    let active: Vec<usize> = (0..z.len()).filter(|i| z[*i] > 0.0).collect();
    let feed = model.ln_fugacity(temperature, pressure, RootSelection::MinimumGibbs)?;
    let mut assessment = StabilityAssessment {
        stable_in_trials: true,
        feed_compressibility: feed.compressibility,
        trials: Vec::new(),
        options,
    };
    // A single active component has only the trivial stationary point once the
    // minimum-Gibbs root is selected.
    if active.len() < 2 {
        return Ok(assessment);
    }
    let d: Vec<f64> = (0..z.len())
        .map(|i| {
            if z[i] > 0.0 {
                z[i].ln() + feed.ln_fugacity_coefficients[i]
            } else {
                0.0
            }
        })
        .collect();
    let k = model.wilson_values(temperature, pressure);
    let mut starts = vec![
        (
            TrialInitialization::WilsonVapor,
            active.iter().map(|&i| z[i] * k[i]).collect::<Vec<_>>(),
        ),
        (
            TrialInitialization::WilsonLiquid,
            active.iter().map(|&i| z[i] / k[i]).collect(),
        ),
    ];
    if options.pure_component_trials {
        for (position, &j) in active.iter().enumerate() {
            let others = (active.len() - 1) as f64;
            starts.push((
                TrialInitialization::PureComponent(j),
                (0..active.len())
                    .map(|p| {
                        if p == position {
                            1.0 - 1e-3
                        } else {
                            1e-3 / others
                        }
                    })
                    .collect(),
            ));
        }
    }
    for (initialization, w0) in starts {
        let trial = run_trial(
            model,
            temperature,
            pressure,
            &active,
            z,
            &d,
            initialization,
            w0,
            options,
            &assessment.trials,
        )?;
        let proves_instability = trial.tangent_plane_distance < options.tm_threshold
            && trial.outcome.locates_stationary_point();
        assessment.trials.push(trial);
        if proves_instability {
            assessment.stable_in_trials = false;
            // Instability is proven; a slower remaining trial must not discard it.
            return Ok(assessment);
        }
    }
    if let Some(slow) = assessment
        .trials
        .iter()
        .find(|t| t.outcome == TrialOutcome::IterationLimitInconclusive)
    {
        return Err(Error::NonConvergence(format!(
            "TPD trial {:?} exhausted {} iterations without a negative TPD; stability undetermined",
            slow.initialization, options.max_iterations
        )));
    }
    Ok(assessment)
}

/// Mole-fraction distance below which a trial iterate is taken to have reached an
/// earlier trial's stationary point (same scale as the trivial-solution test).
const DUPLICATE_DISTANCE: f64 = 1e-6;

#[allow(clippy::too_many_arguments)]
fn run_trial<M: EquilibriumModel>(
    model: &M,
    temperature: Temperature,
    pressure: Pressure,
    active: &[usize],
    z: &[f64],
    d: &[f64],
    initialization: TrialInitialization,
    w0: Vec<f64>,
    options: StabilityOptions,
    earlier: &[TrialResult],
) -> Result<TrialResult, Error> {
    let n = z.len();
    let mut ln_w: Vec<f64> = w0.iter().map(|w| w.ln()).collect();
    let mut best: Option<(f64, Vec<f64>, f64)> = None;
    for iteration in 1..=options.max_iterations {
        let big_w: Vec<f64> = ln_w.iter().map(|l| l.exp()).collect();
        let total: f64 = big_w.iter().sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(Error::NumericalFailure(
                "TPD trial amounts left the representable range".into(),
            ));
        }
        let w = normalized(model_fractions_of(active, &big_w, total, n));
        let fug = model.with_fractions(&w)?.ln_fugacity(
            temperature,
            pressure,
            RootSelection::MinimumGibbs,
        )?;
        let mut tm = 1.0;
        let mut change: f64 = 0.0;
        let mut trivial_distance: f64 = 0.0;
        for (p, &i) in active.iter().enumerate() {
            let phi = fug.ln_fugacity_coefficients[i];
            tm += big_w[p] * (ln_w[p] + phi - d[i] - 1.0);
            let next = d[i] - phi;
            change = change.max((next - ln_w[p]).abs());
            trivial_distance = trivial_distance.max((w[i] - z[i]).abs());
            ln_w[p] = next;
        }
        if !tm.is_finite() {
            return Err(Error::NumericalFailure("nonfinite TPD".into()));
        }
        if best.as_ref().is_none_or(|(b, _, _)| tm < *b) {
            best = Some((tm, w.clone(), fug.compressibility));
        }
        let trivial = trivial_distance <= 1e-6;
        if !trivial && change <= DUPLICATE_DISTANCE {
            let duplicate = earlier.iter().position(|t| {
                t.outcome == TrialOutcome::Converged
                    && active
                        .iter()
                        .all(|&i| (t.composition[i] - w[i]).abs() <= DUPLICATE_DISTANCE)
            });
            if let Some(index) = duplicate {
                return Ok(TrialResult {
                    initialization,
                    outcome: TrialOutcome::DuplicateOf(index),
                    tangent_plane_distance: tm,
                    composition: w,
                    compressibility: fug.compressibility,
                    iterations: iteration,
                });
            }
        }
        if change <= options.tolerance || (trivial && change <= 1e-8) {
            return Ok(TrialResult {
                initialization,
                outcome: if trivial {
                    TrialOutcome::TrivialSolution
                } else {
                    TrialOutcome::Converged
                },
                tangent_plane_distance: tm,
                composition: w,
                compressibility: fug.compressibility,
                iterations: iteration,
            });
        }
    }
    let Some((tm, composition, compressibility)) = best else {
        return Err(Error::NonConvergence(format!(
            "TPD trial {initialization:?} has a zero iteration limit; stability undetermined"
        )));
    };
    Ok(TrialResult {
        initialization,
        outcome: if tm < options.tm_threshold {
            TrialOutcome::IterationLimitAfterNegativeTpd
        } else {
            TrialOutcome::IterationLimitInconclusive
        },
        tangent_plane_distance: tm,
        composition,
        compressibility,
        iterations: options.max_iterations,
    })
}

fn model_fractions_of(active: &[usize], big_w: &[f64], total: f64, n: usize) -> Vec<f64> {
    let mut w = vec![0.0; n];
    for (p, &i) in active.iter().enumerate() {
        w[i] = big_w[p] / total;
    }
    w
}

/// Remove last-bit roundoff so the model's unit-sum check applies to exact input.
pub(crate) fn normalized(mut x: Vec<f64>) -> Vec<f64> {
    let sum: f64 = x.iter().sum();
    x.iter_mut().for_each(|v| *v /= sum);
    x
}
