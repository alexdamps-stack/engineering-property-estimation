//! Bounded PR/PR-BM bubble/dew boundaries; no critical-point or retrograde continuation.
use crate::{
    SATURATION_MAX_REDUCED_TEMPERATURE, StabilityAssessment, StabilityOptions, assess_stability,
    pure_saturation_pressure,
};
use thermo_models::{CubicModel, HomogeneousState, RootSelection};
use thermo_types::{CalculationError as Error, Pressure, Temperature};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubicSaturationKind {
    Bubble,
    Dew,
}

/// Bounds are part of the requested branch definition. Multiple physical roots
/// inside them are rejected rather than silently choosing a retrograde branch.
#[derive(Debug, Clone, Copy)]
pub struct SaturationOptions {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub minimum_pressure_pa: f64,
    pub maximum_pressure_pa: f64,
    pub max_iterations: usize,
}
impl Default for SaturationOptions {
    fn default() -> Self {
        Self {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 350.0,
            minimum_pressure_pa: 10_000.0,
            maximum_pressure_pa: 10_000_000.0,
            max_iterations: 100,
        }
    }
}
impl SaturationOptions {
    pub fn validate(&self) -> Result<(), Error> {
        let values = [
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            self.minimum_pressure_pa,
            self.maximum_pressure_pa,
        ];
        if values.iter().any(|v| !v.is_finite())
            || self.minimum_temperature_k < 200.0
            || self.maximum_temperature_k > 600.0
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || self.minimum_pressure_pa < 0.01
            || self.maximum_pressure_pa > 100_000_000.0
            || self.minimum_pressure_pa >= self.maximum_pressure_pa
            || self.max_iterations == 0
            || self.max_iterations > 2000
        {
            return Err(Error::OutsideValidityDomain("A13 requires ordered finite bounds inside 200..600 K, 0.01..1e8 Pa and iteration budget 1..2000".into()));
        }
        Ok(())
    }
}
#[derive(Debug, Clone)]
pub struct SaturationPoint {
    pub kind: CubicSaturationKind,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub liquid_composition: Vec<f64>,
    pub vapor_composition: Vec<f64>,
    pub liquid: HomogeneousState,
    pub vapor: HomogeneousState,
    pub max_ln_fugacity_residual: f64,
    /// Iterations in the accepted Newton path (other seeds are not included).
    pub iterations: usize,
    /// Finite-trial TPD evidence for liquid and vapor, in that order.
    pub phase_stability: [StabilityAssessment; 2],
}

/// Composition is the specified liquid for bubble, specified vapor for dew.
pub fn saturation_pressure(
    model: &CubicModel,
    temperature: Temperature,
    kind: CubicSaturationKind,
    options: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    solve(model, temperature.value(), true, kind, options)
}
pub fn saturation_temperature(
    model: &CubicModel,
    pressure: Pressure,
    kind: CubicSaturationKind,
    options: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    solve(model, pressure.value(), false, kind, options)
}

fn solve(
    model: &CubicModel,
    fixed: f64,
    pressure_unknown: bool,
    kind: CubicSaturationKind,
    o: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    o.validate()?;
    if model
        .component_ids()
        .iter()
        .zip(model.fractions())
        .any(|(id, x)| {
            *x > 0.0
                && !matches!(
                    *id,
                    "methane" | "ethane" | "propane" | "nitrogen" | "carbon-dioxide"
                )
        })
    {
        return Err(Error::UnsupportedCapability(
            "A13 qualified dry hydrocarbon/N2/CO2 cubic domain only".into(),
        ));
    }
    let (lo, hi) = if pressure_unknown {
        (o.minimum_pressure_pa.ln(), o.maximum_pressure_pa.ln())
    } else {
        (o.minimum_temperature_k.ln(), o.maximum_temperature_k.ln())
    };
    let (flo, fhi) = if pressure_unknown {
        (o.minimum_temperature_k, o.maximum_temperature_k)
    } else {
        (o.minimum_pressure_pa, o.maximum_pressure_pa)
    };
    if !fixed.is_finite() || fixed < flo || fixed > fhi {
        return Err(Error::OutsideValidityDomain(
            "A13 fixed coordinate outside requested supported range".into(),
        ));
    }
    let mut active_components = model
        .fractions()
        .iter()
        .enumerate()
        .filter(|(_, x)| **x > 0.0);
    if let (Some((index, _)), None) = (active_components.next(), active_components.next()) {
        // with_fractions accepts a small unit-sum error without normalizing.
        // Canonicalize the single active component for both pure saturation
        // paths and their returned phase states, retaining inactive slots.
        let mut fractions = vec![0.0; model.fractions().len()];
        fractions[index] = 1.0;
        let pure_model = model.with_fractions(&fractions)?;
        let model = &pure_model;
        let (t, p, it) = if pressure_unknown {
            let r = pure_saturation_pressure(model, Temperature::new(fixed)?, o.max_iterations)?;
            (fixed, r.pressure_pa, r.iterations)
        } else {
            let mut a = o.minimum_temperature_k;
            let mut b = o.maximum_temperature_k.min(
                SATURATION_MAX_REDUCED_TEMPERATURE
                    * model.components()[index].critical_temperature(),
            );
            if b <= a {
                return Err(Error::OutsideValidityDomain(
                    "A13 near-critical diagnostic: no supported pure-temperature interval".into(),
                ));
            }
            let p_a = pure_saturation_pressure(model, Temperature::new(a)?, 200)?.pressure_pa;
            let p_b = pure_saturation_pressure(model, Temperature::new(b)?, 200)?.pressure_pa;
            if fixed < p_a || fixed > p_b {
                return Err(Error::OutsideValidityDomain(
                    "A13 pure saturation not bracketed; near-critical or outside requested range"
                        .into(),
                ));
            }
            let mut count = 0;
            for i in 1..=o.max_iterations {
                count = i;
                let mid = (a + b) / 2.0;
                let r = pure_saturation_pressure(model, Temperature::new(mid)?, 200)?;
                if r.pressure_pa < fixed {
                    a = mid;
                } else {
                    b = mid;
                }
                if (b - a) < 1e-9 {
                    break;
                }
            }
            if b - a >= 1e-9 {
                return Err(Error::NonConvergence(
                    "A13 pure temperature iteration limit".into(),
                ));
            }
            ((a + b) / 2.0, fixed, count)
        };
        if p < o.minimum_pressure_pa || p > o.maximum_pressure_pa {
            return Err(Error::OutsideValidityDomain(
                "A13 saturation pressure outside requested range".into(),
            ));
        }
        return point(
            model,
            t,
            p,
            kind,
            model.fractions().to_vec(),
            model.fractions().to_vec(),
            it,
        );
    }
    let active: Vec<_> = model
        .fractions()
        .iter()
        .enumerate()
        .filter_map(|(i, x)| (*x > 0.0).then_some(i))
        .collect();
    let mut solutions: Vec<SaturationPoint> = Vec::new();
    let mut near_critical = false;
    for seed in 0..=12 {
        let u = lo + (hi - lo) * seed as f64 / 12.0;
        let (t, p) = if pressure_unknown {
            (fixed, u.exp())
        } else {
            (u.exp(), fixed)
        };
        let k = crate::wilson_k(model, Temperature::new(t)?, Pressure::new(p)?);
        let mut q = vec![u];
        q.extend(active.iter().map(|i| k[*i].ln()));
        for iteration in 1..=o.max_iterations {
            let Ok((r, x, y)) = residual(model, fixed, pressure_unknown, kind, &active, &q) else {
                break;
            };
            let norm = r.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
            if norm < 2e-11 {
                let (t, p) = if pressure_unknown {
                    (fixed, q[0].exp())
                } else {
                    (q[0].exp(), fixed)
                };
                match point(model, t, p, kind, x, y, iteration) {
                    Ok(s) => {
                        if !solutions.iter().any(|old| {
                            (old.temperature_k - t).abs() < 1e-5
                                && (old.pressure_pa / p - 1.0).abs() < 1e-7
                        }) {
                            solutions.push(s);
                        }
                    }
                    Err(Error::OutsideValidityDomain(_)) => near_critical = true,
                    _ => {}
                }
                break;
            }
            let n = q.len();
            let mut jac = vec![vec![0.0; n]; n];
            let mut valid = true;
            for j in 0..n {
                let mut shifted = q.clone();
                shifted[j] += 1e-5;
                if let Ok((v, _, _)) =
                    residual(model, fixed, pressure_unknown, kind, &active, &shifted)
                {
                    for i in 0..n {
                        jac[i][j] = (v[i] - r[i]) / 1e-5;
                    }
                } else {
                    valid = false;
                    break;
                }
            }
            if !valid {
                break;
            }
            let Some(step) = linear(jac, r.iter().map(|v| -v).collect()) else {
                near_critical = true;
                break;
            };
            let mut accepted = false;
            let mut scale = 1.0;
            for _ in 0..20 {
                let next: Vec<_> = q.iter().zip(&step).map(|(a, b)| a + scale * b).collect();
                if next[0] >= lo
                    && next[0] <= hi
                    && next.iter().all(|v| v.is_finite() && v.abs() < 100.0)
                    && let Ok((v, _, _)) =
                        residual(model, fixed, pressure_unknown, kind, &active, &next)
                    && v.iter().fold(0.0_f64, |a, v| a.max(v.abs())) < norm
                {
                    q = next;
                    accepted = true;
                    break;
                }
                scale *= 0.5;
            }
            if !accepted {
                break;
            }
        }
    }
    match solutions.len() {
        1=>Ok(solutions.remove(0)),
        0 if near_critical=>Err(Error::OutsideValidityDomain("A13 near-critical diagnostic: trivial/coalescing roots or singular boundary equations; no verified point returned".into())),
        0=>Err(Error::NonConvergence("A13 no verified boundary inside requested bounds; iteration limit, absent bracket or unsupported branch".into())),
        _=>Err(Error::UnsupportedCapability("A13 multiple saturation branches in requested interval; narrow bounds; retrograde branch selection unsupported".into())),
    }
}
type ResidualEvaluation = (Vec<f64>, Vec<f64>, Vec<f64>);

fn residual(
    model: &CubicModel,
    fixed: f64,
    pressure_unknown: bool,
    kind: CubicSaturationKind,
    active: &[usize],
    q: &[f64],
) -> Result<ResidualEvaluation, Error> {
    let (t, p) = if pressure_unknown {
        (fixed, q[0].exp())
    } else {
        (q[0].exp(), fixed)
    };
    let z = model.fractions();
    let mut x = z.to_vec();
    let mut y = z.to_vec();
    for (j, i) in active.iter().enumerate() {
        match kind {
            CubicSaturationKind::Bubble => y[*i] = z[*i] * q[j + 1].exp(),
            CubicSaturationKind::Dew => x[*i] = z[*i] * (-q[j + 1]).exp(),
        }
    }
    let sum = match kind {
        CubicSaturationKind::Bubble => y.iter().sum::<f64>(),
        CubicSaturationKind::Dew => x.iter().sum::<f64>(),
    };
    if !sum.is_finite() || sum <= 0.0 {
        return Err(Error::NumericalFailure(
            "A13 invalid incipient composition".into(),
        ));
    }
    match kind {
        CubicSaturationKind::Bubble => y.iter_mut().for_each(|v| *v /= sum),
        CubicSaturationKind::Dew => x.iter_mut().for_each(|v| *v /= sum),
    };
    let t = Temperature::new(t)?;
    let p = Pressure::new(p)?;
    let l = model
        .with_fractions(&x)?
        .ln_fugacity(t, p, RootSelection::SmallestVolume)?;
    let v = model
        .with_fractions(&y)?
        .ln_fugacity(t, p, RootSelection::LargestVolume)?;
    let mut r = vec![sum.ln()];
    r.extend(
        active.iter().enumerate().map(|(j, i)| {
            q[j + 1] + v.ln_fugacity_coefficients[*i] - l.ln_fugacity_coefficients[*i]
        }),
    );
    Ok((r, x, y))
}
fn point(
    model: &CubicModel,
    t: f64,
    p: f64,
    kind: CubicSaturationKind,
    x: Vec<f64>,
    y: Vec<f64>,
    iterations: usize,
) -> Result<SaturationPoint, Error> {
    let tm = Temperature::new(t)?;
    let pr = Pressure::new(p)?;
    let lm = model.with_fractions(&x)?;
    let vm = model.with_fractions(&y)?;
    let l = lm.evaluate(tm, pr, RootSelection::SmallestVolume)?;
    let v = vm.evaluate(tm, pr, RootSelection::LargestVolume)?;
    // Small absolute composition differences also occur in dilute mixtures,
    // far from criticality. Require distinct physical phase volumes instead;
    // the trivial K=1 solution has coalescing roots and fails this check.
    if v.molar_volume / l.molar_volume < 1.5 {
        return Err(Error::OutsideValidityDomain(
            "A13 near-critical diagnostic: phase volume ratio <1.5".into(),
        ));
    }
    let residual = x
        .iter()
        .zip(&y)
        .enumerate()
        .filter(|(_, (a, b))| **a > 0.0 && **b > 0.0)
        .map(|(i, (a, b))| {
            (a.ln() + l.ln_fugacity_coefficients[i] - b.ln() - v.ln_fugacity_coefficients[i]).abs()
        })
        .fold(0.0, f64::max);
    if residual > 1e-9 {
        return Err(Error::NonConvergence("A13 fugacity equality failed".into()));
    }
    // TPD uses the minimum-Gibbs root. Ensure it assesses the selected phase,
    // rather than silently assessing a different root at the same composition.
    for (phase, composition, selected) in [(&lm, &x, &l), (&vm, &y, &v)] {
        let stable = phase.ln_fugacity(tm, pr, RootSelection::MinimumGibbs)?;
        let excess: f64 = composition
            .iter()
            .zip(&selected.ln_fugacity_coefficients)
            .zip(&stable.ln_fugacity_coefficients)
            .map(|((w, a), b)| w * (a - b))
            .sum();
        if excess > 1e-8 {
            return Err(Error::UnsupportedCapability(
                "A13 selected phase is metastable relative to another root".into(),
            ));
        }
    }
    let phase_stability = [
        assess_stability(&lm, tm, pr, StabilityOptions::default())?,
        assess_stability(&vm, tm, pr, StabilityOptions::default())?,
    ];
    for phase in &phase_stability {
        if !phase.stable_in_trials {
            return Err(Error::UnsupportedCapability(
                "A13 boundary phases unstable in finite TPD trials".into(),
            ));
        }
    }
    Ok(SaturationPoint {
        kind,
        temperature_k: t,
        pressure_pa: p,
        liquid_composition: x,
        vapor_composition: y,
        liquid: l,
        vapor: v,
        max_ln_fugacity_residual: residual,
        iterations,
        phase_stability,
    })
}
fn linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n).max_by(|i, j| a[*i][k].abs().total_cmp(&a[*j][k].abs()))?;
        if a[pivot][k].abs() < 1e-12 {
            return None;
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        let pivot_row = a[k].clone();
        for i in k + 1..n {
            let f = a[i][k] / pivot_row[k];
            for (entry, pivot) in a[i][k..].iter_mut().zip(&pivot_row[k..]) {
                *entry -= f * pivot;
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = (b[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    Some(x)
}

#[derive(Debug, Clone)]
pub enum EnvelopeStop {
    RequestedRangeComplete,
    NearCritical {
        temperature_k: f64,
        diagnostic: String,
    },
    UnverifiedBoundary {
        temperature_k: f64,
        diagnostic: String,
    },
}
#[derive(Debug, Clone)]
pub struct PhaseEnvelope {
    pub kind: CubicSaturationKind,
    pub requested_range: SaturationOptions,
    /// Actual extrema of verified samples; None when the first sample failed.
    pub supported_temperature_range_k: Option<[f64; 2]>,
    pub supported_pressure_range_pa: Option<[f64; 2]>,
    pub points: Vec<SaturationPoint>,
    pub stop: EnvelopeStop,
}
/// Temperature-parameterized bounded sampling of one boundary. Stops at the
/// first failed/near-critical point; never bridges gaps or extrapolates branches.
pub fn trace_phase_envelope(
    model: &CubicModel,
    kind: CubicSaturationKind,
    step_k: f64,
    maximum_points: usize,
    o: SaturationOptions,
) -> Result<PhaseEnvelope, Error> {
    o.validate()?;
    if !step_k.is_finite() || step_k <= 0.0 || maximum_points == 0 || maximum_points > 1000 {
        return Err(Error::UnsupportedCapability(
            "A13 trace requires positive finite step and 1..1000 points".into(),
        ));
    }
    let mut result = PhaseEnvelope {
        kind,
        requested_range: o,
        supported_temperature_range_k: None,
        supported_pressure_range_pa: None,
        points: Vec::new(),
        stop: EnvelopeStop::RequestedRangeComplete,
    };
    let intervals = ((o.maximum_temperature_k - o.minimum_temperature_k) / step_k).ceil();
    if !intervals.is_finite() || intervals >= maximum_points as f64 {
        return Err(Error::UnsupportedCapability(
            "A13 trace point budget cannot cover requested range".into(),
        ));
    }
    let required = intervals as usize + 1;
    for i in 0..required {
        let t = (o.minimum_temperature_k + i as f64 * step_k).min(o.maximum_temperature_k);
        match saturation_pressure(model, Temperature::new(t)?, kind, o) {
            Ok(p) => result.points.push(p),
            Err(e) => {
                let message = e.to_string();
                result.stop = if message.contains("near-critical") {
                    EnvelopeStop::NearCritical {
                        temperature_k: t,
                        diagnostic: message,
                    }
                } else {
                    EnvelopeStop::UnverifiedBoundary {
                        temperature_k: t,
                        diagnostic: message,
                    }
                };
                break;
            }
        }
    }
    if let (Some(first), Some(last)) = (result.points.first(), result.points.last()) {
        result.supported_temperature_range_k = Some([first.temperature_k, last.temperature_k]);
        result.supported_pressure_range_pa = Some([
            result
                .points
                .iter()
                .map(|p| p.pressure_pa)
                .fold(f64::INFINITY, f64::min),
            result
                .points
                .iter()
                .map(|p| p.pressure_pa)
                .fold(0.0, f64::max),
        ]);
    }
    Ok(result)
}
