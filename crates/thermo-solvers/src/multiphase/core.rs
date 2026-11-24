//! Model-independent multiphase stability trials, phase-amount minimisation and
//! the successive-substitution / Newton split.

use thermo_types::{CalculationError as Error, Diagnostic, DiagnosticSeverity, Phase};

use super::adapters::{VlleAdapter, is_hydrocarbon, water_index};
use super::{
    FluidKind, MultiphaseOptions, MultiphasePhase, MultiphaseStability, MultiphaseTpFlash,
    MultiphaseTrial, MultiphaseTrialResult, PhaseRole,
};
use crate::flash::{dominant_eigenvalue_step, newton_step};
use crate::stability::{TrialOutcome, normalized};

/// At most vapor–liquid–liquid.
pub(crate) const MAX_PHASES: usize = 3;
/// Same-kind phases closer than this (max |Δx|) are merged during the split.
const MERGE_DISTANCE: f64 = 1e-5;
/// A nonconverged split with phases closer than this is reported as near a plait point.
pub(crate) const PLAIT_FAILURE_DISTANCE: f64 = 0.05;
/// A converged split with phases closer than this carries a warning.
pub(crate) const PLAIT_WARNING_DISTANCE: f64 = 0.01;
/// Phase-addition rounds after the first split before failing explicitly.
const MAX_ROUNDS: usize = 4;
/// A trial iterate this close to an earlier stationary point is a duplicate.
const DUPLICATE_DISTANCE: f64 = 1e-6;
/// Initial or added phases closer than this to an existing phase are not added.
const DISTINCT_DISTANCE: f64 = 1e-4;
/// Near-pure trial purity.
const TRIAL_IMPURITY: f64 = 1e-3;

pub(crate) fn diagnostic(
    code: &str,
    message: impl Into<String>,
    severity: DiagnosticSeverity,
) -> Diagnostic {
    Diagnostic {
        code: code.into(),
        severity,
        message: message.into(),
    }
}

pub(crate) fn active(z: &[f64]) -> Vec<usize> {
    (0..z.len()).filter(|i| z[*i] > 0.0).collect()
}

pub(crate) fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

/// Σ x_i (ln x_i + ln φ_i) over active components: G/(RT) relative to the pure
/// ideal gases at T and p.
fn reduced_gibbs(x: &[f64], ln_phi: &[f64], active: &[usize]) -> f64 {
    active.iter().map(|&i| x[i] * (x[i].ln() + ln_phi[i])).sum()
}

/// Minimum-Gibbs kind of composition `x`, its ln φ and reduced Gibbs energy.
pub(crate) fn reference<A: VlleAdapter>(
    a: &A,
    x: &[f64],
    active: &[usize],
) -> Result<(FluidKind, Vec<f64>, f64), Error> {
    let mut best: Option<(FluidKind, Vec<f64>, f64)> = None;
    for &kind in a.kinds() {
        let ln_phi = a.ln_phi(x, kind)?;
        let g = reduced_gibbs(x, &ln_phi, active);
        if !g.is_finite() {
            return Err(Error::NumericalFailure(
                "nonfinite reference Gibbs energy".into(),
            ));
        }
        if best.as_ref().is_none_or(|(_, _, b)| g < *b) {
            best = Some((kind, ln_phi, g));
        }
    }
    best.ok_or_else(|| Error::NumericalFailure("adapter declares no phase kind".into()))
}

/// Trial set of the module documentation, as active-component amounts.
fn trial_starts<A: VlleAdapter>(
    a: &A,
    z: &[f64],
    active: &[usize],
) -> Vec<(MultiphaseTrial, FluidKind, Vec<f64>)> {
    let k = a.wilson_k();
    let len = active.len();
    let near_pure = |position: usize| -> Vec<f64> {
        (0..len)
            .map(|q| {
                if q == position {
                    1.0 - TRIAL_IMPURITY
                } else {
                    TRIAL_IMPURITY / (len - 1) as f64
                }
            })
            .collect()
    };
    let mut starts = vec![
        (
            MultiphaseTrial::WilsonVapor,
            a.vapor_kind(),
            active.iter().map(|&i| z[i] * k[i]).collect(),
        ),
        (
            MultiphaseTrial::WilsonLiquid,
            a.liquid_kind(),
            active.iter().map(|&i| z[i] / k[i]).collect(),
        ),
    ];
    let components = a.components();
    let water = water_index(components).and_then(|w| active.iter().position(|&i| i == w));
    if let Some(wpos) = water {
        starts.push((MultiphaseTrial::Aqueous, a.liquid_kind(), near_pure(wpos)));
        let others: Vec<usize> = (0..len).filter(|&q| q != wpos).collect();
        let hydrocarbons: Vec<usize> = others
            .iter()
            .copied()
            .filter(|&q| is_hydrocarbon(&components[active[q]]))
            .collect();
        let group = if hydrocarbons.is_empty() {
            others
        } else {
            hydrocarbons
        };
        let total: f64 = group.iter().map(|&q| z[active[q]]).sum();
        let rest = len - group.len();
        starts.push((
            MultiphaseTrial::Hydrocarbon,
            a.liquid_kind(),
            (0..len)
                .map(|q| {
                    if group.contains(&q) {
                        (1.0 - TRIAL_IMPURITY) * z[active[q]] / total
                    } else {
                        TRIAL_IMPURITY / rest as f64
                    }
                })
                .collect(),
        ));
    }
    for (q, &j) in active.iter().enumerate() {
        if Some(q) != water {
            starts.push((
                MultiphaseTrial::PureComponent(j),
                a.liquid_kind(),
                near_pure(q),
            ));
        }
    }
    starts
}

fn expand(active: &[usize], values: &[f64], n: usize) -> Vec<f64> {
    let mut full = vec![0.0; n];
    for (q, &i) in active.iter().enumerate() {
        full[i] = values[q];
    }
    normalized(full)
}

/// Stability of composition `x` relative to its minimum-Gibbs kind.
pub(crate) fn stability<A: VlleAdapter>(
    a: &A,
    x: &[f64],
    options: MultiphaseOptions,
) -> Result<MultiphaseStability, Error> {
    let act = active(x);
    let (kind, ln_phi, _) = reference(a, x, &act)?;
    stability_with(a, x, kind, &ln_phi, options)
}

/// Stability relative to a reference phase of the given kind and ln φ.
fn stability_with<A: VlleAdapter>(
    a: &A,
    x: &[f64],
    reference_kind: FluidKind,
    ln_phi: &[f64],
    options: MultiphaseOptions,
) -> Result<MultiphaseStability, Error> {
    let act = active(x);
    let mut result = MultiphaseStability {
        stable_in_trials: true,
        reference_kind,
        trials: Vec::new(),
    };
    if act.len() < 2 {
        return Ok(result);
    }
    let d: Vec<f64> = act.iter().map(|&i| x[i].ln() + ln_phi[i]).collect();
    for (trial, kind, w0) in trial_starts(a, x, &act) {
        let r = run_trial(
            a,
            x,
            &act,
            &d,
            reference_kind,
            trial,
            kind,
            &w0,
            options,
            &result.trials,
        )?;
        result.trials.push(r);
    }
    result.stable_in_trials = result.unstable_trials(options.tm_threshold).is_empty();
    if result.stable_in_trials
        && let Some(slow) = result
            .trials
            .iter()
            .find(|t| t.outcome == TrialOutcome::IterationLimitInconclusive)
    {
        return Err(Error::NonConvergence(format!(
            "multiphase TPD trial {:?} exhausted {} iterations without a negative TPD; stability undetermined (trials converge slowly near a plait or critical point)",
            slow.trial, options.trial_max_iterations
        )));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn run_trial<A: VlleAdapter>(
    a: &A,
    x: &[f64],
    act: &[usize],
    d: &[f64],
    reference_kind: FluidKind,
    trial: MultiphaseTrial,
    kind: FluidKind,
    w0: &[f64],
    options: MultiphaseOptions,
    earlier: &[MultiphaseTrialResult],
) -> Result<MultiphaseTrialResult, Error> {
    let n = x.len();
    let mut ln_w: Vec<f64> = w0.iter().map(|w| w.ln()).collect();
    let mut best: Option<(f64, Vec<f64>)> = None;
    let result = |outcome, tm, composition, iterations| MultiphaseTrialResult {
        trial,
        kind,
        outcome,
        tangent_plane_distance: tm,
        composition,
        iterations,
    };
    for iteration in 1..=options.trial_max_iterations {
        let big: Vec<f64> = ln_w.iter().map(|l| l.exp()).collect();
        let total: f64 = big.iter().sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(Error::NumericalFailure(
                "TPD trial amounts left the representable range".into(),
            ));
        }
        let w = expand(act, &big, n);
        let phi = a.ln_phi(&w, kind)?;
        let mut tm = 1.0;
        let mut change: f64 = 0.0;
        let mut trivial_distance: f64 = 0.0;
        for (q, &i) in act.iter().enumerate() {
            tm += big[q] * (ln_w[q] + phi[i] - d[q] - 1.0);
            let next = d[q] - phi[i];
            change = change.max((next - ln_w[q]).abs());
            trivial_distance = trivial_distance.max((w[i] - x[i]).abs());
            ln_w[q] = next;
        }
        if !tm.is_finite() {
            return Err(Error::NumericalFailure("nonfinite multiphase TPD".into()));
        }
        if best.as_ref().is_none_or(|(b, _)| tm < *b) {
            best = Some((tm, w.clone()));
        }
        let trivial = kind == reference_kind && trivial_distance <= 1e-6;
        if !trivial && change <= DUPLICATE_DISTANCE {
            let duplicate = earlier.iter().position(|t| {
                t.outcome == TrialOutcome::Converged
                    && t.kind == kind
                    && distance(&t.composition, &w) <= DUPLICATE_DISTANCE
            });
            if let Some(index) = duplicate {
                return Ok(result(TrialOutcome::DuplicateOf(index), tm, w, iteration));
            }
        }
        if change <= options.trial_tolerance || (trivial && change <= 1e-8) {
            let outcome = if trivial {
                TrialOutcome::TrivialSolution
            } else {
                TrialOutcome::Converged
            };
            return Ok(result(outcome, tm, w, iteration));
        }
    }
    let (tm, w) = best.expect("at least one trial iteration");
    let outcome = if tm < options.tm_threshold {
        TrialOutcome::IterationLimitAfterNegativeTpd
    } else {
        TrialOutcome::IterationLimitInconclusive
    };
    Ok(result(outcome, tm, w, options.trial_max_iterations))
}

/// Dense solve with partial pivoting; `None` when singular.
#[allow(clippy::needless_range_loop)]
pub(crate) fn solve_dense(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for j in 0..n {
        let pivot = (j..n).max_by(|&p, &q| m[p][j].abs().total_cmp(&m[q][j].abs()))?;
        m.swap(j, pivot);
        b.swap(j, pivot);
        if !m[j][j].is_finite() || m[j][j].abs() < 1e-300 {
            return None;
        }
        for i in j + 1..n {
            let f = m[i][j] / m[j][j];
            for c in j..n {
                m[i][c] -= f * m[j][c];
            }
            b[i] -= f * b[j];
        }
    }
    for i in (0..n).rev() {
        for c in i + 1..n {
            b[i] -= m[i][c] * b[c];
        }
        b[i] /= m[i][i];
    }
    b.iter().all(|v| v.is_finite()).then_some(b)
}

/// Shifted weights a_kq = exp(ln_a_kq − max_l ln_a_lq) ∈ (0, 1].
fn shifted(ln_a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let nc = ln_a[0].len();
    let shift: Vec<f64> = (0..nc)
        .map(|q| {
            ln_a.iter()
                .map(|row| row[q])
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .collect();
    ln_a.iter()
        .map(|row| row.iter().zip(&shift).map(|(l, s)| (l - s).exp()).collect())
        .collect()
}

/// Michelsen's multiphase Rachford–Rice: β ≥ 0 minimising
/// Q(β) = Σ_k β_k − Σ_q z_q ln Σ_k β_k a_kq with a_kq = 1/φ_kq = exp(ln_a_kq).
pub(crate) fn phase_amounts(
    z: &[f64],
    ln_a: &[Vec<f64>],
    start: &[f64],
) -> Result<Vec<f64>, Error> {
    let np = ln_a.len();
    if ln_a.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Error::NumericalFailure(
            "nonfinite fugacity coefficient in the phase-amount problem".into(),
        ));
    }
    let a = shifted(ln_a);
    let e = |b: &[f64]| -> Vec<f64> {
        (0..z.len())
            .map(|q| (0..np).map(|k| b[k] * a[k][q]).sum())
            .collect()
    };
    let objective = |b: &[f64]| -> f64 {
        let e = e(b);
        if e.iter().any(|v| v.is_nan() || *v <= 0.0) {
            return f64::INFINITY;
        }
        b.iter().sum::<f64>() - z.iter().zip(&e).map(|(z, e)| z * e.ln()).sum::<f64>()
    };
    let mut beta: Vec<f64> = start.iter().map(|b| b.max(0.0)).collect();
    if beta.len() != np
        || beta.iter().sum::<f64>().is_nan()
        || beta.iter().sum::<f64>() <= 0.0
        || beta.iter().any(|b| !b.is_finite())
    {
        beta = vec![1.0 / np as f64; np];
    }
    // Gradient of Q and the KKT residual (|g_k| for β_k > 0, max(0, −g_k) at β_k = 0).
    let gradient_and_kkt = |b: &[f64]| -> (Vec<f64>, f64, Vec<f64>) {
        let ev = e(b);
        let g: Vec<f64> = (0..np)
            .map(|k| 1.0 - (0..z.len()).map(|q| z[q] * a[k][q] / ev[q]).sum::<f64>())
            .collect();
        let kkt = (0..np)
            .map(|k| {
                if b[k] > 0.0 {
                    g[k].abs()
                } else {
                    (-g[k]).max(0.0)
                }
            })
            .fold(0.0, f64::max);
        (g, kkt, ev)
    };
    let mut previous_kkt = f64::INFINITY;
    for _ in 0..500 {
        let (g, kkt, ev) = gradient_and_kkt(&beta);
        // Converged, or stalled at roundoff (no halving of the KKT residual).
        if kkt <= 1e-13 || (kkt <= 1e-10 && kkt > 0.5 * previous_kkt) {
            return Ok(beta);
        }
        previous_kkt = kkt;
        let mut free: Vec<usize> = (0..np).filter(|&k| beta[k] > 0.0 || g[k] < 0.0).collect();
        // With more free phases than components H is singular: along the null space
        // of the weights E is unchanged and Q is linear, so step to the boundary.
        if free.len() > z.len() {
            if null_space_step(&a, &free, &mut beta) {
                continue;
            }
            free.retain(|&k| beta[k] > 0.0);
        }
        let direction = loop {
            if free.is_empty() {
                return Err(Error::NumericalFailure(
                    "phase-amount problem has no free phase".into(),
                ));
            }
            let mut h: Vec<Vec<f64>> = free
                .iter()
                .map(|&k| {
                    free.iter()
                        .map(|&l| {
                            (0..z.len())
                                .map(|q| z[q] * a[k][q] * a[l][q] / (ev[q] * ev[q]))
                                .sum()
                        })
                        .collect()
                })
                .collect();
            let rhs: Vec<f64> = free.iter().map(|&k| -g[k]).collect();
            let solved = solve_dense(h.clone(), rhs.clone()).or_else(|| {
                // Identical phases make H singular; a tiny ridge keeps the step defined.
                let trace: f64 = (0..h.len()).map(|i| h[i][i]).sum();
                for (i, row) in h.iter_mut().enumerate() {
                    row[i] += 1e-12 * trace.max(1e-300);
                }
                solve_dense(h, rhs)
            });
            // Fall back to the projected gradient when Newton gives no descent.
            let gradient: Vec<f64> = free.iter().map(|&k| -g[k]).collect();
            let d = match solved {
                Some(d) if d.iter().zip(&gradient).map(|(d, m)| -d * m).sum::<f64>() < 0.0 => d,
                _ => gradient,
            };
            if let Some(position) = free
                .iter()
                .zip(&d)
                .position(|(&k, &dk)| beta[k] == 0.0 && dk < 0.0)
            {
                free.remove(position);
                continue;
            }
            let mut full = vec![0.0; np];
            for (&k, dk) in free.iter().zip(d) {
                full[k] = dk;
            }
            break full;
        };
        let mut step_max = f64::INFINITY;
        let mut blocking = None;
        for k in 0..np {
            if direction[k] < 0.0 {
                let s = -beta[k] / direction[k];
                if s < step_max {
                    step_max = s;
                    blocking = Some(k);
                }
            }
        }
        let q0 = objective(&beta);
        let slope: f64 = (0..np).map(|k| g[k] * direction[k]).sum();
        let mut alpha = step_max.min(1.0);
        let mut accepted = None;
        for _ in 0..60 {
            let mut candidate: Vec<f64> = (0..np)
                .map(|k| (beta[k] + alpha * direction[k]).max(0.0))
                .collect();
            if alpha == step_max
                && let Some(k) = blocking
            {
                candidate[k] = 0.0;
            }
            // Near the solution Q changes below its roundoff; a halved KKT residual
            // then accepts the step.
            if objective(&candidate) <= q0 + 1e-4 * alpha * slope
                || (kkt < 1e-6 && gradient_and_kkt(&candidate).1 < 0.5 * kkt)
            {
                accepted = Some(candidate);
                break;
            }
            alpha *= 0.5;
        }
        match accepted {
            Some(candidate) if candidate == beta && kkt <= 1e-10 => return Ok(beta),
            Some(candidate) => beta = candidate,
            // Converged to roundoff: no representable decrease remains.
            None if kkt <= 1e-11 => return Ok(beta),
            None => {
                return Err(Error::NumericalFailure(format!(
                    "phase-amount line search failed (KKT residual {kkt:e})"
                )));
            }
        }
    }
    Err(Error::NonConvergence(
        "multiphase Rachford-Rice exceeded 500 iterations".into(),
    ))
}

/// Moves β along a null vector v of the weights (Σ_k v_k a_kq = 0 over the free
/// phases) in the direction that does not increase Q = Σβ + const, until a phase
/// amount reaches zero. Returns false when no positive step is possible.
fn null_space_step(a: &[Vec<f64>], free: &[usize], beta: &mut [f64]) -> bool {
    let rows: Vec<Vec<f64>> = (0..a[0].len())
        .map(|q| free.iter().map(|&k| a[k][q]).collect())
        .collect();
    let Some(v) = null_vector(rows, free.len()) else {
        return false;
    };
    let s: f64 = v.iter().sum();
    let size = v.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let orientations: &[f64] = if s.abs() > 1e-12 * size {
        if s > 0.0 { &[-1.0] } else { &[1.0] }
    } else {
        &[1.0, -1.0]
    };
    for &sigma in orientations {
        let mut alpha = f64::INFINITY;
        let mut blocking = None;
        for (p, &k) in free.iter().enumerate() {
            let d = sigma * v[p];
            if d < 0.0 && beta[k] / -d < alpha {
                alpha = beta[k] / -d;
                blocking = Some(k);
            }
        }
        if alpha > 0.0 && alpha.is_finite() {
            for (p, &k) in free.iter().enumerate() {
                beta[k] = (beta[k] + alpha * sigma * v[p]).max(0.0);
            }
            if let Some(k) = blocking {
                beta[k] = 0.0;
            }
            return true;
        }
    }
    false
}

/// A nonzero solution of M v = 0 for `rows` × `cols` with cols > rank, by
/// Gaussian elimination with partial pivoting.
fn null_vector(mut m: Vec<Vec<f64>>, cols: usize) -> Option<Vec<f64>> {
    let scale = m.iter().flatten().map(|x| x.abs()).fold(0.0, f64::max);
    let mut pivots: Vec<usize> = Vec::new();
    let mut row = 0;
    for c in 0..cols {
        if row == m.len() {
            break;
        }
        let p = (row..m.len()).max_by(|&i, &j| m[i][c].abs().total_cmp(&m[j][c].abs()))?;
        if m[p][c].abs() <= 1e-13 * scale {
            continue;
        }
        m.swap(row, p);
        let pivot = m[row][c];
        for x in m[row].iter_mut() {
            *x /= pivot;
        }
        for i in 0..m.len() {
            if i != row {
                let f = m[i][c];
                let pivot_row = m[row].clone();
                for (x, y) in m[i].iter_mut().zip(pivot_row) {
                    *x -= f * y;
                }
            }
        }
        pivots.push(c);
        row += 1;
    }
    let free_column = (0..cols).find(|c| !pivots.contains(c))?;
    let mut v = vec![0.0; cols];
    v[free_column] = 1.0;
    for (r, &c) in pivots.iter().enumerate() {
        v[c] = -m[r][free_column];
    }
    v.iter().all(|x| x.is_finite()).then_some(v)
}

/// Unnormalised phase compositions x_kq = z_q a_kq / Σ_l β_l a_lq (active components).
fn compositions(z: &[f64], ln_a: &[Vec<f64>], beta: &[f64]) -> Vec<Vec<f64>> {
    let a = shifted(ln_a);
    let e: Vec<f64> = (0..z.len())
        .map(|q| (0..a.len()).map(|k| beta[k] * a[k][q]).sum())
        .collect();
    a.iter()
        .map(|row| (0..z.len()).map(|q| z[q] * row[q] / e[q]).collect())
        .collect()
}

#[derive(Debug, Clone)]
pub(crate) struct Iterate {
    pub(crate) kind: FluidKind,
    /// Full composition, model order.
    pub(crate) x: Vec<f64>,
    /// Full ln φ, model order.
    pub(crate) ln_phi: Vec<f64>,
}

/// max_q (max_k − min_k) ln(x_kq φ_kq).
fn fugacity_residual(phases: &[Iterate], act: &[usize]) -> f64 {
    act.iter()
        .map(|&i| {
            let mu = phases.iter().map(|p| p.x[i].ln() + p.ln_phi[i]);
            let (lo, hi) = mu.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), m| {
                (lo.min(m), hi.max(m))
            });
            hi - lo
        })
        .fold(0.0, f64::max)
}

/// Smallest composition distance between two phases of the same kind.
pub(crate) fn closest_pair(phases: &[Iterate]) -> Option<(usize, usize, f64)> {
    let mut best: Option<(usize, usize, f64)> = None;
    for k in 0..phases.len() {
        for l in k + 1..phases.len() {
            if phases[k].kind == phases[l].kind {
                let d = distance(&phases[k].x, &phases[l].x);
                if best.is_none_or(|(_, _, b)| d < b) {
                    best = Some((k, l, d));
                }
            }
        }
    }
    best
}

/// First pair of same-kind phases whose ln φ agree within 1e-12 on active components.
fn identical_coefficients(phases: &[Iterate], act: &[usize]) -> Option<(usize, usize)> {
    for k in 0..phases.len() {
        for l in k + 1..phases.len() {
            if phases[k].kind == phases[l].kind
                && act
                    .iter()
                    .all(|&i| (phases[k].ln_phi[i] - phases[l].ln_phi[i]).abs() <= 1e-12)
            {
                return Some((k, l));
            }
        }
    }
    None
}

#[derive(Default)]
struct Counters {
    iterations: usize,
    newton_steps: usize,
    merges: usize,
}

fn plait_note(phases: &[Iterate]) -> String {
    match closest_pair(phases) {
        Some((_, _, d)) if d < PLAIT_FAILURE_DISTANCE => format!(
            "; two phases are {d:.3e} apart in mole fraction: near a plait (liquid-liquid critical) or vapor-liquid critical point the split is not resolved"
        ),
        _ => String::new(),
    }
}

/// Phases and amounts from ln K relative to phase 0 (Newton variables).
fn from_ln_k<A: VlleAdapter>(
    a: &A,
    z_act: &[f64],
    act: &[usize],
    n: usize,
    kinds: &[FluidKind],
    theta: &[f64],
    beta: &[f64],
) -> Result<(Vec<f64>, Vec<Iterate>), Error> {
    let na = act.len();
    let ln_a: Vec<Vec<f64>> = (0..kinds.len())
        .map(|k| {
            if k == 0 {
                vec![0.0; na]
            } else {
                theta[(k - 1) * na..k * na].to_vec()
            }
        })
        .collect();
    let b = phase_amounts(z_act, &ln_a, beta)?;
    if b.iter().any(|v| v.is_nan() || *v <= 0.0) {
        return Err(Error::NumericalFailure(
            "Newton trial removes a phase".into(),
        ));
    }
    let xs = compositions(z_act, &ln_a, &b);
    let phases = kinds
        .iter()
        .zip(xs)
        .map(|(&kind, x)| {
            let x = expand(act, &x, n);
            Ok(Iterate {
                kind,
                ln_phi: a.ln_phi(&x, kind)?,
                x,
            })
        })
        .collect::<Result<Vec<_>, Error>>()?;
    Ok((b, phases))
}

/// One safeguarded Newton step on ln K_kq = ln φ_0q − ln φ_kq (k ≥ 1).
fn newton<A: VlleAdapter>(
    a: &A,
    z_act: &[f64],
    act: &[usize],
    n: usize,
    phases: &[Iterate],
    beta: &[f64],
) -> Option<(Vec<f64>, Vec<Iterate>)> {
    let kinds: Vec<FluidKind> = phases.iter().map(|p| p.kind).collect();
    let na = act.len();
    let theta: Vec<f64> = (1..phases.len())
        .flat_map(|k| {
            act.iter()
                .map(move |&i| phases[0].ln_phi[i] - phases[k].ln_phi[i])
        })
        .collect();
    let evaluate = |th: &[f64]| -> Result<(Vec<f64>, f64), Error> {
        let (b, ph) = from_ln_k(a, z_act, act, n, &kinds, th, beta)?;
        let mut residual = Vec::with_capacity(th.len());
        for k in 1..ph.len() {
            for (q, &i) in act.iter().enumerate() {
                residual.push(th[(k - 1) * na + q] + ph[k].ln_phi[i] - ph[0].ln_phi[i]);
            }
        }
        let g = ph
            .iter()
            .zip(&b)
            .map(|(p, b)| b * reduced_gibbs(&p.x, &p.ln_phi, act))
            .sum();
        Ok((residual, g))
    };
    let (r0, g0) = evaluate(&theta).ok()?;
    let candidate = newton_step(&theta, &r0, g0, evaluate)?;
    from_ln_k(a, z_act, act, n, &kinds, &candidate, beta).ok()
}

/// Successive substitution (with extrapolation and Newton) for a fixed set of
/// candidate phases; vanishing phases are removed and coinciding ones merged.
fn solve_split<A: VlleAdapter>(
    a: &A,
    z: &[f64],
    act: &[usize],
    phases: &mut Vec<Iterate>,
    beta: &mut Vec<f64>,
    options: MultiphaseOptions,
    counters: &mut Counters,
) -> Result<f64, Error> {
    let n = z.len();
    let z_act: Vec<f64> = act.iter().map(|&i| z[i]).collect();
    let mut history: Vec<Vec<f64>> = Vec::new();
    for iteration in 1..=options.max_iterations {
        counters.iterations += 1;
        // Same-kind phases with identical fugacity coefficients (for example two
        // ideal-gas vapors) make the phase-amount problem degenerate: merge them.
        if let Some((k, l)) = identical_coefficients(phases, act) {
            beta[k] += beta[l];
            phases.remove(l);
            beta.remove(l);
            history.clear();
            if phases.len() < 2 {
                return Err(Error::NumericalFailure(
                    "multiphase split reduced to phases with identical fugacity coefficients"
                        .into(),
                ));
            }
        }
        let ln_a: Vec<Vec<f64>> = phases
            .iter()
            .map(|p| act.iter().map(|&i| -p.ln_phi[i]).collect())
            .collect();
        let b = phase_amounts(&z_act, &ln_a, beta)?;
        let xs = compositions(&z_act, &ln_a, &b);
        let mut next_phases = Vec::new();
        let mut next_beta = Vec::new();
        for ((p, x), bk) in phases.iter().zip(xs).zip(b) {
            if bk > 0.0 {
                let x = expand(act, &x, n);
                next_phases.push(Iterate {
                    kind: p.kind,
                    ln_phi: a.ln_phi(&x, p.kind)?,
                    x,
                });
                next_beta.push(bk);
            }
        }
        if next_phases.len() != phases.len() {
            history.clear();
        }
        *phases = next_phases;
        *beta = next_beta;
        if phases.len() < 2 {
            return Err(Error::NumericalFailure(
                "multiphase split collapsed to one phase although a stability trial found a negative tangent-plane distance".into(),
            ));
        }
        let residual = fugacity_residual(phases, act);
        if !residual.is_finite() {
            return Err(Error::NumericalFailure(
                "nonfinite multiphase fugacity residual".into(),
            ));
        }
        if let Some((k, l, d)) = closest_pair(phases)
            && d < MERGE_DISTANCE
        {
            beta[k] += beta[l];
            phases.remove(l);
            beta.remove(l);
            counters.merges += 1;
            history.clear();
            if phases.len() < 2 {
                return Err(Error::NumericalFailure(
                    "VLLE split phases merged into one: the stability trials and the split disagree near a plait (liquid-liquid critical) or critical point".into(),
                ));
            }
            continue;
        }
        if residual <= options.tolerance {
            return Ok(residual);
        }
        if options.newton
            && iteration >= 3
            && residual < 1e-2
            && let Some((b, ph)) = newton(a, &z_act, act, n, phases, beta)
        {
            *beta = b;
            *phases = ph;
            counters.newton_steps += 1;
            history.clear();
            continue;
        }
        history.push(
            phases
                .iter()
                .flat_map(|p| act.iter().map(|&i| p.ln_phi[i]))
                .collect(),
        );
        if history.len() >= 3 && iteration % 5 == 0 {
            if let Some(extrapolated) = dominant_eigenvalue_step(&history) {
                for (k, p) in phases.iter_mut().enumerate() {
                    for (q, &i) in act.iter().enumerate() {
                        p.ln_phi[i] = extrapolated[k * act.len() + q];
                    }
                }
            }
            history.clear();
        }
    }
    Err(Error::NonConvergence(format!(
        "multiphase split exhausted {} iterations{}",
        options.max_iterations,
        plait_note(phases)
    )))
}

/// Multiphase TP flash of composition `z` with the adapter's model.
pub(crate) fn tp_flash<A: VlleAdapter>(
    a: &A,
    z: &[f64],
    options: MultiphaseOptions,
) -> Result<MultiphaseTpFlash<A::State>, Error> {
    let n = z.len();
    let act = active(z);
    let (t, p) = (a.temperature().value(), a.pressure().value());
    let (reference_kind, reference_ln_phi, g_feed) = reference(a, z, &act)?;
    let feed_stability = stability_with(a, z, reference_kind, &reference_ln_phi, options)?;
    let mut diagnostics = vec![diagnostic(
        "VLLE-STABILITY-FINITE-TRIALS",
        format!(
            "Multiphase stability from {} TPD trials (Wilson vapor/liquid, aqueous, hydrocarbon, near-pure components); finite trials do not prove global stability.",
            feed_stability.trials.len()
        ),
        DiagnosticSeverity::Information,
    )];
    if feed_stability.stable_in_trials {
        let state = a.state(z, reference_kind)?;
        let phase = a.label(z, reference_kind, &state)?;
        diagnostics.push(diagnostic(
            "VLLE-SINGLE-PHASE-LABEL",
            "Single-phase label from the model's documented heuristic (cubic) or phase kind (gamma-phi); a label, not a phase-boundary result.",
            DiagnosticSeverity::Information,
        ));
        return Ok(MultiphaseTpFlash {
            temperature_k: t,
            pressure_pa: p,
            phases: vec![MultiphasePhase {
                phase,
                role: if phase == Phase::Liquid {
                    PhaseRole::Liquid
                } else {
                    PhaseRole::Vapor
                },
                kind: reference_kind,
                phase_fraction: 1.0,
                composition: z.to_vec(),
                state,
            }],
            feed_stability,
            final_stability: None,
            iterations: 0,
            max_ln_fugacity_residual: 0.0,
            max_component_balance_residual: 0.0,
            reduced_gibbs_change: 0.0,
            diagnostics,
        });
    }
    // Initial phases: distinct negative-TPD stationary points, then the feed.
    let mut phases: Vec<Iterate> = Vec::new();
    for trial in feed_stability.unstable_trials(options.tm_threshold) {
        if phases.len() == MAX_PHASES {
            break;
        }
        if phases
            .iter()
            .any(|p| p.kind == trial.kind && distance(&p.x, &trial.composition) < DISTINCT_DISTANCE)
        {
            continue;
        }
        phases.push(Iterate {
            kind: trial.kind,
            ln_phi: a.ln_phi(&trial.composition, trial.kind)?,
            x: trial.composition.clone(),
        });
    }
    if phases.len() < MAX_PHASES
        && !phases
            .iter()
            .any(|p| p.kind == reference_kind && distance(&p.x, z) < DISTINCT_DISTANCE)
    {
        phases.push(Iterate {
            kind: reference_kind,
            x: z.to_vec(),
            ln_phi: reference_ln_phi.clone(),
        });
    }
    let mut beta = vec![1.0 / phases.len() as f64; phases.len()];
    let mut counters = Counters::default();
    let mut rounds = 0;
    let (residual, final_stability) = loop {
        let residual = solve_split(a, z, &act, &mut phases, &mut beta, options, &mut counters)?;
        let reference = (0..phases.len())
            .max_by(|&k, &l| beta[k].total_cmp(&beta[l]))
            .expect("at least two phases");
        let check = stability_with(
            a,
            &phases[reference].x,
            phases[reference].kind,
            &phases[reference].ln_phi,
            options,
        )?;
        let added = check
            .unstable_trials(options.tm_threshold)
            .into_iter()
            .find(|t| {
                !phases
                    .iter()
                    .any(|p| p.kind == t.kind && distance(&p.x, &t.composition) < DISTINCT_DISTANCE)
            })
            .cloned();
        let Some(added) = added else {
            break (residual, check);
        };
        if phases.len() >= MAX_PHASES {
            return Err(Error::UnsupportedCapability(
                "post-split stability indicates a fourth phase; at most three phases (VLLE) are supported".into(),
            ));
        }
        rounds += 1;
        if rounds > MAX_ROUNDS {
            return Err(Error::NumericalFailure(format!(
                "VLLE phase set did not settle after {MAX_ROUNDS} phase additions ({} merges){}",
                counters.merges,
                plait_note(&phases)
            )));
        }
        phases.push(Iterate {
            kind: added.kind,
            ln_phi: a.ln_phi(&added.composition, added.kind)?,
            x: added.composition,
        });
        beta.push(0.0);
    };
    let balance = (0..n)
        .map(|i| {
            (phases
                .iter()
                .zip(&beta)
                .map(|(p, b)| b * p.x[i])
                .sum::<f64>()
                - z[i])
                .abs()
        })
        .fold(0.0, f64::max);
    let g_split: f64 = phases
        .iter()
        .zip(&beta)
        .map(|(p, b)| b * reduced_gibbs(&p.x, &p.ln_phi, &act))
        .sum();
    let reduced_gibbs_change = g_split - g_feed;
    if reduced_gibbs_change >= 0.0 {
        return Err(Error::NumericalFailure(format!(
            "converged multiphase split does not lower the Gibbs energy (Δg/RT = {reduced_gibbs_change:e})"
        )));
    }
    if let Some((_, _, d)) = closest_pair(&phases)
        && d < PLAIT_WARNING_DISTANCE
    {
        diagnostics.push(diagnostic(
            "VLLE-NEAR-PLAIT-POINT",
            format!("Two phases are {d:.3e} apart in mole fraction; the split is close to a plait (liquid-liquid critical) or critical point and sensitive to T, p and composition."),
            DiagnosticSeverity::Warning,
        ));
    }
    diagnostics.push(diagnostic(
        "VLLE-SPLIT",
        format!(
            "{} phases after {} split iterations, {} accepted Newton steps, {} merges and {} phase additions.",
            phases.len(),
            counters.iterations,
            counters.newton_steps,
            counters.merges,
            rounds
        ),
        DiagnosticSeverity::Information,
    ));
    diagnostics.push(diagnostic(
        "VLLE-PHASE-LABELS",
        "Phases are ordered vapor first, then liquids by increasing mass density (light, heavy). Cubic: the largest-volume phase is vapor in a three-phase split, and in a two-phase split unless the single-phase heuristic labels it liquid.",
        DiagnosticSeverity::Information,
    ));
    let phases = finish(a, phases, beta)?;
    Ok(MultiphaseTpFlash {
        temperature_k: t,
        pressure_pa: p,
        phases,
        feed_stability,
        final_stability: Some(final_stability),
        iterations: counters.iterations,
        max_ln_fugacity_residual: residual,
        max_component_balance_residual: balance,
        reduced_gibbs_change,
        diagnostics,
    })
}

/// Full states, labels, roles and the reported phase order.
pub(crate) fn finish<A: VlleAdapter>(
    a: &A,
    phases: Vec<Iterate>,
    beta: Vec<f64>,
) -> Result<Vec<MultiphasePhase<A::State>>, Error> {
    let mut out: Vec<MultiphasePhase<A::State>> = Vec::new();
    for (p, b) in phases.into_iter().zip(beta) {
        let state = a.state(&p.x, p.kind)?;
        out.push(MultiphasePhase {
            phase: a.label(&p.x, p.kind, &state)?,
            role: PhaseRole::Liquid,
            kind: p.kind,
            phase_fraction: b,
            composition: p.x,
            state,
        });
    }
    Ok(order::<A>(out))
}

/// Labels, roles and the reported order (vapor first, then liquids by increasing
/// mass density) of phases with their states.
pub(crate) fn order<A: VlleAdapter>(
    mut out: Vec<MultiphasePhase<A::State>>,
) -> Vec<MultiphasePhase<A::State>> {
    out.sort_by(|p, q| A::molar_volume(&q.state).total_cmp(&A::molar_volume(&p.state)));
    if out.iter().all(|p| p.kind == FluidKind::Fluid) {
        let count = out.len();
        for (k, p) in out.iter_mut().enumerate() {
            p.phase = if k == 0 && (count == 3 || p.phase != Phase::Liquid) {
                Phase::Vapor
            } else {
                Phase::Liquid
            };
        }
    }
    out.sort_by(|p, q| {
        let vapor = |x: &MultiphasePhase<A::State>| x.phase != Phase::Liquid;
        vapor(q)
            .cmp(&vapor(p))
            .then(A::mass_density(&p.state).total_cmp(&A::mass_density(&q.state)))
    });
    let liquids = out.iter().filter(|p| p.phase == Phase::Liquid).count();
    let mut seen = 0;
    for p in &mut out {
        p.role = if p.phase != Phase::Liquid {
            PhaseRole::Vapor
        } else {
            seen += 1;
            match (liquids, seen) {
                (1, _) => PhaseRole::Liquid,
                (_, 1) => PhaseRole::LightLiquid,
                (l, s) if s == l => PhaseRole::HeavyLiquid,
                _ => PhaseRole::Liquid,
            }
        };
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_amounts_reproduce_a_two_phase_rachford_rice_split() {
        // K = (2, 0.5) relative to phase 0 for z = (0.5, 0.5): beta_1 = 0.5.
        let z = [0.5, 0.5];
        let ln_a = vec![vec![0.0, 0.0], vec![2f64.ln(), 0.5f64.ln()]];
        let b = phase_amounts(&z, &ln_a, &[0.5, 0.5]).unwrap();
        assert!(
            (b[0] - 0.5).abs() < 1e-12 && (b[1] - 0.5).abs() < 1e-12,
            "{b:?}"
        );
        let x = compositions(&z, &ln_a, &b);
        for row in &x {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn phase_amounts_drop_a_phase_that_cannot_exist() {
        // Third phase with K = 1e-3 everywhere has sum x < 1: its amount is zero.
        let z = [0.5, 0.5];
        let ln_a = vec![
            vec![0.0, 0.0],
            vec![2f64.ln(), 0.5f64.ln()],
            vec![1e-3f64.ln(), 1e-3f64.ln()],
        ];
        let b = phase_amounts(&z, &ln_a, &[0.3, 0.3, 0.4]).unwrap();
        assert_eq!(b[2], 0.0);
        assert!((b[0] + b[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn phase_amounts_give_a_single_phase_outside_the_two_phase_region() {
        let z = [0.9, 0.1];
        // K = (1.05, 0.2): z sum K > 1? 0.945+0.02 < 1 -> phase 1 absent.
        let ln_a = vec![vec![0.0, 0.0], vec![1.05f64.ln(), 0.2f64.ln()]];
        let b = phase_amounts(&z, &ln_a, &[0.5, 0.5]).unwrap();
        assert_eq!(b[1], 0.0);
        assert!((b[0] - 1.0).abs() < 1e-12);
    }
}
