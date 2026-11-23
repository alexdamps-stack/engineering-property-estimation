//! Small dense numerics for the reactors: linear solves, null spaces,
//! constrained Gibbs minimization, bracketed root finding and an adaptive
//! Runge–Kutta integrator. Sizes are the number of package components
//! (tens at most), so dense O(n^3) algebra is adequate.

/// Solves `a x = b` by Gaussian elimination with partial pivoting.
pub(crate) fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    let scale = a
        .iter()
        .flatten()
        .fold(0.0_f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs().partial_cmp(&(1e-14 * scale)) != Some(std::cmp::Ordering::Greater) {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let f = a[row][col] / a[col][col];
            if f != 0.0 {
                let (upper, lower) = a.split_at_mut(row);
                for (target, source) in lower[0][col..].iter_mut().zip(&upper[col][col..]) {
                    *target -= f * source;
                }
                b[row] -= f * b[col];
            }
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - s) / a[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Reduced row echelon form; returns the pivot columns.
fn rref(m: &mut [Vec<f64>], tol: f64) -> Vec<usize> {
    let rows = m.len();
    let cols = m.first().map_or(0, Vec::len);
    let mut pivots = Vec::new();
    let mut r = 0;
    for c in 0..cols {
        if r == rows {
            break;
        }
        let Some(p) = (r..rows).max_by(|&i, &j| m[i][c].abs().total_cmp(&m[j][c].abs())) else {
            break;
        };
        if m[p][c].abs() <= tol {
            continue;
        }
        m.swap(r, p);
        let d = m[r][c];
        for v in m[r].iter_mut() {
            *v /= d;
        }
        for i in 0..rows {
            if i != r && m[i][c] != 0.0 {
                let f = m[i][c];
                let pivot_row = m[r].clone();
                for (target, source) in m[i].iter_mut().zip(&pivot_row) {
                    *target -= f * source;
                }
            }
        }
        pivots.push(c);
        r += 1;
    }
    pivots
}

/// Rank of a matrix given as rows.
pub(crate) fn rank(rows: &[Vec<f64>]) -> usize {
    let scale = rows.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
    rref(&mut rows.to_vec(), 1e-10 * scale.max(1.0)).len()
}

/// Basis of the null space of `rows` (each row has `cols` entries).
pub(crate) fn null_space(rows: &[Vec<f64>], cols: usize) -> Vec<Vec<f64>> {
    let mut m = rows.to_vec();
    let scale = m.iter().flatten().fold(0.0_f64, |s, v| s.max(v.abs()));
    let pivots = rref(&mut m, 1e-10 * scale.max(1.0));
    (0..cols)
        .filter(|c| !pivots.contains(c))
        .map(|free| {
            let mut v = vec![0.0; cols];
            v[free] = 1.0;
            for (r, &p) in pivots.iter().enumerate() {
                v[p] = -m[r][free];
            }
            v
        })
        .collect()
}

/// Linearly independent subset of `rows` (Gram–Schmidt, in order).
fn independent_rows(rows: &[Vec<f64>]) -> Vec<usize> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut keep = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let norm = row.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let mut v = row.clone();
        for q in &basis {
            let d: f64 = v.iter().zip(q).map(|(a, b)| a * b).sum();
            for (a, b) in v.iter_mut().zip(q) {
                *a -= d * b;
            }
        }
        let rest = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if rest > 1e-10 * norm {
            for x in &mut v {
                *x /= rest;
            }
            basis.push(v);
            keep.push(i);
        }
    }
    keep
}

/// Trace amount added to every species before minimization, relative to the
/// total feed. It keeps the feasible set's interior nonempty when a species
/// can neither form nor be consumed (for example a reaction with neither
/// reactants nor products present); such species leave at this trace level.
pub(crate) const TRACE: f64 = 1e-14;

/// Minimizes `G/RT = sum_i n_i (c_i + ln(n_i/N) + ln phi_i(n))` subject to
/// `B n = B n0` for the invariant rows `B`, by Newton iteration on the
/// element potentials (one per independent invariant) and ln N, with the
/// fugacity coefficients updated by successive substitution.
///
/// `c_i` is the standard chemical potential over RT including `ln(p/p0)`.
pub(crate) fn minimize_gibbs(
    c: &[f64],
    invariants: &[Vec<f64>],
    n0: &[f64],
    mut ln_phi: impl FnMut(&[f64]) -> Result<Vec<f64>, String>,
) -> Result<Vec<f64>, String> {
    let s = c.len();
    let total: f64 = n0.iter().sum();
    if !(total.is_finite() && total > 0.0) {
        return Err("equilibrium requires a positive finite feed".into());
    }
    let start: Vec<f64> = n0.iter().map(|n| n + TRACE * total).collect();
    let keep = independent_rows(invariants);
    let b_rows: Vec<&Vec<f64>> = keep.iter().map(|&k| &invariants[k]).collect();
    let b: Vec<f64> = b_rows
        .iter()
        .map(|row| row.iter().zip(&start).map(|(a, n)| a * n).sum())
        .collect();
    let mut phi = vec![0.0; s];
    let mut n = start.clone();
    for _ in 0..200 {
        let cc: Vec<f64> = c.iter().zip(&phi).map(|(c, f)| c + f).collect();
        n = dual_newton(&cc, &b_rows, &b, &n)?;
        let x: Vec<f64> = {
            let sum: f64 = n.iter().sum();
            n.iter().map(|v| v / sum).collect()
        };
        let next = ln_phi(&x)?;
        let change = next
            .iter()
            .zip(&phi)
            .fold(0.0_f64, |d, (a, b)| d.max((a - b).abs()));
        phi = next;
        if change < 1e-11 {
            let cc: Vec<f64> = c.iter().zip(&phi).map(|(c, f)| c + f).collect();
            return dual_newton(&cc, &b_rows, &b, &n);
        }
    }
    Err("equilibrium fugacity-coefficient iteration did not converge in 200 passes".into())
}

/// Ideal-mixture minimization for fixed `c` by Newton on (lambda, ln N).
fn dual_newton(
    c: &[f64],
    rows: &[&Vec<f64>],
    b: &[f64],
    guess: &[f64],
) -> Result<Vec<f64>, String> {
    let s = c.len();
    let m = rows.len();
    let total_guess: f64 = guess.iter().sum();
    // Least-squares start from the guess: ln x_i + c_i = sum_k lambda_k B_ki.
    let target: Vec<f64> = (0..s)
        .map(|i| (guess[i] / total_guess).ln() + c[i])
        .collect();
    let normal: Vec<Vec<f64>> = (0..m)
        .map(|k| {
            (0..m)
                .map(|l| (0..s).map(|i| rows[k][i] * rows[l][i]).sum())
                .collect()
        })
        .collect();
    let rhs: Vec<f64> = (0..m)
        .map(|k| (0..s).map(|i| rows[k][i] * target[i]).sum())
        .collect();
    let mut lambda = if m == 0 {
        Vec::new()
    } else {
        solve_dense(normal, rhs).ok_or("singular equilibrium invariants")?
    };
    let mut ln_n = total_guess.ln();
    let amounts = |lambda: &[f64], ln_n: f64| -> Vec<f64> {
        (0..s)
            .map(|i| {
                let e = ln_n - c[i] + (0..m).map(|k| lambda[k] * rows[k][i]).sum::<f64>();
                e.min(700.0).exp()
            })
            .collect()
    };
    let residual = |n: &[f64], ln_n: f64| -> (Vec<f64>, f64) {
        let mut f = Vec::with_capacity(m + 1);
        let mut worst = 0.0_f64;
        for k in 0..m {
            let sum: f64 = (0..s).map(|i| rows[k][i] * n[i]).sum();
            let scale: f64 = (0..s).map(|i| (rows[k][i] * n[i]).abs()).sum::<f64>() + b[k].abs();
            f.push(sum - b[k]);
            worst = worst.max((sum - b[k]).abs() / scale.max(f64::MIN_POSITIVE));
        }
        let big_n = ln_n.exp();
        let sum: f64 = n.iter().sum();
        f.push(sum - big_n);
        worst = worst.max((sum - big_n).abs() / big_n);
        (f, worst)
    };
    let mut n = amounts(&lambda, ln_n);
    let (mut f, mut worst) = residual(&n, ln_n);
    for _ in 0..200 {
        if worst < 1e-13 {
            return Ok(n);
        }
        let mut jac = vec![vec![0.0; m + 1]; m + 1];
        for k in 0..m {
            for l in 0..m {
                jac[k][l] = (0..s).map(|i| rows[k][i] * rows[l][i] * n[i]).sum();
            }
            let col: f64 = (0..s).map(|i| rows[k][i] * n[i]).sum();
            jac[k][m] = col;
            jac[m][k] = col;
        }
        jac[m][m] = n.iter().sum::<f64>() - ln_n.exp();
        // Row equilibration: an invariant carried only by trace species has
        // tiny entries that would otherwise look singular.
        let mut rhs: Vec<f64> = f.iter().map(|v| -v).collect();
        for (row, r) in jac.iter_mut().zip(rhs.iter_mut()) {
            let scale = row.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
            if scale > 0.0 {
                for v in row.iter_mut() {
                    *v /= scale;
                }
                *r /= scale;
            }
        }
        let step = solve_dense(jac, rhs).ok_or("singular equilibrium Newton matrix")?;
        let largest = step.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        let mut alpha = if largest > 5.0 { 5.0 / largest } else { 1.0 };
        let mut accepted = false;
        for _ in 0..60 {
            let lam: Vec<f64> = (0..m).map(|k| lambda[k] + alpha * step[k]).collect();
            let ln_try = ln_n + alpha * step[m];
            let n_try = amounts(&lam, ln_try);
            let (f_try, w_try) = residual(&n_try, ln_try);
            let norm = |f: &[f64]| f.iter().map(|v| v * v).sum::<f64>();
            if w_try.is_finite() && (norm(&f_try) < norm(&f) || w_try < worst) {
                lambda = lam;
                ln_n = ln_try;
                n = n_try;
                f = f_try;
                worst = w_try;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    if worst < 1e-10 {
        return Ok(n);
    }
    Err(format!(
        "chemical equilibrium did not converge (scaled constraint residual {worst:.3e})"
    ))
}

/// Brent's method on `[lo, hi]`; `f(lo)` and `f(hi)` must differ in sign.
pub(crate) fn brent(
    mut f: impl FnMut(f64) -> Result<f64, String>,
    lo: f64,
    hi: f64,
    xtol: f64,
) -> Result<f64, String> {
    let (mut a, mut b) = (lo, hi);
    let (mut fa, mut fb) = (f(a)?, f(b)?);
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    if fa.signum() == fb.signum() {
        return Err(format!(
            "no solution between {lo:.2} K and {hi:.2} K (energy residual {fa:.4e} and {fb:.4e} W have the same sign)"
        ));
    }
    let (mut c, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    for _ in 0..200 {
        if fb.signum() == fc.signum() {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol = 2.0 * f64::EPSILON * b.abs() + 0.5 * xtol;
        let mid = 0.5 * (c - b);
        if mid.abs() <= tol || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * mid * s;
                q = 1.0 - s;
            } else {
                let qa = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * mid * qa * (qa - r) - (b - a) * (r - 1.0));
                q = (qa - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            if 2.0 * p < (3.0 * mid * q - (tol * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = mid;
                e = d;
            }
        } else {
            d = mid;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol { d } else { tol.copysign(mid) };
        fb = f(b)?;
    }
    Err("temperature iteration did not converge".into())
}

/// Root of `f` near `start`: expands a bracket outward in doubling steps
/// (10, 20, 40, ... K) within `[lo, hi]` until the sign changes, then applies
/// [`brent`]. A side whose evaluation fails (for example outside a property
/// domain) stops expanding. Evaluating near the start avoids the extremes of
/// the temperature interval, where equilibria and flashes are hardest.
pub(crate) fn bracketed_root(
    mut f: impl FnMut(f64) -> Result<f64, String>,
    start: f64,
    lo: f64,
    hi: f64,
    xtol: f64,
) -> Result<f64, String> {
    let start = start.clamp(lo, hi);
    let f0 = f(start)?;
    if f0 == 0.0 {
        return Ok(start);
    }
    let (mut below, mut above) = (Some(start), Some(start));
    let mut step = 10.0;
    let mut last_error = None;
    while below.is_some() || above.is_some() {
        for side in [-1.0, 1.0] {
            let current = if side < 0.0 { &mut below } else { &mut above };
            let Some(from) = *current else { continue };
            let to = (start + side * step).clamp(lo, hi);
            if to == from {
                *current = None;
                continue;
            }
            match f(to) {
                Ok(v) if v.signum() != f0.signum() || v == 0.0 => {
                    let (a, b) = if side < 0.0 { (to, from) } else { (from, to) };
                    return brent(&mut f, a, b, xtol);
                }
                Ok(_) => *current = Some(to),
                Err(e) => {
                    last_error = Some(e);
                    *current = None;
                }
            }
        }
        step *= 2.0;
    }
    Err(match last_error {
        Some(e) => format!("no solution found between {lo:.2} K and {hi:.2} K ({e})"),
        // The residual is outlet minus required enthalpy: negative throughout
        // means even the upper bound is too cold.
        None if f0 < 0.0 => format!(
            "the outlet temperature would be above the package interval ({hi:.2} K); specify an outlet temperature or a cooling duty"
        ),
        None => format!(
            "the outlet temperature would be below the package interval ({lo:.2} K); specify an outlet temperature or a heating duty"
        ),
    })
}

/// Adaptive Dormand–Prince 5(4) integration of `dy/dx = f(x, y)` from `x0`
/// to `x1`. The error per step is held below `rtol |y| + atol` componentwise.
/// Returns the final state and the accepted points `(x, y)`.
pub(crate) fn integrate(
    mut f: impl FnMut(f64, &[f64]) -> Result<Vec<f64>, String>,
    x0: f64,
    x1: f64,
    y0: Vec<f64>,
    rtol: f64,
    atol: &[f64],
    max_steps: usize,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    const C: [f64; 7] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
    const A: [[f64; 6]; 7] = [
        [0.0; 6],
        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
        [
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
            0.0,
            0.0,
        ],
        [
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
            0.0,
        ],
        [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    const B5: [f64; 7] = [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
        0.0,
    ];
    const B4: [f64; 7] = [
        5179.0 / 57600.0,
        0.0,
        7571.0 / 16695.0,
        393.0 / 640.0,
        -92097.0 / 339200.0,
        187.0 / 2100.0,
        1.0 / 40.0,
    ];
    let span = x1 - x0;
    if !(span.is_finite() && span > 0.0) {
        return Err("integration interval must be positive and finite".into());
    }
    let n = y0.len();
    let mut x = x0;
    let mut y = y0;
    let mut points = vec![(x, y.clone())];
    let mut h = span / 100.0;
    let mut k0 = f(x, &y)?;
    for _ in 0..max_steps {
        if x >= x1 {
            return Ok(points);
        }
        h = h.min(x1 - x);
        let mut k = vec![k0.clone()];
        let mut failed = false;
        for stage in 1..7 {
            let ys: Vec<f64> = (0..n)
                .map(|i| y[i] + h * (0..stage).map(|j| A[stage][j] * k[j][i]).sum::<f64>())
                .collect();
            match f(x + C[stage] * h, &ys) {
                Ok(v) => k.push(v),
                Err(_) => {
                    failed = true;
                    break;
                }
            }
        }
        if failed {
            // A trial state outside the property domain: retry smaller.
            h *= 0.25;
            if h < span * 1e-14 {
                return Err("integration step underflow at an invalid state".into());
            }
            continue;
        }
        let y5: Vec<f64> = (0..n)
            .map(|i| y[i] + h * (0..7).map(|j| B5[j] * k[j][i]).sum::<f64>())
            .collect();
        let err = (0..n)
            .map(|i| {
                let e = h * (0..7).map(|j| (B5[j] - B4[j]) * k[j][i]).sum::<f64>();
                e.abs() / (atol[i] + rtol * y[i].abs().max(y5[i].abs()))
            })
            .fold(0.0_f64, f64::max);
        if err.is_finite() && err <= 1.0 {
            x += h;
            if x1 - x < 1e-12 * span {
                x = x1;
            }
            y = y5;
            k0 = k.pop().expect("seven stages");
            points.push((x, y.clone()));
        }
        let factor = if err.is_finite() && err > 0.0 {
            (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
        } else if err == 0.0 {
            5.0
        } else {
            0.2
        };
        h *= factor;
        if h < span * 1e-14 {
            return Err("integration step underflow (stiff or discontinuous rates)".into());
        }
    }
    Err(format!(
        "integration needed more than {max_steps} steps (stiff kinetics are not supported by the explicit integrator)"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_solve_and_null_space() {
        let x = solve_dense(vec![vec![2.0, 1.0], vec![1.0, 3.0]], vec![3.0, 5.0]).unwrap();
        assert!((x[0] - 0.8).abs() < 1e-14 && (x[1] - 1.4).abs() < 1e-14);
        // nu^T of CO + H2O = CO2 + H2 over (CO, H2O, CO2, H2).
        let ns = null_space(&[vec![-1.0, -1.0, 1.0, 1.0]], 4);
        assert_eq!(ns.len(), 3);
        for v in &ns {
            assert!((-v[0] - v[1] + v[2] + v[3]).abs() < 1e-14);
        }
        assert_eq!(rank(&[vec![1.0, 2.0], vec![2.0, 4.0]]), 1);
    }

    #[test]
    fn brent_finds_a_bracketed_root() {
        let r = brent(|x| Ok(x * x - 2.0), 0.0, 2.0, 1e-13).unwrap();
        assert!((r - 2f64.sqrt()).abs() < 1e-12);
        assert!(brent(|x| Ok(x * x + 1.0), 0.0, 2.0, 1e-12).is_err());
    }

    #[test]
    fn integrator_reproduces_exponential_decay() {
        let pts = integrate(
            |_, y| Ok(vec![-2.0 * y[0]]),
            0.0,
            1.0,
            vec![1.0],
            1e-10,
            &[1e-14],
            10_000,
        )
        .unwrap();
        let (x, y) = pts.last().unwrap();
        assert_eq!(*x, 1.0);
        assert!((y[0] - (-2.0f64).exp()).abs() < 1e-9);
    }

    #[test]
    fn ideal_isomerization_equilibrium() {
        // A <-> B with c_B - c_A = -ln 3: x_B / x_A = 3 at equilibrium.
        let n = minimize_gibbs(
            &[0.0, -(3.0f64).ln()],
            &[vec![1.0, 1.0]],
            &[1.0, 0.0],
            |x| Ok(vec![0.0; x.len()]),
        )
        .unwrap();
        assert!((n[1] / n[0] - 3.0).abs() < 1e-11, "{n:?}");
        assert!((n[0] + n[1] - (1.0 + 2.0 * TRACE)).abs() < 1e-12);
    }
}
