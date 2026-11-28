//! Product distillation curves and cut-point specifications of the
//! petroleum fractionator (#190).
//!
//! Each product's TBP curve is built from its composition by
//! [`thermo_data::TbpCurve`] from the normal boiling points and standard
//! liquid volumes (60 °F) of the package's pseudo-components; ASTM D86 follows
//! by the inverse API (Riazi–Daubert) conversion. The curves are on a
//! pseudo-component basis: database components (water from stripping steam,
//! light ends) have no boiling-point characterisation in the package and are
//! left out, and their mole fraction is reported with the curve.
//!
//! A side product may be specified by a cut point, the TBP or D86
//! temperature at a liquid volume percent. Its molar flow then becomes the
//! unknown of an outer quasi-Newton (Broyden) iteration on the logarithms of
//! the cut-point products' molar flows, each evaluation a converged column
//! (including any standard-volume products). The Jacobian starts from
//! forward differences. A cut point that the column cannot reach fails
//! explicitly: a target outside the pseudo-components' TBP range, a column
//! with no solution on the way, a residual that no longer decreases, or no
//! convergence within the iteration limit.

use thermo_data::{CurveComponent, D86_MAX_VOLUME_PERCENT, REPORTED_VOLUME_PERCENTS, TbpCurve};

use super::petroleum::{DrawRate, PetroleumFractionatorSpec, standard_volumes};
use super::{DistillationColumnResult, DistillationCondenser};
use crate::MaterialStream;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};
use crate::units::{CalcContext, UnitOutput};

/// Cut-point residual at which the specification is met, K.
pub const CUT_POINT_TOLERANCE_K: f64 = 1e-6;
const CUT_POINT_ITERATIONS: usize = 25;
/// Forward-difference step in the logarithm of the molar flow.
const DIFFERENCE_STEP: f64 = 1e-4;
/// Largest change of the logarithm of a molar flow in one step.
const MAX_LOG_STEP: f64 = 0.5;
const STEP_HALVINGS: usize = 4;

/// Distillation curve of a cut-point specification.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CurveBasis {
    /// True boiling point curve.
    #[default]
    Tbp,
    /// ASTM D86 curve converted from the TBP curve (0 to 95 vol%).
    AstmD86,
}

/// Cut-point specification of a side product: the temperature of its TBP or
/// D86 curve at a liquid volume percent distilled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutPoint {
    pub basis: CurveBasis,
    /// Liquid volume percent distilled: 0 to 100 (TBP) or 0 to 95 (D86).
    pub volume_percent: f64,
    /// Specified temperature, K.
    pub temperature_k: Option<f64>,
    /// Starting molar flow of the product, mol/s. Without it the iteration
    /// starts from the hydrocarbon feed flow divided by twice the number of
    /// side products plus one.
    pub estimate_mol_s: Option<f64>,
}

impl Default for CutPoint {
    fn default() -> Self {
        Self {
            basis: CurveBasis::Tbp,
            volume_percent: 95.0,
            temperature_k: None,
            estimate_mol_s: None,
        }
    }
}

impl CutPoint {
    pub(super) fn issues(&self, out: &mut Vec<SpecIssue>) {
        let max = match self.basis {
            CurveBasis::Tbp => 100.0,
            CurveBasis::AstmD86 => D86_MAX_VOLUME_PERCENT,
        };
        if !(0.0..=max).contains(&self.volume_percent) {
            out.push(SpecIssue {
                field: "side-product cut-point volume percent",
                problem: SpecProblem::OutOfRange,
            });
        }
        check_optional(
            out,
            "side-product cut-point temperature",
            self.temperature_k,
            Sign::Positive,
        );
        if let Some(estimate) = self.estimate_mol_s
            && !(estimate.is_finite() && estimate > 0.0)
        {
            out.push(SpecIssue {
                field: "side-product cut-point flow estimate",
                problem: SpecProblem::OutOfRange,
            });
        }
    }
}

/// TBP and D86 points of one product at [`REPORTED_VOLUME_PERCENTS`].
#[derive(Debug, Clone, PartialEq)]
pub struct ProductCurve {
    /// Outlet port label (`overhead`, `bottoms`, `side draw k`,
    /// `liquid distillate`).
    pub product: String,
    /// Mole fraction of the product without a boiling-point
    /// characterisation (database components), left out of the curve.
    pub excluded_mole_fraction: f64,
    /// `(volume %, K)` TBP points.
    pub tbp_k: Vec<(f64, f64)>,
    /// `(volume %, K)` ASTM D86 points.
    pub d86_k: Vec<(f64, f64)>,
}

/// Normal boiling point and standard liquid molar volume (m³/mol) of each
/// package component; `None` for a component that is not a pseudo-component.
struct CurveData(Vec<Option<(f64, f64)>>);

impl CurveData {
    fn new(ctx: &CalcContext<'_>) -> Result<Self, String> {
        let volumes = standard_volumes(ctx)?;
        let components = ctx
            .package()
            .formation_components()
            .map_err(|e| e.to_string())?;
        Ok(Self(
            components
                .iter()
                .zip(volumes)
                .map(|(c, v)| {
                    c.pseudo_characterisation()
                        .zip(v)
                        .map(|(p, v)| (p.input.normal_boiling_point_k, v))
                })
                .collect(),
        ))
    }

    /// The pseudo-component TBP curve of a composition and the excluded mole
    /// fraction; `None` without pseudo-components.
    fn curve(&self, x: &[f64]) -> Result<Option<(TbpCurve, f64)>, String> {
        let mut excluded = 0.0;
        let mut present = false;
        let mut components = Vec::new();
        for (data, x) in self.0.iter().zip(x) {
            match data {
                Some((tb, v)) => {
                    present |= *x > 0.0;
                    components.push(CurveComponent {
                        normal_boiling_point_k: *tb,
                        liquid_volume: x.max(0.0) * v,
                    });
                }
                None => excluded += x,
            }
        }
        if !present {
            return Ok(None);
        }
        TbpCurve::from_components(&components)
            .map(|c| Some((c, excluded)))
            .map_err(|e| e.to_string())
    }

    /// TBP range of the pseudo-components, K.
    fn range(&self) -> Result<(f64, f64), String> {
        let unit: Vec<f64> = self.0.iter().map(|d| f64::from(d.is_some())).collect();
        let (curve, _) = self.curve(&unit)?.expect("checked pseudo-components");
        Ok((
            curve.temperature_at(0.0).unwrap(),
            curve.temperature_at(100.0).unwrap(),
        ))
    }
}

/// The material outlet labels of a fractionator, in output order.
fn product_labels(spec: &PetroleumFractionatorSpec) -> Vec<String> {
    ["overhead".to_owned(), "bottoms".to_owned()]
        .into_iter()
        .chain((1..=spec.side_draws).map(|k| format!("side draw {k}")))
        .chain(
            matches!(
                spec.condenser,
                DistillationCondenser::PartialVaporLiquid { .. }
            )
            .then(|| "liquid distillate".to_owned()),
        )
        .collect()
}

/// TBP and D86 points of every product with pseudo-components; a product
/// whose curve cannot be formed is reported as a warning.
pub(super) fn product_curves(
    spec: &PetroleumFractionatorSpec,
    output: &mut UnitOutput,
    ctx: &CalcContext<'_>,
) -> Result<Vec<ProductCurve>, String> {
    let data = CurveData::new(ctx)?;
    let mut curves = Vec::new();
    for (label, outlet) in product_labels(spec).into_iter().zip(&output.outlets) {
        let Some(state) = &outlet.state else { continue };
        let points = data.curve(&state.mole_fractions).and_then(|curve| {
            let Some((curve, excluded)) = curve else {
                return Ok(None);
            };
            let mut tbp = Vec::new();
            let mut d86 = Vec::new();
            for v in REPORTED_VOLUME_PERCENTS {
                tbp.push((v, curve.temperature_at(v).map_err(|e| e.to_string())?));
                d86.push((v, curve.d86_temperature_at(v).map_err(|e| e.to_string())?));
            }
            Ok(Some((tbp, d86, excluded)))
        });
        match points {
            Ok(Some((tbp_k, d86_k, excluded_mole_fraction))) => curves.push(ProductCurve {
                product: label,
                excluded_mole_fraction,
                tbp_k,
                d86_k,
            }),
            Ok(None) => {}
            Err(e) => output
                .warnings
                .push(format!("{label} distillation curve is not reported: {e}")),
        }
    }
    Ok(curves)
}

/// Solves a column calculation with the cut-point products' molar flows.
type Rates = fn(
    &PetroleumFractionatorSpec,
    &[MaterialStream],
    &CalcContext<'_>,
) -> Result<(UnitOutput, DistillationColumnResult), String>;

/// Solves the fractionator with cut-point side products by an outer
/// iteration on their molar flows; `rates` solves the column at fixed flows.
pub(super) fn solve(
    spec: &PetroleumFractionatorSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
    rates: Rates,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    let data = CurveData::new(ctx)?;
    if data.0.iter().all(Option::is_none) {
        return Err(
            "petroleum fractionator needs characterised pseudo-components; the selected package has none"
                .into(),
        );
    }
    let draws = usize::from(spec.side_draws);
    let cuts: Vec<(usize, CutPoint)> = (0..draws)
        .filter_map(|k| match spec.products.get(k).map(|p| p.rate) {
            Some(DrawRate::CutPoint(c)) => Some((k, c)),
            _ => None,
        })
        .collect();
    let (low, high) = data.range()?;
    let mut targets = Vec::new();
    for (k, cut) in &cuts {
        let target = cut.temperature_k.expect("checked cut-point temperature");
        if cut.basis == CurveBasis::Tbp && !(low < target && target < high) {
            return Err(format!(
                "petroleum fractionator side draw {} cut point {target} K is not attainable: the pseudo-components' TBP curve spans {low:.2} to {high:.2} K",
                k + 1
            ));
        }
        targets.push(target);
    }
    let feed: f64 = inlets[..usize::from(spec.feeds)]
        .iter()
        .map(|s| s.molar_flow_mol_s)
        .sum();
    let start = feed / (2.0 * (draws as f64 + 1.0));
    let mut u: Vec<f64> = cuts
        .iter()
        .map(|(_, c)| c.estimate_mol_s.unwrap_or(start).ln())
        .collect();
    let mut column = spec.clone();
    column.products.resize(draws, Default::default());
    let evaluate =
        |u: &[f64]| -> Result<(Vec<f64>, (UnitOutput, DistillationColumnResult)), String> {
            let mut column = column.clone();
            for ((k, _), u) in cuts.iter().zip(u) {
                column.products[*k].rate = DrawRate::Molar(Some(u.exp()));
            }
            let solved = rates(&column, inlets, ctx)?;
            let mut residuals = Vec::new();
            for ((k, cut), target) in cuts.iter().zip(&targets) {
                let state = solved.0.outlets[2 + k]
                    .state
                    .as_ref()
                    .expect("column products carry a state");
                let (curve, _) = data.curve(&state.mole_fractions)?.ok_or_else(|| {
                    format!(
                        "side draw {} carries no pseudo-component for its cut point",
                        k + 1
                    )
                })?;
                let t = match cut.basis {
                    CurveBasis::Tbp => curve.temperature_at(cut.volume_percent),
                    CurveBasis::AstmD86 => curve.d86_temperature_at(cut.volume_percent),
                }
                .map_err(|e| e.to_string())?;
                residuals.push(t - target);
            }
            Ok((residuals, solved))
        };
    let flows = |u: &[f64]| u.iter().map(|u| u.exp()).collect::<Vec<_>>();
    let norm = |r: &[f64]| r.iter().map(|r| r * r).sum::<f64>().sqrt();
    let (mut r, mut solved) = evaluate(&u).map_err(|e| {
        format!(
            "petroleum fractionator cut-point iteration at the starting side-product molar flows {:?} mol/s: {e}",
            flows(&u)
        )
    })?;
    let jacobian = |u: &[f64], r: &[f64]| -> Result<Vec<Vec<f64>>, String> {
        let m = u.len();
        let mut jac = vec![vec![0.0; m]; m];
        for j in 0..m {
            let mut shifted = u.to_vec();
            shifted[j] += DIFFERENCE_STEP;
            let (rj, _) = evaluate(&shifted).map_err(|e| {
                format!(
                    "petroleum fractionator cut-point sensitivity at side-product molar flows {:?} mol/s: {e}",
                    flows(&shifted)
                )
            })?;
            for i in 0..m {
                jac[i][j] = (rj[i] - r[i]) / DIFFERENCE_STEP;
            }
        }
        Ok(jac)
    };
    let mut jac = jacobian(&u, &r)?;
    let mut fresh = true;
    let mut limited = 0;
    for _ in 0..CUT_POINT_ITERATIONS {
        if r.iter().all(|r| r.abs() <= CUT_POINT_TOLERANCE_K) {
            return Ok(solved);
        }
        let step = solve_linear(&jac, &r.iter().map(|r| -r).collect::<Vec<_>>()).ok_or_else(|| {
            format!(
                "petroleum fractionator cut points do not respond independently to the side-product molar flows {:?} mol/s",
                flows(&u)
            )
        })?;
        let largest = step.iter().fold(0.0_f64, |m, s| m.max(s.abs()));
        let mut scale = (MAX_LOG_STEP / largest).min(1.0);
        let mut accepted = None;
        let mut failure = None;
        for _ in 0..=STEP_HALVINGS {
            let trial: Vec<f64> = u.iter().zip(&step).map(|(u, s)| u + scale * s).collect();
            match evaluate(&trial) {
                Ok((rt, st)) if norm(&rt) < norm(&r) => {
                    accepted = Some((trial, rt, st));
                    break;
                }
                Ok(_) => {}
                Err(e) => failure = Some(e),
            }
            scale *= 0.5;
        }
        let Some((trial, rt, st)) = accepted else {
            if !fresh && failure.is_none() {
                // A stale quasi-Newton direction: restart from differences.
                jac = jacobian(&u, &r)?;
                fresh = true;
                continue;
            }
            return Err(unattained(&cuts, &targets, &r, &flows(&u), failure));
        };
        // Steps cut back by column failures twice in a row: the targets lie
        // beyond the flows the column can deliver.
        if failure.is_some() {
            limited += 1;
            if limited >= 2 {
                return Err(unattained(&cuts, &targets, &rt, &flows(&trial), failure));
            }
        } else {
            limited = 0;
        }
        // Broyden update J += (Δr − J Δu) Δuᵀ / (Δuᵀ Δu).
        let du: Vec<f64> = trial.iter().zip(&u).map(|(a, b)| a - b).collect();
        let dd: f64 = du.iter().map(|d| d * d).sum();
        for i in 0..jac.len() {
            let predicted: f64 = jac[i].iter().zip(&du).map(|(j, d)| j * d).sum();
            let correction = (rt[i] - r[i] - predicted) / dd;
            for (j, d) in jac[i].iter_mut().zip(&du) {
                *j += correction * d;
            }
        }
        fresh = false;
        (u, r, solved) = (trial, rt, st);
    }
    if r.iter().all(|r| r.abs() <= CUT_POINT_TOLERANCE_K) {
        return Ok(solved);
    }
    Err(format!(
        "petroleum fractionator cut points did not converge in {CUT_POINT_ITERATIONS} iterations (residuals {r:?} K at side-product molar flows {:?} mol/s)",
        flows(&u)
    ))
}

fn unattained(
    cuts: &[(usize, CutPoint)],
    targets: &[f64],
    residuals: &[f64],
    flows: &[f64],
    failure: Option<String>,
) -> String {
    let parts: Vec<String> = cuts
        .iter()
        .zip(targets)
        .zip(residuals)
        .zip(flows)
        .map(|((((k, cut), target), r), f)| {
            format!(
                "side draw {} {} {} vol% {target} K (closest {:.4} K at {f:.6} mol/s)",
                k + 1,
                match cut.basis {
                    CurveBasis::Tbp => "TBP",
                    CurveBasis::AstmD86 => "D86",
                },
                cut.volume_percent,
                target + r
            )
        })
        .collect();
    // Neither ending proves the target physically unattainable: the column
    // solve may fail numerically, and a stalled residual is a local minimum.
    let reason = failure.map_or_else(
        || "the residual no longer decreases".to_owned(),
        |e| format!("the column solve fails beyond these flows: {e}"),
    );
    format!(
        "petroleum fractionator cut-point specification could not be met: {}; {reason}",
        parts.join(", ")
    )
}

/// Solves `a x = b` by Gaussian elimination with partial pivoting; `None`
/// for a singular matrix.
fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut m: Vec<Vec<f64>> = a
        .iter()
        .zip(b)
        .map(|(row, b)| row.iter().copied().chain([*b]).collect())
        .collect();
    let scale = a
        .iter()
        .flatten()
        .fold(0.0_f64, |s, v| s.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))?;
        // Also refuses NaN pivots.
        if m[pivot][col].abs().is_nan() || m[pivot][col].abs() <= 1e-12 * scale {
            return None;
        }
        m.swap(col, pivot);
        let (top, rest) = m.split_at_mut(col + 1);
        let pivot_row = &top[col];
        for row in rest {
            let factor = row[col] / pivot_row[col];
            for (value, p) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *value -= factor * p;
            }
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let sum: f64 = (row + 1..n).map(|k| m[row][k] * x[k]).sum();
        x[row] = (m[row][n] - sum) / m[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cut_point_issues_follow_the_basis() {
        let mut issues = Vec::new();
        CutPoint {
            basis: CurveBasis::AstmD86,
            volume_percent: 97.0,
            temperature_k: None,
            estimate_mol_s: Some(-1.0),
        }
        .issues(&mut issues);
        let fields: Vec<_> = issues.iter().map(|i| i.field).collect();
        assert_eq!(
            fields,
            [
                "side-product cut-point volume percent",
                "side-product cut-point temperature",
                "side-product cut-point flow estimate"
            ]
        );
        let mut issues = Vec::new();
        CutPoint {
            volume_percent: 97.0,
            temperature_k: Some(500.0),
            ..Default::default()
        }
        .issues(&mut issues);
        assert!(issues.is_empty());
    }

    #[test]
    fn linear_solve_pivots_and_detects_singularity() {
        let x = solve_linear(&[vec![0.0, 2.0], vec![1.0, 1.0]], &[4.0, 3.0]).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-15 && (x[1] - 2.0).abs() < 1e-15);
        assert!(solve_linear(&[vec![1.0, 2.0], vec![2.0, 4.0]], &[1.0, 2.0]).is_none());
    }
}
