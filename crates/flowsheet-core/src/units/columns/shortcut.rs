//! Shortcut distillation column (#185): Fenske–Underwood–Gilliland design
//! estimate with Kirkbride feed location, a total condenser and a partial
//! reboiler at one column pressure.
//!
//! Method (all relative volatilities are with respect to the heavy key):
//!
//! - Feed thermal condition `q = (h_dew − h_F) / (h_dew − h_bub)` from the
//!   saturated-vapour and saturated-liquid enthalpies of the feed composition
//!   at column pressure and the actual feed enthalpy (isenthalpic let-down).
//! - Minimum stages, total reflux: Fenske, M. R. (1932) *Ind. Eng. Chem.* 24,
//!   482–485, <https://doi.org/10.1021/ie50269a003>.
//! - Non-key distribution at total reflux with the same `N_min`
//!   (Hengstebeck–Geddes form of the Fenske equation).
//! - Relative volatilities: geometric mean of the package K-value ratios at
//!   the dew point of the distillate (top) and the bubble point of the
//!   bottoms, iterated with the product split to a fixed point.
//! - Minimum reflux: Underwood, A. J. V. (1948) *Chem. Eng. Prog.* 44,
//!   603–614, for adjacent keys (one root between the key volatilities).
//! - Stages at finite reflux: Gilliland, E. R. (1940) *Ind. Eng. Chem.* 32,
//!   1220–1223, in the fit of Molokanov, Y. K. et al. (1972) *Int. Chem. Eng.*
//!   12, 209–212: `Y = 1 − exp[(1 + 54.4X)/(11 + 117.2X) · (X − 1)/√X]`,
//!   which tends to 1 at minimum reflux and to 0 at total reflux.
//! - Feed location: Kirkbride, C. G. (1944) *Petrol. Refiner* 23(9), 321–336.
//!
//! Stage counts include the partial reboiler and exclude the total
//! condenser. Bubble and dew points use the bounded A13 saturation solver
//! (`CubicPackage::saturation_temperature`), so every component with feed
//! flow must belong to its qualified dry methane/ethane/propane/N2/CO2
//! domain at 200–350 K and 10 kPa–10 MPa.

// Negated comparisons are deliberate: they also reject NaN values.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use thermo_solvers::{CubicSaturationKind, SaturationOptions, SaturationPoint};

use crate::engine::MaterialStream;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};
use crate::units::{CalcContext, UnitOutput};

/// Reflux specification of a shortcut column.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ShortcutReflux {
    /// Operating reflux ratio L/D, greater than the Underwood minimum.
    Ratio(Option<f64>),
    /// Operating reflux as a multiple of the Underwood minimum, greater than 1.
    MinimumMultiple(Option<f64>),
}

impl Default for ShortcutReflux {
    fn default() -> Self {
        Self::MinimumMultiple(None)
    }
}

/// Shortcut column specification (#185). Key identifiers are package
/// component identifiers; recoveries are fractions of the key's feed flow.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ShortcutColumnSpec {
    /// Light key, recovered mainly in the distillate.
    pub light_key: Option<String>,
    /// Heavy key, recovered mainly in the bottoms.
    pub heavy_key: Option<String>,
    /// Fraction of the light-key feed in the distillate, strictly 0 to 1.
    pub light_key_recovery: Option<f64>,
    /// Fraction of the heavy-key feed in the bottoms, strictly 0 to 1.
    pub heavy_key_recovery: Option<f64>,
    /// Absolute column pressure in Pa (condenser and reboiler), at most the
    /// feed pressure.
    pub pressure_pa: Option<f64>,
    /// Operating reflux.
    pub reflux: ShortcutReflux,
}

/// Design estimates of a solved shortcut column.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShortcutColumnResult {
    /// Light-key/heavy-key relative volatility (geometric mean).
    pub relative_volatility: f64,
    /// Feed thermal condition q (1 saturated liquid, 0 saturated vapour).
    pub feed_condition_q: f64,
    /// Fenske minimum theoretical stages at total reflux.
    pub minimum_stages: f64,
    /// Underwood root between the key volatilities.
    pub underwood_root: f64,
    /// Underwood minimum reflux ratio.
    pub minimum_reflux_ratio: f64,
    /// Operating reflux ratio L/D.
    pub reflux_ratio: f64,
    /// Gilliland abscissa X = (R − R_min)/(R + 1).
    pub gilliland_x: f64,
    /// Theoretical stages including the reboiler, excluding the condenser.
    pub theoretical_stages: f64,
    /// Kirkbride stages above the feed; the feed stage from the top is this
    /// value plus one.
    pub rectifying_stages: f64,
    /// Heat removed in the total condenser, W (positive out of the column).
    pub condenser_duty_w: f64,
    /// Heat supplied to the reboiler, W (positive into the column).
    pub reboiler_duty_w: f64,
    /// Dew-point temperature of the overhead vapour entering the condenser, K.
    pub top_vapor_temperature_k: f64,
    /// Relative-volatility fixed-point iterations.
    pub volatility_iterations: usize,
}

impl ShortcutColumnResult {
    /// Every value, for finiteness checks of stored results.
    pub fn numbers(&self) -> [f64; 13] {
        [
            self.relative_volatility,
            self.feed_condition_q,
            self.minimum_stages,
            self.underwood_root,
            self.minimum_reflux_ratio,
            self.reflux_ratio,
            self.gilliland_x,
            self.theoretical_stages,
            self.rectifying_stages,
            self.condenser_duty_w,
            self.reboiler_duty_w,
            self.top_vapor_temperature_k,
            self.volatility_iterations as f64,
        ]
    }
}

fn issue(out: &mut Vec<SpecIssue>, field: &'static str, problem: SpecProblem) {
    out.push(SpecIssue { field, problem });
}

impl ShortcutColumnSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut out = Vec::new();
        let key = |out: &mut Vec<SpecIssue>, field, k: &Option<String>| {
            if k.as_deref().is_none_or(|k| k.trim().is_empty()) {
                issue(out, field, SpecProblem::Missing);
            }
        };
        key(&mut out, "light key", &self.light_key);
        key(&mut out, "heavy key", &self.heavy_key);
        if self.light_key.is_some() && self.light_key == self.heavy_key {
            issue(&mut out, "heavy key", SpecProblem::Duplicate);
        }
        let recovery = |out: &mut Vec<SpecIssue>, field, r: Option<f64>| {
            let before = out.len();
            check_optional(out, field, r, Sign::Fraction);
            if out.len() == before && matches!(r, Some(v) if v == 0.0 || v == 1.0) {
                issue(out, field, SpecProblem::OutOfRange);
            }
        };
        recovery(&mut out, "light key recovery", self.light_key_recovery);
        recovery(&mut out, "heavy key recovery", self.heavy_key_recovery);
        if let (Some(a), Some(b)) = (self.light_key_recovery, self.heavy_key_recovery)
            && a > 0.0
            && a < 1.0
            && b > 0.0
            && b < 1.0
            && a + b <= 1.0
        {
            // Fenske needs a separation: r_LK + r_HK > 1.
            issue(&mut out, "heavy key recovery", SpecProblem::OutOfRange);
        }
        check_optional(
            &mut out,
            "column pressure",
            self.pressure_pa,
            Sign::Positive,
        );
        match self.reflux {
            ShortcutReflux::Ratio(r) => check_optional(&mut out, "reflux ratio", r, Sign::Positive),
            ShortcutReflux::MinimumMultiple(m) => {
                let before = out.len();
                check_optional(&mut out, "reflux multiple", m, Sign::Positive);
                if out.len() == before && m.is_some_and(|m| m <= 1.0) {
                    issue(&mut out, "reflux multiple", SpecProblem::OutOfRange);
                }
            }
        }
        out
    }
}

/// Relative-volatility fixed-point tolerance on max |Δ ln α|.
const ALPHA_TOLERANCE: f64 = 1e-10;
const ALPHA_ITERATIONS: usize = 100;

/// Gilliland's data range in X; results outside it are extrapolations.
const GILLILAND_X_RANGE: [f64; 2] = [0.01, 0.9];

impl CalcContext<'_> {
    /// Bubble or dew point temperature at `p_pa` with the bounded A13 solver.
    fn saturation_t(
        &self,
        x: &[f64],
        p_pa: f64,
        kind: CubicSaturationKind,
        what: &str,
    ) -> Result<SaturationPoint, String> {
        self.cubic_package()?
            .saturation_temperature(x, p_pa, kind, SaturationOptions::default())
            .map_err(|e| format!("shortcut column {what}: {e}"))
    }

    /// Relative volatilities K_i/K_hk at a saturation point.
    fn relative_volatilities(&self, s: &SaturationPoint, hk: usize) -> Result<Vec<f64>, String> {
        let k = self.k_values(
            &s.liquid_composition,
            &s.vapor_composition,
            s.temperature_k,
            s.pressure_pa,
        )?;
        Ok(k.iter().map(|ki| ki / k[hk]).collect())
    }
}

/// Fenske split at total reflux for volatilities `alpha`.
struct Split {
    minimum_stages: f64,
    distillate: Vec<f64>,
    bottoms: Vec<f64>,
}

impl Split {
    fn new(feed: &[f64], alpha: &[f64], lk: usize, r_lk: f64, r_hk: f64) -> Result<Self, String> {
        if !(alpha[lk] > 1.0) {
            return Err(
                "shortcut column light key is not more volatile than the heavy key at column conditions"
                    .into(),
            );
        }
        let ln_separation = (r_lk / (1.0 - r_lk)).ln() + (r_hk / (1.0 - r_hk)).ln();
        let minimum_stages = ln_separation / alpha[lk].ln();
        // ln(d/b) of the heavy key; every other component follows Fenske.
        let ln_hk = ((1.0 - r_hk) / r_hk).ln();
        let (mut distillate, mut bottoms) = (vec![0.0; feed.len()], vec![0.0; feed.len()]);
        for (i, &f) in feed.iter().enumerate() {
            if f > 0.0 {
                let ln_ratio = ln_hk + minimum_stages * alpha[i].ln();
                // Logistic form avoids overflow for very light or heavy components.
                let recovery = if ln_ratio >= 0.0 {
                    1.0 / (1.0 + (-ln_ratio).exp())
                } else {
                    let e = ln_ratio.exp();
                    e / (1.0 + e)
                };
                distillate[i] = f * recovery;
                bottoms[i] = f - distillate[i];
            }
        }
        Ok(Self {
            minimum_stages,
            distillate,
            bottoms,
        })
    }
}

fn fractions(flows: &[f64]) -> (f64, Vec<f64>) {
    let total: f64 = flows.iter().sum();
    (total, flows.iter().map(|f| f / total).collect())
}

/// Underwood root in (1, α_LK) of Σ α_i z_i/(α_i − θ) = 1 − q by bisection
/// (the left side increases monotonically between adjacent poles).
fn underwood_root(z: &[f64], alpha: &[f64], alpha_lk: f64, q: f64) -> Result<f64, String> {
    let f = |theta: f64| -> f64 {
        z.iter()
            .zip(alpha)
            .filter(|(zi, _)| **zi > 0.0)
            .map(|(zi, a)| a * zi / (a - theta))
            .sum::<f64>()
            - (1.0 - q)
    };
    let (mut lo, mut hi) = (1.0_f64, alpha_lk);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let theta = 0.5 * (lo + hi);
    if theta.is_finite() && theta > 1.0 && theta < alpha_lk {
        Ok(theta)
    } else {
        Err("shortcut column Underwood root not bracketed between the key volatilities".into())
    }
}

/// Molokanov et al. (1972) fit of the Gilliland correlation, X in (0, 1).
fn gilliland_y(x: f64) -> f64 {
    1.0 - ((1.0 + 54.4 * x) / (11.0 + 117.2 * x) * (x - 1.0) / x.sqrt()).exp()
}

fn component_index(ids: &[&str], key: &Option<String>, field: &str) -> Result<usize, String> {
    let key = key.as_deref().unwrap_or_default();
    ids.iter()
        .position(|id| *id == key)
        .ok_or_else(|| format!("shortcut column {field} '{key}' is not a package component"))
}

fn finite(name: &str, v: f64) -> Result<f64, String> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(format!("shortcut column {name} is not finite"))
    }
}

pub(super) fn calculate(
    spec: &ShortcutColumnSpec,
    feed: &MaterialStream,
    ctx: &CalcContext,
) -> Result<(UnitOutput, ShortcutColumnResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!("{} {}", issue.field, issue.problem));
    }
    let ids = ctx.component_ids();
    let lk = component_index(&ids, &spec.light_key, "light key")?;
    let hk = component_index(&ids, &spec.heavy_key, "heavy key")?;
    let (r_lk, r_hk) = (
        spec.light_key_recovery.unwrap_or_default(),
        spec.heavy_key_recovery.unwrap_or_default(),
    );
    let p = spec.pressure_pa.unwrap_or_default();

    let flow = feed.molar_flow_mol_s;
    if !flow.is_finite() || flow <= 0.0 {
        return Err("shortcut column requires a positive, finite feed flow".into());
    }
    let state = feed
        .state
        .as_ref()
        .ok_or("shortcut column feed has no thermodynamic state")?;
    if p > state.pressure_pa {
        return Err(format!(
            "shortcut column pressure {p} Pa exceeds the feed pressure {} Pa; pressurize the feed upstream",
            state.pressure_pa
        ));
    }
    let z = &state.mole_fractions;
    for (index, field) in [(lk, "light key"), (hk, "heavy key")] {
        if !(z[index] > 0.0) {
            return Err(format!(
                "shortcut column {field} '{}' is absent from the feed",
                ids[index]
            ));
        }
    }
    let feed_flows: Vec<f64> = z.iter().map(|zi| flow * zi).collect();

    // Feed thermal condition from saturated enthalpies at column pressure.
    let feed_bubble = ctx.saturation_t(z, p, CubicSaturationKind::Bubble, "feed bubble point")?;
    let feed_dew = ctx.saturation_t(z, p, CubicSaturationKind::Dew, "feed dew point")?;
    let latent = feed_dew.vapor.enthalpy - feed_bubble.liquid.enthalpy;
    if !(latent > 0.0) {
        return Err(
            "shortcut column feed has no positive heat of vaporization at column pressure".into(),
        );
    }
    let q = finite(
        "feed condition",
        (feed_dew.vapor.enthalpy - state.enthalpy()) / latent,
    )?;

    // Fixed point of volatilities and Fenske split.
    let mut alpha = ctx.relative_volatilities(&feed_bubble, hk)?;
    let mut iterations = 0;
    let (split, top, bottom) = loop {
        iterations += 1;
        let split = Split::new(&feed_flows, &alpha, lk, r_lk, r_hk)?;
        let (_, x_d) = fractions(&split.distillate);
        let (_, x_b) = fractions(&split.bottoms);
        let top = ctx.saturation_t(&x_d, p, CubicSaturationKind::Dew, "overhead dew point")?;
        let bottom =
            ctx.saturation_t(&x_b, p, CubicSaturationKind::Bubble, "bottoms bubble point")?;
        let a_top = ctx.relative_volatilities(&top, hk)?;
        let a_bottom = ctx.relative_volatilities(&bottom, hk)?;
        let next: Vec<f64> = a_top
            .iter()
            .zip(&a_bottom)
            .map(|(a, b)| (a * b).sqrt())
            .collect();
        let change = next
            .iter()
            .zip(&alpha)
            .zip(z)
            .filter(|(_, zi)| **zi > 0.0)
            .map(|((n, o), _)| (n.ln() - o.ln()).abs())
            .fold(0.0, f64::max);
        if !change.is_finite() {
            return Err("shortcut column relative volatility is not finite".into());
        }
        alpha = next;
        if change < ALPHA_TOLERANCE {
            break (
                Split::new(&feed_flows, &alpha, lk, r_lk, r_hk)?,
                top,
                bottom,
            );
        }
        if iterations >= ALPHA_ITERATIONS {
            return Err(format!(
                "shortcut column relative volatilities did not converge in {ALPHA_ITERATIONS} iterations"
            ));
        }
    };
    let alpha_lk = alpha[lk];
    if let Some(i) = (0..z.len())
        .find(|&i| i != lk && i != hk && z[i] > 0.0 && alpha[i] > 1.0 && alpha[i] < alpha_lk)
    {
        return Err(format!(
            "shortcut column keys are not adjacent: '{}' is between them in volatility",
            ids[i]
        ));
    }
    let (d, x_d) = fractions(&split.distillate);
    let (b, x_b) = fractions(&split.bottoms);

    // Underwood minimum reflux with the converged volatilities.
    let theta = underwood_root(z, &alpha, alpha_lk, q)?;
    let r_min = finite(
        "minimum reflux",
        x_d.iter()
            .zip(&alpha)
            .filter(|(x, _)| **x > 0.0)
            .map(|(x, a)| a * x / (a - theta))
            .sum::<f64>()
            - 1.0,
    )?;
    if !(r_min > 0.0) {
        return Err(format!(
            "shortcut column Underwood minimum reflux {r_min:.6} is not positive; the specification is outside the shortcut method"
        ));
    }
    let reflux = match spec.reflux {
        ShortcutReflux::Ratio(r) => r.unwrap_or_default(),
        ShortcutReflux::MinimumMultiple(m) => m.unwrap_or_default() * r_min,
    };
    if !(reflux > r_min) {
        return Err(format!(
            "shortcut column reflux ratio {reflux:.6} must exceed the minimum reflux ratio {r_min:.6}"
        ));
    }
    let x = (reflux - r_min) / (reflux + 1.0);
    let y = gilliland_y(x);
    let stages = finite(
        "theoretical stage count",
        (split.minimum_stages + y) / (1.0 - y),
    )?;
    let kirkbride = ((b / d) * (z[hk] / z[lk]) * (x_b[lk] / x_d[hk]).powi(2)).powf(0.206);
    let rectifying = finite("feed location", stages * kirkbride / (1.0 + kirkbride))?;

    // Saturated-liquid products and duties.
    let distillate_bubble = ctx.saturation_t(
        &x_d,
        p,
        CubicSaturationKind::Bubble,
        "distillate bubble point",
    )?;
    let distillate_state = ctx.tp(&x_d, distillate_bubble.temperature_k, p)?;
    let bottoms_state = ctx.tp(&x_b, bottom.temperature_k, p)?;
    if !(bottom.temperature_k > top.temperature_k) {
        return Err(
            "shortcut column bottoms bubble point is not above the overhead dew point".into(),
        );
    }
    let condenser = finite(
        "condenser duty",
        (reflux + 1.0) * d * (top.vapor.enthalpy - distillate_state.enthalpy()),
    )?;
    if !(condenser > 0.0) {
        return Err("shortcut column condenser duty is not positive".into());
    }
    let reboiler = finite(
        "reboiler duty",
        d * distillate_state.enthalpy() + b * bottoms_state.enthalpy() + condenser
            - flow * state.enthalpy(),
    )?;
    if !(reboiler > 0.0) {
        return Err(format!(
            "shortcut column reboiler duty {reboiler:.6} W is not positive; the feed carries more heat than the separation removes"
        ));
    }

    let mut warnings = Vec::new();
    if !(GILLILAND_X_RANGE[0]..=GILLILAND_X_RANGE[1]).contains(&x) {
        warnings.push(format!(
            "Gilliland abscissa X = {x:.4} is outside the correlation's data range {}–{}; the stage count is an extrapolation.",
            GILLILAND_X_RANGE[0], GILLILAND_X_RANGE[1]
        ));
    }
    let distributing: Vec<&str> = (0..z.len())
        .filter(|&i| i != lk && i != hk && z[i] > 0.0)
        .filter(|&i| {
            let r = split.distillate[i] / feed_flows[i];
            (0.001..=0.999).contains(&r)
        })
        .map(|i| ids[i])
        .collect();
    if !distributing.is_empty() {
        warnings.push(format!(
            "Non-key components distribute between the products ({}); the Fenske total-reflux split is an approximation.",
            distributing.join(", ")
        ));
    }

    let result = ShortcutColumnResult {
        relative_volatility: alpha_lk,
        feed_condition_q: q,
        minimum_stages: split.minimum_stages,
        underwood_root: theta,
        minimum_reflux_ratio: r_min,
        reflux_ratio: reflux,
        gilliland_x: x,
        theoretical_stages: stages,
        rectifying_stages: rectifying,
        condenser_duty_w: condenser,
        reboiler_duty_w: reboiler,
        top_vapor_temperature_k: top.temperature_k,
        volatility_iterations: iterations,
    };
    let outlets = vec![
        MaterialStream {
            molar_flow_mol_s: d,
            state: Some(distillate_state),
        },
        MaterialStream {
            molar_flow_mol_s: b,
            state: Some(bottoms_state),
        },
    ];
    let mut output = UnitOutput::new(outlets, reboiler - condenser, warnings);
    output.energy = vec![Some(condenser), Some(reboiler)];
    Ok((output, result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use property_packages::CubicPackage;
    use thermo_data::ComponentDatabase;
    use thermo_models::{CubicModel, CubicVariant};
    use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

    fn package(ids: &[&str], z: &[f64]) -> CubicPackage {
        let composition = Composition::new(
            CompositionBasis::Mole,
            ids.iter()
                .zip(z)
                .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
                .collect(),
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let n = ids.len();
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition,
            CubicVariant::PengRobinsonBostonMathias,
            vec![vec![0.0; n]; n],
            "zero interactions for the shortcut column unit test",
        )
        .unwrap();
        CubicPackage::new(model, Default::default())
    }

    fn spec(reflux: ShortcutReflux) -> ShortcutColumnSpec {
        ShortcutColumnSpec {
            light_key: Some("ethane".into()),
            heavy_key: Some("propane".into()),
            light_key_recovery: Some(0.98),
            heavy_key_recovery: Some(0.97),
            pressure_pa: Some(2.0e6),
            reflux,
        }
    }

    fn feed(ctx: &CalcContext, z: &[f64], t: f64, p: f64) -> MaterialStream {
        MaterialStream {
            molar_flow_mol_s: 10.0,
            state: Some(ctx.tp(z, t, p).unwrap()),
        }
    }

    #[test]
    fn molokanov_fit_has_the_gilliland_limits() {
        assert!((gilliland_y(1e-12) - 1.0).abs() < 1e-6);
        assert!(gilliland_y(1.0).abs() < 1e-15);
        // Hand evaluation at X = 0.1: 6.44/22.72 · (−0.9/√0.1) = −0.806714,
        // Y = 1 − e^−0.806714 = 0.55368 (Eduljee 1975 gives 0.5466 nearby).
        assert!((gilliland_y(0.1) - 0.55368).abs() < 1e-4);
    }

    #[test]
    fn binary_column_closes_balances_and_meets_recoveries() {
        let z = [0.4, 0.6];
        let package = package(&["ethane", "propane"], &z);
        let ctx = CalcContext::new(&package, None);
        let inlet = feed(&ctx, &z, 280.0, 2.0e6);
        let (out, r) = calculate(
            &spec(ShortcutReflux::MinimumMultiple(Some(1.3))),
            &inlet,
            &ctx,
        )
        .unwrap();
        let [dist, bot] = out.outlets.as_slice() else {
            panic!("two products")
        };
        let fd = dist.component_flows_mol_s(2);
        let fb = bot.component_flows_mol_s(2);
        assert!((fd[0] / 4.0 - 0.98).abs() < 1e-12);
        assert!((fb[1] / 6.0 - 0.97).abs() < 1e-12);
        for i in 0..2 {
            assert!((fd[i] + fb[i] - 10.0 * z[i]).abs() < 1e-12);
        }
        let energy = dist.enthalpy_flow_w() + bot.enthalpy_flow_w() - inlet.enthalpy_flow_w();
        assert!((energy - out.duty_w).abs() < 1e-6 * r.condenser_duty_w);
        assert_eq!(
            out.energy,
            [Some(r.condenser_duty_w), Some(r.reboiler_duty_w)]
        );
        assert!((r.reflux_ratio / r.minimum_reflux_ratio - 1.3).abs() < 1e-12);
        assert!(r.theoretical_stages > r.minimum_stages && r.minimum_stages > 1.0);
        assert!(r.rectifying_stages > 0.0 && r.rectifying_stages < r.theoretical_stages);
        assert!(r.relative_volatility > 1.0);
        let (td, tb) = (
            dist.state.as_ref().unwrap().temperature_k,
            bot.state.as_ref().unwrap().temperature_k,
        );
        assert!(td < r.top_vapor_temperature_k && r.top_vapor_temperature_k < tb);
        // Binary Fenske with the reported volatility.
        let n_min = ((0.98_f64 / 0.02) * (0.97 / 0.03)).ln() / r.relative_volatility.ln();
        assert!((n_min - r.minimum_stages).abs() < 1e-12);
    }

    #[test]
    fn invalid_specifications_and_states_fail_explicitly() {
        let mut s = spec(ShortcutReflux::Ratio(None));
        s.heavy_key = s.light_key.clone();
        s.light_key_recovery = Some(1.0);
        s.heavy_key_recovery = Some(0.0);
        let problems: Vec<_> = s.issues().iter().map(|i| (i.field, i.problem)).collect();
        assert!(problems.contains(&("heavy key", SpecProblem::Duplicate)));
        assert!(problems.contains(&("light key recovery", SpecProblem::OutOfRange)));
        assert!(problems.contains(&("heavy key recovery", SpecProblem::OutOfRange)));
        assert!(problems.contains(&("reflux ratio", SpecProblem::Missing)));
        let mut s = spec(ShortcutReflux::MinimumMultiple(Some(1.0)));
        s.light_key_recovery = Some(0.5);
        s.heavy_key_recovery = Some(0.5);
        let problems: Vec<_> = s.issues().iter().map(|i| (i.field, i.problem)).collect();
        assert_eq!(
            problems,
            [
                ("heavy key recovery", SpecProblem::OutOfRange),
                ("reflux multiple", SpecProblem::OutOfRange)
            ]
        );

        let z = [0.4, 0.6];
        let package = package(&["ethane", "propane"], &z);
        let ctx = CalcContext::new(&package, None);
        let inlet = feed(&ctx, &z, 280.0, 2.0e6);
        let fail = |s: &ShortcutColumnSpec, f: &MaterialStream| calculate(s, f, &ctx).unwrap_err();
        // Swapped keys.
        let mut swapped = spec(ShortcutReflux::MinimumMultiple(Some(1.3)));
        std::mem::swap(&mut swapped.light_key, &mut swapped.heavy_key);
        assert!(fail(&swapped, &inlet).contains("not more volatile"));
        // Reflux below the minimum.
        assert!(fail(&spec(ShortcutReflux::Ratio(Some(0.01))), &inlet).contains("must exceed"));
        // Column above the feed pressure.
        let mut high = spec(ShortcutReflux::MinimumMultiple(Some(1.3)));
        high.pressure_pa = Some(3.0e6);
        assert!(fail(&high, &inlet).contains("exceeds the feed pressure"));
        // Unknown key, zero flow, missing state.
        let mut unknown = spec(ShortcutReflux::MinimumMultiple(Some(1.3)));
        unknown.light_key = Some("butane".into());
        assert!(fail(&unknown, &inlet).contains("not a package component"));
        let ok = spec(ShortcutReflux::MinimumMultiple(Some(1.3)));
        let empty = MaterialStream {
            molar_flow_mol_s: 0.0,
            state: None,
        };
        assert!(fail(&ok, &empty).contains("positive, finite feed flow"));
        // Near-critical pressure leaves the bounded saturation domain.
        let mut critical = spec(ShortcutReflux::MinimumMultiple(Some(1.3)));
        critical.pressure_pa = Some(5.5e6);
        let hot = feed(&ctx, &z, 340.0, 6.0e6);
        assert!(fail(&critical, &hot).starts_with("shortcut column"));
    }

    #[test]
    fn non_adjacent_keys_are_rejected() {
        let z = [0.2, 0.45, 0.35];
        let package = package(&["methane", "ethane", "propane"], &z);
        let ctx = CalcContext::new(&package, None);
        let inlet = feed(&ctx, &z, 260.0, 2.5e6);
        let mut s = spec(ShortcutReflux::MinimumMultiple(Some(1.3)));
        s.light_key = Some("methane".into());
        s.pressure_pa = Some(2.5e6);
        let err = calculate(&s, &inlet, &ctx).unwrap_err();
        assert!(err.contains("not adjacent"), "{err}");
    }
}
