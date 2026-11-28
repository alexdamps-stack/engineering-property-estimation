//! Three-phase distillation column (#188): equilibrium stages at one
//! pressure, a total condenser with a decanter and a partial reboiler
//! (Aspen Plus RadFrac with a decanter / HYSYS Three Phase Distillation).
//!
//! Supported scope (`docs/a17-multiphase-three-phase-column.md`):
//!
//! - Stage 1 is the total condenser with its decanter at a specified
//!   (subcooled) temperature. The vapor from stage 2 condenses completely and
//!   splits into a light and a heavy liquid (liquid–liquid equilibrium at the
//!   decanter temperature and column pressure).
//! - Stages 2..N−1 are ideal trays (or ideal packed stages) and stage N is the
//!   partial reboiler. These stages carry **one** liquid (vapor–liquid
//!   equilibrium). After convergence every stage is flashed with the package's
//!   three-phase operation; a second liquid on any of them is an explicit
//!   failure. A full VLLE stage model is not implemented.
//! - [`DecanterReflux`] chooses which liquid returns as reflux and how much.
//!   The `distillate` port carries the non-refluxed part of the reflux liquid
//!   (zero for total reflux of that phase); the `decant` port carries the
//!   other liquid. With [`DecanterReflux::Fractions`] the distillate is the
//!   non-refluxed light liquid and the decant the non-refluxed heavy liquid.
//!
//! Stage equations (Seader, Henley and Roper, *Separation Process
//! Principles*, 3rd ed., 2011, §10.1; in − out):
//!
//! - decanter: v_{2,i} − a_i − b_i = 0, ln x^a_i + ln φ^a_i − ln x^b_i −
//!   ln φ^b_i = 0, V_2 h^V_2 − Q_C − A h^a − B h^b = 0 (a light, b heavy
//!   liquid at the decanter temperature);
//! - stage j ≥ 2: l_{j−1,i} + v_{j+1,i} + f_{j,i} − l_{j,i} − v_{j,i} = 0
//!   with l_{1,i} = φ_L a_i + φ_H b_i (reflux fractions),
//!   ln y_ji − ln x_ji − ln K_ji = 0 (K = φ^L/φ^V of the package's liquid and
//!   vapor at the stage state), and the energy balance on stages 2..N−1;
//! - the reboiler replaces its energy balance by the bottoms flow
//!   specification L_N = B; its duty follows from its energy balance.
//!
//! Unknowns: ln component flows of both decanter liquids and the scaled
//! condenser duty; ln liquid and vapor component flows and temperature of
//! every other stage. The solve starts from constant-molar-overflow
//! bubble-point passes with a lagged decanter flash and finishes with a
//! damped Newton iteration on the block-tridiagonal stage Jacobian.

use property_packages::{
    Convergence, EquilibriumPhase, EquilibriumState, PhaseProperties, PhaseRole, RootSelection,
};
use thermo_types::{CalculationError, Phase};

use super::extractor::{block_tridiagonal_solve, thomas};
use crate::MaterialStream;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};
use crate::units::columns::stages::StageFeed;
use crate::units::columns::{ColumnInternals, MAX_COLUMN_STAGES, packing::PackedStages};
use crate::units::{CalcContext, UnitOutput};

const TOLERANCE: f64 = 1e-11;
const MAX_NEWTON: usize = 80;
const MAX_PASSES: usize = 60;
const MAX_LN_STEP: f64 = 2.0;
const MAX_T_STEP: f64 = 10.0;
const MIN_VAPOR_LIQUID_VOLUME_RATIO: f64 = 1.05;
/// Decanter liquids closer than this (max mole-fraction difference) are a
/// plait point or the trivial solution: explicit failure.
pub const DECANTER_PLAIT_DISTANCE: f64 = 1e-3;
/// Largest deviation of the solved decanter and stage phases from the
/// independent check flash.
const FLASH_AGREEMENT: f64 = 1e-6;

/// How much of one decanter liquid returns to the column.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PhaseReflux {
    /// The whole liquid returns; the distillate is zero.
    Total,
    /// Reflux divided by distillate of that liquid (molar), > 0.
    Ratio(Option<f64>),
}

/// Which decanter liquid is refluxed (mutually exclusive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DecanterReflux {
    /// The light liquid is refluxed; the heavy liquid leaves as decant.
    Light(PhaseReflux),
    /// The heavy liquid is refluxed; the light liquid leaves as decant.
    Heavy(PhaseReflux),
    /// Specified fractions (0..=1) of both liquids are refluxed; the rest of
    /// the light liquid is the distillate, the rest of the heavy liquid the
    /// decant.
    Fractions {
        light: Option<f64>,
        heavy: Option<f64>,
    },
}

impl Default for DecanterReflux {
    fn default() -> Self {
        Self::Light(PhaseReflux::Total)
    }
}

impl DecanterReflux {
    /// Refluxed fractions of the light and heavy liquid.
    fn fractions(self) -> Option<[f64; 2]> {
        let phase = |r: PhaseReflux| match r {
            PhaseReflux::Total => Some(1.0),
            PhaseReflux::Ratio(ratio) => ratio.map(|r| r / (1.0 + r)),
        };
        match self {
            Self::Light(r) => Some([phase(r)?, 0.0]),
            Self::Heavy(r) => Some([0.0, phase(r)?]),
            Self::Fractions { light, heavy } => Some([light?, heavy?]),
        }
    }

    fn issues(self, out: &mut Vec<SpecIssue>) {
        match self {
            Self::Light(PhaseReflux::Ratio(r)) | Self::Heavy(PhaseReflux::Ratio(r)) => {
                check_optional(out, "reflux ratio", r, Sign::Positive)
            }
            Self::Fractions { light, heavy } => {
                check_optional(out, "light reflux fraction", light, Sign::Fraction);
                check_optional(out, "heavy reflux fraction", heavy, Sign::Fraction);
                if let (Some(l), Some(h)) = (light, heavy)
                    && ((l == 0.0 && h == 0.0) || (l == 1.0 && h == 1.0))
                {
                    out.push(SpecIssue {
                        field: "reflux fractions",
                        problem: SpecProblem::OutOfRange,
                    });
                }
            }
            _ => {}
        }
    }
}

/// Three-phase distillation specification (#188).
#[derive(Debug, Clone, PartialEq)]
pub struct ThreePhaseColumnSpec {
    /// Ideal trays or one ideal packed bed between decanter and reboiler.
    pub internals: ColumnInternals,
    /// One-based feed stage: 2 (top tray) to N (reboiler).
    pub feed_stage: u16,
    /// Uniform absolute column pressure, Pa.
    pub pressure_pa: Option<f64>,
    /// Decanter (condensate) temperature, K, below the condensate bubble point.
    pub decanter_temperature_k: Option<f64>,
    /// Reflux phase and amount.
    pub reflux: DecanterReflux,
    /// Fraction of the feed withdrawn as bottoms, strictly between 0 and 1.
    pub bottoms_fraction: Option<f64>,
}

impl Default for ThreePhaseColumnSpec {
    fn default() -> Self {
        Self {
            internals: ColumnInternals::Trays(5),
            feed_stage: 2,
            pressure_pa: None,
            decanter_temperature_k: None,
            reflux: DecanterReflux::default(),
            bottoms_fraction: None,
        }
    }
}

/// Internal stage count and bed (same rules as the rigorous column).
fn resolve(internals: ColumnInternals) -> Result<(usize, Option<PackedStages>), String> {
    let (count, bed) = match internals {
        ColumnInternals::Trays(count) => (count, None),
        ColumnInternals::Packed(spec) => {
            let packed = spec.resolve()?;
            if !packed.packed_height_m.is_finite() {
                return Err("calculated packed height is not finite".into());
            }
            (packed.stages, Some(packed))
        }
    };
    if !(1..=MAX_COLUMN_STAGES).contains(&count) {
        return Err(format!(
            "three-phase column needs 1 to {MAX_COLUMN_STAGES} contacting stages, got {count}"
        ));
    }
    Ok((usize::from(count), bed))
}

impl ThreePhaseColumnSpec {
    /// Total stage count N including decanter (stage 1) and reboiler (stage N).
    pub fn stage_count(&self) -> Result<usize, String> {
        resolve(self.internals).map(|(n, _)| n + 2)
    }

    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        check_optional(
            &mut issues,
            "column pressure",
            self.pressure_pa,
            Sign::Positive,
        );
        check_optional(
            &mut issues,
            "decanter temperature",
            self.decanter_temperature_k,
            Sign::Positive,
        );
        self.reflux.issues(&mut issues);
        check_optional(
            &mut issues,
            "bottoms fraction",
            self.bottoms_fraction,
            Sign::Fraction,
        );
        if matches!(self.bottoms_fraction, Some(0.0 | 1.0)) {
            issues.push(SpecIssue {
                field: "bottoms fraction",
                problem: SpecProblem::OutOfRange,
            });
        }
        if let ColumnInternals::Packed(packed) = self.internals {
            issues.extend(packed.issues());
        }
        match self.stage_count() {
            Ok(n) => {
                if !(2..=n).contains(&usize::from(self.feed_stage)) {
                    issues.push(SpecIssue {
                        field: "feed stage",
                        problem: SpecProblem::OutOfRange,
                    });
                }
            }
            Err(_) if issues.iter().all(|i| i.field != "internals") => issues.push(SpecIssue {
                field: "internals",
                problem: SpecProblem::CountOutOfRange,
            }),
            Err(_) => {}
        }
        issues
    }
}

/// Solved profile of a three-phase column. Stage vectors run from the
/// decanter (stage 1) to the reboiler (stage N).
#[derive(Debug, Clone, PartialEq)]
pub struct ThreePhaseColumnResult {
    /// Stage temperatures; stage 1 is the decanter temperature.
    pub stage_temperatures_k: Vec<f64>,
    /// Liquid leaving each stage downward, mol/s; stage 1: the reflux.
    pub liquid_flows_mol_s: Vec<f64>,
    /// Vapor leaving each stage upward, mol/s; stage 1: zero.
    pub vapor_flows_mol_s: Vec<f64>,
    /// Liquid mole fractions; stage 1: the reflux composition.
    pub liquid_compositions: Vec<Vec<f64>>,
    /// Vapor mole fractions; stage 1: the condensate (vapor from stage 2).
    pub vapor_compositions: Vec<Vec<f64>>,
    /// Light decanter liquid: flow (mol/s) and mole fractions.
    pub light_liquid_mol_s: f64,
    pub light_liquid_composition: Vec<f64>,
    /// Heavy decanter liquid: flow (mol/s) and mole fractions.
    pub heavy_liquid_mol_s: f64,
    pub heavy_liquid_composition: Vec<f64>,
    /// Whether the light liquid is (part of) the reflux; false when only the
    /// heavy liquid returns.
    pub light_refluxed: bool,
    /// Heat removed by the condenser, W.
    pub condenser_duty_w: f64,
    /// Heat supplied to the reboiler, W.
    pub reboiler_duty_w: f64,
    /// Bed height supplied or derived from the specified HETP.
    pub packed_height_m: Option<f64>,
    /// Accepted Newton iterations.
    pub iterations: usize,
}

/// Properties of one block: the two decanter liquids, or a stage's liquid
/// and vapor.
#[derive(Debug, Clone)]
struct BlockProperties {
    first: PhaseProperties,
    second: PhaseProperties,
}

/// Decoded block: decanter (`l` light, `v` heavy, `t` = scaled condenser
/// duty) or stage (`l` liquid, `v` vapor component flows, `t` temperature).
#[derive(Debug, Clone)]
struct Block {
    l: Vec<f64>,
    v: Vec<f64>,
    t: f64,
}

fn total(flows: &[f64]) -> f64 {
    flows.iter().sum()
}

fn normalized(flows: &[f64]) -> Vec<f64> {
    let sum = total(flows);
    flows.iter().map(|f| f / sum).collect()
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|x| x.abs()).fold(0.0, f64::max)
}

fn sum_squares(values: &[f64]) -> f64 {
    values.iter().map(|x| x * x).sum()
}

fn max_difference(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

struct Column<'a, 'p> {
    ctx: &'a CalcContext<'p>,
    pressure: f64,
    decanter_t: f64,
    range: [f64; 2],
    /// Blocks: decanter, contacting stages, reboiler.
    blocks: usize,
    components: usize,
    active: Vec<usize>,
    /// Feed component flows and enthalpy flow per block.
    feed_flows: Vec<Vec<f64>>,
    feed_enthalpy: Vec<f64>,
    /// Per active component: feed flow.
    scale: Vec<f64>,
    flow: f64,
    energy_scale: f64,
    bottoms: f64,
    /// Refluxed fractions of the light (a) and heavy (b) liquid.
    reflux: [f64; 2],
}

impl Column<'_, '_> {
    fn width(&self) -> usize {
        2 * self.active.len() + 1
    }

    fn decode(&self, v: &[f64], j: usize) -> Block {
        let (w, m) = (self.width(), self.active.len());
        let block = &v[j * w..(j + 1) * w];
        let mut l = vec![0.0; self.components];
        let mut u = vec![0.0; self.components];
        for (k, &i) in self.active.iter().enumerate() {
            l[i] = block[k].exp();
            u[i] = block[m + k].exp();
        }
        Block {
            l,
            v: u,
            t: block[2 * m],
        }
    }

    fn encode(&self, blocks: &[Block]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.blocks * self.width());
        for b in blocks {
            out.extend(self.active.iter().map(|&i| b.l[i].ln()));
            out.extend(self.active.iter().map(|&i| b.v[i].ln()));
            out.push(b.t);
        }
        out
    }

    fn properties(&self, v: &[f64], j: usize) -> Result<BlockProperties, String> {
        let b = self.decode(v, j);
        if b.l.iter().chain(&b.v).any(|f| !f.is_finite()) || !b.t.is_finite() {
            return Err(format!("stage {} flows are not finite", j + 1));
        }
        let liquid = |flows: &[f64], t| {
            self.ctx.homogeneous(
                &normalized(flows),
                t,
                self.pressure,
                RootSelection::SmallestVolume,
            )
        };
        if j == 0 {
            return Ok(BlockProperties {
                first: liquid(&b.l, self.decanter_t)?,
                second: liquid(&b.v, self.decanter_t)?,
            });
        }
        let [lo, hi] = self.range;
        if !(lo..=hi).contains(&b.t) {
            return Err(format!(
                "stage {} temperature {:.3} K leaves the column interval {lo:.2}..={hi:.2} K",
                j + 1,
                b.t
            ));
        }
        let first = liquid(&b.l, b.t)?;
        let second = self.ctx.homogeneous(
            &normalized(&b.v),
            b.t,
            self.pressure,
            RootSelection::LargestVolume,
        )?;
        if second.molar_volume <= MIN_VAPOR_LIQUID_VOLUME_RATIO * first.molar_volume {
            return Err(format!(
                "stage {} lacks distinct liquid and vapor roots",
                j + 1
            ));
        }
        Ok(BlockProperties { first, second })
    }

    fn all_properties(&self, v: &[f64]) -> Result<Vec<BlockProperties>, String> {
        (0..self.blocks).map(|j| self.properties(v, j)).collect()
    }

    /// Reflux component flows and enthalpy flow.
    fn reflux_of(&self, decanter: &Block, props: &BlockProperties) -> (Vec<f64>, f64) {
        let [fa, fb] = self.reflux;
        let flows = decanter
            .l
            .iter()
            .zip(&decanter.v)
            .map(|(a, b)| fa * a + fb * b)
            .collect();
        let h = fa * total(&decanter.l) * props.first.enthalpy
            + fb * total(&decanter.v) * props.second.enthalpy;
        (flows, h)
    }

    fn residual(&self, v: &[f64], props: &[BlockProperties]) -> Vec<f64> {
        let n = self.blocks;
        let blocks: Vec<Block> = (0..n).map(|j| self.decode(v, j)).collect();
        let mut out = Vec::with_capacity(n * self.width());
        let top = &blocks[1];
        // Decanter.
        {
            let d = &blocks[0];
            let p = &props[0];
            for (k, &i) in self.active.iter().enumerate() {
                out.push((top.v[i] - d.l[i] - d.v[i]) / self.scale[k]);
            }
            let (ln_a, ln_b) = (total(&d.l).ln(), total(&d.v).ln());
            for &i in &self.active {
                out.push(
                    d.l[i].ln() - ln_a + p.first.ln_fugacity_coefficients[i] - d.v[i].ln() + ln_b
                        - p.second.ln_fugacity_coefficients[i],
                );
            }
            let h_in = total(&top.v) * props[1].second.enthalpy;
            let h_out = d.t * self.energy_scale
                + total(&d.l) * p.first.enthalpy
                + total(&d.v) * p.second.enthalpy;
            out.push((h_in - h_out) / self.energy_scale);
        }
        let (reflux, reflux_h) = self.reflux_of(&blocks[0], &props[0]);
        for j in 1..n {
            let s = &blocks[j];
            let p = &props[j];
            let (l_in, h_l_in) = if j == 1 {
                (reflux.clone(), reflux_h)
            } else {
                (
                    blocks[j - 1].l.clone(),
                    total(&blocks[j - 1].l) * props[j - 1].first.enthalpy,
                )
            };
            let zero = vec![0.0; self.components];
            let (v_in, h_v_in) = if j + 1 < n {
                (
                    &blocks[j + 1].v,
                    total(&blocks[j + 1].v) * props[j + 1].second.enthalpy,
                )
            } else {
                (&zero, 0.0)
            };
            for (k, &i) in self.active.iter().enumerate() {
                out.push(
                    (l_in[i] + v_in[i] + self.feed_flows[j][i] - s.l[i] - s.v[i]) / self.scale[k],
                );
            }
            let (ln_l, ln_v) = (total(&s.l).ln(), total(&s.v).ln());
            for &i in &self.active {
                out.push(
                    s.v[i].ln() - ln_v - s.l[i].ln() + ln_l
                        - (p.first.ln_fugacity_coefficients[i]
                            - p.second.ln_fugacity_coefficients[i]),
                );
            }
            if j + 1 < n {
                let h_in = h_l_in + h_v_in + self.feed_enthalpy[j];
                let h_out = total(&s.l) * p.first.enthalpy + total(&s.v) * p.second.enthalpy;
                out.push((h_in - h_out) / self.energy_scale);
            } else {
                out.push((total(&s.l) - self.bottoms) / self.flow);
            }
        }
        out
    }

    /// Reboiler duty from the reboiler energy balance, W.
    fn reboiler_duty(&self, v: &[f64], props: &[BlockProperties]) -> f64 {
        let n = self.blocks;
        let last = self.decode(v, n - 1);
        let h_out = total(&last.l) * props[n - 1].first.enthalpy
            + total(&last.v) * props[n - 1].second.enthalpy;
        let h_in = if n == 2 {
            self.reflux_of(&self.decode(v, 0), &props[0]).1
        } else {
            total(&self.decode(v, n - 2).l) * props[n - 2].first.enthalpy
        };
        h_out - h_in - self.feed_enthalpy[n - 1]
    }

    /// Block-tridiagonal forward-difference Jacobian with three colors.
    #[allow(clippy::type_complexity)]
    fn jacobian(
        &self,
        v: &[f64],
        props: &[BlockProperties],
        r: &[f64],
    ) -> Result<(Vec<Matrix>, Vec<Matrix>, Vec<Matrix>), String> {
        let (n, w, m) = (self.blocks, self.width(), self.active.len());
        let zero = vec![vec![0.0; w]; w];
        let mut lower = vec![zero.clone(); n];
        let mut diag = vec![zero.clone(); n];
        let mut upper = vec![zero; n];
        for color in 0..3.min(n) {
            let group: Vec<usize> = (color..n).step_by(3).collect();
            for k in 0..w {
                let base = if k < 2 * m { 1e-7 } else { 1e-6 };
                let mut attempt = None;
                for h in [base, -base] {
                    let mut trial = v.to_vec();
                    for &j in &group {
                        trial[j * w + k] += h;
                    }
                    let mut perturbed = props.to_vec();
                    let ok = group.iter().try_for_each(|&j| {
                        perturbed[j] = self.properties(&trial, j)?;
                        Ok::<(), String>(())
                    });
                    if ok.is_ok() {
                        attempt = Some((h, self.residual(&trial, &perturbed)));
                        break;
                    }
                }
                let (h, shifted) =
                    attempt.ok_or("three-phase column Jacobian evaluation failed")?;
                for &j in &group {
                    for row_block in j.saturating_sub(1)..=(j + 1).min(n - 1) {
                        for row in 0..w {
                            let index = row_block * w + row;
                            let d = (shifted[index] - r[index]) / h;
                            if row_block + 1 == j {
                                upper[row_block][row][k] = d;
                            } else if row_block == j {
                                diag[j][row][k] = d;
                            } else {
                                lower[row_block][row][k] = d;
                            }
                        }
                    }
                }
            }
        }
        Ok((lower, diag, upper))
    }

    fn newton(&self, mut v: Vec<f64>) -> Result<(Vec<f64>, usize), String> {
        let (w, m) = (self.width(), self.active.len());
        let mut props = self
            .all_properties(&v)
            .map_err(|e| format!("three-phase column initialization: {e}"))?;
        let mut r = self.residual(&v, &props);
        for iteration in 0..=MAX_NEWTON {
            let error = max_abs(&r);
            if !error.is_finite() {
                return Err("three-phase column stage residuals are not finite".into());
            }
            if error <= TOLERANCE {
                return Ok((v, iteration));
            }
            if iteration == MAX_NEWTON {
                break;
            }
            let (lower, diag, upper) = self.jacobian(&v, &props, &r)?;
            let step = block_tridiagonal_solve(lower, diag, upper, r.iter().map(|x| -x).collect())
                .ok_or("three-phase column stage Jacobian is singular (vanishing phase or decanter plait point)")?;
            let mut limit = 1.0_f64;
            for (index, d) in step.iter().enumerate() {
                let cap = match (index / w, index % w) {
                    (_, k) if k < 2 * m => MAX_LN_STEP,
                    (0, _) => f64::INFINITY,
                    _ => MAX_T_STEP,
                };
                if d.abs() > cap {
                    limit = limit.min(cap / d.abs());
                }
            }
            let merit = sum_squares(&r);
            let mut accepted = None;
            let mut lambda = limit;
            for _ in 0..30 {
                let trial: Vec<f64> = v.iter().zip(&step).map(|(a, d)| a + lambda * d).collect();
                if let Ok(p) = self.all_properties(&trial) {
                    let rt = self.residual(&trial, &p);
                    if sum_squares(&rt) < merit {
                        accepted = Some((trial, p, rt));
                        break;
                    }
                }
                lambda *= 0.5;
            }
            let Some((trial, p, rt)) = accepted else {
                return Err(format!(
                    "three-phase column stage equations stalled (scaled residual {error:.3e}); change the reflux or bottoms specification"
                ));
            };
            v = trial;
            props = p;
            r = rt;
        }
        Err(format!(
            "three-phase column did not converge in {MAX_NEWTON} Newton iterations (scaled residual {:.3e}); check that the bottoms and reflux specification is feasible for the decanter liquids",
            max_abs(&r)
        ))
    }
}

type Matrix = Vec<Vec<f64>>;

/// K-values and the bubble function ln Σ K_i x_i of a single liquid `x` at T,
/// with the incipient vapor `y` (updated in place).
fn bubble_function(
    ctx: &CalcContext<'_>,
    x: &[f64],
    y: &mut Vec<f64>,
    t: f64,
    p: f64,
) -> Result<(f64, Vec<f64>), String> {
    let liquid = ctx.homogeneous(x, t, p, RootSelection::SmallestVolume)?;
    let mut k = Vec::new();
    for _ in 0..30 {
        let vapor = ctx.homogeneous(y, t, p, RootSelection::LargestVolume)?;
        k = liquid
            .ln_fugacity_coefficients
            .iter()
            .zip(&vapor.ln_fugacity_coefficients)
            .map(|(l, v)| (l - v).exp())
            .collect();
        let next: Vec<f64> = k.iter().zip(x).map(|(k, x)| k * x).collect();
        let next = normalized(&next);
        let change = max_difference(&next, y);
        *y = next;
        if change < 1e-12 {
            break;
        }
    }
    let sum: f64 = k.iter().zip(x).map(|(k, x)| k * x).sum();
    Ok((sum.ln(), k))
}

/// Bubble temperature of `x` treated as one liquid (no stability test): a
/// scan for a sign change of ln Σ K x followed by bisection.
fn single_liquid_bubble(
    ctx: &CalcContext<'_>,
    x: &[f64],
    p: f64,
    [lo, hi]: [f64; 2],
) -> Result<(f64, Vec<f64>), String> {
    let mut y = x.to_vec();
    let steps = 60;
    let mut previous: Option<(f64, f64)> = None;
    for s in 0..=steps {
        let t = lo + (hi - lo) * f64::from(s) / f64::from(steps);
        let Ok((g, _)) = bubble_function(ctx, x, &mut y, t, p) else {
            continue;
        };
        if let Some((t0, g0)) = previous
            && g0 < 0.0
            && g >= 0.0
        {
            let (mut a, mut b) = (t0, t);
            for _ in 0..60 {
                let mid = 0.5 * (a + b);
                let (gm, _) = bubble_function(ctx, x, &mut y, mid, p)?;
                if gm < 0.0 { a = mid } else { b = mid }
            }
            let t = 0.5 * (a + b);
            bubble_function(ctx, x, &mut y, t, p)?;
            return Ok((t, y));
        }
        previous = Some((t, g));
    }
    Err(format!(
        "three-phase column initialization: no bubble point of the feed in {lo:.2}..{hi:.2} K at {p:.0} Pa"
    ))
}

/// The light and heavy liquid of a VLLE flash of the condensate at the
/// decanter: explicit failures for vapor and for a single liquid.
#[allow(clippy::type_complexity)]
fn decanter_split(
    ctx: &CalcContext<'_>,
    z: &[f64],
    t: f64,
    p: f64,
) -> Result<Result<[(Vec<f64>, f64, f64); 2], String>, String> {
    let flash = ctx
        .vlle(z, t, p)
        .map_err(|e| format!("three-phase column needs the VLLE/LLE flash (decanter): {e}"))?;
    if flash.roles.contains(&PhaseRole::Vapor) {
        return Ok(Err(format!(
            "three-phase column: vapor forms in the decanter at {t:.3} K and {p:.0} Pa; the decanter temperature must be below the bubble point of the condensate"
        )));
    }
    let pick = |role| {
        flash.phase_with_role(role).map(|phase| {
            (
                phase.composition.clone(),
                phase.phase_fraction,
                phase.state.mass_density,
            )
        })
    };
    match (pick(PhaseRole::LightLiquid), pick(PhaseRole::HeavyLiquid)) {
        (Some(light), Some(heavy)) => Ok(Ok([light, heavy])),
        _ => Ok(Err(format!(
            "three-phase column: no liquid-liquid split in the decanter at {t:.3} K (the condensate is one liquid); use the distillation column"
        ))),
    }
}

fn product(
    ctx: &CalcContext<'_>,
    flow: f64,
    composition: &[f64],
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<MaterialStream, String> {
    let composition = composition.to_vec();
    let state = ctx.homogeneous(
        &composition,
        temperature_k,
        pressure_pa,
        RootSelection::SmallestVolume,
    )?;
    Ok(MaterialStream {
        molar_flow_mol_s: flow,
        state: Some(EquilibriumState {
            temperature_k,
            pressure_pa,
            mole_fractions: composition.clone(),
            phases: vec![EquilibriumPhase {
                phase: Phase::Liquid,
                phase_fraction: 1.0,
                composition,
                state,
            }],
            iterations: 0,
            convergence: Convergence {
                max_ln_fugacity_residual: 0.0,
                max_component_balance_residual: None,
                scaled_energy_residual: None,
                scaled_entropy_residual: None,
            },
            diagnostics: Vec::new(),
        }),
    })
}

/// Column temperature interval: the package interval intersected with its
/// bubble/dew interval.
fn temperature_range(ctx: &CalcContext<'_>) -> Result<[f64; 2], String> {
    let [vle_lo, vle_hi] = ctx
        .package()
        .saturation_temperature_interval()
        .map_err(|e| e.to_string())?;
    let [lo, hi] = ctx.package().temperature_interval();
    let (lo, hi) = (lo.max(vle_lo), hi.min(vle_hi));
    if lo >= hi {
        Err("package temperature interval does not overlap the column VLE interval".into())
    } else {
        Ok([lo, hi])
    }
}

pub(super) fn calculate(
    spec: &ThreePhaseColumnSpec,
    inputs: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, ThreePhaseColumnResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!(
            "three-phase column {} {}",
            issue.field, issue.problem
        ));
    }
    let (internal, bed) = resolve(spec.internals)?;
    let blocks = internal + 2;
    let pressure = spec.pressure_pa.expect("checked column pressure");
    let decanter_t = spec.decanter_temperature_k.expect("checked decanter T");
    let reflux = spec.reflux.fractions().expect("checked reflux");
    let components = ctx.component_count();
    let stream = &inputs[0];
    if !stream.molar_flow_mol_s.is_finite() || stream.molar_flow_mol_s <= 0.0 {
        return Err("three-phase column feed needs a positive finite flow".into());
    }
    let feed_block = usize::from(spec.feed_stage) - 1;
    let feed = StageFeed::from_stream(feed_block + 1, stream, components)?;
    if stream
        .state
        .as_ref()
        .is_some_and(|s| s.pressure_pa < pressure)
    {
        return Err(
            "three-phase column pressure exceeds the feed pressure; pressurize the feed upstream"
                .into(),
        );
    }
    let [plo, phi] = ctx.package().temperature_interval();
    if !(plo..=phi).contains(&decanter_t) {
        return Err(format!(
            "three-phase column decanter temperature {decanter_t} K is outside the package interval {plo:.2}..{phi:.2} K"
        ));
    }
    let range = temperature_range(ctx)?;
    let flow = feed.molar_flow_mol_s;
    let z = feed.z.clone();
    let active: Vec<usize> = (0..components).filter(|&i| z[i] > 0.0).collect();
    if active.len() < 2 {
        return Err("three-phase column needs at least two components in the feed".into());
    }
    // The decanter needs the package's three-phase flash; fail before any
    // column work when the package lacks it (GERG-2008, EOS-CG-2021).
    if let Err(CalculationError::UnsupportedCapability(message)) =
        ctx.package().vlle_tp(&z, decanter_t, pressure)
    {
        return Err(format!(
            "three-phase column needs the VLLE/LLE flash for its decanter: {message}"
        ));
    }
    let mut feed_flows = vec![vec![0.0; components]; blocks];
    let mut feed_enthalpy = vec![0.0; blocks];
    feed_flows[feed_block] = z.iter().map(|zi| zi * flow).collect();
    feed_enthalpy[feed_block] = flow * feed.molar_enthalpy_j_mol;
    let bottoms = flow * spec.bottoms_fraction.expect("checked bottoms fraction");
    let overhead = flow - bottoms;
    let column = Column {
        ctx,
        pressure,
        decanter_t,
        range,
        blocks,
        components,
        scale: active.iter().map(|&i| z[i] * flow).collect(),
        active,
        feed_flows,
        feed_enthalpy,
        flow,
        energy_scale: stream.thermal_scale_w().max(flow),
        bottoms,
        reflux,
    };

    // Initialization: constant molar overflow, bubble-point temperatures and
    // a lagged decanter flash of the top vapor.
    let (t0, y0) = single_liquid_bubble(ctx, &z, pressure, range)?;
    let stages = blocks - 1;
    let mut temperatures = vec![t0; stages];
    let mut x = vec![z.clone(); stages];
    let mut y = vec![y0; stages];
    let mut split = None;
    let mut failure = String::new();
    let mut vapor = 2.0 * overhead + bottoms;
    let mut reflux_flows: Vec<f64> = y[0].iter().map(|yi| yi * (vapor - overhead)).collect();
    for _ in 0..MAX_PASSES {
        // Decanter of the current top vapor.
        match decanter_split(ctx, &y[0], decanter_t, pressure)? {
            Ok([light, heavy]) => {
                let product_share = (1.0 - reflux[0]) * light.1 + (1.0 - reflux[1]) * heavy.1;
                vapor = (overhead / product_share.max(1e-3)).clamp(overhead * 1.05, 50.0 * flow);
                reflux_flows = (0..components)
                    .map(|i| {
                        vapor
                            * (reflux[0] * light.1 * light.0[i] + reflux[1] * heavy.1 * heavy.0[i])
                    })
                    .collect();
                split = Some([light, heavy]);
            }
            Err(reason) => {
                failure = reason;
                split = None;
                reflux_flows = y[0].iter().map(|yi| yi * (vapor - overhead)).collect();
            }
        }
        let reflux_total = total(&reflux_flows);
        let liquid: Vec<f64> = (0..stages)
            .map(|s| {
                if s + 1 == stages {
                    bottoms
                } else if s + 1 >= feed_block {
                    reflux_total + flow
                } else {
                    reflux_total
                }
            })
            .collect();
        // With a liquid feed, constant molar overflow gives V = R + D on
        // every stage, including the reboiler.
        let vapor_flows = vec![reflux_total + overhead; stages];
        let mut k_values = Vec::with_capacity(stages);
        for s in 0..stages {
            let (_, k) = bubble_function(ctx, &x[s], &mut y[s].clone(), temperatures[s], pressure)
                .map_err(|e| format!("three-phase column initialization: {e}"))?;
            k_values.push(k);
        }
        let mut next_x = vec![vec![0.0; components]; stages];
        for &i in &column.active {
            let factor: Vec<f64> = (0..stages)
                .map(|s| k_values[s][i] * vapor_flows[s] / liquid[s])
                .collect();
            let diag: Vec<f64> = factor.iter().map(|f| 1.0 + f).collect();
            let mut rhs: Vec<f64> = (0..stages).map(|s| column.feed_flows[s + 1][i]).collect();
            rhs[0] += reflux_flows[i];
            let l = thomas(&diag, &factor, &rhs)
                .ok_or("three-phase column initialization: singular component balance")?;
            for s in 0..stages {
                next_x[s][i] = l[s].max(1e-30);
            }
        }
        let mut change: f64 = 0.0;
        for s in 0..stages {
            let xs = normalized(&next_x[s]);
            // One damped Newton step on the bubble function per pass.
            let mut ys = y[s].clone();
            let (g, _) = bubble_function(ctx, &xs, &mut ys, temperatures[s], pressure)
                .map_err(|e| format!("three-phase column initialization: {e}"))?;
            let dt = 0.01;
            let mut y2 = ys.clone();
            let (g2, _) = bubble_function(ctx, &xs, &mut y2, temperatures[s] + dt, pressure)
                .map_err(|e| format!("three-phase column initialization: {e}"))?;
            let slope = (g2 - g) / dt;
            let step = if slope > 0.0 {
                (-g / slope).clamp(-10.0, 10.0)
            } else {
                0.0
            };
            let t = (temperatures[s] + step).clamp(range[0], range[1]);
            change = change
                .max((t - temperatures[s]).abs() / 10.0)
                .max(max_difference(&xs, &x[s]));
            temperatures[s] = t;
            x[s] = xs;
            y[s] = ys;
        }
        if change < 1e-6 && split.is_some() {
            break;
        }
    }
    let Some([light, heavy]) = split else {
        return Err(format!("{failure} (initial column profile)"));
    };

    // Newton vector from the initial profile.
    let floor = 1e-14 * flow;
    let keep = |f: f64, i: usize| if z[i] > 0.0 { f.max(floor) } else { 0.0 };
    let reflux_total = total(&reflux_flows);
    let mut initial = Vec::with_capacity(blocks);
    initial.push(Block {
        l: (0..components)
            .map(|i| keep(vapor * light.1 * light.0[i], i))
            .collect(),
        v: (0..components)
            .map(|i| keep(vapor * heavy.1 * heavy.0[i], i))
            .collect(),
        t: 0.0,
    });
    for s in 0..stages {
        let liquid = if s + 1 == stages {
            bottoms
        } else if s + 1 >= feed_block {
            reflux_total + flow
        } else {
            reflux_total
        };
        let vapor_flow = reflux_total + overhead;
        initial.push(Block {
            l: (0..components).map(|i| keep(liquid * x[s][i], i)).collect(),
            v: (0..components)
                .map(|i| keep(vapor_flow * y[s][i], i))
                .collect(),
            t: temperatures[s],
        });
    }
    // Condenser duty from the decanter energy balance of the initial profile.
    let mut v0 = column.encode(&initial);
    if let Ok(props) = column.all_properties(&v0) {
        let w = column.width();
        let r = column.residual(&v0, &props);
        v0[w - 1] = r[w - 1];
    }
    let (solved, iterations) = column.newton(v0)?;
    let props = column.all_properties(&solved)?;
    let decoded: Vec<Block> = (0..blocks).map(|j| column.decode(&solved, j)).collect();

    // Verification: decanter against an independent LLE flash.
    let mut warnings = Vec::new();
    let decanter = &decoded[0];
    let condensate = normalized(&decoded[1].v);
    let (xa, xb) = (normalized(&decanter.l), normalized(&decanter.v));
    let distance = max_difference(&xa, &xb);
    if distance < DECANTER_PLAIT_DISTANCE {
        return Err(format!(
            "three-phase column: the decanter liquids coincide (max mole-fraction difference {distance:.2e}); the condensate reaches a plait point or the trivial solution"
        ));
    }
    let [check_light, check_heavy] = decanter_split(ctx, &condensate, decanter_t, pressure)??;
    let deviation = max_difference(&check_light.0, &xa).max(max_difference(&check_heavy.0, &xb));
    if deviation > FLASH_AGREEMENT {
        let swapped = max_difference(&check_light.0, &xb).max(max_difference(&check_heavy.0, &xa));
        return Err(if swapped <= FLASH_AGREEMENT {
            "three-phase column: the decanter liquids exchanged their density order during the solve; specify the other reflux phase".into()
        } else {
            format!(
                "three-phase column: the solved decanter split is not the stable liquid-liquid split of the condensate (deviation {deviation:.2e} from the VLLE flash)"
            )
        });
    }
    // Stages: exactly one vapor and one liquid, matching the solved tie line.
    for (j, b) in decoded.iter().enumerate().skip(1) {
        let stage = j + 1;
        let (xs, ys) = (normalized(&b.l), normalized(&b.v));
        let midpoint: Vec<f64> = xs.iter().zip(&ys).map(|(x, y)| 0.5 * (x + y)).collect();
        let flash = ctx
            .vlle(&midpoint, b.t, pressure)
            .map_err(|e| format!("three-phase column stage {stage} phase check: {e}"))?;
        if flash
            .roles
            .iter()
            .any(|r| matches!(r, PhaseRole::LightLiquid | PhaseRole::HeavyLiquid))
        {
            return Err(format!(
                "three-phase column: stage {stage} liquid splits into two liquids; this column allows two liquids only in the decanter (internal three-phase stages are not supported)"
            ));
        }
        let matches = match (
            flash.phase_with_role(PhaseRole::Vapor),
            flash.phase_with_role(PhaseRole::Liquid),
        ) {
            (Some(v), Some(l)) => {
                (v.phase_fraction - 0.5).abs() <= 2e-5
                    && max_difference(&v.composition, &ys) <= 2e-6
                    && max_difference(&l.composition, &xs) <= 2e-6
            }
            _ => false,
        };
        if !matches {
            return Err(format!(
                "three-phase column: stage {stage} did not reproduce a stable vapor-liquid split"
            ));
        }
    }
    let light_flow = total(&decanter.l);
    let heavy_flow = total(&decanter.v);
    if props[0].first.mass_density > props[0].second.mass_density {
        return Err(
            "three-phase column: the decanter liquids exchanged their density order during the solve; specify the other reflux phase".into(),
        );
    }
    if distance < 0.02 {
        warnings.push(format!(
            "TPC-NEAR-PLAIT-POINT: decanter liquids differ by only {distance:.4} in mole fraction"
        ));
    }
    let condenser = solved[column.width() - 1] * column.energy_scale;
    let reboiler = column.reboiler_duty(&solved, &props);
    if condenser <= 0.0 || reboiler <= 0.0 {
        return Err(format!(
            "three-phase column requires heat removal at the condenser and heat input at the reboiler (condenser {condenser:.3e} W removed, reboiler {reboiler:.3e} W supplied); the specification is infeasible"
        ));
    }
    let [fa, fb] = reflux;
    let (distillate, decant) = match spec.reflux {
        DecanterReflux::Light(_) => (((1.0 - fa) * light_flow, &xa), (heavy_flow, &xb)),
        DecanterReflux::Heavy(_) => (((1.0 - fb) * heavy_flow, &xb), (light_flow, &xa)),
        DecanterReflux::Fractions { .. } => (
            ((1.0 - fa) * light_flow, &xa),
            ((1.0 - fb) * heavy_flow, &xb),
        ),
    };
    let last = &decoded[blocks - 1];
    let outlets = vec![
        product(ctx, distillate.0, distillate.1, decanter_t, pressure)?,
        product(ctx, decant.0, decant.1, decanter_t, pressure)?,
        product(ctx, total(&last.l), &normalized(&last.l), last.t, pressure)?,
    ];
    let (reflux_flows, _) = column.reflux_of(decanter, &props[0]);
    let mut liquid_flows = vec![total(&reflux_flows)];
    let mut vapor_flows = vec![0.0];
    let mut liquid_compositions = vec![normalized(&reflux_flows)];
    let mut vapor_compositions = vec![condensate];
    let mut temperatures = vec![decanter_t];
    for b in &decoded[1..] {
        liquid_flows.push(total(&b.l));
        vapor_flows.push(total(&b.v));
        liquid_compositions.push(normalized(&b.l));
        vapor_compositions.push(normalized(&b.v));
        temperatures.push(b.t);
    }
    let mut output = UnitOutput::new(outlets, reboiler - condenser, warnings);
    output.energy = vec![Some(condenser), Some(reboiler)];
    Ok((
        output,
        ThreePhaseColumnResult {
            stage_temperatures_k: temperatures,
            liquid_flows_mol_s: liquid_flows,
            vapor_flows_mol_s: vapor_flows,
            liquid_compositions,
            vapor_compositions,
            light_liquid_mol_s: light_flow,
            light_liquid_composition: xa,
            heavy_liquid_mol_s: heavy_flow,
            heavy_liquid_composition: xb,
            light_refluxed: fa > 0.0,
            condenser_duty_w: condenser,
            reboiler_duty_w: reboiler,
            packed_height_m: bed.map(|b| b.packed_height_m),
            iterations,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflux_choices_map_to_fractions() {
        assert_eq!(
            DecanterReflux::Light(PhaseReflux::Total).fractions(),
            Some([1.0, 0.0])
        );
        assert_eq!(
            DecanterReflux::Heavy(PhaseReflux::Ratio(Some(3.0))).fractions(),
            Some([0.0, 0.75])
        );
        assert_eq!(
            DecanterReflux::Fractions {
                light: Some(0.5),
                heavy: Some(0.2)
            }
            .fractions(),
            Some([0.5, 0.2])
        );
        assert_eq!(
            DecanterReflux::Light(PhaseReflux::Ratio(None)).fractions(),
            None
        );
    }

    #[test]
    fn specification_issues_are_explicit() {
        let fields: Vec<_> = ThreePhaseColumnSpec::default()
            .issues()
            .iter()
            .map(|i| i.field)
            .collect();
        assert_eq!(
            fields,
            [
                "column pressure",
                "decanter temperature",
                "bottoms fraction"
            ]
        );
        let spec = ThreePhaseColumnSpec {
            internals: ColumnInternals::Trays(3),
            feed_stage: 6,
            pressure_pa: Some(1e5),
            decanter_temperature_k: Some(320.0),
            reflux: DecanterReflux::Fractions {
                light: Some(0.0),
                heavy: Some(0.0),
            },
            bottoms_fraction: Some(1.0),
        };
        let fields: Vec<_> = spec.issues().iter().map(|i| i.field).collect();
        assert_eq!(
            fields,
            ["reflux fractions", "bottoms fraction", "feed stage"]
        );
        let spec = ThreePhaseColumnSpec {
            internals: ColumnInternals::Trays(0),
            ..Default::default()
        };
        assert!(
            spec.issues()
                .iter()
                .any(|i| i.field == "internals" && i.problem == SpecProblem::CountOutOfRange)
        );
        let spec = ThreePhaseColumnSpec {
            reflux: DecanterReflux::Light(PhaseReflux::Ratio(Some(-1.0))),
            ..Default::default()
        };
        assert!(spec.issues().iter().any(|i| i.field == "reflux ratio"));
    }
}
