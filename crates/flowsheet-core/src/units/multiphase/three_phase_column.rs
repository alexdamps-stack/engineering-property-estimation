//! Three-phase distillation column (#188): equilibrium stages at one
//! pressure, a total condenser with a decanter and a partial reboiler
//! (Aspen Plus RadFrac with a decanter and free-water/three-phase stages /
//! HYSYS Three Phase Distillation).
//!
//! Supported scope (`docs/a17-multiphase-three-phase-column.md`):
//!
//! - Stage 1 is the total condenser with its decanter at a specified
//!   (subcooled) temperature. The vapor from stage 2 condenses completely and
//!   splits into a light and a heavy liquid (liquid–liquid equilibrium at the
//!   decanter temperature and column pressure).
//! - Stages 2..N−1 are ideal trays (or ideal packed stages) and stage N is the
//!   partial reboiler. Each of them carries one liquid (vapor–liquid
//!   equilibrium) or two liquids (vapor–liquid–liquid equilibrium). The set
//!   of three-phase stages is found automatically: after each converged pass
//!   every single-liquid stage is flashed with the package's three-phase
//!   operation (a split adds its second liquid), every three-phase stage is
//!   checked against the package's VLLE PH flash, and a second liquid that
//!   vanishes or merges with the first during the iteration is removed. A
//!   phase set that recurs is an explicit phase-count oscillation failure.
//! - [`DecanterReflux`] chooses which liquid returns as reflux and how much.
//!   The `distillate` port carries the non-refluxed part of the reflux liquid
//!   (zero for total reflux of that phase); the `decant` port carries the
//!   other liquid. With [`DecanterReflux::Fractions`] the distillate is the
//!   non-refluxed light liquid and the decant the non-refluxed heavy liquid.
//!
//! Stage equations (Seader, Henley and Roper, *Separation Process
//! Principles*, 3rd ed., 2011, §10.1; three-phase stages as in Cairns and
//! Furzer, Ind. Eng. Chem. Res. 29 (1990) 1349, in a Naphtali–Sandholm
//! (AIChE J. 17 (1971) 148) simultaneous-correction form; in − out):
//!
//! - decanter: v_{2,i} − a_i − b_i = 0, ln x^a_i + ln φ^a_i − ln x^b_i −
//!   ln φ^b_i = 0, V_2 h^V_2 − Q_C − A h^a − B h^b = 0 (a light, b heavy
//!   liquid at the decanter temperature);
//! - stage j ≥ 2: l_{j−1,i} + v_{j+1,i} + f_{j,i} − l_{j,i} − l'_{j,i} −
//!   v_{j,i} = 0 with l_{1,i} = φ_L a_i + φ_H b_i (reflux fractions) and
//!   l_{j−1} the total liquid of the stage above, ln y_ji − ln x_ji −
//!   ln K_ji = 0 (K = φ^L/φ^V of the package's liquid and vapor at the stage
//!   state), on a three-phase stage also ln x_ji + ln φ_ji − ln x'_ji −
//!   ln φ'_ji = 0 for its second liquid l', and the energy balance on stages
//!   2..N−1;
//! - the reboiler replaces its energy balance by the bottoms flow
//!   specification L_N + L'_N = B; its duty follows from its energy balance.
//!
//! Unknowns: ln component flows of both decanter liquids and the scaled
//! condenser duty; ln liquid, vapor and second-liquid component flows and
//! temperature of every other stage (the second-liquid variables of a
//! single-liquid stage are inert). The solve starts from
//! constant-molar-overflow bubble-point passes with a lagged decanter flash
//! and finishes with a damped Newton iteration on the block-tridiagonal stage
//! Jacobian.

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
/// Largest deviation of a three-phase stage from the package VLLE PH flash
/// of its contents (temperature in K, mole fractions, phase fractions / 10).
const THREE_PHASE_AGREEMENT: f64 = 1e-5;
/// A second stage liquid below this share of the stage liquid vanishes.
const VANISHING_LIQUID_SHARE: f64 = 1e-9;
/// Two stage liquids closer than this (max mole-fraction difference) merge.
const MERGE_DISTANCE: f64 = 1e-4;
/// Phase-set rounds before an explicit failure.
const MAX_PHASE_ROUNDS: usize = 12;

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

/// Two liquids on one stage. The light liquid is the less dense one.
#[derive(Debug, Clone, PartialEq)]
pub struct StageLiquidSplit {
    pub light_mol_s: f64,
    pub light_composition: Vec<f64>,
    pub heavy_mol_s: f64,
    pub heavy_composition: Vec<f64>,
}

/// Solved profile of a three-phase column. Stage vectors run from the
/// decanter (stage 1) to the reboiler (stage N).
#[derive(Debug, Clone, PartialEq)]
pub struct ThreePhaseColumnResult {
    /// Stage temperatures; stage 1 is the decanter temperature.
    pub stage_temperatures_k: Vec<f64>,
    /// Liquid leaving each stage downward (both liquids on a three-phase
    /// stage), mol/s; stage 1: the reflux.
    pub liquid_flows_mol_s: Vec<f64>,
    /// Vapor leaving each stage upward, mol/s; stage 1: zero.
    pub vapor_flows_mol_s: Vec<f64>,
    /// Liquid mole fractions (total liquid on a three-phase stage); stage 1:
    /// the reflux composition.
    pub liquid_compositions: Vec<Vec<f64>>,
    /// Vapor mole fractions; stage 1: the condensate (vapor from stage 2).
    pub vapor_compositions: Vec<Vec<f64>>,
    /// Phases on each stage: 2 on the decanter (two liquids), 2 on a
    /// vapor–liquid stage, 3 on a vapor–liquid–liquid stage.
    pub stage_phase_counts: Vec<u8>,
    /// Both liquids of every stage that carries two: stage 1 the decanter,
    /// `None` on a single-liquid stage.
    pub stage_liquid_splits: Vec<Option<StageLiquidSplit>>,
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
    /// Accepted Newton iterations (all phase-set rounds).
    pub iterations: usize,
}

/// Properties of one block: the two decanter liquids, or a stage's liquid
/// and vapor plus its second liquid on a three-phase stage.
#[derive(Debug, Clone)]
struct BlockProperties {
    first: PhaseProperties,
    second: PhaseProperties,
    third: Option<PhaseProperties>,
}

/// Decoded block: decanter (`l` light, `v` heavy, `t` = scaled condenser
/// duty) or stage (`l` liquid in vapor–liquid equilibrium, `v` vapor,
/// `l2` second liquid (zero on a single-liquid stage), `t` temperature).
#[derive(Debug, Clone)]
struct Block {
    l: Vec<f64>,
    v: Vec<f64>,
    l2: Vec<f64>,
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

/// Why a Newton pass ended without a solution.
enum NewtonFailure {
    Failed(String),
    /// A second liquid vanished or merged with the first on this block: the
    /// stage returns to one liquid and the column is solved again.
    Merge {
        block: usize,
        v: Vec<f64>,
        iterations: usize,
    },
}

/// Phase-set change of one stage found after a converged pass.
enum PhaseChange {
    /// The stage liquid splits: light and heavy liquid of the check flash.
    Split(usize, Vec<f64>, Vec<f64>),
    Merge(usize),
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
    /// Blocks carrying a second liquid (never the decanter).
    split: Vec<bool>,
}

impl Column<'_, '_> {
    fn width(&self) -> usize {
        3 * self.active.len() + 1
    }

    fn decode(&self, v: &[f64], j: usize) -> Block {
        let (w, m) = (self.width(), self.active.len());
        let block = &v[j * w..(j + 1) * w];
        let mut l = vec![0.0; self.components];
        let mut u = vec![0.0; self.components];
        let mut l2 = vec![0.0; self.components];
        for (k, &i) in self.active.iter().enumerate() {
            l[i] = block[k].exp();
            u[i] = block[m + k].exp();
            if self.split[j] {
                l2[i] = block[2 * m + k].exp();
            }
        }
        Block {
            l,
            v: u,
            l2,
            t: block[3 * m],
        }
    }

    fn encode(&self, blocks: &[Block]) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.blocks * self.width());
        for (j, b) in blocks.iter().enumerate() {
            out.extend(self.active.iter().map(|&i| b.l[i].ln()));
            out.extend(self.active.iter().map(|&i| b.v[i].ln()));
            if self.split[j] {
                out.extend(self.active.iter().map(|&i| b.l2[i].ln()));
            } else {
                out.extend(self.active.iter().map(|_| 0.0));
            }
            out.push(b.t);
        }
        out
    }

    fn decode_all(&self, v: &[f64]) -> Vec<Block> {
        (0..self.blocks).map(|j| self.decode(v, j)).collect()
    }

    fn properties(&self, v: &[f64], j: usize) -> Result<BlockProperties, String> {
        let b = self.decode(v, j);
        if b.l.iter().chain(&b.v).chain(&b.l2).any(|f| !f.is_finite()) || !b.t.is_finite() {
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
                third: None,
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
        let third = if self.split[j] {
            Some(liquid(&b.l2, b.t)?)
        } else {
            None
        };
        Ok(BlockProperties {
            first,
            second,
            third,
        })
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

    /// Liquid leaving a stage downward (both liquids): component flows and
    /// enthalpy flow.
    fn liquid_out(b: &Block, p: &BlockProperties) -> (Vec<f64>, f64) {
        let flows = b.l.iter().zip(&b.l2).map(|(a, c)| a + c).collect();
        let h = total(&b.l) * p.first.enthalpy
            + p.third
                .as_ref()
                .map_or(0.0, |third| total(&b.l2) * third.enthalpy);
        (flows, h)
    }

    fn residual(&self, v: &[f64], props: &[BlockProperties]) -> Vec<f64> {
        let (n, w, m) = (self.blocks, self.width(), self.active.len());
        let blocks = self.decode_all(v);
        let mut out = Vec::with_capacity(n * w);
        // Rows of an unused second liquid: its (zero) variables themselves.
        let dummy = |out: &mut Vec<f64>, j: usize| {
            out.extend((0..m).map(|k| v[j * w + 2 * m + k]));
        };
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
            dummy(&mut out, 0);
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
                Self::liquid_out(&blocks[j - 1], &props[j - 1])
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
                    (l_in[i] + v_in[i] + self.feed_flows[j][i] - s.l[i] - s.l2[i] - s.v[i])
                        / self.scale[k],
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
            match &p.third {
                Some(third) => {
                    let ln_l2 = total(&s.l2).ln();
                    for &i in &self.active {
                        out.push(
                            s.l[i].ln() - ln_l + p.first.ln_fugacity_coefficients[i] - s.l2[i].ln()
                                + ln_l2
                                - third.ln_fugacity_coefficients[i],
                        );
                    }
                }
                None => dummy(&mut out, j),
            }
            if j + 1 < n {
                let h_in = h_l_in + h_v_in + self.feed_enthalpy[j];
                let h_out = Self::liquid_out(s, p).1 + total(&s.v) * p.second.enthalpy;
                out.push((h_in - h_out) / self.energy_scale);
            } else {
                out.push((total(&s.l) + total(&s.l2) - self.bottoms) / self.flow);
            }
        }
        out
    }

    /// Reboiler duty from the reboiler energy balance, W.
    fn reboiler_duty(&self, v: &[f64], props: &[BlockProperties]) -> f64 {
        let n = self.blocks;
        let last = self.decode(v, n - 1);
        let h_out = Self::liquid_out(&last, &props[n - 1]).1
            + total(&last.v) * props[n - 1].second.enthalpy;
        let h_in = if n == 2 {
            self.reflux_of(&self.decode(v, 0), &props[0]).1
        } else {
            Self::liquid_out(&self.decode(v, n - 2), &props[n - 2]).1
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
                let base = if k < 3 * m { 1e-7 } else { 1e-6 };
                // A second-liquid variable changes the properties only on a
                // three-phase stage.
                let second_liquid = (2 * m..3 * m).contains(&k);
                let mut attempt = None;
                for h in [base, -base] {
                    let mut trial = v.to_vec();
                    for &j in &group {
                        trial[j * w + k] += h;
                    }
                    let mut perturbed = props.to_vec();
                    let ok = group.iter().try_for_each(|&j| {
                        if !second_liquid || self.split[j] {
                            perturbed[j] = self.properties(&trial, j)?;
                        }
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

    /// A three-phase block whose second liquid vanished or coincides with
    /// the first.
    fn degenerate_split(&self, b: &Block) -> bool {
        let (a, c) = (total(&b.l), total(&b.l2));
        a.min(c) < VANISHING_LIQUID_SHARE * (a + c)
            || max_difference(&normalized(&b.l), &normalized(&b.l2)) < MERGE_DISTANCE
    }

    fn newton(&self, mut v: Vec<f64>) -> Result<(Vec<f64>, usize), NewtonFailure> {
        let failed = NewtonFailure::Failed;
        let (w, m) = (self.width(), self.active.len());
        let mut props = self
            .all_properties(&v)
            .map_err(|e| failed(format!("three-phase column initialization: {e}")))?;
        let mut r = self.residual(&v, &props);
        for iteration in 0..=MAX_NEWTON {
            let error = max_abs(&r);
            if !error.is_finite() {
                return Err(failed(
                    "three-phase column stage residuals are not finite".into(),
                ));
            }
            if error <= TOLERANCE {
                return Ok((v, iteration));
            }
            if iteration == MAX_NEWTON {
                break;
            }
            let (lower, diag, upper) = self.jacobian(&v, &props, &r).map_err(failed)?;
            let step = block_tridiagonal_solve(lower, diag, upper, r.iter().map(|x| -x).collect())
                .ok_or_else(|| failed("three-phase column stage Jacobian is singular (vanishing phase or decanter plait point)".into()))?;
            let mut limit = 1.0_f64;
            for (index, d) in step.iter().enumerate() {
                let cap = match (index / w, index % w) {
                    (_, k) if k < 3 * m => MAX_LN_STEP,
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
                return Err(failed(format!(
                    "three-phase column stage equations stalled (scaled residual {error:.3e}); change the reflux or bottoms specification"
                )));
            };
            v = trial;
            props = p;
            r = rt;
            if let Some(block) = (1..self.blocks)
                .find(|&j| self.split[j] && self.degenerate_split(&self.decode(&v, j)))
            {
                return Err(NewtonFailure::Merge {
                    block,
                    v,
                    iterations: iteration + 1,
                });
            }
        }
        Err(failed(format!(
            "three-phase column did not converge in {MAX_NEWTON} Newton iterations (scaled residual {:.3e}); check that the bottoms and reflux specification is feasible for the decanter liquids",
            max_abs(&r)
        )))
    }

    /// Returns block `j` to one liquid (both liquids combined).
    fn merge(&mut self, v: &[f64], j: usize) -> Vec<f64> {
        let mut blocks = self.decode_all(v);
        let b = &mut blocks[j];
        b.l = b.l.iter().zip(&b.l2).map(|(a, c)| a + c).collect();
        b.l2 = vec![0.0; self.components];
        self.split[j] = false;
        self.encode(&blocks)
    }

    /// Splits the liquid of block `j` onto the tie line `light`–`heavy` of a
    /// check flash (lever rule of the current liquid, kept inside 5..95 %).
    fn split_block(&mut self, v: &[f64], j: usize, light: &[f64], heavy: &[f64]) -> Vec<f64> {
        let mut blocks = self.decode_all(v);
        let floor = 1e-14 * self.flow;
        let b = &mut blocks[j];
        let liquid = total(&b.l);
        let x = normalized(&b.l);
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..self.components {
            let d = heavy[i] - light[i];
            num += (x[i] - light[i]) * d;
            den += d * d;
        }
        let beta = if den > 0.0 {
            (num / den).clamp(0.05, 0.95)
        } else {
            0.5
        };
        for &i in &self.active {
            b.l[i] = (liquid * (1.0 - beta) * light[i]).max(floor);
            b.l2[i] = (liquid * beta * heavy[i]).max(floor);
        }
        self.split[j] = true;
        self.encode(&blocks)
    }

    /// Light and heavy liquid of a VLLE flash of `z` at the stage state, if
    /// it has both; flash errors count as no split.
    fn liquid_split(&self, z: &[f64], t: f64) -> Option<(Vec<f64>, Vec<f64>)> {
        let flash = self.ctx.vlle(z, t, self.pressure).ok()?;
        let light = flash.phase_with_role(PhaseRole::LightLiquid)?;
        let heavy = flash.phase_with_role(PhaseRole::HeavyLiquid)?;
        Some((light.composition.clone(), heavy.composition.clone()))
    }

    /// Liquid–liquid tangent-plane test of the stage liquid `x` (Michelsen,
    /// Fluid Phase Equilib. 9 (1982) 1) with the package liquid fugacity
    /// coefficients, from a near-pure trial of every component: the
    /// incipient liquid of the deepest negative stationary point. Used when
    /// the package flash of the stage fails (close to a binary three-phase
    /// temperature).
    fn incipient_liquid(&self, x: &[f64], t: f64) -> Option<Vec<f64>> {
        let ln_phi = |w: &[f64]| {
            self.ctx
                .homogeneous(w, t, self.pressure, RootSelection::SmallestVolume)
                .ok()
                .map(|p| p.ln_fugacity_coefficients)
        };
        let parent = ln_phi(x)?;
        let d: Vec<f64> = (0..self.components)
            .map(|i| {
                if x[i] > 0.0 {
                    x[i].ln() + parent[i]
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect();
        let mut best: Option<(f64, Vec<f64>)> = None;
        for &dominant in &self.active {
            let mut w: Vec<f64> = (0..self.components)
                .map(|i| match (i == dominant, x[i] > 0.0) {
                    (true, _) => 1.0,
                    (false, true) => 1e-3,
                    (false, false) => 0.0,
                })
                .collect();
            w = normalized(&w);
            let mut tm = 0.0;
            for _ in 0..300 {
                let Some(phi) = ln_phi(&w) else { break };
                let big: Vec<f64> = (0..self.components)
                    .map(|i| {
                        if x[i] > 0.0 {
                            (d[i] - phi[i]).exp()
                        } else {
                            0.0
                        }
                    })
                    .collect();
                tm = 1.0 - total(&big);
                let next = normalized(&big);
                let change = max_difference(&next, &w);
                w = next;
                if change < 1e-12 {
                    break;
                }
            }
            if tm < -1e-8
                && max_difference(&w, x) > 1e-3
                && best.as_ref().is_none_or(|(b, _)| tm < *b)
            {
                best = Some((tm, w));
            }
        }
        best.map(|(_, w)| w)
    }

    /// Stages whose liquid splits in a VLLE flash of the tie-line midpoint
    /// or of the stage liquid (light and heavy liquid of that flash).
    fn splitting_stages(&self, v: &[f64]) -> Vec<(usize, Vec<f64>, Vec<f64>)> {
        let mut out = Vec::new();
        for j in 1..self.blocks {
            let b = self.decode(v, j);
            let (xs, ys) = (normalized(&b.l), normalized(&b.v));
            let midpoint: Vec<f64> = xs.iter().zip(&ys).map(|(x, y)| 0.5 * (x + y)).collect();
            if let Some((light, heavy)) = self
                .liquid_split(&midpoint, b.t)
                .or_else(|| self.liquid_split(&xs, b.t))
                .or_else(|| self.incipient_liquid(&xs, b.t).map(|w| (xs.clone(), w)))
            {
                out.push((j, light, heavy));
            }
        }
        out
    }

    /// Phase-set check of a converged pass: single-liquid stages against a
    /// VLLE flash of their tie-line midpoint, three-phase stages against a
    /// VLLE PH flash of their contents.
    fn phase_changes(
        &self,
        v: &[f64],
        props: &[BlockProperties],
    ) -> Result<Vec<PhaseChange>, String> {
        let mut changes = Vec::new();
        for (j, p) in props.iter().enumerate().skip(1) {
            let stage = j + 1;
            let b = self.decode(v, j);
            let (xs, ys) = (normalized(&b.l), normalized(&b.v));
            if self.split[j] {
                if self.degenerate_split(&b) {
                    changes.push(PhaseChange::Merge(j));
                    continue;
                }
                let third = p.third.as_ref().expect("three-phase stage properties");
                let (la, lb, lv) = (total(&b.l), total(&b.l2), total(&b.v));
                let f = la + lb + lv;
                let z: Vec<f64> = (0..self.components)
                    .map(|i| (b.l[i] + b.l2[i] + b.v[i]) / f)
                    .collect();
                let h = (la * p.first.enthalpy + lb * third.enthalpy + lv * p.second.enthalpy) / f;
                let x2 = normalized(&b.l2);
                let flash = self
                    .ctx
                    .vlle_ph(&z, self.pressure, h)
                    .map_err(|e| format!("three-phase column stage {stage} phase check: {e}"))?;
                let phases = (
                    flash.phase_with_role(PhaseRole::Vapor),
                    flash.phase_with_role(PhaseRole::LightLiquid),
                    flash.phase_with_role(PhaseRole::HeavyLiquid),
                );
                let (Some(vap), Some(light), Some(heavy)) = phases else {
                    // The package finds fewer phases: the second liquid is
                    // not stable, solve the stage with one liquid.
                    changes.push(PhaseChange::Merge(j));
                    continue;
                };
                let liquids = |a: &[f64], c: &[f64], fa: f64, fc: f64| {
                    max_difference(&light.composition, a)
                        .max(max_difference(&heavy.composition, c))
                        .max((light.phase_fraction - fa).abs() / 10.0)
                        .max((heavy.phase_fraction - fc).abs() / 10.0)
                };
                let deviation = (flash.state.temperature_k - b.t)
                    .abs()
                    .max(max_difference(&vap.composition, &ys))
                    .max((vap.phase_fraction - lv / f).abs() / 10.0)
                    .max(liquids(&xs, &x2, la / f, lb / f).min(liquids(&x2, &xs, lb / f, la / f)));
                if deviation > THREE_PHASE_AGREEMENT {
                    return Err(format!(
                        "three-phase column: stage {stage} did not reproduce a stable vapor-liquid-liquid split (deviation {deviation:.2e} from the VLLE PH flash)"
                    ));
                }
                continue;
            }
            let midpoint: Vec<f64> = xs.iter().zip(&ys).map(|(x, y)| 0.5 * (x + y)).collect();
            let flash = match self.ctx.vlle(&midpoint, b.t, self.pressure) {
                Ok(flash) => flash,
                Err(e) => {
                    // Close to a binary three-phase temperature the midpoint
                    // flash can fail; a second liquid is then still found
                    // from the stage liquid alone.
                    if let Some((light, heavy)) = self.liquid_split(&xs, b.t) {
                        changes.push(PhaseChange::Split(j, light, heavy));
                        continue;
                    }
                    if let Some(w) = self.incipient_liquid(&xs, b.t) {
                        changes.push(PhaseChange::Split(j, xs, w));
                        continue;
                    }
                    return Err(format!("three-phase column stage {stage} phase check: {e}"));
                }
            };
            if let (Some(a), Some(c)) = (
                flash.phase_with_role(PhaseRole::LightLiquid),
                flash.phase_with_role(PhaseRole::HeavyLiquid),
            ) {
                changes.push(PhaseChange::Split(
                    j,
                    a.composition.clone(),
                    c.composition.clone(),
                ));
                continue;
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
        Ok(changes)
    }

    fn three_phase_stages(&self) -> String {
        let stages: Vec<String> = (1..self.blocks)
            .filter(|&j| self.split[j])
            .map(|j| (j + 1).to_string())
            .collect();
        if stages.is_empty() {
            String::new()
        } else {
            format!(" (two liquids on stages {})", stages.join(", "))
        }
    }
}

type Matrix = Vec<Vec<f64>>;

/// Records a phase set; a set seen before is a phase-count oscillation.
fn remember(seen: &mut Vec<Vec<bool>>, split: &[bool]) -> Result<(), String> {
    if seen.iter().any(|s| s == split) {
        let stages: Vec<String> = split
            .iter()
            .enumerate()
            .filter(|(_, s)| **s)
            .map(|(j, _)| (j + 1).to_string())
            .collect();
        return Err(format!(
            "three-phase column: the stage phase counts oscillate (the set with two liquids on stages [{}] recurs); change the reflux or bottoms specification",
            stages.join(", ")
        ));
    }
    seen.push(split.to_vec());
    Ok(())
}

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

/// A liquid product: one liquid, or (bottoms of a three-phase reboiler) two
/// liquids given as (fraction, composition, properties).
fn product(
    ctx: &CalcContext<'_>,
    flow: f64,
    composition: &[f64],
    temperature_k: f64,
    pressure_pa: f64,
    liquids: Option<[(f64, Vec<f64>, PhaseProperties); 2]>,
) -> Result<MaterialStream, String> {
    let composition = composition.to_vec();
    let phases = match liquids {
        Some(pair) => pair
            .into_iter()
            .map(|(phase_fraction, composition, state)| EquilibriumPhase {
                phase: Phase::Liquid,
                phase_fraction,
                composition,
                state,
            })
            .collect(),
        None => vec![EquilibriumPhase {
            phase: Phase::Liquid,
            phase_fraction: 1.0,
            composition: composition.clone(),
            state: ctx.homogeneous(
                &composition,
                temperature_k,
                pressure_pa,
                RootSelection::SmallestVolume,
            )?,
        }],
    };
    Ok(MaterialStream {
        molar_flow_mol_s: flow,
        state: Some(EquilibriumState {
            temperature_k,
            pressure_pa,
            mole_fractions: composition,
            phases,
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
    let mut column = Column {
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
        split: vec![false; blocks],
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
        l2: vec![0.0; components],
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
            l2: vec![0.0; components],
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

    // Phase-set rounds: solve with the current set of three-phase stages,
    // then split stages whose liquid is unstable and merge stages whose
    // second liquid vanished, until the set no longer changes.
    let mut seen = vec![column.split.clone()];
    let mut start = v0.clone();
    let mut iterations = 0;
    let mut round = 0;
    let solved = loop {
        round += 1;
        if round > MAX_PHASE_ROUNDS {
            return Err(format!(
                "three-phase column: the stage phase counts did not settle in {MAX_PHASE_ROUNDS} rounds{}",
                column.three_phase_stages()
            ));
        }
        let (v, used) = match column.newton(start) {
            Ok(solution) => solution,
            Err(NewtonFailure::Merge {
                block,
                v,
                iterations: used,
            }) => {
                iterations += used;
                start = column.merge(&v, block);
                remember(&mut seen, &column.split)?;
                continue;
            }
            Err(NewtonFailure::Failed(message)) => {
                // Single-liquid stages failed: retry with two liquids where
                // the initial profile already splits.
                let splits = if round == 1 {
                    column.splitting_stages(&v0)
                } else {
                    Vec::new()
                };
                if splits.is_empty() {
                    return Err(format!("{message}{}", column.three_phase_stages()));
                }
                start = v0.clone();
                for (j, light, heavy) in splits {
                    start = column.split_block(&start, j, &light, &heavy);
                }
                remember(&mut seen, &column.split)?;
                continue;
            }
        };
        iterations += used;
        let props = column.all_properties(&v)?;
        let changes = column.phase_changes(&v, &props)?;
        if changes.is_empty() {
            break v;
        }
        start = v;
        for change in changes {
            start = match change {
                PhaseChange::Split(j, light, heavy) => {
                    column.split_block(&start, j, &light, &heavy)
                }
                PhaseChange::Merge(j) => column.merge(&start, j),
            };
        }
        remember(&mut seen, &column.split)?;
    };
    let props = column.all_properties(&solved)?;
    let decoded = column.decode_all(&solved);

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
    // Both liquids of every three-phase stage, ordered by mass density.
    let mut splits = vec![Some(StageLiquidSplit {
        light_mol_s: light_flow,
        light_composition: xa.clone(),
        heavy_mol_s: heavy_flow,
        heavy_composition: xb.clone(),
    })];
    let mut phase_counts = vec![2_u8];
    for (j, b) in decoded.iter().enumerate().skip(1) {
        let Some(third) = &props[j].third else {
            splits.push(None);
            phase_counts.push(2);
            continue;
        };
        let first = (total(&b.l), normalized(&b.l));
        let second = (total(&b.l2), normalized(&b.l2));
        let (light, heavy) = if props[j].first.mass_density <= third.mass_density {
            (first, second)
        } else {
            (second, first)
        };
        let distance = max_difference(&light.1, &heavy.1);
        if distance < 0.02 {
            warnings.push(format!(
                "TPC-NEAR-PLAIT-POINT: stage {} liquids differ by only {distance:.4} in mole fraction",
                j + 1
            ));
        }
        splits.push(Some(StageLiquidSplit {
            light_mol_s: light.0,
            light_composition: light.1,
            heavy_mol_s: heavy.0,
            heavy_composition: heavy.1,
        }));
        phase_counts.push(3);
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
    let (bottoms_flows, _) = Column::liquid_out(last, &props[blocks - 1]);
    let bottoms_flow = total(&bottoms_flows);
    // A three-phase reboiler delivers both liquids, the light one first.
    let bottoms_liquids = props[blocks - 1].third.clone().map(|third| {
        let first = (
            total(&last.l) / bottoms_flow,
            normalized(&last.l),
            props[blocks - 1].first.clone(),
        );
        let second = (total(&last.l2) / bottoms_flow, normalized(&last.l2), third);
        if first.2.mass_density <= second.2.mass_density {
            [first, second]
        } else {
            [second, first]
        }
    });
    let outlets = vec![
        product(ctx, distillate.0, distillate.1, decanter_t, pressure, None)?,
        product(ctx, decant.0, decant.1, decanter_t, pressure, None)?,
        product(
            ctx,
            bottoms_flow,
            &normalized(&bottoms_flows),
            last.t,
            pressure,
            bottoms_liquids,
        )?,
    ];
    let (reflux_flows, _) = column.reflux_of(decanter, &props[0]);
    let mut liquid_flows = vec![total(&reflux_flows)];
    let mut vapor_flows = vec![0.0];
    let mut liquid_compositions = vec![normalized(&reflux_flows)];
    let mut vapor_compositions = vec![condensate];
    let mut temperatures = vec![decanter_t];
    for (b, p) in decoded[1..].iter().zip(&props[1..]) {
        let (liquid, _) = Column::liquid_out(b, p);
        liquid_flows.push(total(&liquid));
        vapor_flows.push(total(&b.v));
        liquid_compositions.push(normalized(&liquid));
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
            stage_phase_counts: phase_counts,
            stage_liquid_splits: splits,
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
    fn a_recurring_phase_set_is_an_oscillation_failure() {
        let mut seen = vec![vec![false; 4]];
        remember(&mut seen, &[false, true, false, false]).unwrap();
        remember(&mut seen, &[false, true, true, false]).unwrap();
        let message = remember(&mut seen, &[false, true, false, false]).unwrap_err();
        assert!(
            message.contains("phase counts oscillate") && message.contains("stages [2]"),
            "{message}"
        );
        assert_eq!(seen.len(), 3);
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
