//! Countercurrent liquid–liquid extractor (#189): ideal equilibrium stages
//! (or an ideal packed section via HETP) at one column pressure, isothermal
//! or adiabatic.
//!
//! Stages are numbered 1..N from the **feed end**. The feed enters stage 1 as
//! the raffinate phase and the solvent enters stage N as the extract phase;
//! the extract leaves stage 1 and the raffinate leaves stage N (Seader,
//! Henley and Roper, *Separation Process Principles*, 3rd ed., 2011, ch. 8
//! and §10.5; Treybal, *Mass-Transfer Operations*, 3rd ed., 1980, ch. 10).
//! Whether the feed end is the top or the bottom of the column follows from
//! the phase densities ([`FeedPhase`]).
//!
//! Stage equations (in − out, stage j, R_0 = feed, E_{N+1} = solvent):
//!
//! - component: r_{i,j−1} + e_{i,j+1} − r_ij − e_ij = 0;
//! - equilibrium: ln x^E_ij + ln φ^E_ij − ln x^R_ij − ln φ^R_ij = 0, with the
//!   package's liquid fugacity coefficients (NRTL-IG: activity coefficient
//!   and pure-liquid standard state; cubic: liquid root);
//! - energy (adiabatic only): H^R_{j−1} + H^E_{j+1} − H^R_j − H^E_j = 0.
//!
//! The unknowns are ln e_ij, ln r_ij of the components that enter the column
//! and, when adiabatic, T_j. The solve starts from the two-liquid flash of the
//! combined inlets, runs isothermal sum-rates passes (Tsuboka–Katayama type
//! successive substitution with a tridiagonal component balance) and finishes
//! with a damped Newton iteration on the block-tridiagonal stage Jacobian.
//! Every solved stage is then checked against an independent VLLE flash of its
//! overall composition (no vapor, exactly two liquids, same tie line).

use property_packages::{
    Convergence, EquilibriumPhase, EquilibriumState, PhaseProperties, PhaseRole, RootSelection,
};
use thermo_types::Phase;

use crate::MaterialStream;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};
use crate::units::columns::packing::{PackedSectionSpec, PackedStages};
use crate::units::{CalcContext, UnitOutput};

/// Largest number of ideal stages of one extractor.
pub const MAX_EXTRACTOR_STAGES: u16 = 50;

/// Newton tolerance on the scaled stage residuals (max norm).
const TOLERANCE: f64 = 1e-12;
const MAX_NEWTON: usize = 60;
const MAX_SUBSTITUTION: usize = 40;
/// Largest Newton change of a log flow and of a stage temperature (K).
const MAX_LN_STEP: f64 = 2.0;
const MAX_T_STEP: f64 = 15.0;
/// Stage phases closer than this (max |x^E − x^R|) are a plait point or the
/// trivial solution: explicit failure.
pub const PLAIT_FAILURE_DISTANCE: f64 = 1e-3;
/// Stage phases closer than this give a near-plait-point warning.
pub const PLAIT_WARNING_DISTANCE: f64 = 0.02;
/// Largest deviation of the solved stage tie line from the VLLE check flash.
const FLASH_AGREEMENT: f64 = 1e-6;

/// Ideal stages or a packed section of the extractor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExtractorInternals {
    /// Number of ideal equilibrium stages, 1 to [`MAX_EXTRACTOR_STAGES`].
    Stages(u16),
    /// One packed section; its ideal stages follow from height / HETP.
    Packed(PackedSectionSpec),
}

impl Default for ExtractorInternals {
    fn default() -> Self {
        Self::Stages(5)
    }
}

impl ExtractorInternals {
    /// Ideal stage count and the resolved packed section, if any.
    pub fn resolve(self) -> Result<(usize, Option<PackedStages>), String> {
        let (count, bed) = match self {
            Self::Stages(count) => (count, None),
            Self::Packed(spec) => {
                let packed = spec.resolve()?;
                (packed.stages, Some(packed))
            }
        };
        if !(1..=MAX_EXTRACTOR_STAGES).contains(&count) {
            return Err(format!(
                "liquid-liquid extractor needs 1 to {MAX_EXTRACTOR_STAGES} ideal stages, got {count}"
            ));
        }
        Ok((usize::from(count), bed))
    }

    fn issues(self) -> Vec<SpecIssue> {
        match self {
            Self::Stages(count) if !(1..=MAX_EXTRACTOR_STAGES).contains(&count) => {
                vec![SpecIssue {
                    field: "stages",
                    problem: SpecProblem::CountOutOfRange,
                }]
            }
            Self::Stages(_) => Vec::new(),
            Self::Packed(spec) => {
                let issues = spec.issues();
                if issues.is_empty() && self.resolve().is_err() {
                    vec![SpecIssue {
                        field: "internals",
                        problem: SpecProblem::CountOutOfRange,
                    }]
                } else {
                    issues
                }
            }
        }
    }
}

/// Thermal mode of the extractor (mutually exclusive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExtractorThermal {
    /// Every stage at one specified temperature; the duty that holds it is
    /// calculated.
    Isothermal {
        /// Stage temperature, K.
        temperature_k: Option<f64>,
    },
    /// No heat exchange; stage temperatures follow from the energy balances
    /// with the inlet enthalpies.
    Adiabatic,
}

impl Default for ExtractorThermal {
    fn default() -> Self {
        Self::Isothermal {
            temperature_k: None,
        }
    }
}

/// Which liquid the feed (raffinate) phase is, which fixes the column
/// orientation: a heavy feed enters at the top, a light feed at the bottom.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum FeedPhase {
    /// Determined from the stage phase densities.
    #[default]
    Detect,
    /// The raffinate is the heavier liquid on every stage (feed at the top).
    Heavy,
    /// The raffinate is the lighter liquid on every stage (feed at the bottom).
    Light,
}

/// Liquid–liquid extractor specification (#189).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LiquidExtractorSpec {
    /// Ideal stages or packed section.
    pub internals: ExtractorInternals,
    /// Column pressure in Pa (absolute), at most both inlet pressures.
    pub pressure_pa: Option<f64>,
    /// Isothermal or adiabatic.
    pub thermal: ExtractorThermal,
    /// Feed phase and orientation.
    pub feed_phase: FeedPhase,
}

impl LiquidExtractorSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        check_optional(
            &mut issues,
            "column pressure",
            self.pressure_pa,
            Sign::Positive,
        );
        if let ExtractorThermal::Isothermal { temperature_k } = self.thermal {
            check_optional(&mut issues, "temperature", temperature_k, Sign::Positive);
        }
        issues.extend(self.internals.issues());
        issues
    }
}

/// Solved stage profile of a liquid–liquid extractor. Every vector is in
/// stage order from the feed end (stage 1: feed in, extract out) to the
/// solvent end (stage N: solvent in, raffinate out).
#[derive(Debug, Clone, PartialEq)]
pub struct LiquidExtractorResult {
    pub stage_temperatures_k: Vec<f64>,
    /// Extract-phase flow leaving each stage, mol/s.
    pub extract_flows_mol_s: Vec<f64>,
    /// Raffinate-phase flow leaving each stage, mol/s.
    pub raffinate_flows_mol_s: Vec<f64>,
    /// Extract-phase mole fractions per stage, package order.
    pub extract_compositions: Vec<Vec<f64>>,
    /// Raffinate-phase mole fractions per stage, package order.
    pub raffinate_compositions: Vec<Vec<f64>>,
    /// Whether the raffinate (feed) phase is the heavier liquid, so that the
    /// feed enters at the top.
    pub feed_is_heavy: bool,
    /// Duty holding an isothermal column, W, positive into the unit (zero
    /// when adiabatic).
    pub duty_w: f64,
    /// Specified or derived packed height, m.
    pub packed_height_m: Option<f64>,
    /// Accepted Newton iterations.
    pub iterations: usize,
}

/// One inlet as component flows.
struct Inlet {
    flows: Vec<f64>,
    enthalpy_flow_w: f64,
    thermal_scale_w: f64,
}

fn inlet(stream: &MaterialStream, pressure: f64, port: &str) -> Result<Inlet, String> {
    let state = stream
        .state
        .as_ref()
        .ok_or_else(|| format!("liquid-liquid extractor {port} is missing or has no state"))?;
    if !stream.molar_flow_mol_s.is_finite() || stream.molar_flow_mol_s <= 0.0 {
        return Err(format!(
            "liquid-liquid extractor {port} needs a positive finite flow"
        ));
    }
    if state.pressure_pa < pressure {
        return Err(format!(
            "liquid-liquid extractor pressure exceeds the {port} pressure"
        ));
    }
    Ok(Inlet {
        flows: state
            .mole_fractions
            .iter()
            .map(|z| z * stream.molar_flow_mol_s)
            .collect(),
        enthalpy_flow_w: stream.enthalpy_flow_w(),
        thermal_scale_w: stream.thermal_scale_w(),
    })
}

/// Liquid properties of the two phases of one stage.
#[derive(Debug, Clone)]
struct StageProperties {
    extract: PhaseProperties,
    raffinate: PhaseProperties,
}

/// Decoded stage: component flows (package order) and temperature.
struct StageFlows {
    e: Vec<f64>,
    r: Vec<f64>,
    t: f64,
}

impl StageFlows {
    fn extract_total(&self) -> f64 {
        self.e.iter().sum()
    }
    fn raffinate_total(&self) -> f64 {
        self.r.iter().sum()
    }
}

fn normalized(flows: &[f64]) -> Vec<f64> {
    let total: f64 = flows.iter().sum();
    flows.iter().map(|f| f / total).collect()
}

struct Cascade<'a, 'p> {
    ctx: &'a CalcContext<'p>,
    pressure: f64,
    stages: usize,
    components: usize,
    active: Vec<usize>,
    feed: Inlet,
    solvent: Inlet,
    /// Per active component: total inlet flow.
    scale: Vec<f64>,
    energy_scale: f64,
    /// `Some(T)` isothermal, `None` adiabatic.
    fixed_temperature: Option<f64>,
    temperature_range: [f64; 2],
}

impl Cascade<'_, '_> {
    fn width(&self) -> usize {
        2 * self.active.len() + usize::from(self.fixed_temperature.is_none())
    }

    fn decode(&self, v: &[f64], j: usize) -> StageFlows {
        let (w, m) = (self.width(), self.active.len());
        let block = &v[j * w..(j + 1) * w];
        let mut e = vec![0.0; self.components];
        let mut r = vec![0.0; self.components];
        for (k, &i) in self.active.iter().enumerate() {
            e[i] = block[k].exp();
            r[i] = block[m + k].exp();
        }
        let t = self.fixed_temperature.unwrap_or_else(|| block[2 * m]);
        StageFlows { e, r, t }
    }

    fn encode(&self, stages: &[StageFlows]) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.stages * self.width());
        for s in stages {
            v.extend(self.active.iter().map(|&i| s.e[i].ln()));
            v.extend(self.active.iter().map(|&i| s.r[i].ln()));
            if self.fixed_temperature.is_none() {
                v.push(s.t);
            }
        }
        v
    }

    fn properties(&self, v: &[f64], j: usize) -> Result<StageProperties, String> {
        let s = self.decode(v, j);
        if !s.t.is_finite()
            || !(self.temperature_range[0]..=self.temperature_range[1]).contains(&s.t)
        {
            return Err(format!(
                "stage {} temperature {:.3} K leaves the package interval",
                j + 1,
                s.t
            ));
        }
        if s.e.iter().chain(&s.r).any(|f| !f.is_finite()) {
            return Err(format!("stage {} flows are not finite", j + 1));
        }
        let liquid = |flows: &[f64]| {
            self.ctx.homogeneous(
                &normalized(flows),
                s.t,
                self.pressure,
                RootSelection::SmallestVolume,
            )
        };
        Ok(StageProperties {
            extract: liquid(&s.e)?,
            raffinate: liquid(&s.r)?,
        })
    }

    fn all_properties(&self, v: &[f64]) -> Result<Vec<StageProperties>, String> {
        (0..self.stages).map(|j| self.properties(v, j)).collect()
    }

    fn residual(&self, v: &[f64], props: &[StageProperties]) -> Vec<f64> {
        let n = self.stages;
        let flows: Vec<StageFlows> = (0..n).map(|j| self.decode(v, j)).collect();
        let mut out = Vec::with_capacity(n * self.width());
        for j in 0..n {
            let s = &flows[j];
            let r_in = if j == 0 {
                &self.feed.flows
            } else {
                &flows[j - 1].r
            };
            let e_in = if j + 1 == n {
                &self.solvent.flows
            } else {
                &flows[j + 1].e
            };
            for (k, &i) in self.active.iter().enumerate() {
                out.push((r_in[i] + e_in[i] - s.r[i] - s.e[i]) / self.scale[k]);
            }
            let (ln_e, ln_r) = (s.extract_total().ln(), s.raffinate_total().ln());
            let p = &props[j];
            for &i in &self.active {
                out.push(
                    s.e[i].ln() - ln_e + p.extract.ln_fugacity_coefficients[i] - s.r[i].ln() + ln_r
                        - p.raffinate.ln_fugacity_coefficients[i],
                );
            }
            if self.fixed_temperature.is_none() {
                let h_in = if j == 0 {
                    self.feed.enthalpy_flow_w
                } else {
                    flows[j - 1].raffinate_total() * props[j - 1].raffinate.enthalpy
                } + if j + 1 == n {
                    self.solvent.enthalpy_flow_w
                } else {
                    flows[j + 1].extract_total() * props[j + 1].extract.enthalpy
                };
                let h_out = s.raffinate_total() * p.raffinate.enthalpy
                    + s.extract_total() * p.extract.enthalpy;
                out.push((h_in - h_out) / self.energy_scale);
            }
        }
        out
    }

    /// Isothermal sum-rates passes at the current stage temperatures: fixed
    /// K_ij = φ^R/φ^E and phase flows give a tridiagonal component balance
    /// per component (Thomas algorithm). Stops early on any failure; Newton
    /// continues from the last good profile.
    fn substitution(&self, v: &mut Vec<f64>) {
        let n = self.stages;
        for _ in 0..MAX_SUBSTITUTION {
            let Ok(props) = self.all_properties(v) else {
                return;
            };
            let flows: Vec<StageFlows> = (0..n).map(|j| self.decode(v, j)).collect();
            let mut next: Vec<StageFlows> = flows
                .iter()
                .map(|s| StageFlows {
                    e: vec![0.0; self.components],
                    r: vec![0.0; self.components],
                    t: s.t,
                })
                .collect();
            for &i in &self.active {
                // Extraction factors s_j = K_ij E_j / R_j.
                let factor: Vec<f64> = (0..n)
                    .map(|j| {
                        (props[j].raffinate.ln_fugacity_coefficients[i]
                            - props[j].extract.ln_fugacity_coefficients[i])
                            .exp()
                            * flows[j].extract_total()
                            / flows[j].raffinate_total()
                    })
                    .collect();
                // (1 + s_j) r_j − r_{j−1} − s_{j+1} r_{j+1} = F δ_j0 + S δ_j,N−1.
                let diag: Vec<f64> = factor.iter().map(|s| 1.0 + s).collect();
                let mut rhs = vec![0.0; n];
                rhs[0] += self.feed.flows[i];
                rhs[n - 1] += self.solvent.flows[i];
                let Some(r) = thomas(&diag, &factor, &rhs) else {
                    return;
                };
                for j in 0..n {
                    if !(r[j].is_finite() && r[j] > 0.0 && factor[j].is_finite()) {
                        return;
                    }
                    next[j].r[i] = r[j];
                    next[j].e[i] = factor[j] * r[j];
                    if next[j].e[i] <= 0.0 {
                        return;
                    }
                }
            }
            let candidate = self.encode(&next);
            let change = candidate
                .iter()
                .zip(v.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            *v = candidate;
            if change < 1e-9 {
                return;
            }
        }
    }

    /// Block-tridiagonal forward-difference Jacobian: stage j's residuals
    /// depend on stages j−1, j, j+1 only, so stages three apart are perturbed
    /// together.
    #[allow(clippy::type_complexity)]
    fn jacobian(
        &self,
        v: &[f64],
        props: &[StageProperties],
        r: &[f64],
    ) -> Result<(Vec<Block>, Vec<Block>, Vec<Block>), String> {
        let (n, w, m) = (self.stages, self.width(), self.active.len());
        let mut lower = vec![zero_block(w); n];
        let mut diag = vec![zero_block(w); n];
        let mut upper = vec![zero_block(w); n];
        for color in 0..3.min(n) {
            let group: Vec<usize> = (color..n).step_by(3).collect();
            for k in 0..w {
                let base = if k < 2 * m { 1e-7 } else { 1e-5 };
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
                    attempt.ok_or("liquid-liquid extractor Jacobian evaluation failed")?;
                for &j in &group {
                    for row_stage in j.saturating_sub(1)..=(j + 1).min(n - 1) {
                        for row in 0..w {
                            let index = row_stage * w + row;
                            let d = (shifted[index] - r[index]) / h;
                            if row_stage + 1 == j {
                                upper[row_stage][row][k] = d;
                            } else if row_stage == j {
                                diag[j][row][k] = d;
                            } else {
                                lower[row_stage][row][k] = d;
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
            .map_err(|e| format!("liquid-liquid extractor initialization: {e}"))?;
        let mut r = self.residual(&v, &props);
        for iteration in 0..=MAX_NEWTON {
            let error = max_abs(&r);
            if !error.is_finite() {
                return Err("liquid-liquid extractor stage residuals are not finite".into());
            }
            if error <= TOLERANCE {
                return Ok((v, iteration));
            }
            if iteration == MAX_NEWTON {
                break;
            }
            let (lower, diag, upper) = self.jacobian(&v, &props, &r)?;
            let step = block_tridiagonal_solve(lower, diag, upper, r.iter().map(|x| -x).collect())
                .ok_or("liquid-liquid extractor stage Jacobian is singular (plait point or vanishing phase)")?;
            let mut limit = 1.0_f64;
            for (index, d) in step.iter().enumerate() {
                let cap = if index % w < 2 * m {
                    MAX_LN_STEP
                } else {
                    MAX_T_STEP
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
                    "liquid-liquid extractor stage equations stalled (scaled residual {error:.3e}); the column may approach a plait point or lose its two-liquid region"
                ));
            };
            v = trial;
            props = p;
            r = rt;
        }
        Err(format!(
            "liquid-liquid extractor did not converge in {MAX_NEWTON} Newton iterations (scaled residual {:.3e})",
            max_abs(&r)
        ))
    }
}

type Block = Vec<Vec<f64>>;

fn zero_block(w: usize) -> Block {
    vec![vec![0.0; w]; w]
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|x| x.abs()).fold(0.0, f64::max)
}

fn sum_squares(values: &[f64]) -> f64 {
    values.iter().map(|x| x * x).sum()
}

/// Thomas algorithm for diag_j r_j − r_{j−1} − upper_{j+1} r_{j+1} = rhs_j
/// (sub-diagonal −1, super-diagonal −factor_{j+1}).
#[allow(clippy::needless_range_loop)]
fn thomas(diag: &[f64], factor: &[f64], rhs: &[f64]) -> Option<Vec<f64>> {
    let n = diag.len();
    let mut c = vec![0.0; n];
    let mut d = vec![0.0; n];
    let mut b = diag[0];
    if b == 0.0 {
        return None;
    }
    c[0] = if n > 1 { -factor[1] / b } else { 0.0 };
    d[0] = rhs[0] / b;
    for j in 1..n {
        b = diag[j] + c[j - 1];
        if b == 0.0 || !b.is_finite() {
            return None;
        }
        c[j] = if j + 1 < n { -factor[j + 1] / b } else { 0.0 };
        d[j] = (rhs[j] + d[j - 1]) / b;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1];
    for j in (0..n - 1).rev() {
        x[j] = d[j] - c[j] * x[j + 1];
    }
    Some(x)
}

/// Solves `a X = b` (b with any number of columns) by Gaussian elimination
/// with partial pivoting.
#[allow(clippy::needless_range_loop)] // dense elimination reads clearest with indices
fn solve_dense(mut a: Block, mut b: Block) -> Option<Block> {
    let n = a.len();
    let scale = a
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);
    for col in 0..n {
        let pivot = (col..n).max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))?;
        if a[pivot][col].abs() <= 1e-300 * scale || !a[pivot][col].is_finite() {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let f = a[row][col] / a[col][col];
            if f != 0.0 {
                for k in col..n {
                    a[row][k] -= f * a[col][k];
                }
                for k in 0..b[row].len() {
                    b[row][k] -= f * b[col][k];
                }
            }
        }
    }
    for row in (0..n).rev() {
        for k in 0..b[row].len() {
            let mut value = b[row][k];
            for j in row + 1..n {
                value -= a[row][j] * b[j][k];
            }
            b[row][k] = value / a[row][row];
        }
    }
    Some(b)
}

/// Block Thomas algorithm for L_j x_{j−1} + D_j x_j + U_j x_{j+1} = rhs_j.
#[allow(clippy::needless_range_loop)]
fn block_tridiagonal_solve(
    lower: Vec<Block>,
    mut diag: Vec<Block>,
    upper: Vec<Block>,
    rhs: Vec<f64>,
) -> Option<Vec<f64>> {
    let n = diag.len();
    let w = diag[0].len();
    let mut g: Vec<Vec<f64>> = rhs.chunks(w).map(<[f64]>::to_vec).collect();
    // Z_j = D'_j^{-1} U_j, Y_j = D'_j^{-1} g_j.
    let mut z: Vec<Block> = Vec::with_capacity(n);
    let mut y: Vec<Vec<f64>> = Vec::with_capacity(n);
    for j in 0..n {
        if j > 0 {
            let l = &lower[j];
            for row in 0..w {
                for col in 0..w {
                    let mut sum = 0.0;
                    for k in 0..w {
                        sum += l[row][k] * z[j - 1][k][col];
                    }
                    diag[j][row][col] -= sum;
                }
                let mut sum = 0.0;
                for k in 0..w {
                    sum += l[row][k] * y[j - 1][k];
                }
                g[j][row] -= sum;
            }
        }
        let columns = if j + 1 < n { w + 1 } else { 1 };
        let b: Block = (0..w)
            .map(|row| {
                let mut line = Vec::with_capacity(columns);
                if j + 1 < n {
                    line.extend_from_slice(&upper[j][row]);
                }
                line.push(g[j][row]);
                line
            })
            .collect();
        let x = solve_dense(diag[j].clone(), b)?;
        y.push(x.iter().map(|line| line[columns - 1]).collect());
        z.push(x.iter().map(|line| line[..columns - 1].to_vec()).collect());
    }
    let mut out = vec![vec![0.0; w]; n];
    out[n - 1] = y[n - 1].clone();
    for j in (0..n - 1).rev() {
        for row in 0..w {
            let mut value = y[j][row];
            for k in 0..w {
                value -= z[j][row][k] * out[j + 1][k];
            }
            out[j][row] = value;
        }
    }
    let flat: Vec<f64> = out.into_iter().flatten().collect();
    flat.iter().all(|x| x.is_finite()).then_some(flat)
}

fn product(
    ctx: &CalcContext<'_>,
    flows: &[f64],
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<MaterialStream, String> {
    let composition = normalized(flows);
    let state = ctx.homogeneous(
        &composition,
        temperature_k,
        pressure_pa,
        RootSelection::SmallestVolume,
    )?;
    Ok(MaterialStream {
        molar_flow_mol_s: flows.iter().sum(),
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

fn max_difference(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// The two liquids of a VLLE flash, explicit failures for vapor and for a
/// single liquid.
fn two_liquids(
    ctx: &CalcContext<'_>,
    z: &[f64],
    t: f64,
    p: f64,
    place: &str,
) -> Result<[(Vec<f64>, f64); 2], String> {
    let flash = ctx
        .vlle(z, t, p)
        .map_err(|e| format!("liquid-liquid extractor needs the VLLE/LLE flash ({place}): {e}"))?;
    if flash.roles.contains(&PhaseRole::Vapor) {
        return Err(format!(
            "liquid-liquid extractor: vapor forms {place} at {t:.3} K and {p:.0} Pa; raise the column pressure"
        ));
    }
    let liquids: Vec<(Vec<f64>, f64)> = flash
        .roles
        .iter()
        .zip(&flash.state.phases)
        .filter(|(role, _)| matches!(role, PhaseRole::LightLiquid | PhaseRole::HeavyLiquid))
        .map(|(_, phase)| (phase.composition.clone(), phase.phase_fraction))
        .collect();
    match <[(Vec<f64>, f64); 2]>::try_from(liquids) {
        Ok(pair) => Ok(pair),
        Err(_) => Err(format!(
            "liquid-liquid extractor: no liquid-liquid split {place} at {t:.3} K (one liquid phase; miscible region or beyond the plait point)"
        )),
    }
}

pub(super) fn calculate(
    spec: &LiquidExtractorSpec,
    inputs: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, LiquidExtractorResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!(
            "liquid-liquid extractor {} {}",
            issue.field, issue.problem
        ));
    }
    let (stages, bed) = spec.internals.resolve()?;
    let pressure = spec.pressure_pa.expect("checked extractor pressure");
    let feed = inlet(&inputs[0], pressure, "feed")?;
    let solvent = inlet(&inputs[1], pressure, "solvent")?;
    let components = ctx.component_count();
    let total: Vec<f64> = feed
        .flows
        .iter()
        .zip(&solvent.flows)
        .map(|(a, b)| a + b)
        .collect();
    let total_flow: f64 = total.iter().sum();
    let active: Vec<usize> = (0..components).filter(|&i| total[i] > 0.0).collect();
    if active.len() < 2 {
        return Err("liquid-liquid extractor needs at least two components in its inlets".into());
    }
    let z: Vec<f64> = total.iter().map(|f| f / total_flow).collect();
    let temperature_range = ctx.package().temperature_interval();
    let fixed_temperature = match spec.thermal {
        ExtractorThermal::Isothermal { temperature_k } => {
            let t = temperature_k.expect("checked extractor temperature");
            if !(temperature_range[0]..=temperature_range[1]).contains(&t) {
                return Err(format!(
                    "liquid-liquid extractor temperature {t} K is outside the package interval {:.2}..{:.2} K",
                    temperature_range[0], temperature_range[1]
                ));
            }
            Some(t)
        }
        ExtractorThermal::Adiabatic => None,
    };
    // Starting temperature: the specified one, or the adiabatic mixing
    // temperature of both inlets at column pressure.
    let t0 = match fixed_temperature {
        Some(t) => t,
        None => {
            let h = (feed.enthalpy_flow_w + solvent.enthalpy_flow_w) / total_flow;
            ctx.vlle_ph(&z, pressure, h)
                .map_err(|e| {
                    format!(
                        "liquid-liquid extractor needs the VLLE/LLE flash (adiabatic mixing of the inlets): {e}"
                    )
                })?
                .state
                .temperature_k
        }
    };
    let [first, second] = two_liquids(ctx, &z, t0, pressure, "in the mixed inlets")?;
    // The extract is the liquid carrying more of the solvent's material.
    let solvent_share: Vec<f64> = total
        .iter()
        .zip(&solvent.flows)
        .map(|(t, s)| if *t > 0.0 { s / t } else { 0.0 })
        .collect();
    let weight = |x: &[f64]| -> f64 { x.iter().zip(&solvent_share).map(|(a, b)| a * b).sum() };
    let (extract0, raffinate0) = if weight(&first.0) >= weight(&second.0) {
        (first, second)
    } else {
        (second, first)
    };
    let floor = 1e-14 * total_flow;
    let initial: Vec<StageFlows> = (0..stages)
        .map(|_| {
            let flows = |x: &[f64], beta: f64| -> Vec<f64> {
                x.iter()
                    .enumerate()
                    .map(|(i, xi)| {
                        if total[i] > 0.0 {
                            (xi * beta * total_flow).max(floor)
                        } else {
                            0.0
                        }
                    })
                    .collect()
            };
            StageFlows {
                e: flows(&extract0.0, extract0.1),
                r: flows(&raffinate0.0, raffinate0.1),
                t: t0,
            }
        })
        .collect();
    let scale: Vec<f64> = active.iter().map(|&i| total[i]).collect();
    let energy_scale = (feed.thermal_scale_w + solvent.thermal_scale_w).max(total_flow);
    let cascade = Cascade {
        ctx,
        pressure,
        stages,
        components,
        active,
        feed,
        solvent,
        scale,
        energy_scale,
        fixed_temperature,
        temperature_range,
    };
    let mut v = cascade.encode(&initial);
    cascade.substitution(&mut v);
    let (solved, iterations) = cascade.newton(v)?;
    let props = cascade.all_properties(&solved)?;
    let flows: Vec<StageFlows> = (0..stages).map(|j| cascade.decode(&solved, j)).collect();

    // Verification of every stage: distinct phases, no vapor, and the
    // stable liquid-liquid split of the stage's overall composition.
    let mut warnings = Vec::new();
    let mut heavier = Vec::with_capacity(stages);
    for (j, s) in flows.iter().enumerate() {
        let stage = j + 1;
        let (xe, xr) = (normalized(&s.e), normalized(&s.r));
        let distance = max_difference(&xe, &xr);
        if distance < PLAIT_FAILURE_DISTANCE {
            return Err(format!(
                "liquid-liquid extractor stage {stage}: the two liquids coincide (max mole-fraction difference {distance:.2e}); the column reaches a plait point or the trivial solution"
            ));
        }
        if distance < PLAIT_WARNING_DISTANCE {
            warnings.push(format!(
                "LLX-NEAR-PLAIT-POINT: stage {stage} liquids differ by only {distance:.4} in mole fraction"
            ));
        }
        let stage_total = s.extract_total() + s.raffinate_total();
        let zj: Vec<f64> =
            s.e.iter()
                .zip(&s.r)
                .map(|(a, b)| (a + b) / stage_total)
                .collect();
        let [a, b] = two_liquids(ctx, &zj, s.t, pressure, &format!("at stage {stage}"))?;
        let (fe, fr) = if max_difference(&a.0, &xe) <= max_difference(&b.0, &xe) {
            (a, b)
        } else {
            (b, a)
        };
        let deviation = max_difference(&fe.0, &xe).max(max_difference(&fr.0, &xr));
        if deviation > FLASH_AGREEMENT {
            return Err(format!(
                "liquid-liquid extractor stage {stage}: the solved tie line is not the stable liquid-liquid split of the stage (deviation {deviation:.2e} from the VLLE flash)"
            ));
        }
        heavier.push(props[j].raffinate.mass_density > props[j].extract.mass_density);
    }
    let feed_is_heavy = heavier[0];
    if let Some(j) = heavier.iter().position(|h| *h != feed_is_heavy) {
        return Err(format!(
            "liquid-liquid extractor: the raffinate is the {} liquid at stage 1 but the {} at stage {}; the phases invert their density order and countercurrent gravity flow is not defined",
            if feed_is_heavy { "heavier" } else { "lighter" },
            if feed_is_heavy { "lighter" } else { "heavier" },
            j + 1
        ));
    }
    match (spec.feed_phase, feed_is_heavy) {
        (FeedPhase::Heavy, false) => {
            return Err(
                "liquid-liquid extractor: the feed phase is specified heavy, but the raffinate is the lighter liquid".into(),
            );
        }
        (FeedPhase::Light, true) => {
            return Err(
                "liquid-liquid extractor: the feed phase is specified light, but the raffinate is the heavier liquid".into(),
            );
        }
        _ => {}
    }
    let last = &flows[stages - 1];
    let outlets = vec![
        product(ctx, &flows[0].e, flows[0].t, pressure)?,
        product(ctx, &last.r, last.t, pressure)?,
    ];
    let duty_w = match fixed_temperature {
        Some(_) => {
            outlets
                .iter()
                .map(MaterialStream::enthalpy_flow_w)
                .sum::<f64>()
                - cascade.feed.enthalpy_flow_w
                - cascade.solvent.enthalpy_flow_w
        }
        None => 0.0,
    };
    let result = LiquidExtractorResult {
        stage_temperatures_k: flows.iter().map(|s| s.t).collect(),
        extract_flows_mol_s: flows.iter().map(StageFlows::extract_total).collect(),
        raffinate_flows_mol_s: flows.iter().map(StageFlows::raffinate_total).collect(),
        extract_compositions: flows.iter().map(|s| normalized(&s.e)).collect(),
        raffinate_compositions: flows.iter().map(|s| normalized(&s.r)).collect(),
        feed_is_heavy,
        duty_w,
        packed_height_m: bed.map(|b| b.packed_height_m),
        iterations,
    };
    Ok((UnitOutput::new(outlets, duty_w, warnings), result))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn block_thomas_matches_a_dense_solve() {
        // Three 2×2 blocks of a diagonally dominant system.
        let lower = vec![
            zero_block(2),
            vec![vec![0.3, -0.1], vec![0.2, 0.4]],
            vec![vec![-0.2, 0.1], vec![0.05, 0.3]],
        ];
        let diag = vec![
            vec![vec![4.0, 1.0], vec![0.5, 3.0]],
            vec![vec![5.0, -1.0], vec![1.0, 4.0]],
            vec![vec![3.0, 0.2], vec![-0.4, 2.5]],
        ];
        let upper = vec![
            vec![vec![0.1, 0.2], vec![-0.3, 0.1]],
            vec![vec![0.4, 0.0], vec![0.1, -0.2]],
            zero_block(2),
        ];
        let rhs = vec![1.0, -2.0, 0.5, 3.0, -1.0, 2.0];
        let mut dense = vec![vec![0.0; 6]; 6];
        for j in 0..3 {
            for r in 0..2 {
                for c in 0..2 {
                    dense[2 * j + r][2 * j + c] = diag[j][r][c];
                    if j > 0 {
                        dense[2 * j + r][2 * (j - 1) + c] = lower[j][r][c];
                    }
                    if j < 2 {
                        dense[2 * j + r][2 * (j + 1) + c] = upper[j][r][c];
                    }
                }
            }
        }
        let expected = solve_dense(dense, rhs.iter().map(|x| vec![*x]).collect()).unwrap();
        let x = block_tridiagonal_solve(lower, diag, upper, rhs).unwrap();
        for (a, b) in x.iter().zip(expected) {
            assert!((a - b[0]).abs() < 1e-13);
        }
    }

    #[test]
    fn thomas_solves_the_component_balance() {
        let factor = [0.5, 2.0, 1.5];
        let diag: Vec<f64> = factor.iter().map(|s| 1.0 + s).collect();
        let rhs = [1.0, 0.0, 0.25];
        let r = thomas(&diag, &factor, &rhs).unwrap();
        // Residuals of (1 + s_j) r_j − r_{j−1} − s_{j+1} r_{j+1} = rhs_j.
        for j in 0..3 {
            let prev = if j > 0 { r[j - 1] } else { 0.0 };
            let next = if j < 2 { factor[j + 1] * r[j + 1] } else { 0.0 };
            assert!((diag[j] * r[j] - prev - next - rhs[j]).abs() < 1e-14);
        }
        // Overall balance: raffinate out of the last stage plus extract out
        // of the first equals the inlets.
        let out = r[2] + factor[0] * r[0];
        assert!((out - 1.25).abs() < 1e-14);
    }

    #[test]
    fn specification_issues_are_explicit() {
        let spec = LiquidExtractorSpec::default();
        let fields: Vec<_> = spec.issues().iter().map(|i| i.field).collect();
        assert_eq!(fields, ["column pressure", "temperature"]);
        let spec = LiquidExtractorSpec {
            internals: ExtractorInternals::Stages(51),
            pressure_pa: Some(1e5),
            thermal: ExtractorThermal::Adiabatic,
            feed_phase: FeedPhase::Heavy,
        };
        assert_eq!(spec.issues()[0].problem, SpecProblem::CountOutOfRange);
        let packed = ExtractorInternals::Packed(PackedSectionSpec {
            hetp_m: Some(0.5),
            basis: crate::units::columns::packing::PackedSectionBasis::TheoreticalStages(60),
        });
        assert_eq!(packed.issues()[0].problem, SpecProblem::CountOutOfRange);
        let packed = ExtractorInternals::Packed(PackedSectionSpec {
            hetp_m: Some(0.5),
            basis: crate::units::columns::packing::PackedSectionBasis::PackedHeight {
                height_m: Some(2.0),
                rounding: Default::default(),
            },
        });
        assert!(packed.issues().is_empty());
        assert_eq!(packed.resolve().unwrap().0, 4);
    }
}
