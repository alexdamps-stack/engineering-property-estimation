//! Equilibrium-stage foundation for rigorous distillation (#186) and the
//! absorber / stripper (#187): stage state, stream conventions, MESH
//! residuals, overall balances, the tridiagonal component balance and the
//! condenser / reboiler boundary contract. The rigorous distillation
//! solver uses these pieces to calculate a column in the flowsheet.
//!
//! # Conventions
//!
//! Stages are numbered `1..=N` from the top; slices use index `j - 1`.
//! On stage `j`
//!
//! - `L_j` (`liquid_mol_s`) leaves downward to stage `j + 1`; `L_N` is the
//!   bottom liquid product;
//! - `V_j` (`vapor_mol_s`) leaves upward to stage `j - 1`; `V_1` is the top
//!   vapor product;
//! - `U_j` / `W_j` are liquid / vapor side draws leaving at the stage
//!   compositions `x_j` / `y_j`;
//! - `F_j` is the sum of the feeds assigned to the stage, each with its own
//!   overall composition and molar enthalpy (the feed enters as a whole);
//! - `Q_j` is the stage duty in W, **positive into the stage** like every unit
//!   duty in the project.
//!
//! The residuals follow Seader, Henley and Roper, *Separation Process
//! Principles*, 3rd ed., Eqs. 10-1 to 10-4 (Wang–Henke MESH), with the duty
//! sign flipped to the project convention:
//!
//! ```text
//! M_ij = L_{j-1} x_{i,j-1} + V_{j+1} y_{i,j+1} + F_j z_ij − (L_j + U_j) x_ij − (V_j + W_j) y_ij
//! E_ij = η_j K_ij x_ij + (1 − η_j) y_{i,j+1} − y_ij        (Murphree vapor efficiency η_j)
//! S_j  = Σ x_ij − 1,  Σ y_ij − 1
//! H_j  = L_{j-1} h^L_{j-1} + V_{j+1} h^V_{j+1} + F_j h^F_j + Q_j − (L_j + U_j) h^L_j − (V_j + W_j) h^V_j
//! ```
//!
//! Residuals are "in minus out": positive means accumulation. With η = 1 the
//! stage is an equilibrium stage; the bottom stage has no vapor from below,
//! so it is always treated as an equilibrium stage.
//!
//! # Condenser and reboiler
//!
//! A condenser, when present, is stage 1 and a reboiler is stage N, matching
//! the HYSYS / Aspen Plus numbering. A total condenser has `V_1 = 0` and
//! draws the distillate as the liquid side draw `U_1`; the reflux is `L_1`. A
//! partial condenser draws vapor distillate `V_1` and optionally `U_1`. The
//! reboiler's bottom product is `L_N`. The end duties map onto the existing
//! column energy ports (`condenser duty` outlet, `reboiler duty` inlet) with
//! the project's energy-stream sign: positive in the port direction, so the
//! condenser port reports `−Q_1` and the reboiler port `Q_N`
//! ([`ColumnEnds::port_duties`]).

use super::packing::PackedStages;
use crate::engine::MaterialStream;
use crate::units::CalcContext;
use thermo_models::RootSelection;

/// Relative tolerance on the conditions of [`ColumnEnds::check`].
const END_TOLERANCE: f64 = 1e-12;

/// Variables of one stage: conditions, leaving flows, phase compositions and
/// duty. Whether a value is an unknown or a specification is decided by the
/// column solver.
#[derive(Debug, Clone, PartialEq)]
pub struct Stage {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// `L_j`, liquid to the stage below (bottom product on the last stage).
    pub liquid_mol_s: f64,
    /// `V_j`, vapor to the stage above (top product on the first stage).
    pub vapor_mol_s: f64,
    /// `U_j`, liquid side draw (distillate of a total condenser on stage 1).
    pub liquid_draw_mol_s: f64,
    /// `W_j`, vapor side draw.
    pub vapor_draw_mol_s: f64,
    /// Liquid mole fractions `x_j`, package order.
    pub x: Vec<f64>,
    /// Vapor mole fractions `y_j`, package order.
    pub y: Vec<f64>,
    /// `Q_j` in W, positive into the stage.
    pub duty_w: f64,
    /// Murphree vapor efficiency `η_j` in (0, 1]; 1 is an equilibrium stage.
    pub efficiency: f64,
}

impl Stage {
    /// An equilibrium stage without draws or duty.
    pub fn equilibrium(
        temperature_k: f64,
        pressure_pa: f64,
        liquid_mol_s: f64,
        vapor_mol_s: f64,
        x: Vec<f64>,
        y: Vec<f64>,
    ) -> Self {
        Self {
            temperature_k,
            pressure_pa,
            liquid_mol_s,
            vapor_mol_s,
            liquid_draw_mol_s: 0.0,
            vapor_draw_mol_s: 0.0,
            x,
            y,
            duty_w: 0.0,
            efficiency: 1.0,
        }
    }

    /// Total liquid leaving the stage, `L_j + U_j`.
    pub fn liquid_out_mol_s(&self) -> f64 {
        self.liquid_mol_s + self.liquid_draw_mol_s
    }

    /// Total vapor leaving the stage, `V_j + W_j`.
    pub fn vapor_out_mol_s(&self) -> f64 {
        self.vapor_mol_s + self.vapor_draw_mol_s
    }
}

/// A feed entering one stage as a whole (both phases of a two-phase feed).
#[derive(Debug, Clone, PartialEq)]
pub struct StageFeed {
    /// One-based stage number.
    pub stage: usize,
    pub molar_flow_mol_s: f64,
    /// Overall mole fractions, package order.
    pub z: Vec<f64>,
    /// Overall molar enthalpy in J/mol on the package caloric reference.
    pub molar_enthalpy_j_mol: f64,
}

impl StageFeed {
    /// The feed carried by a flowsheet stream. A zero-flow stream without a
    /// state (an unconnected optional inlet) becomes an empty feed.
    pub fn from_stream(
        stage: usize,
        stream: &MaterialStream,
        component_count: usize,
    ) -> Result<Self, String> {
        match &stream.state {
            Some(s) => Ok(Self {
                stage,
                molar_flow_mol_s: stream.molar_flow_mol_s,
                z: s.mole_fractions.clone(),
                molar_enthalpy_j_mol: s.enthalpy(),
            }),
            None if stream.molar_flow_mol_s == 0.0 => Ok(Self {
                stage,
                molar_flow_mol_s: 0.0,
                z: vec![0.0; component_count],
                molar_enthalpy_j_mol: 0.0,
            }),
            None => Err(format!(
                "feed to stage {stage} has flow but no thermodynamic state"
            )),
        }
    }
}

/// Phase properties of one stage evaluated at its T, p, x and y.
#[derive(Debug, Clone, PartialEq)]
pub struct StageProperties {
    /// Equilibrium ratios `K_ij = φ_i^L / φ_i^V`.
    pub k: Vec<f64>,
    /// Molar enthalpy of the stage liquid, J/mol.
    pub liquid_enthalpy_j_mol: f64,
    /// Molar enthalpy of the stage vapor, J/mol.
    pub vapor_enthalpy_j_mol: f64,
}

impl CalcContext<'_> {
    /// [`StageProperties`] of `stage` from the cubic package: K-values from
    /// the fugacity coefficients of `x` (liquid root) and `y` (vapor root)
    /// and the homogeneous phase enthalpies on the same roots. No stability
    /// check; the column solver owns phase existence.
    pub fn stage_properties(&self, stage: &Stage) -> Result<StageProperties, String> {
        let (t, p) = (stage.temperature_k, stage.pressure_pa);
        let k = self.k_values(&stage.x, &stage.y, t, p)?;
        let liquid = self.homogeneous(&stage.x, t, p, RootSelection::SmallestVolume)?;
        let vapor = self.homogeneous(&stage.y, t, p, RootSelection::LargestVolume)?;
        Ok(StageProperties {
            k,
            liquid_enthalpy_j_mol: liquid.enthalpy,
            vapor_enthalpy_j_mol: vapor.enthalpy,
        })
    }
}

/// MESH residuals of one stage (see the module conventions).
#[derive(Debug, Clone, PartialEq)]
pub struct StageResiduals {
    /// `M_ij`, mol/s, in minus out.
    pub component_mol_s: Vec<f64>,
    /// `E_ij`, mole fraction.
    pub equilibrium: Vec<f64>,
    /// `Σ x_ij − 1`.
    pub liquid_summation: f64,
    /// `Σ y_ij − 1`.
    pub vapor_summation: f64,
    /// `H_j`, W, in minus out including the stage duty.
    pub energy_w: f64,
}

/// Overall component and energy balance of a cascade: feeds plus duties
/// minus all products (top vapor, bottom liquid, side draws).
#[derive(Debug, Clone, PartialEq)]
pub struct CascadeBalance {
    /// Component flow into all products, mol/s.
    pub product_component_mol_s: Vec<f64>,
    /// Feed minus product component flows, mol/s.
    pub component_residuals_mol_s: Vec<f64>,
    /// `Σ Q_j`, the net external duty in W (positive into the column).
    pub duty_w: f64,
    /// Feed enthalpy plus duty minus product enthalpy, W.
    pub energy_residual_w: f64,
}

/// A countercurrent stage cascade with its feeds.
#[derive(Debug, Clone, PartialEq)]
pub struct StageCascade {
    /// Stages 1..=N from the top.
    pub stages: Vec<Stage>,
    pub feeds: Vec<StageFeed>,
}

impl StageCascade {
    /// Number of stages.
    pub fn len(&self) -> usize {
        self.stages.len()
    }

    /// Whether the cascade has no stages.
    pub fn is_empty(&self) -> bool {
        self.stages.is_empty()
    }

    /// Checks dimensions, finiteness, nonnegative flows, positive T and p,
    /// efficiencies in (0, 1] and feed stage numbers. Compositions need not
    /// sum to one (that is the summation residual).
    pub fn validate(&self, component_count: usize) -> Result<(), String> {
        if self.stages.is_empty() {
            return Err("a stage cascade needs at least one stage".into());
        }
        for (j, s) in self.stages.iter().enumerate().map(|(j, s)| (j + 1, s)) {
            if !(s.temperature_k.is_finite() && s.temperature_k > 0.0) {
                return Err(format!("stage {j} temperature must be positive"));
            }
            if !(s.pressure_pa.is_finite() && s.pressure_pa > 0.0) {
                return Err(format!("stage {j} pressure must be positive"));
            }
            for (name, v) in [
                ("liquid flow", s.liquid_mol_s),
                ("vapor flow", s.vapor_mol_s),
                ("liquid draw", s.liquid_draw_mol_s),
                ("vapor draw", s.vapor_draw_mol_s),
            ] {
                if !(v.is_finite() && v >= 0.0) {
                    return Err(format!("stage {j} {name} must be finite and nonnegative"));
                }
            }
            if !s.duty_w.is_finite() {
                return Err(format!("stage {j} duty must be finite"));
            }
            if !(s.efficiency > 0.0 && s.efficiency <= 1.0) {
                return Err(format!("stage {j} efficiency must be in (0, 1]"));
            }
            for (name, v) in [("x", &s.x), ("y", &s.y)] {
                if v.len() != component_count || !v.iter().all(|c| c.is_finite() && *c >= 0.0) {
                    return Err(format!(
                        "stage {j} {name} needs {component_count} finite nonnegative fractions"
                    ));
                }
            }
        }
        for f in &self.feeds {
            if !(1..=self.len()).contains(&f.stage) {
                return Err(format!(
                    "feed stage {} is outside stages 1 to {}",
                    f.stage,
                    self.len()
                ));
            }
            if !(f.molar_flow_mol_s.is_finite() && f.molar_flow_mol_s >= 0.0)
                || !f.molar_enthalpy_j_mol.is_finite()
                || f.z.len() != component_count
                || !f.z.iter().all(|c| c.is_finite() && *c >= 0.0)
            {
                return Err(format!("feed to stage {} is invalid", f.stage));
            }
        }
        Ok(())
    }

    /// Feeds of stage `j` (one-based): total flow, component flows and
    /// enthalpy flow.
    fn stage_feed(&self, j: usize, components: usize) -> (f64, Vec<f64>, f64) {
        let mut comp = vec![0.0; components];
        let (mut flow, mut enthalpy) = (0.0, 0.0);
        for f in self.feeds.iter().filter(|f| f.stage == j) {
            flow += f.molar_flow_mol_s;
            enthalpy += f.molar_flow_mol_s * f.molar_enthalpy_j_mol;
            for (c, z) in comp.iter_mut().zip(&f.z) {
                *c += f.molar_flow_mol_s * z;
            }
        }
        (flow, comp, enthalpy)
    }

    /// Total-flow residual of every stage (in minus out), mol/s. Zero for
    /// any stage flows consistent with the feeds and draws, independently
    /// of compositions.
    pub fn total_flow_residuals(&self) -> Vec<f64> {
        let n = self.len();
        (0..n)
            .map(|j| {
                let s = &self.stages[j];
                let above = if j > 0 {
                    self.stages[j - 1].liquid_mol_s
                } else {
                    0.0
                };
                let below = if j + 1 < n {
                    self.stages[j + 1].vapor_mol_s
                } else {
                    0.0
                };
                let feed: f64 = self
                    .feeds
                    .iter()
                    .filter(|f| f.stage == j + 1)
                    .map(|f| f.molar_flow_mol_s)
                    .sum();
                above + below + feed - s.liquid_out_mol_s() - s.vapor_out_mol_s()
            })
            .collect()
    }

    /// MESH residuals of every stage with `props[j]` evaluated at stage `j`.
    pub fn residuals(&self, props: &[StageProperties]) -> Result<Vec<StageResiduals>, String> {
        let c = self.components()?;
        self.validate(c)?;
        self.check_properties(props, c)?;
        let n = self.len();
        Ok((0..n)
            .map(|j| {
                let s = &self.stages[j];
                let (_, feed, feed_h) = self.stage_feed(j + 1, c);
                let above = (j > 0).then(|| (&self.stages[j - 1], &props[j - 1]));
                let below = (j + 1 < n).then(|| (&self.stages[j + 1], &props[j + 1]));
                let component_mol_s = (0..c)
                    .map(|i| {
                        above.map_or(0.0, |(a, _)| a.liquid_mol_s * a.x[i])
                            + below.map_or(0.0, |(b, _)| b.vapor_mol_s * b.y[i])
                            + feed[i]
                            - s.liquid_out_mol_s() * s.x[i]
                            - s.vapor_out_mol_s() * s.y[i]
                    })
                    .collect();
                let eta = if below.is_some() { s.efficiency } else { 1.0 };
                let equilibrium = (0..c)
                    .map(|i| {
                        let entering = below.map_or(0.0, |(b, _)| b.y[i]);
                        eta * props[j].k[i] * s.x[i] + (1.0 - eta) * entering - s.y[i]
                    })
                    .collect();
                let energy_w = above.map_or(0.0, |(a, p)| a.liquid_mol_s * p.liquid_enthalpy_j_mol)
                    + below.map_or(0.0, |(b, p)| b.vapor_mol_s * p.vapor_enthalpy_j_mol)
                    + feed_h
                    + s.duty_w
                    - s.liquid_out_mol_s() * props[j].liquid_enthalpy_j_mol
                    - s.vapor_out_mol_s() * props[j].vapor_enthalpy_j_mol;
                StageResiduals {
                    component_mol_s,
                    equilibrium,
                    liquid_summation: s.x.iter().sum::<f64>() - 1.0,
                    vapor_summation: s.y.iter().sum::<f64>() - 1.0,
                    energy_w,
                }
            })
            .collect())
    }

    /// Overall balance around the cascade. Its residuals equal the sums of
    /// the stage component and energy residuals, because internal streams
    /// cancel.
    pub fn balance(&self, props: &[StageProperties]) -> Result<CascadeBalance, String> {
        let c = self.components()?;
        self.validate(c)?;
        self.check_properties(props, c)?;
        let n = self.len();
        let (top, bottom) = (&self.stages[0], &self.stages[n - 1]);
        let mut product = vec![0.0; c];
        let mut product_h = top.vapor_mol_s * props[0].vapor_enthalpy_j_mol
            + bottom.liquid_mol_s * props[n - 1].liquid_enthalpy_j_mol;
        for (i, p) in product.iter_mut().enumerate() {
            *p = top.vapor_mol_s * top.y[i] + bottom.liquid_mol_s * bottom.x[i];
        }
        for (s, pr) in self.stages.iter().zip(props) {
            product_h += s.liquid_draw_mol_s * pr.liquid_enthalpy_j_mol
                + s.vapor_draw_mol_s * pr.vapor_enthalpy_j_mol;
            for (i, p) in product.iter_mut().enumerate() {
                *p += s.liquid_draw_mol_s * s.x[i] + s.vapor_draw_mol_s * s.y[i];
            }
        }
        let mut feed = vec![0.0; c];
        let mut feed_h = 0.0;
        for j in 1..=n {
            let (_, comp, h) = self.stage_feed(j, c);
            feed_h += h;
            feed.iter_mut().zip(comp).for_each(|(f, v)| *f += v);
        }
        let duty_w: f64 = self.stages.iter().map(|s| s.duty_w).sum();
        Ok(CascadeBalance {
            component_residuals_mol_s: feed.iter().zip(&product).map(|(f, p)| f - p).collect(),
            product_component_mol_s: product,
            duty_w,
            energy_residual_w: feed_h + duty_w - product_h,
        })
    }

    fn components(&self) -> Result<usize, String> {
        self.stages
            .first()
            .map(|s| s.x.len())
            .ok_or_else(|| "a stage cascade needs at least one stage".into())
    }

    fn check_properties(&self, props: &[StageProperties], c: usize) -> Result<(), String> {
        if props.len() != self.len() {
            return Err(format!(
                "{} stage property sets for {} stages",
                props.len(),
                self.len()
            ));
        }
        for (j, p) in props.iter().enumerate() {
            if p.k.len() != c
                || !p.k.iter().all(|k| k.is_finite() && *k > 0.0)
                || !p.liquid_enthalpy_j_mol.is_finite()
                || !p.vapor_enthalpy_j_mol.is_finite()
            {
                return Err(format!("stage {} properties are invalid", j + 1));
            }
        }
        Ok(())
    }

    /// Solves the stage component balances for the liquid mole fractions of
    /// every component with the stage flows and the equilibrium ratios
    /// `k[j][i]` held fixed and `y_ij = K_ij x_ij` (equilibrium stages), by
    /// the Thomas algorithm (Wang–Henke; Seader, Henley and Roper Eqs. 10-6
    /// to 10-13). Returns `x[j][i]`, generally not summing to one: the
    /// bubble-point or sum-rates outer loop normalizes and updates
    /// temperatures and flows. Stage compositions, temperatures and duties
    /// of `self` are not used.
    pub fn solve_liquid_fractions(&self, k: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
        let c = self.components()?;
        self.validate(c)?;
        let n = self.len();
        if k.len() != n
            || k.iter()
                .any(|r| r.len() != c || !r.iter().all(|v| v.is_finite() && *v > 0.0))
        {
            return Err("equilibrium ratios need one positive finite row per stage".into());
        }
        if self.stages.iter().any(|s| s.efficiency != 1.0) {
            return Err("the tridiagonal component balance assumes equilibrium stages".into());
        }
        let feed: Vec<Vec<f64>> = (1..=n).map(|j| self.stage_feed(j, c).1).collect();
        let mut x = vec![vec![0.0; c]; n];
        let (mut p, mut q) = (vec![0.0; n], vec![0.0; n]);
        for i in 0..c {
            // A_j x_{j-1} + B_j x_j + C_j x_{j+1} = D_j, forward elimination.
            for j in 0..n {
                let s = &self.stages[j];
                let a = if j > 0 {
                    self.stages[j - 1].liquid_mol_s
                } else {
                    0.0
                };
                let b = -(s.liquid_out_mol_s() + s.vapor_out_mol_s() * k[j][i]);
                let cc = if j + 1 < n {
                    self.stages[j + 1].vapor_mol_s * k[j + 1][i]
                } else {
                    0.0
                };
                let d = -feed[j][i];
                let (p_prev, q_prev) = if j > 0 {
                    (p[j - 1], q[j - 1])
                } else {
                    (0.0, 0.0)
                };
                let pivot = b - a * p_prev;
                if !pivot.is_finite() || pivot.abs() <= f64::MIN_POSITIVE {
                    return Err(format!(
                        "singular component balance on stage {} (no liquid or vapor leaves it)",
                        j + 1
                    ));
                }
                p[j] = cc / pivot;
                q[j] = (d - a * q_prev) / pivot;
            }
            x[n - 1][i] = q[n - 1];
            for j in (0..n - 1).rev() {
                x[j][i] = q[j] - p[j] * x[j + 1][i];
            }
        }
        if x.iter().flatten().all(|v| v.is_finite()) {
            Ok(x)
        } else {
            Err("nonfinite liquid composition in the component balance".into())
        }
    }
}

/// Condenser at the top stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Condenser {
    /// All overhead vapor condensed: `V_1 = 0`, liquid distillate `U_1`.
    Total,
    /// Vapor distillate `V_1`, optional liquid distillate `U_1`.
    Partial,
}

/// Which ends of a cascade are a condenser (stage 1) and a reboiler
/// (stage N). An absorber or stripper has neither and exchanges no heat
/// through energy ports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ColumnEnds {
    pub condenser: Option<Condenser>,
    pub reboiler: bool,
}

impl ColumnEnds {
    /// Rigorous distillation column (#186) with both ends.
    pub const DISTILLATION: Self = Self {
        condenser: Some(Condenser::Total),
        reboiler: true,
    };

    /// Absorber / stripper (#187) without condenser or reboiler.
    pub const NONE: Self = Self {
        condenser: None,
        reboiler: false,
    };

    /// Checks the cascade against the ends: at least two stages with both
    /// ends, `V_1 = 0` for a total condenser, no vapor draw from a
    /// condenser or reboiler, and zero duty on every stage that has no
    /// energy port (intermediate stages and absent ends), so the net duty
    /// is exactly the sum of the port duties.
    pub fn check(&self, cascade: &StageCascade) -> Result<(), String> {
        let n = cascade.len();
        if n == 0 {
            return Err("a stage cascade needs at least one stage".into());
        }
        if self.condenser.is_some() && self.reboiler && n < 2 {
            return Err("a column with condenser and reboiler needs at least two stages".into());
        }
        let top = &cascade.stages[0];
        let bottom = &cascade.stages[n - 1];
        let scale = |s: &Stage| s.liquid_out_mol_s() + s.vapor_out_mol_s();
        if self.condenser == Some(Condenser::Total) && top.vapor_mol_s > END_TOLERANCE * scale(top)
        {
            return Err("a total condenser has no vapor distillate".into());
        }
        if self.condenser.is_some() && top.vapor_draw_mol_s != 0.0 {
            return Err("the condenser vapor product is the top vapor flow, not a draw".into());
        }
        if self.reboiler && bottom.vapor_draw_mol_s != 0.0 {
            return Err("a reboiler vapor draw is not supported".into());
        }
        for (j, s) in cascade.stages.iter().enumerate() {
            let has_port = (j == 0 && self.condenser.is_some()) || (j == n - 1 && self.reboiler);
            if !has_port && s.duty_w != 0.0 {
                return Err(format!(
                    "stage {} has a duty but no energy port; side heaters and coolers are not available",
                    j + 1
                ));
            }
        }
        Ok(())
    }

    /// Values for the column's energy ports in port order (`condenser duty`
    /// outlet, then `reboiler duty` inlet, as present), positive in the port
    /// direction: `−Q_1` for the condenser and `Q_N` for the reboiler. This
    /// is `UnitOutput::energy`; the unit's `duty_w` is `Q_1 + Q_N`.
    pub fn port_duties(&self, cascade: &StageCascade) -> Result<Vec<Option<f64>>, String> {
        self.check(cascade)?;
        let mut out = Vec::new();
        if self.condenser.is_some() {
            out.push(Some(-cascade.stages[0].duty_w));
        }
        if self.reboiler {
            out.push(Some(cascade.stages[cascade.len() - 1].duty_w));
        }
        Ok(out)
    }
}

/// A packed section placed in the cascade: its ideal theoretical stages
/// ([`PackedStages`]) occupy the consecutive column stages
/// `top_stage..top_stage + stages` (one-based, numbered from the top like
/// every [`Stage`]). The packing module counts stages from zero within the
/// section; this is the only translation between the two numberings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlacedPackedSection {
    /// One-based column stage at the top of the bed.
    pub top_stage: usize,
    pub packed: PackedStages,
}

impl PlacedPackedSection {
    /// One-based column stages of the section.
    pub fn column_stages(&self) -> std::ops::RangeInclusive<usize> {
        self.top_stage..=self.top_stage + usize::from(self.packed.stages) - 1
    }

    /// One-based column stage receiving a feed or draw at `depth_m` below
    /// the top of the bed. A feed at the top face enters the first section
    /// stage, one at the bottom face the last; a point on an internal stage
    /// boundary belongs to the lower stage ([`PackedStages::stage_at_depth`]).
    /// The feed enters that stage as a whole ([`StageFeed`]).
    pub fn column_stage_at_depth(&self, depth_m: f64) -> Result<usize, String> {
        Ok(self.top_stage + usize::from(self.packed.stage_at_depth(depth_m)?))
    }

    /// Checks the section against a cascade: within its stages, clear of a
    /// condenser (stage 1) and reboiler (stage N), which are separate
    /// equilibrium stages, and with Murphree efficiency
    /// [`PackedStages::stage_efficiency`] (1) on every packed stage, since
    /// the HETP already accounts for non-ideal contacting.
    pub fn check(&self, cascade: &StageCascade, ends: ColumnEnds) -> Result<(), String> {
        let n = cascade.len();
        let range = self.column_stages();
        let first = if ends.condenser.is_some() { 2 } else { 1 };
        let last = if ends.reboiler {
            n.saturating_sub(1)
        } else {
            n
        };
        if self.top_stage < first || *range.end() > last {
            return Err(format!(
                "packed stages {} to {} lie outside the column stages {first} to {last} \
                 available to packing",
                range.start(),
                range.end()
            ));
        }
        let eta = self.packed.stage_efficiency();
        for j in range {
            if cascade.stages[j - 1].efficiency != eta {
                return Err(format!(
                    "packed stage {j} has efficiency {}; HETP stages are ideal",
                    cascade.stages[j - 1].efficiency
                ));
            }
        }
        Ok(())
    }
}

/// A fixed terminal temperature approach at a condenser or reboiler against
/// a utility inlet temperature, as for the exchangers: the condenser outlet
/// is the hot outlet facing the coolant inlet (`T_1 − T_c,in = ΔT`), the
/// reboiler liquid is the cold outlet facing the heating-medium inlet
/// (`T_h,in − T_N = ΔT`). The approach becomes an end-stage temperature
/// specification that replaces a duty or temperature specification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EndApproach {
    Condenser {
        coolant_inlet_k: f64,
        approach_k: f64,
    },
    Reboiler {
        heating_inlet_k: f64,
        approach_k: f64,
    },
}

impl EndApproach {
    /// The end-stage temperature in K the approach prescribes.
    pub fn stage_temperature_k(&self) -> Result<f64, String> {
        let (utility, dt, t) = match *self {
            Self::Condenser {
                coolant_inlet_k,
                approach_k,
            } => (coolant_inlet_k, approach_k, coolant_inlet_k + approach_k),
            Self::Reboiler {
                heating_inlet_k,
                approach_k,
            } => (heating_inlet_k, approach_k, heating_inlet_k - approach_k),
        };
        if !(utility.is_finite() && utility > 0.0) {
            return Err("utility inlet temperature must be positive".into());
        }
        if !(dt.is_finite() && dt >= 0.0) {
            return Err("terminal approach must be finite and nonnegative".into());
        }
        if t > 0.0 {
            Ok(t)
        } else {
            Err("terminal approach leaves no positive stage temperature".into())
        }
    }

    /// The attained approach of a solved cascade, hot minus cold at the
    /// terminal.
    pub fn attained_k(&self, cascade: &StageCascade) -> Option<f64> {
        match *self {
            Self::Condenser {
                coolant_inlet_k, ..
            } => cascade
                .stages
                .first()
                .map(|s| s.temperature_k - coolant_inlet_k),
            Self::Reboiler {
                heating_inlet_k, ..
            } => cascade
                .stages
                .last()
                .map(|s| heating_inlet_k - s.temperature_k),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::columns::packing::{PackedSectionBasis, PackedSectionSpec, StageRounding};

    fn cascade_of(n: usize) -> StageCascade {
        let s = Stage::equilibrium(350.0, 1e5, 1.0, 1.0, vec![0.5, 0.5], vec![0.5, 0.5]);
        StageCascade {
            stages: vec![s; n],
            feeds: vec![],
        }
    }

    /// Packed numbering (zero-based within the bed) meets column numbering
    /// (one-based from the top), condenser/reboiler exclusion and ideal
    /// packed stages.
    #[test]
    fn packed_sections_place_ideal_stages_in_column_numbering() {
        // 4 m of packing at HETP 1 m: four ideal stages.
        let packed = PackedSectionSpec {
            hetp_m: Some(1.0),
            basis: PackedSectionBasis::PackedHeight {
                height_m: Some(4.0),
                rounding: StageRounding::Nearest,
            },
        }
        .resolve()
        .unwrap();
        // Distillation: condenser 1, bed on stages 2..=5, reboiler 6.
        let mut column = cascade_of(6);
        let bed = PlacedPackedSection {
            top_stage: 2,
            packed,
        };
        assert_eq!(bed.column_stages(), 2..=5);
        assert_eq!(bed.column_stage_at_depth(0.0), Ok(2));
        assert_eq!(bed.column_stage_at_depth(2.0), Ok(4));
        assert_eq!(bed.column_stage_at_depth(4.0), Ok(5));
        assert!(bed.column_stage_at_depth(4.5).is_err());
        bed.check(&column, ColumnEnds::DISTILLATION).unwrap();
        for top_stage in [1, 3] {
            let off = PlacedPackedSection { top_stage, packed };
            assert!(
                off.check(&column, ColumnEnds::DISTILLATION)
                    .unwrap_err()
                    .contains("outside")
            );
        }
        column.stages[3].efficiency = 0.7;
        assert!(
            bed.check(&column, ColumnEnds::DISTILLATION)
                .unwrap_err()
                .contains("ideal")
        );
        // Absorber: the whole cascade is the bed; lean liquid at the top
        // face enters stage 1 and rich gas at the bottom face stage N.
        let absorber = PlacedPackedSection {
            top_stage: 1,
            packed,
        };
        absorber.check(&cascade_of(4), ColumnEnds::NONE).unwrap();
        assert_eq!(absorber.column_stage_at_depth(0.0), Ok(1));
        assert_eq!(absorber.column_stage_at_depth(4.0), Ok(4));
    }

    fn props(k: &[f64], hl: f64, hv: f64) -> StageProperties {
        StageProperties {
            k: k.to_vec(),
            liquid_enthalpy_j_mol: hl,
            vapor_enthalpy_j_mol: hv,
        }
    }

    /// One equilibrium stage with a single feed is an isothermal flash:
    /// with constant K the Rachford–Rice solution must give zero residuals.
    #[test]
    fn single_stage_flash_closes_with_rachford_rice_split() {
        let (z, k) = ([0.5, 0.5], [2.0, 0.5]);
        // Rachford–Rice: Σ z (K − 1) / (1 + β (K − 1)) = 0 gives β = 0.5.
        let beta = 0.5;
        let x: Vec<f64> = (0..2).map(|i| z[i] / (1.0 + beta * (k[i] - 1.0))).collect();
        let y: Vec<f64> = (0..2).map(|i| k[i] * x[i]).collect();
        let cascade = StageCascade {
            stages: vec![Stage {
                duty_w: 100.0 * (0.5 * -1000.0 + 0.5 * 30_000.0 - 2_000.0),
                ..Stage::equilibrium(350.0, 1e5, 50.0, 50.0, x.clone(), y.clone())
            }],
            feeds: vec![StageFeed {
                stage: 1,
                molar_flow_mol_s: 100.0,
                z: z.to_vec(),
                molar_enthalpy_j_mol: 2_000.0,
            }],
        };
        let r = &cascade.residuals(&[props(&k, -1000.0, 30_000.0)]).unwrap()[0];
        for v in r.component_mol_s.iter().chain(&r.equilibrium) {
            assert!(v.abs() < 1e-12, "{r:?}");
        }
        assert!(r.liquid_summation.abs() < 1e-15 && r.vapor_summation.abs() < 1e-15);
        assert!(r.energy_w.abs() < 1e-9, "{r:?}");
        assert_eq!(cascade.total_flow_residuals(), [0.0]);
        let solved = cascade.solve_liquid_fractions(&[k.to_vec()]).unwrap();
        for i in 0..2 {
            assert!((solved[0][i] - x[i]).abs() < 1e-14);
        }
    }

    #[test]
    fn murphree_efficiency_mixes_equilibrium_and_entering_vapor() {
        let mut top = Stage::equilibrium(350.0, 1e5, 1.0, 1.0, vec![0.4, 0.6], vec![0.0; 2]);
        top.efficiency = 0.6;
        let mut below = Stage::equilibrium(355.0, 1e5, 1.0, 1.0, vec![0.3, 0.7], vec![0.45, 0.55]);
        below.efficiency = 0.5;
        let k = [2.0, 0.5];
        // y = 0.6 K x + 0.4 y_below.
        top.y = vec![0.6 * 0.8 + 0.4 * 0.45, 0.6 * 0.3 + 0.4 * 0.55];
        let cascade = StageCascade {
            stages: vec![top, below],
            feeds: vec![],
        };
        let r = cascade
            .residuals(&[props(&k, 0.0, 0.0), props(&k, 0.0, 0.0)])
            .unwrap();
        assert!(r[0].equilibrium.iter().all(|e| e.abs() < 1e-15));
        // The bottom stage ignores its efficiency (no vapor from below).
        assert!((r[1].equilibrium[0] - (2.0 * 0.3 - 0.45)).abs() < 1e-15);
        assert!(
            cascade
                .solve_liquid_fractions(&[k.to_vec(), k.to_vec()])
                .is_err()
        );
    }

    #[test]
    fn invalid_cascades_fail_explicitly() {
        let stage = Stage::equilibrium(350.0, 1e5, 1.0, 1.0, vec![0.5, 0.5], vec![0.5, 0.5]);
        let p = props(&[1.0, 1.0], 0.0, 0.0);
        let empty = StageCascade {
            stages: vec![],
            feeds: vec![],
        };
        assert!(empty.residuals(&[]).is_err());
        let mut bad = StageCascade {
            stages: vec![stage.clone()],
            feeds: vec![StageFeed {
                stage: 2,
                molar_flow_mol_s: 1.0,
                z: vec![0.5, 0.5],
                molar_enthalpy_j_mol: 0.0,
            }],
        };
        assert!(
            bad.residuals(std::slice::from_ref(&p))
                .unwrap_err()
                .contains("feed stage 2")
        );
        bad.feeds.clear();
        assert!(bad.residuals(&[]).unwrap_err().contains("property sets"));
        bad.stages[0].vapor_mol_s = -1.0;
        assert!(bad.validate(2).unwrap_err().contains("vapor flow"));
        bad.stages[0] = Stage {
            efficiency: 0.0,
            ..stage.clone()
        };
        assert!(bad.validate(2).unwrap_err().contains("efficiency"));
        let stopped = StageCascade {
            stages: vec![Stage::equilibrium(
                350.0,
                1e5,
                0.0,
                0.0,
                vec![0.0; 2],
                vec![0.0; 2],
            )],
            feeds: vec![],
        };
        assert!(
            stopped
                .solve_liquid_fractions(&[vec![1.0, 1.0]])
                .unwrap_err()
                .contains("singular")
        );
        let flowing_without_state = MaterialStream {
            molar_flow_mol_s: 1.0,
            state: None,
        };
        assert!(StageFeed::from_stream(1, &flowing_without_state, 2).is_err());
        let unconnected = MaterialStream {
            molar_flow_mol_s: 0.0,
            state: None,
        };
        assert_eq!(
            StageFeed::from_stream(1, &unconnected, 2).unwrap().z,
            [0.0; 2]
        );
    }

    fn two_stage(top_vapor: f64, duties: [f64; 2]) -> StageCascade {
        let s = |v, q| Stage {
            duty_w: q,
            ..Stage::equilibrium(350.0, 1e5, 1.0, v, vec![0.5, 0.5], vec![0.5, 0.5])
        };
        StageCascade {
            stages: vec![s(top_vapor, duties[0]), s(1.0, duties[1])],
            feeds: vec![],
        }
    }

    #[test]
    fn end_duties_follow_the_energy_port_convention() {
        let column = two_stage(0.0, [-5_000.0, 6_000.0]);
        assert_eq!(
            ColumnEnds::DISTILLATION.port_duties(&column).unwrap(),
            [Some(5_000.0), Some(6_000.0)]
        );
        let condenser_only = ColumnEnds {
            condenser: Some(Condenser::Partial),
            reboiler: false,
        };
        let partial = two_stage(0.3, [-5_000.0, 0.0]);
        assert_eq!(
            condenser_only.port_duties(&partial).unwrap(),
            [Some(5_000.0)]
        );
        assert!(
            ColumnEnds::DISTILLATION
                .check(&partial)
                .unwrap_err()
                .contains("total condenser")
        );
        assert!(
            ColumnEnds::NONE
                .check(&column)
                .unwrap_err()
                .contains("no energy port")
        );
        assert_eq!(
            ColumnEnds::NONE
                .port_duties(&two_stage(1.0, [0.0; 2]))
                .unwrap(),
            []
        );
    }

    #[test]
    fn fixed_end_approaches_map_to_stage_temperatures() {
        let condenser = EndApproach::Condenser {
            coolant_inlet_k: 300.0,
            approach_k: 10.0,
        };
        let reboiler = EndApproach::Reboiler {
            heating_inlet_k: 420.0,
            approach_k: 15.0,
        };
        assert_eq!(condenser.stage_temperature_k().unwrap(), 310.0);
        assert_eq!(reboiler.stage_temperature_k().unwrap(), 405.0);
        let mut column = two_stage(0.0, [0.0; 2]);
        column.stages[0].temperature_k = 310.0;
        column.stages[1].temperature_k = 405.0;
        assert_eq!(condenser.attained_k(&column), Some(10.0));
        assert_eq!(reboiler.attained_k(&column), Some(15.0));
        let negative = EndApproach::Condenser {
            coolant_inlet_k: 300.0,
            approach_k: -1.0,
        };
        assert!(negative.stage_temperature_k().is_err());
        let impossible = EndApproach::Reboiler {
            heating_inlet_k: 300.0,
            approach_k: 300.0,
        };
        assert!(impossible.stage_temperature_k().is_err());
    }
}
