//! MESH solve of a rigorous distillation column (#186): total or partial
//! condenser, partial reboiler, condenser and stage pressure drops, liquid
//! and vapor side draws, and two column specifications.
//!
//! The specification set follows the Aspen Plus RadFrac / HYSYS column
//! convention of two degrees of freedom once the feeds, draws, pressures and
//! stages are fixed: one [`TopSpecification`] and one
//! [`BottomSpecification`]. Specifications that fix the distillate and
//! reflux flows directly are substituted into the flows; every other
//! specification adds one residual equation and frees the distillate or the
//! reflux as an additional Newton unknown.

use property_packages::SaturationKind;
use thermo_types::Phase;

use super::mesh::{self, TEMPERATURE_SCALE};
use super::side_sections::{Cooling, SideSections, SideStripperProfile, Stripping};
use super::stages::{
    ColumnEnds, Condenser, PlacedPackedSection, Stage, StageCascade, StageFeed, StageProperties,
};
use super::{ColumnInternals, count_issue};
use crate::MaterialStream;
use crate::numerics::solve_dense;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional, check_value};
use crate::units::{CalcContext, UnitOutput};

/// Same enthalpy scale as the MESH energy residuals.
const ENTHALPY_SCALE: f64 = 10_000.0;
/// Scaled specification residual at which the outer flow iteration hands
/// over to the simultaneous Newton solve.
const OUTER_TOLERANCE: f64 = 1e-8;
const OUTER_ITERATIONS: usize = 40;
/// Starting liquid fraction, relative to its feed fraction, of a component
/// that is absent from the seed inlets (stripping steam).
const ABSENT_LIQUID_SEED: f64 = 0.1;

/// Liquid side draw from an internal stage, in port order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LiquidSideDraw {
    /// One-based stage number, including the condenser as stage 1.
    pub stage: u16,
    /// Prescribed molar withdrawal in mol/s.
    pub molar_flow_mol_s: Option<f64>,
}

impl Default for LiquidSideDraw {
    fn default() -> Self {
        Self {
            stage: 2,
            molar_flow_mol_s: None,
        }
    }
}

/// Vapor side draw from an internal stage (port `vapor draw k`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VaporSideDraw {
    /// One-based stage number, including the condenser as stage 1. Only
    /// internal stages (not the condenser or reboiler) may carry a vapor draw.
    pub stage: u16,
    /// Prescribed molar withdrawal in mol/s.
    pub molar_flow_mol_s: Option<f64>,
}

impl Default for VaporSideDraw {
    fn default() -> Self {
        Self {
            stage: 2,
            molar_flow_mol_s: None,
        }
    }
}

/// Condenser type (Aspen Plus RadFrac naming).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub enum DistillationCondenser {
    /// All overhead vapor condensed to its bubble point; liquid distillate on
    /// the `distillate` port.
    #[default]
    Total,
    /// Vapor distillate on the `distillate` port, all condensate refluxed.
    PartialVapor,
    /// Vapor distillate on the `distillate` port and a liquid distillate at a
    /// specified molar rate on the `liquid distillate` port.
    PartialVaporLiquid {
        /// Liquid distillate, mol/s; positive.
        liquid_distillate_mol_s: Option<f64>,
    },
}

impl DistillationCondenser {
    fn stage_condenser(self) -> Condenser {
        match self {
            Self::Total => Condenser::Total,
            _ => Condenser::Partial,
        }
    }
}

/// Pressure drop below the condenser's first tray.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColumnPressureDrop {
    /// Drop per stage, Pa, between consecutive stages from the first tray to
    /// the reboiler.
    PerStage(f64),
    /// Drop from the first tray to the reboiler, Pa, distributed uniformly.
    Column(f64),
}

impl Default for ColumnPressureDrop {
    fn default() -> Self {
        Self::PerStage(0.0)
    }
}

/// Mole fraction or recovery target of one component.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ComponentTarget {
    /// Package component identifier.
    pub component: Option<String>,
    /// Mole fraction in the product, or fraction of the component's total
    /// feed recovered in the product; strictly between 0 and 1.
    pub value: Option<f64>,
}

/// First column specification (Aspen "top"/condenser specification).
#[derive(Debug, Clone, PartialEq)]
pub enum TopSpecification {
    /// Molar reflux over total molar distillate, `L_1 / D`.
    RefluxRatio(Option<f64>),
    /// Molar reflux `L_1`, mol/s.
    RefluxRate(Option<f64>),
    /// Total molar distillate `D` (vapor plus liquid), mol/s.
    DistillateRate(Option<f64>),
    /// Heat removed by the condenser, W; positive.
    CondenserDuty(Option<f64>),
    /// Mole fraction of a component in the total distillate.
    DistillatePurity(ComponentTarget),
    /// Fraction of a component's feed recovered in the total distillate.
    DistillateRecovery(ComponentTarget),
}

impl Default for TopSpecification {
    fn default() -> Self {
        Self::RefluxRatio(None)
    }
}

/// Second column specification (Aspen "bottom"/reboiler specification).
#[derive(Debug, Clone, PartialEq)]
pub enum BottomSpecification {
    /// Total distillate over total feed, strictly between 0 and 1.
    DistillateFraction(Option<f64>),
    /// Bottoms molar flow `L_N`, mol/s.
    BottomsRate(Option<f64>),
    /// Heat supplied to the reboiler, W; positive.
    ReboilerDuty(Option<f64>),
    /// Boilup over bottoms, `V_N / L_N`.
    BoilupRatio(Option<f64>),
    /// Mole fraction of a component in the bottoms.
    BottomsPurity(ComponentTarget),
    /// Fraction of a component's feed recovered in the bottoms.
    BottomsRecovery(ComponentTarget),
}

impl Default for BottomSpecification {
    fn default() -> Self {
        Self::DistillateFraction(None)
    }
}

/// Rigorous distillation column specification (#186).
#[derive(Debug, Clone, PartialEq)]
pub struct DistillationColumnSpec {
    /// Number of material feed ports.
    pub feeds: u16,
    /// Number of liquid side-draw ports (`side draw k`).
    pub side_draws: u16,
    /// One-based feed stage per port. Unset entries default to stage 2.
    pub feed_stages: Vec<u16>,
    /// Specification of each side-draw port. Unset entries have no flow specification.
    pub liquid_draws: Vec<LiquidSideDraw>,
    /// Number of vapor side-draw ports (`vapor draw k`).
    pub vapor_side_draws: u16,
    /// Specification of each vapor side-draw port.
    pub vapor_draws: Vec<VaporSideDraw>,
    /// Ideal trays or ideal packed stages.
    pub internals: ColumnInternals,
    /// Condenser (top stage) absolute pressure, Pa. Without pressure drops
    /// this is the uniform column pressure.
    pub pressure_pa: Option<f64>,
    /// Pressure drop from the condenser to the first tray, Pa.
    pub condenser_pressure_drop_pa: f64,
    /// Pressure drop over the trays or packed stages and the reboiler.
    pub pressure_drop: ColumnPressureDrop,
    /// Condenser type.
    pub condenser: DistillationCondenser,
    /// First specification.
    pub top: TopSpecification,
    /// Second specification.
    pub bottom: BottomSpecification,
}

impl Default for DistillationColumnSpec {
    fn default() -> Self {
        Self {
            feeds: 1,
            side_draws: 0,
            feed_stages: Vec::new(),
            liquid_draws: Vec::new(),
            vapor_side_draws: 0,
            vapor_draws: Vec::new(),
            internals: ColumnInternals::Trays(5),
            pressure_pa: None,
            condenser_pressure_drop_pa: 0.0,
            pressure_drop: ColumnPressureDrop::default(),
            condenser: DistillationCondenser::Total,
            top: TopSpecification::default(),
            bottom: BottomSpecification::default(),
        }
    }
}

fn target_issues(out: &mut Vec<SpecIssue>, field: &'static str, target: &ComponentTarget) {
    if target.component.as_deref().is_none_or(str::is_empty) {
        out.push(SpecIssue {
            field,
            problem: SpecProblem::Missing,
        });
    }
    open_fraction(out, field, target.value);
}

fn open_fraction(out: &mut Vec<SpecIssue>, field: &'static str, value: Option<f64>) {
    check_optional(out, field, value, Sign::Fraction);
    if matches!(value, Some(0.0 | 1.0)) {
        out.push(SpecIssue {
            field,
            problem: SpecProblem::OutOfRange,
        });
    }
}

impl DistillationColumnSpec {
    /// Total column stages and optional resolved bed, including both ends.
    pub(super) fn stage_layout(&self) -> Result<(usize, Option<PlacedPackedSection>), String> {
        let (internal, packed) = self.internals.resolve()?;
        Ok((
            internal + 2,
            packed.map(|packed| PlacedPackedSection {
                top_stage: 2,
                packed,
            }),
        ))
    }

    /// Absolute pressure of every stage, condenser first.
    pub fn stage_pressures_pa(&self, stages: usize) -> Option<Vec<f64>> {
        let top = self.pressure_pa?;
        let first_tray = top + self.condenser_pressure_drop_pa;
        let per_stage = match self.pressure_drop {
            ColumnPressureDrop::PerStage(dp) => dp,
            ColumnPressureDrop::Column(dp) => dp / stages.saturating_sub(2).max(1) as f64,
        };
        Some(
            (0..stages)
                .map(|j| match j {
                    0 => top,
                    _ => first_tray + (j - 1) as f64 * per_stage,
                })
                .collect(),
        )
    }

    /// Whether the two specifications both fix the distillate flow.
    fn overspecified(&self) -> bool {
        matches!(self.top, TopSpecification::DistillateRate(_))
            && matches!(
                self.bottom,
                BottomSpecification::DistillateFraction(_) | BottomSpecification::BottomsRate(_)
            )
    }

    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        self.issues_with(BottomEnd::Reboiler)
    }

    /// Issues of the specification; an adiabatic bottom ignores `bottom`.
    pub(super) fn issues_with(&self, end: BottomEnd) -> Vec<SpecIssue> {
        let reboiler = end == BottomEnd::Reboiler;
        let mut issues = Vec::new();
        count_issue(&mut issues, "feeds", self.feeds, 1);
        count_issue(&mut issues, "side draws", self.side_draws, 0);
        count_issue(&mut issues, "vapor side draws", self.vapor_side_draws, 0);
        if self.feed_stages.len() > usize::from(self.feeds) {
            issues.push(SpecIssue {
                field: "feed stages",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        if self.liquid_draws.len() > usize::from(self.side_draws) {
            issues.push(SpecIssue {
                field: "side-draw specifications",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        if self.vapor_draws.len() > usize::from(self.vapor_side_draws) {
            issues.push(SpecIssue {
                field: "vapor-draw specifications",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        check_optional(
            &mut issues,
            "column pressure",
            self.pressure_pa,
            Sign::Positive,
        );
        check_value(
            &mut issues,
            "condenser pressure drop",
            self.condenser_pressure_drop_pa,
            Sign::NonNegative,
        );
        let (ColumnPressureDrop::PerStage(dp) | ColumnPressureDrop::Column(dp)) =
            self.pressure_drop;
        check_value(&mut issues, "stage pressure drop", dp, Sign::NonNegative);
        if let DistillationCondenser::PartialVaporLiquid {
            liquid_distillate_mol_s,
        } = self.condenser
        {
            check_optional(
                &mut issues,
                "liquid distillate rate",
                liquid_distillate_mol_s,
                Sign::Positive,
            );
        }
        match &self.top {
            TopSpecification::RefluxRatio(v) => {
                check_optional(&mut issues, "reflux ratio", *v, Sign::Positive)
            }
            TopSpecification::RefluxRate(v) => {
                check_optional(&mut issues, "reflux rate", *v, Sign::Positive)
            }
            TopSpecification::DistillateRate(v) => {
                check_optional(&mut issues, "distillate rate", *v, Sign::Positive)
            }
            TopSpecification::CondenserDuty(v) => {
                check_optional(&mut issues, "condenser duty", *v, Sign::Positive)
            }
            TopSpecification::DistillatePurity(t) => {
                target_issues(&mut issues, "distillate mole fraction", t)
            }
            TopSpecification::DistillateRecovery(t) => {
                target_issues(&mut issues, "distillate recovery", t)
            }
        }
        match &self.bottom {
            _ if !reboiler => {}
            BottomSpecification::DistillateFraction(v) => {
                open_fraction(&mut issues, "distillate fraction", *v)
            }
            BottomSpecification::BottomsRate(v) => {
                check_optional(&mut issues, "bottoms rate", *v, Sign::Positive)
            }
            BottomSpecification::ReboilerDuty(v) => {
                check_optional(&mut issues, "reboiler duty", *v, Sign::Positive)
            }
            BottomSpecification::BoilupRatio(v) => {
                check_optional(&mut issues, "boilup ratio", *v, Sign::Positive)
            }
            BottomSpecification::BottomsPurity(t) => {
                target_issues(&mut issues, "bottoms mole fraction", t)
            }
            BottomSpecification::BottomsRecovery(t) => {
                target_issues(&mut issues, "bottoms recovery", t)
            }
        }
        if reboiler && self.overspecified() {
            issues.push(SpecIssue {
                field: "bottom specification",
                problem: SpecProblem::Duplicate,
            });
        }
        if let Ok((count, _)) = self.stage_layout() {
            if self
                .feed_stages
                .iter()
                .take(usize::from(self.feeds))
                .any(|&s| !(2..=count).contains(&usize::from(s)))
            {
                issues.push(SpecIssue {
                    field: "feed stage",
                    problem: SpecProblem::OutOfRange,
                });
            }
            if self
                .liquid_draws
                .iter()
                .take(usize::from(self.side_draws))
                .any(|draw| !(2..count).contains(&usize::from(draw.stage)))
            {
                issues.push(SpecIssue {
                    field: "side-draw stage",
                    problem: SpecProblem::OutOfRange,
                });
            }
            if self
                .vapor_draws
                .iter()
                .take(usize::from(self.vapor_side_draws))
                .any(|draw| !(2..count).contains(&usize::from(draw.stage)))
            {
                issues.push(SpecIssue {
                    field: "vapor-draw stage",
                    problem: SpecProblem::OutOfRange,
                });
            }
        }
        issues.extend(self.internals.issues());
        if self.side_draws <= super::MAX_COLUMN_PORTS {
            for i in 0..usize::from(self.side_draws) {
                check_optional(
                    &mut issues,
                    "side-draw molar flow",
                    self.liquid_draws
                        .get(i)
                        .and_then(|draw| draw.molar_flow_mol_s),
                    Sign::Positive,
                );
            }
        }
        if self.vapor_side_draws <= super::MAX_COLUMN_PORTS {
            for i in 0..usize::from(self.vapor_side_draws) {
                check_optional(
                    &mut issues,
                    "vapor-draw molar flow",
                    self.vapor_draws
                        .get(i)
                        .and_then(|draw| draw.molar_flow_mol_s),
                    Sign::Positive,
                );
            }
        }
        issues
    }
}

/// Accepted solution evidence for a rigorous distillation column.
#[derive(Debug, Clone, PartialEq)]
pub struct DistillationColumnResult {
    /// Stage temperatures from the condenser to the reboiler.
    pub stage_temperatures_k: Vec<f64>,
    /// Stage pressures from the condenser to the reboiler, Pa.
    pub stage_pressures_pa: Vec<f64>,
    /// Heat removed by the condenser, W.
    pub condenser_duty_w: f64,
    /// Heat supplied to the reboiler, W.
    pub reboiler_duty_w: f64,
    /// Total distillate (vapor plus liquid), mol/s.
    pub distillate_mol_s: f64,
    /// Reflux `L_1`, mol/s.
    pub reflux_mol_s: f64,
    /// Bottoms `L_N`, mol/s.
    pub bottoms_mol_s: f64,
    /// Boilup `V_N`, mol/s.
    pub boilup_mol_s: f64,
    /// Bed height supplied or derived from the specified HETP.
    pub packed_height_m: Option<f64>,
    /// Accepted Newton iterations of the final solve.
    pub iterations: usize,
    /// TBP/D86 curves of the products (petroleum fractionator, #190);
    /// empty for other columns.
    pub product_curves: Vec<super::ProductCurve>,
    /// Heat removed by each pump-around of a petroleum fractionator (#190),
    /// W; empty without pump-arounds.
    pub pump_around_duties_w: Vec<f64>,
    /// Profile of each side stripper of a petroleum fractionator (#190);
    /// empty without side strippers.
    pub side_strippers: Vec<SideStripperProfile>,
}

/// Bottom end of a column solved by [`calculate_with`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) enum BottomEnd {
    /// Partial reboiler with the specification's bottom specification.
    #[default]
    Reboiler,
    /// Adiabatic bottom stage without a reboiler, for example stripped by
    /// steam entering as a bottom feed (#190). The bottom specification is
    /// not used: the zero bottom duty takes its place, so the distillate
    /// becomes a Newton unknown.
    Adiabatic,
}

/// Options of [`calculate_with`] beyond the specification.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct ColumnOptions {
    pub(super) bottom: BottomEnd,
    /// Number of leading inlets whose blend sets the bubble-point start of
    /// the profile; all inlets when `None`. A stripping-steam inlet is left
    /// out so that the start is the hydrocarbon feed's bubble point.
    pub(super) seed_inlets: Option<usize>,
}

/// How the solve obtains the total distillate.
#[derive(Debug, Clone, Copy)]
enum DistillateFlow {
    Fixed(f64),
    /// Newton unknown `u`, `D = A / (1 + e^−u)` with `A` the feed less draws.
    Free,
}

/// How the solve obtains the reflux.
#[derive(Debug, Clone, Copy)]
enum RefluxFlow {
    Fixed(f64),
    Ratio(f64),
    /// Newton unknown `v`, `L_1 = F e^v`.
    Free,
}

/// A specification that is a residual equation rather than a flow.
#[derive(Debug, Clone, Copy)]
enum SpecEquation {
    CondenserDuty(f64),
    ReboilerDuty(f64),
    BoilupRatio(f64),
    DistillatePurity(usize, f64),
    DistillateRecovery(usize, f64),
    BottomsPurity(usize, f64),
    BottomsRecovery(usize, f64),
}

struct Decoded {
    cascade: StageCascade,
    distillate: f64,
    reflux: f64,
    bottoms: f64,
    /// Side-stripper cascades, in specification order.
    strippers: Vec<StageCascade>,
}

/// A decoded column with its stage properties and every coupling feed's
/// enthalpy set.
struct Evaluated {
    decoded: Decoded,
    properties: Vec<StageProperties>,
    stripper_properties: Vec<Vec<StageProperties>>,
}

#[derive(Clone)]
struct Column<'a, 'p> {
    ctx: &'a CalcContext<'p>,
    pressures: Vec<f64>,
    temperature_range: [f64; 2],
    feeds: Vec<StageFeed>,
    /// Total feed per stage.
    feed_by_stage: Vec<f64>,
    /// Component feed flows, mol/s.
    feed_components: Vec<f64>,
    liquid_draws: Vec<f64>,
    vapor_draws: Vec<f64>,
    active_components: Vec<usize>,
    flow: f64,
    /// Feed less side draws: the upper bound of the distillate.
    available: f64,
    condenser: Condenser,
    liquid_distillate: f64,
    distillate: DistillateFlow,
    reflux: RefluxFlow,
    equations: Vec<SpecEquation>,
    /// Pump-arounds and side strippers solved with the column (#190).
    sections: SideSections,
}

impl Column<'_, '_> {
    fn ends(&self) -> ColumnEnds {
        ColumnEnds {
            condenser: Some(self.condenser),
            reboiler: true,
        }
    }

    fn partial(&self) -> bool {
        self.condenser == Condenser::Partial
    }

    /// Number of profile variables, before any free flow: the main column,
    /// then each side stripper.
    fn profile_len(&self) -> usize {
        let m = self.active_components.len() - 1;
        self.main_profile_len()
            + self
                .sections
                .strippers
                .iter()
                .map(|s| s.variable_len(m))
                .sum::<usize>()
    }

    /// Number of main-column profile variables.
    fn main_profile_len(&self) -> usize {
        let n = self.pressures.len();
        let m = self.active_components.len() - 1;
        let condenser = 1 + m + if self.partial() { m } else { 0 };
        condenser + (n - 1) * (1 + 2 * m) + (n - 2)
    }

    fn flows(&self, extra: &[f64]) -> (f64, f64) {
        let mut extra = extra.iter();
        let distillate = match self.distillate {
            DistillateFlow::Fixed(d) => d,
            DistillateFlow::Free => {
                let u = *extra.next().expect("free distillate variable");
                self.available / (1.0 + (-u).exp())
            }
        };
        let reflux = match self.reflux {
            RefluxFlow::Fixed(l) => l,
            RefluxFlow::Ratio(r) => distillate * r,
            RefluxFlow::Free => self.flow * extra.next().expect("free reflux variable").exp(),
        };
        (distillate, reflux)
    }

    /// Free-flow variables for the given flows.
    fn flow_variables(&self, distillate: f64, reflux: f64) -> Vec<f64> {
        let mut out = Vec::new();
        if matches!(self.distillate, DistillateFlow::Free) {
            let s = distillate / self.available;
            out.push((s / (1.0 - s)).ln());
        }
        if matches!(self.reflux, RefluxFlow::Free) {
            out.push((reflux / self.flow).ln());
        }
        out
    }

    /// The starting column of a column with side strippers: each stripper
    /// replaced by its bottoms drawn at its draw stage and its steam fed to
    /// its return stage, which keeps the external balance.
    fn precursor(&self) -> Self {
        let mut precursor = self.clone();
        precursor.sections.strippers.clear();
        for stripper in &self.sections.strippers {
            precursor.liquid_draws[stripper.draw] += stripper.bottoms;
            if let Stripping::Steam(steam) = &stripper.bottom {
                precursor.feeds.push(StageFeed {
                    stage: stripper.ret + 1,
                    ..steam.clone()
                });
                precursor.feed_by_stage[stripper.ret] += steam.molar_flow_mol_s;
            }
        }
        precursor
    }

    /// The same column with both flows fixed.
    fn with_fixed_flows(&self, distillate: f64, reflux: f64) -> Self {
        let mut fixed = self.clone();
        fixed.distillate = DistillateFlow::Fixed(distillate);
        fixed.reflux = RefluxFlow::Fixed(reflux);
        fixed.equations.clear();
        fixed
    }

    fn decode(&self, values: &[f64]) -> Result<Decoded, String> {
        let n = self.pressures.len();
        let c = self.ctx.component_count();
        let m = self.active_components.len() - 1;
        let (profile, extra) = values.split_at(self.profile_len());
        let (profile, mut section_values) = profile.split_at(self.main_profile_len());
        let mut stripper_values = Vec::with_capacity(self.sections.strippers.len());
        for stripper in &self.sections.strippers {
            let (own, rest) = section_values.split_at(stripper.variable_len(m));
            stripper_values.push(own);
            section_values = rest;
        }
        let (distillate, reflux) = self.flows(extra);
        if !(distillate.is_finite() && reflux.is_finite() && distillate > 0.0 && reflux > 0.0) {
            return Err("column distillate and reflux must be positive and finite".into());
        }
        let vapor_distillate = distillate - self.liquid_distillate;
        if self.partial() && vapor_distillate <= 0.0 {
            return Err(
                "column vapor distillate must be positive; the liquid distillate takes the whole distillate"
                    .into(),
            );
        }
        // Pump-around and side-stripper streams are internal draws and feeds.
        let mut feed_by_stage = self.feed_by_stage.clone();
        let mut liquid_draws = self.liquid_draws.clone();
        for pump in &self.sections.pump_arounds {
            liquid_draws[pump.draw] += pump.flow;
            feed_by_stage[pump.ret] += pump.flow;
        }
        let mut tops = Vec::with_capacity(self.sections.strippers.len());
        for (k, (stripper, own)) in self
            .sections
            .strippers
            .iter()
            .zip(&stripper_values)
            .enumerate()
        {
            let (vapor, y) = stripper.top_vapor(own, &self.active_components, c, self.flow);
            let draw = stripper.draw_flow(vapor);
            if !(draw.is_finite() && draw > 0.0) {
                return Err(format!(
                    "side stripper {} draws no liquid from the column: its top vapor and bottoms are below its steam",
                    k + 1
                ));
            }
            liquid_draws[stripper.draw] += draw;
            feed_by_stage[stripper.ret] += vapor;
            tops.push((vapor, y, draw));
        }
        let mut remaining = -distillate;
        let liquid_offset: Vec<_> = (0..n)
            .map(|j| {
                remaining += feed_by_stage[j] - liquid_draws[j] - self.vapor_draws[j];
                remaining
            })
            .collect();
        let bottoms = liquid_offset[n - 1];
        if !bottoms.is_finite() || bottoms <= 0.0 {
            return Err("column withdrawals must leave positive bottoms flow".into());
        }
        let expand = |v| mesh::fractions(v, &self.active_components, c);
        let mut tail = profile;
        let mut liquid = vec![0.0; n];
        let mut stages: Vec<Stage> = Vec::with_capacity(n);
        liquid[0] = reflux;
        liquid[n - 1] = bottoms;
        for j in 0..n {
            let t = tail[0] * TEMPERATURE_SCALE;
            tail = &tail[1..];
            let x = if j == 0 && !self.partial() {
                Vec::new()
            } else {
                let (x, rest) = tail.split_at(m);
                tail = rest;
                expand(x)
            };
            let (y, rest) = tail.split_at(m);
            tail = rest;
            let vapor = match j {
                0 if self.partial() => vapor_distillate,
                0 => 0.0,
                1 => reflux + distillate,
                _ => {
                    let extra = self.flow * tail[0].exp();
                    tail = &tail[1..];
                    let offset = liquid_offset[j - 1];
                    liquid[j - 1] = extra + offset.max(0.0);
                    extra + (-offset).max(0.0)
                }
            };
            let mut stage = Stage::equilibrium(t, self.pressures[j], 0.0, vapor, x, expand(y));
            stage.liquid_draw_mol_s = match j {
                0 if self.partial() => self.liquid_distillate,
                0 => distillate,
                _ => liquid_draws[j],
            };
            stage.vapor_draw_mol_s = self.vapor_draws[j];
            stages.push(stage);
        }
        if !self.partial() {
            stages[0].x = stages[1].y.clone();
        }
        for (s, l) in stages.iter_mut().zip(liquid) {
            s.liquid_mol_s = l;
        }
        // Coupling feeds carry placeholder enthalpies until `evaluate`.
        let mut feeds = self.feeds.clone();
        for pump in &self.sections.pump_arounds {
            feeds.push(StageFeed {
                stage: pump.ret + 1,
                molar_flow_mol_s: pump.flow,
                z: stages[pump.draw].x.clone(),
                molar_enthalpy_j_mol: 0.0,
            });
        }
        let mut strippers = Vec::with_capacity(tops.len());
        for ((stripper, own), (vapor, y, draw)) in self
            .sections
            .strippers
            .iter()
            .zip(stripper_values)
            .zip(tops)
        {
            feeds.push(StageFeed {
                stage: stripper.ret + 1,
                molar_flow_mol_s: vapor,
                z: y,
                molar_enthalpy_j_mol: 0.0,
            });
            strippers.push(stripper.decode(
                own,
                &self.active_components,
                c,
                self.flow,
                self.pressures[stripper.draw],
                draw,
                stages[stripper.draw].x.clone(),
            )?);
        }
        Ok(Decoded {
            cascade: StageCascade { stages, feeds },
            distillate,
            reflux,
            bottoms,
            strippers,
        })
    }

    /// Decodes `values` and evaluates the stage properties of the column and
    /// its side strippers, setting the enthalpy of every coupling feed: a
    /// side stripper's liquid feed at its draw-stage liquid, its vapor return
    /// at its top-stage vapor and a pump-around return from its cooling.
    fn evaluate(&self, values: &[f64]) -> Result<Evaluated, String> {
        let mut decoded = self.decode(values)?;
        let properties = mesh::properties(
            self.ctx,
            &decoded.cascade,
            self.ends(),
            self.temperature_range,
        )?;
        let mut stripper_properties = Vec::with_capacity(decoded.strippers.len());
        for (stripper, cascade) in self.sections.strippers.iter().zip(&mut decoded.strippers) {
            cascade.feeds[0].molar_enthalpy_j_mol = properties[stripper.draw].liquid_enthalpy_j_mol;
            stripper_properties.push(mesh::properties(
                self.ctx,
                cascade,
                stripper.ends(),
                self.temperature_range,
            )?);
        }
        let external = self.feeds.len();
        for (k, pump) in self.sections.pump_arounds.iter().enumerate() {
            let h = pump.return_enthalpy(self.ctx, &decoded.cascade, &properties)?;
            decoded.cascade.feeds[external + k].molar_enthalpy_j_mol = h;
        }
        let returns = external + self.sections.pump_arounds.len();
        for (k, props) in stripper_properties.iter().enumerate() {
            decoded.cascade.feeds[returns + k].molar_enthalpy_j_mol = props[0].vapor_enthalpy_j_mol;
        }
        Ok(Evaluated {
            decoded,
            properties,
            stripper_properties,
        })
    }

    /// Component `i` mole fraction of the total distillate.
    fn distillate_fraction(&self, decoded: &Decoded, i: usize) -> f64 {
        let top = &decoded.cascade.stages[0];
        if self.partial() {
            (top.vapor_mol_s * top.y[i] + top.liquid_draw_mol_s * top.x[i]) / decoded.distillate
        } else {
            top.x[i]
        }
    }

    fn residual(&self, values: &[f64]) -> Result<Vec<f64>, String> {
        let Evaluated {
            decoded,
            properties,
            stripper_properties,
        } = self.evaluate(values)?;
        let cascade = &decoded.cascade;
        let ends = self.ends();
        let mut out = mesh::residuals(
            cascade,
            &properties,
            ends,
            &self.active_components,
            self.flow,
        )?;
        for ((stripper, own), props) in self
            .sections
            .strippers
            .iter()
            .zip(&decoded.strippers)
            .zip(&stripper_properties)
        {
            out.extend(stripper.residuals(own, props, &self.active_components, self.flow)?);
        }
        if self.equations.is_empty() {
            return Ok(out);
        }
        let balances = cascade.residuals(&properties)?;
        let n = cascade.len();
        let (f, b) = (self.flow, decoded.bottoms);
        for equation in &self.equations {
            out.push(match *equation {
                // The stage duty that closes the condenser energy balance is
                // minus its residual without duty; the condenser removes it.
                SpecEquation::CondenserDuty(q) => (balances[0].energy_w - q) / f / ENTHALPY_SCALE,
                SpecEquation::ReboilerDuty(q) => {
                    (-balances[n - 1].energy_w - q) / f / ENTHALPY_SCALE
                }
                SpecEquation::BoilupRatio(r) => (cascade.stages[n - 1].vapor_mol_s - r * b) / f,
                SpecEquation::DistillatePurity(i, v) => self.distillate_fraction(&decoded, i) - v,
                SpecEquation::DistillateRecovery(i, v) => {
                    (decoded.distillate * self.distillate_fraction(&decoded, i)
                        - v * self.feed_components[i])
                        / f
                }
                SpecEquation::BottomsPurity(i, v) => cascade.stages[n - 1].x[i] - v,
                SpecEquation::BottomsRecovery(i, v) => {
                    (b * cascade.stages[n - 1].x[i] - v * self.feed_components[i]) / f
                }
            });
        }
        if out.iter().any(|value| !value.is_finite()) {
            return Err("column has a nonfinite specification residual".into());
        }
        Ok(out)
    }

    /// Specification residuals of a column solved at fixed flows.
    fn spec_residuals(&self, profile: &[f64], flows: &[f64]) -> Result<Vec<f64>, String> {
        let values: Vec<f64> = profile.iter().chain(flows).copied().collect();
        let r = self.residual(&values)?;
        Ok(r[r.len() - self.equations.len()..].to_vec())
    }

    /// Solves the profile at the flows of `flows`, starting from `profile`.
    fn solve_profile(&self, profile: &[f64], flows: &[f64]) -> Result<Vec<f64>, String> {
        let (distillate, reflux) = self.flows(flows);
        let fixed = self.with_fixed_flows(distillate, reflux);
        fixed.decode(profile)?;
        Ok(mesh::solve(profile.to_vec(), |v| fixed.residual(v))?.0)
    }

    /// Outer Newton iteration on the free flows: every evaluation solves the
    /// profile at fixed flows, so the specification residuals become a
    /// function of the one or two free flows alone. It supplies the start of
    /// the simultaneous Newton solve.
    fn converge_flows(
        &self,
        mut profile: Vec<f64>,
        mut flows: Vec<f64>,
    ) -> Result<(Vec<f64>, Vec<f64>), String> {
        let norm = |r: &[f64]| r.iter().map(|v| v.abs()).fold(0.0, f64::max);
        profile = self.solve_profile(&profile, &flows)?;
        let mut g = self.spec_residuals(&profile, &flows)?;
        for _ in 0..OUTER_ITERATIONS {
            let error = norm(&g);
            if error <= OUTER_TOLERANCE {
                return Ok((profile, flows));
            }
            let k = flows.len();
            let mut jacobian = vec![vec![0.0; k]; k];
            for col in 0..k {
                let mut shifted = flows.clone();
                let h = 1e-5 * flows[col].abs().max(1.0);
                shifted[col] += h;
                let solved = self.solve_profile(&profile, &shifted)?;
                let gs = self.spec_residuals(&solved, &shifted)?;
                for (row, (new, old)) in jacobian.iter_mut().zip(gs.iter().zip(&g)) {
                    row[col] = (new - old) / h;
                }
            }
            let step = solve_dense(jacobian, g.iter().map(|v| -v).collect())
                .ok_or("the specifications do not change with the free flows")?;
            let mut accepted = None;
            for backtrack in 0..20 {
                // Limit each step to a factor e in the free flows.
                let scale = 0.5_f64.powi(backtrack)
                    * (1.0 / step.iter().map(|s| s.abs()).fold(0.0, f64::max)).min(1.0);
                let trial: Vec<f64> = flows
                    .iter()
                    .zip(&step)
                    .map(|(f, s)| f + scale * s)
                    .collect();
                let Ok(solved) = self.solve_profile(&profile, &trial) else {
                    continue;
                };
                if let Ok(gt) = self.spec_residuals(&solved, &trial)
                    && norm(&gt) < error
                {
                    accepted = Some((solved, trial, gt));
                    break;
                }
            }
            let Some((solved, trial, gt)) = accepted else {
                return Err(format!(
                    "no column profile approaches the specifications (scaled residual {error:.3e})"
                ));
            };
            profile = solved;
            flows = trial;
            g = gt;
        }
        Err(format!(
            "specification iteration did not converge in {OUTER_ITERATIONS} iterations \
             (scaled residual {:.3e})",
            norm(&g)
        ))
    }
}

fn component(
    ctx: &CalcContext<'_>,
    active: &[usize],
    target: &ComponentTarget,
    field: &str,
) -> Result<(usize, f64), String> {
    let id = target.component.as_deref().unwrap_or_default();
    let index = ctx
        .component_ids()
        .iter()
        .position(|c| *c == id)
        .ok_or_else(|| format!("column {field} component '{id}' is not a package component"))?;
    if !active.contains(&index) {
        return Err(format!(
            "column {field} component '{id}' is not in the column feed"
        ));
    }
    Ok((index, target.value.expect("checked component target")))
}

pub(super) fn calculate(
    spec: &DistillationColumnSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    calculate_with(spec, inlets, ctx, ColumnOptions::default())
}

/// [`calculate`] with a choice of bottom end and initialization. With an
/// adiabatic bottom the energy result holds the condenser port only and the
/// unit duty is the condenser duty.
pub(super) fn calculate_with(
    spec: &DistillationColumnSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
    options: ColumnOptions,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    calculate_sections(spec, inlets, ctx, options, SideSections::default())
}

/// [`calculate_with`] with pump-arounds and side strippers solved
/// simultaneously with the column (#190). `inlets` are the main-column
/// feeds; side-stripper steam is part of `sections`. The outlets are those
/// of [`calculate_with`] followed by each side stripper's bottoms, and the
/// energy values those of [`calculate_with`] followed by each pump-around
/// duty (outlet) and each side-stripper reboiler duty (inlet).
pub(super) fn calculate_sections(
    spec: &DistillationColumnSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
    options: ColumnOptions,
    sections: SideSections,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    let reboiled = options.bottom == BottomEnd::Reboiler;
    if reboiled && spec.overspecified() {
        return Err(
            "column distillate rate and the bottom specification both fix the \
             distillate; choose a reboiler duty, boilup ratio, purity or recovery"
                .into(),
        );
    }
    if let Some(issue) = spec.issues_with(options.bottom).first() {
        return Err(format!("column {} {}", issue.field, issue.problem));
    }
    let (n, packed) = spec.stage_layout()?;
    let pressures = spec.stage_pressures_pa(n).expect("checked column pressure");
    let temperature_range = mesh::temperature_range(ctx)?;
    let mut feed_by_stage = vec![0.0; n];
    let mut feed_components = vec![0.0; ctx.component_count()];
    let mut feeds = Vec::with_capacity(inlets.len());
    let mut seed_components = vec![0.0; ctx.component_count()];
    for (i, inlet) in inlets.iter().enumerate() {
        let stage = usize::from(spec.feed_stages.get(i).copied().unwrap_or(2));
        let feed = StageFeed::from_stream(stage, inlet, ctx.component_count())?;
        if feed.molar_flow_mol_s > 0.0
            && inlet
                .state
                .as_ref()
                .is_some_and(|state| state.pressure_pa < pressures[stage - 1])
        {
            return Err(format!(
                "column pressure exceeds feed {} pressure; pressurize the feed upstream",
                i + 1
            ));
        }
        feed_by_stage[stage - 1] += feed.molar_flow_mol_s;
        for (total, z) in feed_components.iter_mut().zip(&feed.z) {
            *total += feed.molar_flow_mol_s * z;
        }
        if options.seed_inlets.is_none_or(|count| i < count) {
            for (total, z) in seed_components.iter_mut().zip(&feed.z) {
                *total += feed.molar_flow_mol_s * z;
            }
        }
        feeds.push(feed);
    }
    let flow: f64 = feed_by_stage.iter().sum();
    if !flow.is_finite() || flow <= 0.0 {
        return Err("column requires a positive, finite total feed flow".into());
    }
    let feed_composition: Vec<f64> = feed_components.iter().map(|z| z / flow).collect();
    // Side-stripper steam can bring components the main feeds lack.
    let mut steam_components = feed_components.clone();
    let mut steam_flow = flow;
    for stripper in &sections.strippers {
        if let Stripping::Steam(steam) = &stripper.bottom {
            steam_flow += steam.molar_flow_mol_s;
            for (total, z) in steam_components.iter_mut().zip(&steam.z) {
                *total += steam.molar_flow_mol_s * z;
            }
        }
    }
    let all_composition: Vec<f64> = if steam_flow > flow {
        steam_components.iter().map(|z| z / steam_flow).collect()
    } else {
        feed_composition.clone()
    };
    let active_components = mesh::active_components(&all_composition)?;
    let liquid = &spec.liquid_draws[..usize::from(spec.side_draws)];
    let vapor = &spec.vapor_draws[..usize::from(spec.vapor_side_draws)];
    let mut liquid_draws = vec![0.0; n];
    for draw in liquid {
        liquid_draws[usize::from(draw.stage) - 1] +=
            draw.molar_flow_mol_s.expect("checked side-draw flow");
    }
    let mut vapor_draws = vec![0.0; n];
    for draw in vapor {
        vapor_draws[usize::from(draw.stage) - 1] +=
            draw.molar_flow_mol_s.expect("checked vapor-draw flow");
    }
    let mut available = flow - liquid_draws.iter().sum::<f64>() - vapor_draws.iter().sum::<f64>();
    // Each side stripper withdraws its bottoms net of its steam.
    for stripper in &sections.strippers {
        available -= stripper.bottoms - stripper.steam_flow();
    }
    if !available.is_finite() || available <= 0.0 {
        return Err("column side draws must leave positive bottoms flow and distillate".into());
    }
    let liquid_distillate = match spec.condenser {
        DistillationCondenser::PartialVaporLiquid {
            liquid_distillate_mol_s,
        } => liquid_distillate_mol_s.expect("checked liquid distillate"),
        _ => 0.0,
    };

    // Resolve the two specifications into fixed flows and equations.
    let mut equations = Vec::new();
    let mut distillate = match &spec.bottom {
        // The zero duty of an adiabatic bottom stage replaces the bottom
        // specification; its residual is the reboiler-duty residual at zero.
        _ if !reboiled => {
            equations.push(SpecEquation::ReboilerDuty(0.0));
            None
        }
        BottomSpecification::DistillateFraction(f) => Some(flow * f.expect("checked")),
        BottomSpecification::BottomsRate(b) => Some(available - b.expect("checked")),
        BottomSpecification::ReboilerDuty(q) => {
            equations.push(SpecEquation::ReboilerDuty(q.expect("checked")));
            None
        }
        BottomSpecification::BoilupRatio(r) => {
            equations.push(SpecEquation::BoilupRatio(r.expect("checked")));
            None
        }
        BottomSpecification::BottomsPurity(t) => {
            let (i, v) = component(ctx, &active_components, t, "bottoms mole fraction")?;
            equations.push(SpecEquation::BottomsPurity(i, v));
            None
        }
        BottomSpecification::BottomsRecovery(t) => {
            let (i, v) = component(ctx, &active_components, t, "bottoms recovery")?;
            equations.push(SpecEquation::BottomsRecovery(i, v));
            None
        }
    };
    let reflux = match &spec.top {
        TopSpecification::RefluxRatio(r) => RefluxFlow::Ratio(r.expect("checked")),
        TopSpecification::RefluxRate(l) => RefluxFlow::Fixed(l.expect("checked")),
        TopSpecification::DistillateRate(d) => {
            distillate = Some(d.expect("checked"));
            RefluxFlow::Free
        }
        TopSpecification::CondenserDuty(q) => {
            equations.insert(0, SpecEquation::CondenserDuty(q.expect("checked")));
            RefluxFlow::Free
        }
        TopSpecification::DistillatePurity(t) => {
            let (i, v) = component(ctx, &active_components, t, "distillate mole fraction")?;
            equations.insert(0, SpecEquation::DistillatePurity(i, v));
            RefluxFlow::Free
        }
        TopSpecification::DistillateRecovery(t) => {
            let (i, v) = component(ctx, &active_components, t, "distillate recovery")?;
            equations.insert(0, SpecEquation::DistillateRecovery(i, v));
            RefluxFlow::Free
        }
    };
    if let Some(d) = distillate {
        if !d.is_finite() || d <= 0.0 || d >= available {
            return Err(
                "column withdrawals must leave positive bottoms flow and distillate".into(),
            );
        }
        if d <= liquid_distillate {
            return Err(
                "column liquid distillate must be smaller than the total distillate".into(),
            );
        }
        if let RefluxFlow::Ratio(r) = reflux
            && !(d * r).is_finite()
        {
            return Err("column reflux must be finite".into());
        }
    }

    let seed_pressure = pressures[1];
    let seed_total: f64 = seed_components.iter().sum();
    let seed_composition: Vec<f64> = if options.seed_inlets.is_some() && seed_total > 0.0 {
        seed_components.iter().map(|z| z / seed_total).collect()
    } else {
        feed_composition.clone()
    };
    let seed = ctx
        .saturation_temperature(&seed_composition, seed_pressure, SaturationKind::Bubble)
        .map_err(|e| format!("column feed bubble-point initialization: {e}"))?;
    // The condenser liquid starts at the composition of the vapor leaving
    // the first tray. Its own bubble point can be much colder than the tray.
    // When that liquid has no single-liquid bubble point (it splits into two
    // liquids), the condenser starts at the feed bubble point instead; the
    // converged column is then rejected by the stage liquid–liquid check.
    let condenser_seed = match ctx.saturation_temperature(
        &seed.vapor_composition,
        pressures[0],
        SaturationKind::Bubble,
    ) {
        Ok(state) => state,
        Err(e) if is_liquid_split(&e) => seed.clone(),
        Err(e) => return Err(format!("column condenser bubble-point initialization: {e}")),
    };
    let column = Column {
        ctx,
        pressures: pressures.clone(),
        temperature_range,
        feeds,
        feed_by_stage,
        feed_components,
        liquid_draws,
        vapor_draws,
        active_components,
        flow,
        available,
        condenser: spec.condenser.stage_condenser(),
        liquid_distillate,
        distillate: distillate.map_or(DistillateFlow::Free, DistillateFlow::Fixed),
        reflux,
        equations,
        sections,
    };
    // Starting flows: a free distillate starts at half the available
    // product (or above the liquid distillate), a free reflux at L/D = 2.
    let d0 = distillate.unwrap_or_else(|| {
        let half = 0.5 * available;
        if half > liquid_distillate {
            half
        } else {
            0.5 * (liquid_distillate + available)
        }
    });
    let l0 = match reflux {
        RefluxFlow::Fixed(l) => l,
        RefluxFlow::Ratio(r) => r * d0,
        RefluxFlow::Free => 2.0 * d0,
    };
    // A component absent from the seed inlets (stripping steam) starts at
    // its overall feed fraction in the vapor and at a tenth of it in the
    // liquid; with every inlet in the seed this changes nothing. (Started
    // at the full fraction in the liquid, water/pseudo-component columns
    // under PR either stalled or converged to a condenser liquid that splits.)
    let fill = |phase: &[f64], scale: f64| {
        let phase: Vec<f64> = phase
            .iter()
            .zip(&all_composition)
            .map(|(&p, &z)| if p > 0.0 { p } else { z * scale })
            .collect();
        mesh::log_ratios(&phase, &column.active_components)
    };
    let ratios = |phase: &[f64]| fill(phase, 1.0);
    let x = fill(&seed.liquid_composition, ABSENT_LIQUID_SEED);
    let y = ratios(&seed.vapor_composition);
    let mut initial = vec![condenser_seed.temperature_k / TEMPERATURE_SCALE];
    if column.partial() {
        initial.extend(ratios(&seed.vapor_composition));
    }
    initial.extend(ratios(&condenser_seed.vapor_composition));
    for j in 1..n {
        initial.push(seed.temperature_k / TEMPERATURE_SCALE);
        initial.extend(&x);
        initial.extend(&y);
        if j > 1 {
            initial.push(((l0 + d0) / flow).ln());
        }
    }
    let (mut d_start, mut l_start) = (d0, l0);
    let mut warm = false;
    if !column.sections.strippers.is_empty() {
        // Side strippers start from the converged precursor column: each
        // stripper replaced by its bottoms drawn at the draw stage and its
        // steam fed to the return stage (the same external balance). Its
        // stages start at the draw-stage state. The coupled solve below is
        // the simultaneous Newton solve; a precursor that does not converge
        // leaves the plain start.
        let precursor = column.precursor();
        if let Ok((solved, _)) = column_solve(&precursor, initial.clone(), d0, l0)
            && let Ok(decoded) = precursor.decode(&solved)
        {
            initial = solved[..column.main_profile_len()].to_vec();
            for stripper in &column.sections.strippers {
                let stage = &decoded.cascade.stages[stripper.draw];
                initial.extend(stripper.initial(
                    stage.temperature_k,
                    &mesh::log_ratios(&stage.x, &column.active_components),
                    &mesh::log_ratios(&stage.y, &column.active_components),
                    flow,
                ));
            }
            (d_start, l_start) = (decoded.distillate, decoded.reflux);
            warm = true;
        } else {
            for stripper in &column.sections.strippers {
                initial.extend(stripper.initial(seed.temperature_k, &x, &y, flow));
            }
        }
    }
    // From the precursor the simultaneous Newton solve starts directly
    // (free flows included); the outer flow iteration is the fallback.
    let warm_solve = warm
        .then(|| {
            let start: Vec<f64> = initial
                .iter()
                .copied()
                .chain(column.flow_variables(d_start, l_start))
                .collect();
            mesh::solve(start, |values| column.residual(values)).ok()
        })
        .flatten();
    let (solved, iterations) = match warm_solve {
        Some(solved) => solved,
        None => column_solve(&column, initial, d_start, l_start)?,
    };
    let Evaluated {
        decoded,
        properties,
        stripper_properties,
    } = column.evaluate(&solved)?;
    let ends = column.ends();
    let mut cascade = decoded.cascade;
    let mut warnings = mesh::check_liquid_liquid(ctx, &cascade)?;
    mesh::verify_phases(ctx, &cascade)?;
    for (k, stripper) in decoded.strippers.iter().enumerate() {
        let context = |e: String| format!("side stripper {}: {e}", k + 1);
        for warning in mesh::check_liquid_liquid(ctx, stripper).map_err(context)? {
            if !warnings.contains(&warning) {
                warnings.push(warning);
            }
        }
        mesh::verify_phases(ctx, stripper).map_err(context)?;
    }
    let residuals = cascade.residuals(&properties)?;
    cascade.stages[0].duty_w = -residuals[0].energy_w;
    // An adiabatic bottom stage has no duty; its energy residual is within
    // the solve tolerance and stays in the unit's energy balance.
    cascade.stages[n - 1].duty_w = if reboiled {
        -residuals[n - 1].energy_w
    } else {
        0.0
    };
    if let Some(section) = &packed {
        section.check(&cascade, ends)?;
    }
    let mut energy = ends.port_duties(&cascade)?;
    let condenser = -cascade.stages[0].duty_w;
    let reboiler = cascade.stages[n - 1].duty_w;
    if !reboiled {
        energy.truncate(1);
        if condenser <= 0.0 {
            return Err("column requires heat removal at the condenser".into());
        }
    } else if condenser <= 0.0 || reboiler <= 0.0 {
        return Err(
            "column requires heat removal at the condenser and heat input at the reboiler".into(),
        );
    }
    let mut duty = reboiler - condenser;
    // Pump-arounds: a cooled liquid returned at or above its draw.
    let external = column.feeds.len();
    let mut pump_around_duties = Vec::with_capacity(column.sections.pump_arounds.len());
    for (k, pump) in column.sections.pump_arounds.iter().enumerate() {
        let returned = &cascade.feeds[external + k];
        let draw = &cascade.stages[pump.draw];
        let removed = pump.flow
            * (properties[pump.draw].liquid_enthalpy_j_mol - returned.molar_enthalpy_j_mol);
        if let Cooling::ReturnTemperature(t) = pump.cooling {
            if t >= draw.temperature_k {
                return Err(format!(
                    "pump-around {} return temperature {t} K is not below its draw-stage temperature {:.3} K",
                    k + 1,
                    draw.temperature_k
                ));
            }
            let pressure = cascade.stages[pump.ret].pressure_pa;
            let state = ctx
                .tp(&returned.z, t, pressure)
                .map_err(|e| format!("pump-around {} return state: {e}", k + 1))?;
            if !matches!(state.phases.as_slice(), [p] if p.phase == Phase::Liquid) {
                return Err(format!(
                    "pump-around {} return at {t} K and {pressure} Pa is not a single liquid",
                    k + 1
                ));
            }
        }
        if !(removed.is_finite() && removed > 0.0) {
            return Err(format!("pump-around {} must remove heat", k + 1));
        }
        duty -= removed;
        pump_around_duties.push(removed);
        energy.push(Some(removed));
    }
    // Side strippers: bottoms products and reboiler duties.
    let returns = external + column.sections.pump_arounds.len();
    let mut profiles = Vec::with_capacity(decoded.strippers.len());
    let mut stripper_products = Vec::with_capacity(decoded.strippers.len());
    for (k, ((stripper, own), props)) in column
        .sections
        .strippers
        .iter()
        .zip(&decoded.strippers)
        .zip(&stripper_properties)
        .enumerate()
    {
        let last = own.stages.len() - 1;
        let reboiler_duty = if stripper.reboiled() {
            let q = -own.residuals(props)?[last].energy_w;
            if !(q.is_finite() && q > 0.0) {
                return Err(format!(
                    "side stripper {} requires heat input at its reboiler",
                    k + 1
                ));
            }
            duty += q;
            energy.push(Some(q));
            q
        } else {
            0.0
        };
        stripper_products.push(mesh::product(
            ctx,
            &own.stages[last],
            stripper.bottoms,
            Phase::Liquid,
        )?);
        profiles.push(SideStripperProfile {
            stage_temperatures_k: own.stages.iter().map(|s| s.temperature_k).collect(),
            draw_mol_s: own.feeds[0].molar_flow_mol_s,
            vapor_return_mol_s: cascade.feeds[returns + k].molar_flow_mol_s,
            reboiler_duty_w: reboiler_duty,
        });
    }
    let top = &cascade.stages[0];
    let mut outlets = vec![
        if column.partial() {
            mesh::product(ctx, top, top.vapor_mol_s, Phase::Vapor)?
        } else {
            mesh::product(ctx, top, decoded.distillate, Phase::Liquid)?
        },
        mesh::product(ctx, &cascade.stages[n - 1], decoded.bottoms, Phase::Liquid)?,
    ];
    for draw in liquid {
        outlets.push(mesh::product(
            ctx,
            &cascade.stages[usize::from(draw.stage) - 1],
            draw.molar_flow_mol_s.expect("checked side-draw flow"),
            Phase::Liquid,
        )?);
    }
    for draw in vapor {
        outlets.push(mesh::product(
            ctx,
            &cascade.stages[usize::from(draw.stage) - 1],
            draw.molar_flow_mol_s.expect("checked vapor-draw flow"),
            Phase::Vapor,
        )?);
    }
    if liquid_distillate > 0.0 {
        outlets.push(mesh::product(ctx, top, liquid_distillate, Phase::Liquid)?);
    }
    outlets.extend(stripper_products);
    let mut output = UnitOutput::new(outlets, duty, warnings);
    output.energy = energy;
    Ok((
        output,
        DistillationColumnResult {
            stage_temperatures_k: cascade.stages.iter().map(|s| s.temperature_k).collect(),
            stage_pressures_pa: pressures,
            condenser_duty_w: condenser,
            reboiler_duty_w: reboiler,
            distillate_mol_s: decoded.distillate,
            reflux_mol_s: decoded.reflux,
            bottoms_mol_s: decoded.bottoms,
            boilup_mol_s: cascade.stages[n - 1].vapor_mol_s,
            packed_height_m: packed.map(|section| section.packed.packed_height_m),
            iterations,
            product_curves: Vec::new(),
            pump_around_duties_w: pump_around_duties,
            side_strippers: profiles,
        },
    ))
}

/// Solves `column` from the profile `initial` and the starting distillate
/// and reflux: directly, or with free flows by the outer flow iteration
/// followed by the simultaneous Newton solve.
fn column_solve(
    column: &Column<'_, '_>,
    initial: Vec<f64>,
    distillate: f64,
    reflux: f64,
) -> Result<(Vec<f64>, usize), String> {
    if column.equations.is_empty() {
        return mesh::solve(initial, |values| column.residual(values));
    }
    let flows = column.flow_variables(distillate, reflux);
    let (profile, flows) = column
        .converge_flows(initial, flows)
        .map_err(|e| format!("column specifications cannot be met: {e}"))?;
    let start: Vec<f64> = profile.into_iter().chain(flows).collect();
    mesh::solve(start, |values| column.residual(values))
        .map_err(|e| format!("column specifications cannot be met: {e}"))
}

/// Saturation failures caused by a liquid that splits into two liquids: a
/// three-phase bubble point (cubic) or a saturation liquid that fails the
/// NRTL tangent-plane trials.
fn is_liquid_split(message: &str) -> bool {
    message.contains("splits into two liquids") || message.contains("liquid-liquid split")
}
