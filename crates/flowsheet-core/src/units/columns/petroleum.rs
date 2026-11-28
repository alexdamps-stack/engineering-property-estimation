//! Petroleum fractionator (#190): a single main column for feeds that carry
//! characterised petroleum pseudo-components (#229), with a total or partial
//! condenser, a reboiler or stripping steam, and liquid side products
//! specified by molar or standard liquid volume flow.
//!
//! The column is composed on the rigorous MESH solve of the distillation
//! column (#186): the fractionator specification maps to a
//! [`DistillationColumnSpec`] with the same stage numbering, feeds, liquid
//! side draws, pressure profile, condenser and top specification. Stripping
//! steam is one more feed to the bottom stage, which then has no reboiler
//! and is adiabatic ([`BottomEnd::Adiabatic`]). A side product specified by
//! standard liquid volume is converted to its molar flow by successive
//! substitution on the product's own composition, each pass a converged
//! column. Pump-arounds and side strippers ([`super::side_sections`]) are
//! solved simultaneously with the main column in one Newton system; a
//! stripped side product is the stripper's bottoms. Product TBP/D86 curves
//! and cut-point side products are in [`super::petroleum_curves`]. General
//! coupled multi-column systems are not modelled
//! (`docs/a17-columns-petroleum.md`).

use super::distillation::{self, BottomEnd, ColumnOptions};
use super::petroleum_curves::{self, CutPoint};
use super::side_sections::{Cooling, PumpAroundLoop, SideSections, StripperSection, Stripping};
use super::stages::StageFeed;
use super::{
    BottomSpecification, ColumnInternals, ColumnPressureDrop, DistillationColumnResult,
    DistillationColumnSpec, DistillationCondenser, LiquidSideDraw, MAX_COLUMN_PORTS,
    MAX_COLUMN_STAGES, TopSpecification, count_issue,
};
use crate::MaterialStream;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};
use crate::units::{CalcContext, UnitOutput};

/// Density of water at 60 °F (288.7056 K) and 101.325 kPa, kg/m³: IAPWS-95
/// as evaluated by CoolProp 7.2.0. A pseudo-component's standard liquid
/// volume is `M / (SG · ρ_w)` with its specific gravity SG at 60 °F/60 °F.
pub const WATER_DENSITY_60F_KG_M3: f64 = 999.017_082_407_830_6;

/// Relative change of a volume-specified draw's molar flow at which the
/// successive substitution is converged.
const VOLUME_TOLERANCE: f64 = 1e-11;
const VOLUME_ITERATIONS: usize = 30;

/// Flow specification of a liquid side product.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DrawRate {
    /// Molar withdrawal, mol/s.
    Molar(Option<f64>),
    /// Standard liquid volume flow at 60 °F, m³/s, from the ideal mixing of
    /// the pseudo-components' standard liquid volumes. Every component of
    /// the column feed must then be a pseudo-component.
    StandardLiquidVolume(Option<f64>),
    /// TBP or ASTM D86 temperature at a liquid volume percent; the molar
    /// flow is found by an outer iteration (`petroleum_curves`).
    CutPoint(CutPoint),
}

impl Default for DrawRate {
    fn default() -> Self {
        Self::Molar(None)
    }
}

/// One liquid side product (port `side draw k`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SideProduct {
    /// One-based stage number, the condenser being stage 1. Only internal
    /// stages (not the condenser or the bottom stage) may carry a draw.
    pub stage: u16,
    pub rate: DrawRate,
}

impl Default for SideProduct {
    fn default() -> Self {
        Self {
            stage: 2,
            rate: DrawRate::default(),
        }
    }
}

/// Cooling of a pump-around.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PumpAroundCooling {
    /// Temperature of the returned liquid, K; below the draw-stage
    /// temperature, and the return must be a single liquid at the
    /// return-stage pressure.
    ReturnTemperature(Option<f64>),
    /// Heat removed, W; positive.
    Duty(Option<f64>),
}

impl Default for PumpAroundCooling {
    fn default() -> Self {
        Self::ReturnTemperature(None)
    }
}

/// A pump-around (port `pump-around duty k`): liquid drawn from an internal
/// stage, cooled and returned to the same or a higher stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PumpAround {
    /// One-based internal stage the liquid is drawn from.
    pub draw_stage: u16,
    /// One-based stage the cooled liquid returns to, `2..=draw_stage`.
    pub return_stage: u16,
    /// Circulated molar flow, mol/s.
    pub molar_flow_mol_s: Option<f64>,
    pub cooling: PumpAroundCooling,
}

impl Default for PumpAround {
    fn default() -> Self {
        Self {
            draw_stage: 3,
            return_stage: 2,
            molar_flow_mol_s: None,
            cooling: PumpAroundCooling::default(),
        }
    }
}

/// Stripping agent of a side stripper.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum SideStripperBottom {
    /// Steam (port `stripper steam k`) to the adiabatic bottom stage.
    #[default]
    Steam,
    /// Partial reboiler (port `stripper reboiler duty k`) with its duty, W.
    ReboilerDuty(Option<f64>),
    /// Partial reboiler (port `stripper reboiler duty k`) with the boilup
    /// ratio, vapor from the reboiler over the product.
    BoilupRatio(Option<f64>),
}

/// A side stripper on one side product: the product's liquid draw feeds its
/// top stage, its top vapor returns to the main column and its bottoms
/// leave on the product's `side draw` port at the product's flow
/// specification. It operates at the pressure of the draw stage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SideStripper {
    /// One-based side product (`side draw k`) the stripper is attached to.
    pub product: u16,
    /// Number of equilibrium stages, including a reboiler.
    pub stages: u16,
    /// One-based main-column stage receiving the top vapor, `2..=` the
    /// product's draw stage.
    pub return_stage: u16,
    pub bottom: SideStripperBottom,
}

impl Default for SideStripper {
    fn default() -> Self {
        Self {
            product: 1,
            stages: 3,
            return_stage: 2,
            bottom: SideStripperBottom::default(),
        }
    }
}

/// Bottom end of the fractionator.
#[derive(Debug, Clone, PartialEq)]
pub enum FractionatorBottom {
    /// Partial reboiler (energy port `reboiler duty`) with its column
    /// specification.
    Reboiler(BottomSpecification),
    /// Stripping steam (port `stripping steam`) enters the adiabatic bottom
    /// stage; the bottom has no reboiler and no second specification.
    StrippingSteam,
}

impl Default for FractionatorBottom {
    fn default() -> Self {
        Self::Reboiler(BottomSpecification::default())
    }
}

/// Petroleum fractionator specification (#190).
///
/// Ports beyond the fixed `feed`, `overhead`, `bottoms` and `condenser duty`
/// follow in this order: `side draw 1..k`, `feed 2..n`, `liquid distillate`
/// (partial vapor-and-liquid condenser), `stripping steam` (steam bottom)
/// or the `reboiler duty` energy inlet (reboiled bottom), `pump-around
/// duty 1..p` energy outlets and, per side stripper in list order, `stripper
/// steam k` or the `stripper reboiler duty k` energy inlet. Existing side
/// draw connections of a drawn-only fractionator keep their port indices.
#[derive(Debug, Clone, PartialEq)]
pub struct PetroleumFractionatorSpec {
    /// Number of liquid side products, 0 to [`MAX_COLUMN_PORTS`].
    pub side_draws: u16,
    /// Number of feeds, including the fixed `feed` port; at least 1.
    pub feeds: u16,
    /// One-based feed stage per feed. Unset entries default to stage 2.
    pub feed_stages: Vec<u16>,
    /// Stage and flow of each side product. Unset entries have no flow.
    pub products: Vec<SideProduct>,
    /// Ideal trays or one ideal packed bed between condenser and bottom.
    pub internals: ColumnInternals,
    /// Condenser (top stage) absolute pressure, Pa.
    pub pressure_pa: Option<f64>,
    /// Pressure drop from the condenser to the first tray, Pa.
    pub condenser_pressure_drop_pa: f64,
    /// Pressure drop over the trays and the bottom stage.
    pub pressure_drop: ColumnPressureDrop,
    /// Condenser type; the vapor or liquid distillate is the `overhead`.
    pub condenser: DistillationCondenser,
    /// Top specification (reflux ratio, rates, condenser duty, purity or
    /// recovery).
    pub top: TopSpecification,
    /// Reboiler with its specification, or stripping steam.
    pub bottom: FractionatorBottom,
    /// Pump-arounds, at most [`MAX_COLUMN_PORTS`].
    pub pump_arounds: Vec<PumpAround>,
    /// Side strippers, at most one per side product.
    pub side_strippers: Vec<SideStripper>,
}

impl Default for PetroleumFractionatorSpec {
    fn default() -> Self {
        Self {
            side_draws: 0,
            feeds: 1,
            feed_stages: Vec::new(),
            products: Vec::new(),
            internals: ColumnInternals::Trays(10),
            pressure_pa: None,
            condenser_pressure_drop_pa: 0.0,
            pressure_drop: ColumnPressureDrop::default(),
            condenser: DistillationCondenser::Total,
            top: TopSpecification::default(),
            bottom: FractionatorBottom::default(),
            pump_arounds: Vec::new(),
            side_strippers: Vec::new(),
        }
    }
}

impl PetroleumFractionatorSpec {
    fn bottom_end(&self) -> BottomEnd {
        match self.bottom {
            FractionatorBottom::Reboiler(_) => BottomEnd::Reboiler,
            FractionatorBottom::StrippingSteam => BottomEnd::Adiabatic,
        }
    }

    /// Whether stripping steam enters the bottom stage.
    pub fn steam(&self) -> bool {
        self.bottom == FractionatorBottom::StrippingSteam
    }

    fn product(&self, k: usize) -> SideProduct {
        self.products.get(k).copied().unwrap_or_default()
    }

    /// Index in `side_strippers` of the stripper on side product `k`
    /// (zero-based).
    fn stripper_of(&self, k: usize) -> Option<usize> {
        self.side_strippers
            .iter()
            .position(|s| usize::from(s.product) == k + 1)
    }

    /// Whether the column has pump-arounds or side strippers.
    fn has_sections(&self) -> bool {
        !(self.pump_arounds.is_empty() && self.side_strippers.is_empty())
    }

    /// The equivalent distillation column. `molar` gives the molar flow of
    /// each side product, and the column then draws only the products
    /// without a side stripper (a stripped product is the stripper's
    /// bottoms; its main-column draw is part of the side sections). `None`
    /// keeps every product as a draw at its molar specification (and marks
    /// a volume-specified product as set, its value is checked separately),
    /// for the specification checks.
    fn column(&self, molar: Option<&[f64]>) -> DistillationColumnSpec {
        let draws: Vec<usize> = (0..usize::from(self.side_draws.min(MAX_COLUMN_PORTS)))
            .filter(|&k| molar.is_none() || self.stripper_of(k).is_none())
            .collect();
        let side_draws = if molar.is_none() {
            self.side_draws
        } else {
            u16::try_from(draws.len()).expect("at most MAX_COLUMN_PORTS draws")
        };
        let mut feed_stages: Vec<u16> = (0..usize::from(self.feeds.min(MAX_COLUMN_PORTS)))
            .map(|i| self.feed_stages.get(i).copied().unwrap_or(2))
            .collect();
        let mut feeds = self.feeds;
        if self.steam() {
            let stages = self.internals.resolve().map_or(2, |(count, _)| count + 2);
            feed_stages.push(u16::try_from(stages).unwrap_or(u16::MAX));
            feeds = feeds.saturating_add(1);
        }
        DistillationColumnSpec {
            feeds,
            side_draws,
            feed_stages,
            liquid_draws: draws
                .into_iter()
                .map(|k| {
                    let product = self.product(k);
                    LiquidSideDraw {
                        stage: product.stage,
                        molar_flow_mol_s: match (molar, product.rate) {
                            (Some(flows), _) => Some(flows[k]),
                            (None, DrawRate::Molar(v)) => v,
                            (None, DrawRate::StandardLiquidVolume(_) | DrawRate::CutPoint(_)) => {
                                Some(1.0)
                            }
                        },
                    }
                })
                .collect(),
            vapor_side_draws: 0,
            vapor_draws: Vec::new(),
            internals: self.internals,
            pressure_pa: self.pressure_pa,
            condenser_pressure_drop_pa: self.condenser_pressure_drop_pa,
            pressure_drop: self.pressure_drop,
            condenser: self.condenser,
            top: self.top.clone(),
            bottom: match &self.bottom {
                FractionatorBottom::Reboiler(b) => b.clone(),
                FractionatorBottom::StrippingSteam => BottomSpecification::default(),
            },
        }
    }

    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        count_issue(&mut issues, "feeds", self.feeds, 1);
        count_issue(&mut issues, "side draws", self.side_draws, 0);
        if self.products.len() > usize::from(self.side_draws) {
            issues.push(SpecIssue {
                field: "side-product specifications",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        if self.steam() && self.feeds >= MAX_COLUMN_PORTS {
            issues.push(SpecIssue {
                field: "feeds",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        if !issues.is_empty() {
            return issues;
        }
        issues.extend(self.column(None).issues_with(self.bottom_end()));
        for k in 0..usize::from(self.side_draws) {
            match self.product(k).rate {
                DrawRate::StandardLiquidVolume(v) => check_optional(
                    &mut issues,
                    "side-product standard liquid volume flow",
                    v,
                    Sign::Positive,
                ),
                DrawRate::CutPoint(cut) => cut.issues(&mut issues),
                DrawRate::Molar(_) => {}
            }
        }
        self.section_issues(&mut issues);
        issues
    }

    /// Issues of the pump-arounds and side strippers.
    fn section_issues(&self, issues: &mut Vec<SpecIssue>) {
        let mut push = |field, problem| issues.push(SpecIssue { field, problem });
        if self.pump_arounds.len() > usize::from(MAX_COLUMN_PORTS) {
            push("pump-arounds", SpecProblem::CountOutOfRange);
        }
        if self.side_strippers.len() > usize::from(self.side_draws) {
            push("side strippers", SpecProblem::CountOutOfRange);
        }
        if !issues.is_empty() {
            return;
        }
        let stages = self.internals.resolve().ok().map(|(count, _)| count + 2);
        for pump in &self.pump_arounds {
            check_optional(
                issues,
                "pump-around molar flow",
                pump.molar_flow_mol_s,
                Sign::Positive,
            );
            match pump.cooling {
                PumpAroundCooling::ReturnTemperature(t) => {
                    check_optional(issues, "pump-around return temperature", t, Sign::Positive)
                }
                PumpAroundCooling::Duty(q) => {
                    check_optional(issues, "pump-around duty", q, Sign::Positive)
                }
            }
            let draw = usize::from(pump.draw_stage);
            if stages.is_some_and(|n| {
                !(2..n).contains(&draw) || !(2..=draw).contains(&usize::from(pump.return_stage))
            }) {
                issues.push(SpecIssue {
                    field: "pump-around stage",
                    problem: SpecProblem::OutOfRange,
                });
            }
        }
        for (i, stripper) in self.side_strippers.iter().enumerate() {
            let product = usize::from(stripper.product);
            let mut push = |field, problem| issues.push(SpecIssue { field, problem });
            if !(1..=usize::from(self.side_draws)).contains(&product) {
                push("side-stripper product", SpecProblem::OutOfRange);
            } else {
                if self.side_strippers[..i]
                    .iter()
                    .any(|s| s.product == stripper.product)
                {
                    push("side-stripper product", SpecProblem::Duplicate);
                }
                let draw = self.product(product - 1).stage;
                if !(2..=draw).contains(&stripper.return_stage) {
                    push("side-stripper return stage", SpecProblem::OutOfRange);
                }
            }
            if !(1..=MAX_COLUMN_STAGES).contains(&stripper.stages) {
                push("side-stripper stages", SpecProblem::OutOfRange);
            }
            match stripper.bottom {
                SideStripperBottom::Steam => {}
                SideStripperBottom::ReboilerDuty(q) => {
                    check_optional(issues, "side-stripper reboiler duty", q, Sign::Positive)
                }
                SideStripperBottom::BoilupRatio(r) => {
                    check_optional(issues, "side-stripper boilup ratio", r, Sign::Positive)
                }
            }
        }
    }

    /// Pump-arounds and side strippers for the calculation; `molar` is the
    /// molar flow of every side product and `steam` the `stripper steam`
    /// inlets in order.
    fn sections(
        &self,
        molar: &[f64],
        steam: &[MaterialStream],
        ctx: &CalcContext<'_>,
    ) -> Result<SideSections, String> {
        let column = self.column(None);
        let (stages, _) = column.stage_layout()?;
        let pressures = column
            .stage_pressures_pa(stages)
            .expect("checked column pressure");
        let pump_arounds = self
            .pump_arounds
            .iter()
            .map(|pump| PumpAroundLoop {
                draw: usize::from(pump.draw_stage) - 1,
                ret: usize::from(pump.return_stage) - 1,
                flow: pump.molar_flow_mol_s.expect("checked pump-around flow"),
                cooling: match pump.cooling {
                    PumpAroundCooling::ReturnTemperature(t) => {
                        Cooling::ReturnTemperature(t.expect("checked return temperature"))
                    }
                    PumpAroundCooling::Duty(q) => Cooling::Duty(q.expect("checked duty")),
                },
            })
            .collect();
        let mut steam = steam.iter();
        let mut strippers = Vec::with_capacity(self.side_strippers.len());
        for (i, stripper) in self.side_strippers.iter().enumerate() {
            let k = usize::from(stripper.product) - 1;
            let draw = usize::from(self.product(k).stage) - 1;
            let count = usize::from(stripper.stages);
            let bottom = match stripper.bottom {
                SideStripperBottom::Steam => {
                    let stream = steam
                        .next()
                        .expect("one inlet per steam-stripped side stripper");
                    let feed = StageFeed::from_stream(count, stream, ctx.component_count())?;
                    if !(feed.molar_flow_mol_s.is_finite() && feed.molar_flow_mol_s > 0.0) {
                        return Err(format!(
                            "petroleum fractionator side stripper {} needs a steam flow on stripper steam {}",
                            i + 1,
                            i + 1
                        ));
                    }
                    if stream
                        .state
                        .as_ref()
                        .is_some_and(|s| s.pressure_pa < pressures[draw])
                    {
                        return Err(format!(
                            "petroleum fractionator side stripper {} steam pressure is below the stripper (draw-stage) pressure",
                            i + 1
                        ));
                    }
                    Stripping::Steam(feed)
                }
                SideStripperBottom::ReboilerDuty(q) => {
                    Stripping::ReboilerDuty(q.expect("checked reboiler duty"))
                }
                SideStripperBottom::BoilupRatio(r) => {
                    Stripping::BoilupRatio(r.expect("checked boilup ratio"))
                }
            };
            strippers.push(StripperSection {
                draw,
                ret: usize::from(stripper.return_stage) - 1,
                stages: count,
                bottoms: molar[k],
                bottom,
            });
        }
        Ok(SideSections {
            pump_arounds,
            strippers,
        })
    }
}

/// Solves the column at the molar flow of every side product, with its
/// pump-arounds and side strippers, and returns the outlets in port order
/// (`overhead`, `bottoms`, `side draw 1..k`, `liquid distillate`).
fn solve(
    spec: &PetroleumFractionatorSpec,
    molar: &[f64],
    main: &[MaterialStream],
    steam: &[MaterialStream],
    ctx: &CalcContext<'_>,
    options: ColumnOptions,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    let column = spec.column(Some(molar));
    if !spec.has_sections() {
        return distillation::calculate_with(&column, main, ctx, options);
    }
    let sections = spec.sections(molar, steam, ctx)?;
    let stripped = sections.strippers.len();
    let (mut output, result) =
        distillation::calculate_sections(&column, main, ctx, options, sections)?;
    // Column outlets: overhead, bottoms, unstripped draws, liquid
    // distillate, then the side-stripper bottoms.
    let mut rest = output.outlets.split_off(2);
    let bottoms = rest.split_off(rest.len() - stripped);
    let liquid_distillate = rest.split_off(usize::from(column.side_draws));
    let mut unstripped = rest.into_iter();
    for k in 0..usize::from(spec.side_draws) {
        output.outlets.push(match spec.stripper_of(k) {
            Some(i) => bottoms[i].clone(),
            None => unstripped.next().expect("one outlet per unstripped draw"),
        });
    }
    output.outlets.extend(liquid_distillate);
    Ok((output, result))
}

/// Standard liquid molar volume of every package component at 60 °F, m³/mol;
/// `None` for a component without a specific gravity (not a
/// pseudo-component).
pub(super) fn standard_volumes(ctx: &CalcContext<'_>) -> Result<Vec<Option<f64>>, String> {
    let components = ctx.package().formation_components().map_err(|e| {
        format!(
            "petroleum fractionator needs a cubic package with characterised pseudo-components: {e}"
        )
    })?;
    Ok(components
        .iter()
        .map(|c| {
            c.pseudo_characterisation().map(|p| {
                p.molar_mass_kg_per_mol / (p.input.specific_gravity * WATER_DENSITY_60F_KG_M3)
            })
        })
        .collect())
}

/// Standard liquid volume per mole of a mixture, m³/mol.
fn mixture_volume(volumes: &[Option<f64>], x: &[f64]) -> Option<f64> {
    x.iter()
        .zip(volumes)
        .filter(|(x, _)| **x > 0.0)
        .map(|(x, v)| v.map(|v| x * v))
        .sum()
}

/// Calculates the fractionator and its product curves. Inlets are `feed`,
/// `feed 2..n`, with a steam bottom `stripping steam`, and the
/// `stripper steam k` inlets, in port order.
pub(super) fn calculate(
    spec: &PetroleumFractionatorSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!(
            "petroleum fractionator {} {}",
            issue.field, issue.problem
        ));
    }
    let cut_points = (0..usize::from(spec.side_draws))
        .any(|k| matches!(spec.product(k).rate, DrawRate::CutPoint(_)));
    let (mut output, mut result) = if cut_points {
        petroleum_curves::solve(spec, inlets, ctx, calculate_rates)?
    } else {
        calculate_rates(spec, inlets, ctx)?
    };
    result.product_curves = petroleum_curves::product_curves(spec, &mut output, ctx)?;
    Ok((output, result))
}

/// Calculates the fractionator at fixed molar or standard-volume products.
fn calculate_rates(
    spec: &PetroleumFractionatorSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!(
            "petroleum fractionator {} {}",
            issue.field, issue.problem
        ));
    }
    let volumes = standard_volumes(ctx)?;
    if volumes.iter().all(Option::is_none) {
        return Err(
            "petroleum fractionator needs characterised pseudo-components; the selected package has none"
                .into(),
        );
    }
    let feeds = usize::from(spec.feeds);
    let mut hydrocarbon = vec![0.0; ctx.component_count()];
    for inlet in &inlets[..feeds] {
        if let Some(state) = &inlet.state {
            for (total, z) in hydrocarbon.iter_mut().zip(&state.mole_fractions) {
                *total += inlet.molar_flow_mol_s * z;
            }
        }
    }
    if !hydrocarbon
        .iter()
        .zip(&volumes)
        .any(|(flow, volume)| *flow > 0.0 && volume.is_some())
    {
        return Err(
            "petroleum fractionator feed carries no characterised pseudo-component; use the distillation column"
                .into(),
        );
    }
    let draws = usize::from(spec.side_draws);
    let targets: Vec<Option<f64>> = (0..draws)
        .map(|k| match spec.product(k).rate {
            DrawRate::Molar(_) | DrawRate::CutPoint(_) => None,
            DrawRate::StandardLiquidVolume(v) => v,
        })
        .collect();
    let options = ColumnOptions {
        bottom: spec.bottom_end(),
        seed_inlets: spec.steam().then_some(feeds),
    };
    // Side-stripper steam follows the main-column inlets.
    let (main, steam) = inlets.split_at(feeds + usize::from(spec.steam()));
    if targets.iter().all(Option::is_none) {
        if !spec.has_sections() {
            return distillation::calculate_with(&spec.column(None), inlets, ctx, options);
        }
        let molar: Vec<f64> = spec
            .column(None)
            .liquid_draws
            .iter()
            .map(|d| d.molar_flow_mol_s.expect("checked side-product flow"))
            .collect();
        return solve(spec, &molar, main, steam, ctx, options);
    }

    // Every component that can reach a draw needs a standard liquid volume.
    for inlet in inlets {
        if let Some(state) = &inlet.state
            && inlet.molar_flow_mol_s > 0.0
        {
            for (i, z) in state.mole_fractions.iter().enumerate() {
                if *z > 0.0 && volumes[i].is_none() {
                    return Err(format!(
                        "petroleum fractionator side products specified by standard liquid volume need pseudo-components only; {} has no specific gravity",
                        ctx.component_ids()[i]
                    ));
                }
            }
        }
    }
    // Start from the feed's mean standard volume.
    let feed_total: f64 = hydrocarbon.iter().sum();
    let feed_x: Vec<f64> = hydrocarbon.iter().map(|f| f / feed_total).collect();
    let feed_volume = mixture_volume(&volumes, &feed_x).expect("checked standard volumes");
    let mut molar: Vec<f64> = (0..draws)
        .map(|k| match (spec.product(k).rate, targets[k]) {
            (DrawRate::StandardLiquidVolume(_), target) => {
                target.expect("checked side-product volume") / feed_volume
            }
            (DrawRate::Molar(v), _) => v.expect("checked side-product flow"),
            (DrawRate::CutPoint(_), _) => unreachable!("cut points are solved at molar flows"),
        })
        .collect();
    for _ in 0..VOLUME_ITERATIONS {
        let (output, result) = solve(spec, &molar, main, steam, ctx, options).map_err(|e| {
            format!("petroleum fractionator at side-product molar flows {molar:?} mol/s: {e}")
        })?;
        let mut converged = true;
        for (k, target) in targets.iter().enumerate() {
            let Some(target) = target else { continue };
            let state = output.outlets[2 + k]
                .state
                .as_ref()
                .expect("column products carry a state");
            let volume =
                mixture_volume(&volumes, &state.mole_fractions).expect("checked standard volumes");
            let next = target / volume;
            if (next - molar[k]).abs() > VOLUME_TOLERANCE * molar[k] {
                converged = false;
            }
            molar[k] = next;
        }
        if converged {
            return Ok((output, result));
        }
    }
    Err(format!(
        "petroleum fractionator side-product volume flows did not converge in {VOLUME_ITERATIONS} column solves"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steam_feeds_the_bottom_stage_and_drops_the_bottom_specification() {
        let spec = PetroleumFractionatorSpec {
            feeds: 2,
            feed_stages: vec![8, 9],
            internals: ColumnInternals::Trays(8),
            bottom: FractionatorBottom::StrippingSteam,
            pressure_pa: Some(2e5),
            top: TopSpecification::RefluxRatio(Some(2.0)),
            ..Default::default()
        };
        let column = spec.column(None);
        assert_eq!(column.feeds, 3);
        assert_eq!(column.feed_stages, [8, 9, 10]);
        assert!(spec.issues().is_empty(), "{:?}", spec.issues());
        let reboiled = PetroleumFractionatorSpec {
            bottom: FractionatorBottom::default(),
            ..spec
        };
        assert_eq!(reboiled.issues()[0].field, "distillate fraction");
    }

    #[test]
    fn volume_products_are_checked_separately() {
        let spec = PetroleumFractionatorSpec {
            side_draws: 2,
            products: vec![
                SideProduct {
                    stage: 4,
                    rate: DrawRate::StandardLiquidVolume(Some(-1.0)),
                },
                SideProduct {
                    stage: 7,
                    rate: DrawRate::Molar(None),
                },
            ],
            pressure_pa: Some(2e5),
            top: TopSpecification::RefluxRatio(Some(2.0)),
            bottom: FractionatorBottom::Reboiler(BottomSpecification::DistillateFraction(Some(
                0.3,
            ))),
            ..Default::default()
        };
        let fields: Vec<_> = spec.issues().iter().map(|i| i.field).collect();
        assert_eq!(
            fields,
            [
                "side-draw molar flow",
                "side-product standard liquid volume flow"
            ]
        );
        let x = [0.5, 0.5];
        assert_eq!(mixture_volume(&[Some(1.0), None], &x), None);
        assert_eq!(mixture_volume(&[Some(1.0), Some(3.0)], &x), Some(2.0));
    }
}
