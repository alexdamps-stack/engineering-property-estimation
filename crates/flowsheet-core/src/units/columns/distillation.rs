//! Constant-pressure MESH solve for a total condenser and partial reboiler.

use property_packages::SaturationKind;
use thermo_types::Phase;

use super::mesh::{self, TEMPERATURE_SCALE};
use super::stages::{ColumnEnds, PlacedPackedSection, Stage, StageCascade, StageFeed};
use super::{ColumnInternals, count_issue};
use crate::MaterialStream;
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional};
use crate::units::{CalcContext, UnitOutput};

/// Liquid side draw from an internal stage, in port order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LiquidSideDraw {
    /// One-based stage number, including the total condenser as stage 1.
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

/// Rigorous distillation column specification (#186).
#[derive(Debug, Clone, PartialEq)]
pub struct DistillationColumnSpec {
    /// Number of material feed ports.
    pub feeds: u16,
    /// Number of material side-draw ports.
    pub side_draws: u16,
    /// One-based feed stage per port. Unset entries default to stage 2.
    pub feed_stages: Vec<u16>,
    /// Specification of each side-draw port. Unset entries have no flow specification.
    pub liquid_draws: Vec<LiquidSideDraw>,
    /// Ideal trays or ideal packed stages.
    pub internals: ColumnInternals,
    /// Uniform absolute pressure, Pa.
    pub pressure_pa: Option<f64>,
    /// Molar reflux divided by molar distillate flow.
    pub reflux_ratio: Option<f64>,
    /// Fraction of the total feed withdrawn as distillate, strictly between 0 and 1.
    pub distillate_fraction: Option<f64>,
}

impl Default for DistillationColumnSpec {
    fn default() -> Self {
        Self {
            feeds: 1,
            side_draws: 0,
            feed_stages: Vec::new(),
            liquid_draws: Vec::new(),
            internals: ColumnInternals::Trays(5),
            pressure_pa: None,
            reflux_ratio: None,
            distillate_fraction: None,
        }
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

    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        count_issue(&mut issues, "feeds", self.feeds, 1);
        count_issue(&mut issues, "side draws", self.side_draws, 0);
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
        check_optional(
            &mut issues,
            "column pressure",
            self.pressure_pa,
            Sign::Positive,
        );
        check_optional(
            &mut issues,
            "reflux ratio",
            self.reflux_ratio,
            Sign::Positive,
        );
        check_optional(
            &mut issues,
            "distillate fraction",
            self.distillate_fraction,
            Sign::Fraction,
        );
        if matches!(self.distillate_fraction, Some(0.0 | 1.0)) {
            issues.push(SpecIssue {
                field: "distillate fraction",
                problem: SpecProblem::OutOfRange,
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
        issues
    }
}

/// Accepted solution evidence for a rigorous distillation column.
#[derive(Debug, Clone, PartialEq)]
pub struct DistillationColumnResult {
    /// Stage temperatures from the total condenser to the reboiler.
    pub stage_temperatures_k: Vec<f64>,
    /// Heat removed by the condenser, W.
    pub condenser_duty_w: f64,
    /// Heat supplied to the reboiler, W.
    pub reboiler_duty_w: f64,
    /// Bed height supplied or derived from the specified HETP.
    pub packed_height_m: Option<f64>,
    /// Accepted Newton iterations.
    pub iterations: usize,
}

struct Column<'a, 'p> {
    ctx: &'a CalcContext<'p>,
    pressure: f64,
    temperature_range: [f64; 2],
    feeds: Vec<StageFeed>,
    draws_by_stage: Vec<f64>,
    active_components: Vec<usize>,
    flow: f64,
    distillate: f64,
    reflux: f64,
    /// Cumulative feed less draw and distillate through each stage.
    liquid_offset: Vec<f64>,
}

impl Column<'_, '_> {
    fn decode(&self, values: &[f64]) -> StageCascade {
        let n = self.liquid_offset.len();
        let c = self.ctx.component_count();
        let m = self.active_components.len() - 1;
        let expand = |v| mesh::fractions(v, &self.active_components, c);
        let mut tail = values;
        let mut liquid = vec![0.0; n];
        let mut stages: Vec<Stage> = Vec::with_capacity(n);
        liquid[0] = self.reflux;
        liquid[n - 1] = self.liquid_offset[n - 1];
        for j in 0..n {
            let t = tail[0] * TEMPERATURE_SCALE;
            tail = &tail[1..];
            let x = if j == 0 {
                Vec::new()
            } else {
                let (x, rest) = tail.split_at(m);
                tail = rest;
                expand(x)
            };
            let (y, rest) = tail.split_at(m);
            tail = rest;
            let vapor = match j {
                0 => 0.0,
                1 => self.reflux + self.distillate,
                _ => {
                    let extra = self.flow * tail[0].exp();
                    tail = &tail[1..];
                    let offset = self.liquid_offset[j - 1];
                    liquid[j - 1] = extra + offset.max(0.0);
                    extra + (-offset).max(0.0)
                }
            };
            let mut stage = Stage::equilibrium(t, self.pressure, 0.0, vapor, x, expand(y));
            stage.liquid_draw_mol_s = if j == 0 {
                self.distillate
            } else {
                self.draws_by_stage[j]
            };
            stages.push(stage);
        }
        stages[0].x = stages[1].y.clone();
        for (s, l) in stages.iter_mut().zip(liquid) {
            s.liquid_mol_s = l;
        }
        StageCascade {
            stages,
            feeds: self.feeds.clone(),
        }
    }

    fn residual(&self, values: &[f64]) -> Result<Vec<f64>, String> {
        let cascade = self.decode(values);
        let ends = ColumnEnds::DISTILLATION;
        let properties = mesh::properties(self.ctx, &cascade, ends, self.temperature_range)?;
        mesh::residuals(
            &cascade,
            &properties,
            ends,
            &self.active_components,
            self.flow,
        )
    }
}

pub(super) fn calculate(
    spec: &DistillationColumnSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, DistillationColumnResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!("column {} {}", issue.field, issue.problem));
    }
    let (n, packed) = spec.stage_layout()?;
    let pressure = spec.pressure_pa.expect("checked column pressure");
    let temperature_range = mesh::temperature_range(ctx)?;
    let mut totals = vec![0.0; n];
    let mut feed_composition = vec![0.0; ctx.component_count()];
    let mut feeds = Vec::with_capacity(inlets.len());
    for (i, inlet) in inlets.iter().enumerate() {
        let stage = usize::from(spec.feed_stages.get(i).copied().unwrap_or(2));
        let feed = StageFeed::from_stream(stage, inlet, ctx.component_count())?;
        if feed.molar_flow_mol_s > 0.0
            && inlet
                .state
                .as_ref()
                .is_some_and(|state| state.pressure_pa < pressure)
        {
            return Err(format!(
                "column pressure exceeds feed {} pressure; pressurize the feed upstream",
                i + 1
            ));
        }
        totals[stage - 1] += feed.molar_flow_mol_s;
        for (total, z) in feed_composition.iter_mut().zip(&feed.z) {
            *total += feed.molar_flow_mol_s * z;
        }
        feeds.push(feed);
    }
    let flow: f64 = totals.iter().sum();
    if !flow.is_finite() || flow <= 0.0 {
        return Err("column requires a positive, finite total feed flow".into());
    }
    feed_composition.iter_mut().for_each(|z| *z /= flow);
    let active_components = mesh::active_components(&feed_composition)?;
    let draws = &spec.liquid_draws[..usize::from(spec.side_draws)];
    let mut draws_by_stage = vec![0.0; n];
    for draw in draws {
        draws_by_stage[usize::from(draw.stage) - 1] +=
            draw.molar_flow_mol_s.expect("checked side-draw flow");
    }
    let distillate = flow
        * spec
            .distillate_fraction
            .expect("checked distillate fraction");
    let reflux = distillate * spec.reflux_ratio.expect("checked reflux ratio");
    let mut remaining = -distillate;
    let liquid_offset: Vec<_> = totals
        .iter()
        .zip(&draws_by_stage)
        .map(|(f, d)| {
            remaining += f - d;
            remaining
        })
        .collect();
    if !remaining.is_finite() || remaining <= 0.0 || !reflux.is_finite() || reflux <= 0.0 {
        return Err("column withdrawals must leave positive bottoms flow and finite reflux".into());
    }
    let seed = ctx
        .saturation_temperature(&feed_composition, pressure, SaturationKind::Bubble)
        .map_err(|e| format!("column feed bubble-point initialization: {e}"))?;
    // The condenser liquid starts at the composition of the vapor leaving
    // the first tray. Its own bubble point can be much colder than the tray.
    let condenser_seed = ctx
        .saturation_temperature(&seed.vapor_composition, pressure, SaturationKind::Bubble)
        .map_err(|e| format!("column condenser bubble-point initialization: {e}"))?;
    let column = Column {
        ctx,
        pressure,
        temperature_range,
        feeds,
        draws_by_stage,
        active_components,
        flow,
        distillate,
        reflux,
        liquid_offset,
    };
    let ratios = |phase: &[f64]| mesh::log_ratios(phase, &column.active_components);
    let x = ratios(&seed.liquid_composition);
    let y = ratios(&seed.vapor_composition);
    let mut initial = vec![condenser_seed.temperature_k / TEMPERATURE_SCALE];
    initial.extend(ratios(&condenser_seed.vapor_composition));
    for j in 1..n {
        initial.push(seed.temperature_k / TEMPERATURE_SCALE);
        initial.extend(&x);
        initial.extend(&y);
        if j > 1 {
            initial.push(((reflux + distillate) / flow).ln());
        }
    }
    let (solved, iterations) = mesh::solve(initial, |values| column.residual(values))?;
    let mut cascade = column.decode(&solved);
    let properties = mesh::properties(ctx, &cascade, ColumnEnds::DISTILLATION, temperature_range)?;
    mesh::verify_phases(ctx, &cascade)?;
    let residuals = cascade.residuals(&properties)?;
    cascade.stages[0].duty_w = -residuals[0].energy_w;
    cascade.stages[n - 1].duty_w = -residuals[n - 1].energy_w;
    if let Some(section) = &packed {
        section.check(&cascade, ColumnEnds::DISTILLATION)?;
    }
    let energy = ColumnEnds::DISTILLATION.port_duties(&cascade)?;
    let condenser = -cascade.stages[0].duty_w;
    let reboiler = cascade.stages[n - 1].duty_w;
    if condenser <= 0.0 || reboiler <= 0.0 {
        return Err(
            "column requires heat removal at the condenser and heat input at the reboiler".into(),
        );
    }
    let mut outlets = vec![
        mesh::product(ctx, &cascade.stages[0], distillate, Phase::Liquid)?,
        mesh::product(ctx, &cascade.stages[n - 1], remaining, Phase::Liquid)?,
    ];
    for draw in draws {
        outlets.push(mesh::product(
            ctx,
            &cascade.stages[usize::from(draw.stage) - 1],
            draw.molar_flow_mol_s.expect("checked side-draw flow"),
            Phase::Liquid,
        )?);
    }
    let mut output = UnitOutput::new(outlets, reboiler - condenser, Vec::new());
    output.energy = energy;
    Ok((
        output,
        DistillationColumnResult {
            stage_temperatures_k: cascade.stages.iter().map(|s| s.temperature_k).collect(),
            condenser_duty_w: condenser,
            reboiler_duty_w: reboiler,
            packed_height_m: packed.map(|section| section.packed.packed_height_m),
            iterations,
        },
    ))
}
