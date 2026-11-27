//! Adiabatic countercurrent equilibrium absorber or stripper (#187).

use thermo_types::Phase;

use super::ColumnInternals;
use super::mesh::{self, TEMPERATURE_SCALE};
use super::stages::{ColumnEnds, PlacedPackedSection, Stage, StageCascade, StageFeed};
use crate::MaterialStream;
use crate::unit::{Sign, SpecIssue, check_optional};
use crate::units::{CalcContext, UnitOutput};

/// Countercurrent contacting at one absolute pressure. Liquid
/// enters the top stage and the gas enters the bottom stage.
#[derive(Debug, Clone, PartialEq)]
pub struct AbsorberSpec {
    /// Column pressure in Pa, at most the lower inlet pressure.
    pub pressure_pa: Option<f64>,
    /// Ideal trays or a packed bed sized by HETP.
    pub internals: ColumnInternals,
}

impl Default for AbsorberSpec {
    fn default() -> Self {
        Self {
            pressure_pa: None,
            internals: ColumnInternals::Trays(5),
        }
    }
}

impl AbsorberSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        check_optional(
            &mut issues,
            "column pressure",
            self.pressure_pa,
            Sign::Positive,
        );
        issues.extend(self.internals.issues());
        issues
    }
}

/// Temperature profile of a solved adiabatic absorber or stripper.
#[derive(Debug, Clone, PartialEq)]
pub struct AbsorberResult {
    /// Equilibrium-stage temperatures, from the liquid inlet to the gas inlet.
    pub stage_temperatures_k: Vec<f64>,
    /// Specified or derived packed height, m.
    pub packed_height_m: Option<f64>,
    /// Accepted Newton iterations.
    pub iterations: usize,
}

struct Absorber<'a, 'p> {
    ctx: &'a CalcContext<'p>,
    pressure: f64,
    temperature_range: [f64; 2],
    stage_count: usize,
    gas: StageFeed,
    liquid: StageFeed,
    active: Vec<usize>,
    total_flow: f64,
}

fn logistic(value: f64) -> f64 {
    if value >= 0.0 {
        1.0 / (1.0 + (-value).exp())
    } else {
        let ratio = value.exp();
        ratio / (1.0 + ratio)
    }
}

impl Absorber<'_, '_> {
    fn decode(&self, values: &[f64]) -> StageCascade {
        let n = self.stage_count;
        let ratios = self.active.len() - 1;
        let expand = |v| mesh::fractions(v, &self.active, self.ctx.component_count());
        let gas_out = self.total_flow * logistic(values[0]);
        // Every internal cut has L_j - V_(j+1) = L_feed - V_out.
        // The positive free term preserves both internal flows and the net.
        let offset = self.liquid.molar_flow_mol_s - gas_out;
        let pairs: Vec<_> = values[1..n]
            .iter()
            .map(|value| {
                let free = self.total_flow * value.exp();
                (free + offset.max(0.0), free + (-offset).max(0.0))
            })
            .collect();
        let mut stages = Vec::with_capacity(n);
        for (j, vars) in values[n..].chunks_exact(1 + 2 * ratios).enumerate() {
            let liquid = pairs
                .get(j)
                .map_or(self.total_flow - gas_out, |pair| pair.0);
            let vapor = if j == 0 { gas_out } else { pairs[j - 1].1 };
            stages.push(Stage::equilibrium(
                vars[0] * TEMPERATURE_SCALE,
                self.pressure,
                liquid,
                vapor,
                expand(&vars[1..=ratios]),
                expand(&vars[ratios + 1..]),
            ));
        }
        StageCascade {
            stages,
            feeds: vec![self.liquid.clone(), self.gas.clone()],
        }
    }

    fn residual(&self, values: &[f64]) -> Result<Vec<f64>, String> {
        let cascade = self.decode(values);
        let properties =
            mesh::properties(self.ctx, &cascade, ColumnEnds::NONE, self.temperature_range)?;
        mesh::residuals(
            &cascade,
            &properties,
            ColumnEnds::NONE,
            &self.active,
            self.total_flow,
        )
    }

    fn initial(&self) -> Result<Vec<f64>, String> {
        // Flash a blend of both feeds at the total incoming enthalpy.
        // It is the exact one-stage solution and seeds each contact stage.
        let l = self.liquid.molar_flow_mol_s / self.total_flow;
        let v = self.gas.molar_flow_mol_s / self.total_flow;
        let z: Vec<_> = self
            .liquid
            .z
            .iter()
            .zip(&self.gas.z)
            .map(|(x, y)| l * x + v * y)
            .collect();
        let h = l * self.liquid.molar_enthalpy_j_mol + v * self.gas.molar_enthalpy_j_mol;
        let seed = self
            .ctx
            .ph(&z, self.pressure, h)
            .map_err(|e| format!("absorber feed-blend PH initialization: {e}"))?;
        let [gas, liquid] = seed.phases.as_slice() else {
            return Err("absorber feed-blend PH initialization requires two phases".into());
        };
        if gas.phase != Phase::Vapor
            || liquid.phase != Phase::Liquid
            || gas.phase_fraction <= 0.0
            || liquid.phase_fraction <= 0.0
        {
            return Err("absorber feed-blend PH initialization requires positive liquid and vapor fractions".into());
        }
        let mut initial = vec![(gas.phase_fraction / liquid.phase_fraction).ln()];
        // When the contacts have equal compositions initially, the two
        // inlet flows give the scale for the free part of each cut.
        initial.resize(self.stage_count, l.min(v).ln());
        let mut stage = vec![seed.temperature_k / TEMPERATURE_SCALE];
        stage.extend(mesh::log_ratios(&liquid.composition, &self.active));
        stage.extend(mesh::log_ratios(&gas.composition, &self.active));
        for _ in 0..self.stage_count {
            initial.extend(&stage);
        }
        Ok(initial)
    }
}

fn inlet(
    stream: &MaterialStream,
    stage: usize,
    count: usize,
    pressure: f64,
    required: Phase,
    port: &str,
) -> Result<StageFeed, String> {
    let state = stream
        .state
        .as_ref()
        .ok_or_else(|| format!("absorber {port} needs a flowing feed"))?;
    if !stream.molar_flow_mol_s.is_finite() || stream.molar_flow_mol_s <= 0.0 {
        return Err(format!("absorber {port} needs a positive finite flow"));
    }
    if state.pressure_pa < pressure {
        return Err(format!("absorber pressure exceeds {port} pressure"));
    }
    if !matches!(state.phases.as_slice(), [phase] if phase.phase == required) {
        return Err(format!(
            "absorber {port} requires a single {required:?} phase at its inlet state"
        ));
    }
    StageFeed::from_stream(stage, stream, count)
}

pub(super) fn calculate(
    spec: &AbsorberSpec,
    inputs: &[MaterialStream],
    ctx: &CalcContext<'_>,
) -> Result<(UnitOutput, AbsorberResult), String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!("absorber {} {}", issue.field, issue.problem));
    }
    let (stage_count, bed) = spec.internals.resolve()?;
    let pressure = spec.pressure_pa.expect("checked absorber pressure");
    // The column's thermodynamic domain is the package's: stage temperatures
    // stay inside its bubble/dew interval, and the stage equations and checks
    // use only its PH, homogeneous-phase, TP and VLLE operations, which reject
    // states outside their own domains.
    let temperature_range = mesh::temperature_range(ctx)?;
    let c = ctx.component_count();
    let gas = inlet(
        &inputs[0],
        stage_count,
        c,
        pressure,
        Phase::Vapor,
        "gas inlet",
    )?;
    let liquid = inlet(&inputs[1], 1, c, pressure, Phase::Liquid, "liquid inlet")?;
    let active: Vec<f64> = gas.z.iter().zip(&liquid.z).map(|(y, x)| y + x).collect();
    let active = mesh::active_components(&active)?;
    let total_flow = gas.molar_flow_mol_s + liquid.molar_flow_mol_s;
    if !total_flow.is_finite() {
        return Err("absorber total feed flow is not finite".into());
    }
    let column = Absorber {
        ctx,
        pressure,
        temperature_range,
        stage_count,
        gas,
        liquid,
        active,
        total_flow,
    };
    let (solved, iterations) = mesh::solve(column.initial()?, |v| column.residual(v))?;
    let cascade = column.decode(&solved);
    let warnings = mesh::check_liquid_liquid(ctx, &cascade)?;
    mesh::verify_phases(ctx, &cascade)?;
    ColumnEnds::NONE.check(&cascade)?;
    if let Some(packed) = bed {
        PlacedPackedSection {
            top_stage: 1,
            packed,
        }
        .check(&cascade, ColumnEnds::NONE)?;
    }
    let top = &cascade.stages[0];
    let bottom = &cascade.stages[stage_count - 1];
    let outlets = vec![
        mesh::product(ctx, top, top.vapor_mol_s, Phase::Vapor)?,
        mesh::product(ctx, bottom, bottom.liquid_mol_s, Phase::Liquid)?,
    ];
    Ok((
        UnitOutput::new(outlets, 0.0, warnings),
        AbsorberResult {
            stage_temperatures_k: cascade
                .stages
                .iter()
                .map(|stage| stage.temperature_k)
                .collect(),
            packed_height_m: bed.map(|packed| packed.packed_height_m),
            iterations,
        },
    ))
}
