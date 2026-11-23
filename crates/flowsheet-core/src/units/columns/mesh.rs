//! Shared thermodynamic validation and Newton iteration for equilibrium cascades.

use property_packages::{Convergence, EquilibriumState};
use thermo_models::RootSelection;
use thermo_solvers::{FlashPhase, SaturationOptions};
use thermo_types::Phase;

use super::stages::{ColumnEnds, Condenser, Stage, StageCascade, StageProperties};
use crate::MaterialStream;
use crate::numerics::solve_dense;
use crate::units::CalcContext;

pub(super) const TEMPERATURE_SCALE: f64 = 100.0;
const ENTHALPY_SCALE: f64 = 10_000.0;
const MIN_VAPOR_LIQUID_VOLUME_RATIO: f64 = 1.05;
const MAX_ITERATIONS: usize = 40;
const RESIDUAL_TOLERANCE: f64 = 2e-11;

pub(super) fn temperature_range(ctx: &CalcContext<'_>) -> Result<[f64; 2], String> {
    let options = SaturationOptions::default();
    let [lo, hi] = ctx.package().temperature_interval();
    let lo = lo.max(options.minimum_temperature_k);
    let hi = hi.min(options.maximum_temperature_k);
    if lo >= hi {
        Err("package temperature interval does not overlap the column VLE interval".into())
    } else {
        Ok([lo, hi])
    }
}

pub(super) fn active_components(z: &[f64]) -> Result<Vec<usize>, String> {
    let active: Vec<_> = (0..z.len()).filter(|&i| z[i] > 0.0).collect();
    if active.len() < 2 {
        Err("column requires at least two components in the feed".into())
    } else {
        Ok(active)
    }
}

pub(super) fn log_ratios(phase: &[f64], active: &[usize]) -> Vec<f64> {
    let (&reference, others) = active.split_last().expect("checked active components");
    others
        .iter()
        .map(|&i| phase[i].ln() - phase[reference].ln())
        .collect()
}

/// Normalized active fractions use log ratios to the last active component.
/// Components absent from the feeds stay exactly zero.
pub(super) fn fractions(log_ratios: &[f64], active: &[usize], count: usize) -> Vec<f64> {
    let max = log_ratios.iter().copied().fold(0.0, f64::max);
    let weights: Vec<_> = log_ratios
        .iter()
        .copied()
        .chain([0.0])
        .map(|x| (x - max).exp())
        .collect();
    let total: f64 = weights.iter().sum();
    let mut result = vec![0.0; count];
    for (&index, weight) in active.iter().zip(weights) {
        result[index] = weight / total;
    }
    result
}

pub(super) fn properties(
    ctx: &CalcContext<'_>,
    cascade: &StageCascade,
    ends: ColumnEnds,
    [lo, hi]: [f64; 2],
) -> Result<Vec<StageProperties>, String> {
    cascade.validate(ctx.component_count())?;
    cascade
        .stages
        .iter()
        .enumerate()
        .map(|(j, stage)| {
            let t = stage.temperature_k;
            let p = stage.pressure_pa;
            let total_condenser = j == 0 && ends.condenser == Some(Condenser::Total);
            if !(lo..=hi).contains(&t) {
                return Err(format!(
                    "stage {} temperature is outside {lo}..={hi} K",
                    j + 1
                ));
            }
            if stage.liquid_mol_s <= 0.0 || (!total_condenser && stage.vapor_mol_s <= 0.0) {
                return Err(format!(
                    "stage {} has no internal liquid or vapor flow",
                    j + 1
                ));
            }
            let liquid = ctx.homogeneous(&stage.x, t, p, RootSelection::SmallestVolume)?;
            let vapor = ctx.homogeneous(&stage.y, t, p, RootSelection::LargestVolume)?;
            if vapor.molar_volume <= MIN_VAPOR_LIQUID_VOLUME_RATIO * liquid.molar_volume {
                return Err(format!(
                    "stage {} lacks distinct liquid and vapor roots",
                    j + 1
                ));
            }
            Ok(StageProperties {
                k: liquid
                    .ln_fugacity_coefficients
                    .iter()
                    .zip(&vapor.ln_fugacity_coefficients)
                    .map(|(l, v)| (l - v).exp())
                    .collect(),
                liquid_enthalpy_j_mol: liquid.enthalpy,
                vapor_enthalpy_j_mol: vapor.enthalpy,
            })
        })
        .collect()
}

pub(super) fn residuals(
    cascade: &StageCascade,
    properties: &[StageProperties],
    ends: ColumnEnds,
    active: &[usize],
    flow: f64,
) -> Result<Vec<f64>, String> {
    let balances = cascade.residuals(properties)?;
    let mut out = Vec::new();
    for (j, ((stage, property), balance)) in cascade
        .stages
        .iter()
        .zip(properties)
        .zip(balances)
        .enumerate()
    {
        let total_condenser = j == 0 && ends.condenser == Some(Condenser::Total);
        let end_duty =
            (j == 0 && ends.condenser.is_some()) || (j + 1 == cascade.len() && ends.reboiler);
        out.extend(
            active
                .iter()
                .map(|&i| stage.y[i].ln() - stage.x[i].ln() - property.k[i].ln()),
        );
        // The total condenser inherits the incoming vapor's composition.
        if !total_condenser {
            out.extend(
                active[..active.len() - 1]
                    .iter()
                    .map(|&i| balance.component_mol_s[i] / flow),
            );
        }
        if !end_duty {
            out.push(balance.energy_w / flow / ENTHALPY_SCALE);
        }
    }
    if out.iter().any(|value| !value.is_finite()) {
        return Err("column has a nonfinite MESH residual".into());
    }
    Ok(out)
}

pub(super) fn verify_phases(ctx: &CalcContext<'_>, cascade: &StageCascade) -> Result<(), String> {
    for (j, stage) in cascade.stages.iter().enumerate() {
        // The midpoint of each tie line must independently flash to both
        // solved phases. Forced cubic roots alone do not certify VLE.
        let midpoint: Vec<_> = stage
            .x
            .iter()
            .zip(&stage.y)
            .map(|(x, y)| 0.5 * (x + y))
            .collect();
        let flashed = ctx
            .tp(&midpoint, stage.temperature_k, stage.pressure_pa)
            .map_err(|e| format!("stage {} phase check: {e}", j + 1))?;
        let matches = |phase: &FlashPhase, expected: &[f64]| {
            (phase.phase_fraction - 0.5).abs() <= 2e-5
                && phase
                    .composition
                    .iter()
                    .zip(expected)
                    .all(|(actual, given)| (actual - given).abs() <= 2e-6)
        };
        if !matches!(flashed.phases.as_slice(), [v, l]
            if v.phase == Phase::Vapor && l.phase == Phase::Liquid
                && matches(v, &stage.y) && matches(l, &stage.x))
        {
            return Err(format!(
                "stage {} did not reproduce a stable vapor–liquid split",
                j + 1
            ));
        }
    }
    Ok(())
}

pub(super) fn product(
    ctx: &CalcContext<'_>,
    stage: &Stage,
    flow: f64,
    phase: Phase,
) -> Result<MaterialStream, String> {
    let (composition, root) = match phase {
        Phase::Liquid => (&stage.x, RootSelection::SmallestVolume),
        Phase::Vapor => (&stage.y, RootSelection::LargestVolume),
        _ => unreachable!("column product must be liquid or vapor"),
    };
    let state = ctx.homogeneous(composition, stage.temperature_k, stage.pressure_pa, root)?;
    Ok(MaterialStream {
        molar_flow_mol_s: flow,
        state: Some(EquilibriumState {
            temperature_k: stage.temperature_k,
            pressure_pa: stage.pressure_pa,
            mole_fractions: composition.clone(),
            phases: vec![FlashPhase {
                phase,
                phase_fraction: 1.0,
                composition: composition.clone(),
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

fn norm(values: &[f64]) -> f64 {
    values.iter().map(|x| x.abs()).fold(0.0, f64::max)
}

pub(super) fn solve(
    mut variables: Vec<f64>,
    residual: impl Fn(&[f64]) -> Result<Vec<f64>, String>,
) -> Result<(Vec<f64>, usize), String> {
    let mut r = residual(&variables).map_err(|e| format!("column initialization: {e}"))?;
    for iteration in 0..=MAX_ITERATIONS {
        let error = norm(&r);
        if error <= RESIDUAL_TOLERANCE {
            return Ok((variables, iteration));
        }
        if iteration == MAX_ITERATIONS {
            break;
        }
        let n = variables.len();
        let mut jacobian = vec![vec![0.0; n]; n];
        let mut trial = variables.clone();
        for (k, &value) in variables.iter().enumerate() {
            let h = 1e-6 * value.abs().max(1.0);
            trial[k] = value + h;
            let shifted = match residual(&trial) {
                Ok(shifted) => shifted,
                Err(_) => {
                    trial[k] = value - h;
                    residual(&trial).map_err(|e| format!("column Jacobian: {e}"))?
                }
            };
            let delta = trial[k] - value;
            for (row, (&new, old)) in jacobian.iter_mut().zip(shifted.iter().zip(&r)) {
                row[k] = (new - old) / delta;
            }
            trial[k] = value;
        }
        let step = solve_dense(jacobian, r.iter().map(|v| -v).collect())
            .ok_or("column MESH Jacobian is singular")?;
        let mut accepted = None;
        for backtrack in 0..24 {
            let scale = 0.5_f64.powi(backtrack);
            for ((trial, value), delta) in trial.iter_mut().zip(&variables).zip(&step) {
                *trial = value + scale * delta;
            }
            if let Ok(rt) = residual(&trial)
                && norm(&rt) < error
            {
                accepted = Some(rt);
                break;
            }
        }
        let Some(next) = accepted else {
            return Err(format!(
                "column MESH line search stalled (scaled residual {error:.3e})"
            ));
        };
        variables = trial;
        r = next;
    }
    Err(format!(
        "column MESH did not converge in {MAX_ITERATIONS} iterations (scaled residual {:.3e})",
        norm(&r)
    ))
}
