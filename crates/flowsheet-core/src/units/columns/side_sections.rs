//! Pump-arounds and side strippers attached to the rigorous column (#190).
//!
//! Both are solved simultaneously with the main column in one Newton system
//! ([`super::distillation`]): the unknowns of every side stripper follow the
//! main-column profile, its MESH residuals follow the main-column residuals,
//! and the streams that couple the cascades are feeds and draws whose flows,
//! compositions and enthalpies are functions of the unknowns.
//!
//! - A **pump-around** draws a fixed molar flow of liquid from main stage
//!   `draw` (at the stage liquid state), cools it and returns it as a feed to
//!   stage `ret ≤ draw` at the return-stage pressure. Its cooling is either a
//!   return temperature (the returned liquid's enthalpy on the liquid root at
//!   that temperature) or a duty `Q` removed (`h_return = h_draw − Q / P`).
//! - A **side stripper** is a cascade of `stages` equilibrium stages at the
//!   pressure of its draw stage. Liquid from main stage `draw` enters its
//!   top stage, its top vapor returns as a feed to main stage `ret ≤ draw`,
//!   and its bottom liquid is the side product, with a fixed molar flow `B`.
//!   It is stripped by steam entering its adiabatic bottom stage, or by a
//!   reboiler with a duty or boilup-ratio specification. The draw from the
//!   main column follows from the stripper's total balance,
//!   `S = B + V_1 − F_steam`, so the stripper adds no flow unknown beyond its
//!   stage vapor flows.

use super::mesh::{self, TEMPERATURE_SCALE};
use super::stages::{ColumnEnds, Stage, StageCascade, StageFeed, StageProperties};
use crate::units::CalcContext;
use thermo_models::RootSelection;

/// Same enthalpy scale as the MESH energy residuals.
const ENTHALPY_SCALE: f64 = 10_000.0;

/// How a pump-around is cooled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Cooling {
    /// Return temperature, K.
    ReturnTemperature(f64),
    /// Heat removed, W; positive.
    Duty(f64),
}

/// A liquid pump-around of the main column; stages are zero-based.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct PumpAroundLoop {
    pub(super) draw: usize,
    pub(super) ret: usize,
    pub(super) flow: f64,
    pub(super) cooling: Cooling,
}

impl PumpAroundLoop {
    /// Molar enthalpy of the returned liquid, J/mol.
    pub(super) fn return_enthalpy(
        &self,
        ctx: &CalcContext<'_>,
        cascade: &StageCascade,
        properties: &[StageProperties],
    ) -> Result<f64, String> {
        match self.cooling {
            Cooling::ReturnTemperature(t) => Ok(ctx
                .homogeneous(
                    &cascade.stages[self.draw].x,
                    t,
                    cascade.stages[self.ret].pressure_pa,
                    RootSelection::SmallestVolume,
                )?
                .enthalpy),
            Cooling::Duty(q) => Ok(properties[self.draw].liquid_enthalpy_j_mol - q / self.flow),
        }
    }
}

/// Stripping agent and bottom specification of a side stripper.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum Stripping {
    /// Steam (any vapor stream) fed to the adiabatic bottom stage.
    Steam(StageFeed),
    /// Partial reboiler with its duty, W.
    ReboilerDuty(f64),
    /// Partial reboiler with the boilup ratio `V_N / B`.
    BoilupRatio(f64),
}

/// A side stripper on a main-column draw; stages are zero-based.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct StripperSection {
    pub(super) draw: usize,
    pub(super) ret: usize,
    pub(super) stages: usize,
    /// Bottom product `B`, mol/s.
    pub(super) bottoms: f64,
    pub(super) bottom: Stripping,
}

impl StripperSection {
    pub(super) fn steam_flow(&self) -> f64 {
        match &self.bottom {
            Stripping::Steam(feed) => feed.molar_flow_mol_s,
            _ => 0.0,
        }
    }

    pub(super) fn reboiled(&self) -> bool {
        !matches!(self.bottom, Stripping::Steam(_))
    }

    pub(super) fn ends(&self) -> ColumnEnds {
        ColumnEnds {
            condenser: None,
            reboiler: self.reboiled(),
        }
    }

    /// Number of unknowns: per stage the temperature, `m` liquid and `m`
    /// vapor log ratios and the log vapor flow.
    pub(super) fn variable_len(&self, m: usize) -> usize {
        self.stages * (2 * m + 2)
    }

    /// Top vapor flow and composition.
    pub(super) fn top_vapor(
        &self,
        values: &[f64],
        active: &[usize],
        components: usize,
        flow: f64,
    ) -> (f64, Vec<f64>) {
        let m = active.len() - 1;
        (
            flow * values[1 + 2 * m].exp(),
            mesh::fractions(&values[1 + m..1 + 2 * m], active, components),
        )
    }

    /// Liquid drawn from the main column for a top vapor flow.
    pub(super) fn draw_flow(&self, top_vapor: f64) -> f64 {
        self.bottoms + top_vapor - self.steam_flow()
    }

    /// The stripper cascade. Its liquid feed (the main-column draw at
    /// composition `x_draw`) carries a placeholder enthalpy, set from the
    /// main-column properties before the residuals.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn decode(
        &self,
        values: &[f64],
        active: &[usize],
        components: usize,
        flow: f64,
        pressure: f64,
        draw: f64,
        x_draw: Vec<f64>,
    ) -> Result<StageCascade, String> {
        let m = active.len() - 1;
        let width = 2 * m + 2;
        let vapor: Vec<f64> = (0..self.stages)
            .map(|j| flow * values[j * width + 1 + 2 * m].exp())
            .collect();
        let lift = self.bottoms - self.steam_flow();
        let mut stages = Vec::with_capacity(self.stages);
        for (j, v) in vapor.iter().enumerate() {
            let chunk = &values[j * width..(j + 1) * width];
            let liquid = if j + 1 == self.stages {
                self.bottoms
            } else {
                lift + vapor[j + 1]
            };
            if !(liquid.is_finite() && liquid > 0.0 && v.is_finite()) {
                return Err("side stripper stage liquid flow must be positive".into());
            }
            stages.push(Stage::equilibrium(
                chunk[0] * TEMPERATURE_SCALE,
                pressure,
                liquid,
                *v,
                mesh::fractions(&chunk[1..1 + m], active, components),
                mesh::fractions(&chunk[1 + m..1 + 2 * m], active, components),
            ));
        }
        let mut feeds = vec![StageFeed {
            stage: 1,
            molar_flow_mol_s: draw,
            z: x_draw,
            molar_enthalpy_j_mol: 0.0,
        }];
        if let Stripping::Steam(steam) = &self.bottom {
            feeds.push(StageFeed {
                stage: self.stages,
                ..steam.clone()
            });
        }
        Ok(StageCascade { stages, feeds })
    }

    /// MESH residuals of the stripper and, with a reboiler, its bottom
    /// specification.
    pub(super) fn residuals(
        &self,
        cascade: &StageCascade,
        properties: &[StageProperties],
        active: &[usize],
        flow: f64,
    ) -> Result<Vec<f64>, String> {
        let mut out = mesh::residuals(cascade, properties, self.ends(), active, flow)?;
        let last = self.stages - 1;
        match self.bottom {
            Stripping::Steam(_) => {}
            Stripping::ReboilerDuty(q) => {
                let balances = cascade.residuals(properties)?;
                out.push((-balances[last].energy_w - q) / flow / ENTHALPY_SCALE);
            }
            Stripping::BoilupRatio(r) => {
                out.push((cascade.stages[last].vapor_mol_s - r * self.bottoms) / flow);
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err("side stripper has a nonfinite residual".into());
        }
        Ok(out)
    }

    /// Starting unknowns: every stage at the main column's starting tray
    /// state (temperature and log ratios) and a uniform vapor flow.
    pub(super) fn initial(&self, temperature: f64, x: &[f64], y: &[f64], flow: f64) -> Vec<f64> {
        let vapor = match &self.bottom {
            Stripping::Steam(steam) => steam.molar_flow_mol_s + 0.25 * self.bottoms,
            Stripping::ReboilerDuty(_) => 0.5 * self.bottoms,
            Stripping::BoilupRatio(r) => r * self.bottoms,
        };
        let mut out = Vec::with_capacity(self.variable_len(x.len()));
        for _ in 0..self.stages {
            out.push(temperature / TEMPERATURE_SCALE);
            out.extend(x);
            out.extend(y);
            out.push((vapor / flow).ln());
        }
        out
    }
}

/// Pump-arounds and side strippers of one column.
#[derive(Debug, Clone, Default, PartialEq)]
pub(super) struct SideSections {
    pub(super) pump_arounds: Vec<PumpAroundLoop>,
    pub(super) strippers: Vec<StripperSection>,
}

/// Converged profile of one side stripper (#190).
#[derive(Debug, Clone, PartialEq)]
pub struct SideStripperProfile {
    /// Stage temperatures, top first, K.
    pub stage_temperatures_k: Vec<f64>,
    /// Liquid drawn from the main column, mol/s.
    pub draw_mol_s: f64,
    /// Top vapor returned to the main column, mol/s.
    pub vapor_return_mol_s: f64,
    /// Reboiler duty, W; zero when steam-stripped.
    pub reboiler_duty_w: f64,
}
