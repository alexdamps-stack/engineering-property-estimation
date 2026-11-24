//! Equilibrium, Gibbs and kinetic (CSTR, plug-flow, batch) reactors.
//!
//! All reactors share the formation-basis energy balance of
//! [`crate::reaction`]: outlet states come from [`crate::engine::react`],
//! which adds the standard formation enthalpy of the generation.
//!
//! * Equilibrium and Gibbs reactors minimize the Gibbs energy of a single
//!   reacting fluid phase ([`crate::numerics::minimize_gibbs`]): subject to
//!   element conservation (Gibbs) or to the invariants of the listed
//!   reactions (equilibrium). The equilibrium composition is evaluated at the
//!   outlet temperature plus the temperature approach.
//! * Kinetic reactors evaluate rate laws on the reacting phase of a TP or PH
//!   equilibrium state at local conditions ([`rates`]).
//!
//! Adiabatic and duty specifications of the equilibrium, Gibbs and CSTR
//! reactors find the outlet temperature by Brent's method over the package
//! temperature interval, then close the energy balance exactly with a PH
//! flash at the converged composition.

use crate::engine::{MaterialStream, react};
use crate::numerics::{bracketed_root, integrate, minimize_gibbs, null_space, rank, solve_dense};
use crate::reaction::{Chemistry, ReactionTerms};
use crate::units::UnitOutput;
use crate::{
    BatchThermal, ConcentrationBasis, EquilibriumConstant, EquilibriumOptions, EquilibriumReaction,
    FugacityBasis, KineticReaction, Kinetics, PlugFlowThermal, RateBasis, RateLaw, ReactingPhase,
    SeparatorThermal,
};
use property_packages::{EquilibriumState, FlowsheetPackage, RootSelection};
use thermo_types::{MOLAR_GAS_CONSTANT as R, Phase, Temperature};

const P0: f64 = 1.0e5;

/// What an equilibrium calculation conserves.
enum Constraint<'a> {
    /// Element balances; flagged components are inert.
    Elements { inert: &'a [bool] },
    /// Invariants of the given reactions (columns of `nu`).
    Reactions {
        nu: &'a [Vec<f64>],
        constants: &'a [EquilibriumConstant],
    },
}

/// Standard Gibbs energies over RT of all package components at `t`.
fn gibbs_over_rt(pkg: &dyn FlowsheetPackage, t: f64) -> Result<Vec<f64>, String> {
    let temperature = Temperature::new(t).map_err(|e| e.to_string())?;
    pkg.formation_components()
        .map_err(|e| e.to_string())?
        .iter()
        .map(|c| {
            c.standard_gibbs_energy(temperature)
                .map(|g| g / (R * t))
                .map_err(|e| e.to_string())
        })
        .collect()
}

/// Least-squares extents `xi` with `nu xi = g` for independent reactions.
fn extents_of(columns: &[Vec<f64>], generation: &[f64]) -> Vec<f64> {
    let r = columns.len();
    let normal: Vec<Vec<f64>> = (0..r)
        .map(|a| (0..r).map(|b| dot(&columns[a], &columns[b])).collect())
        .collect();
    let rhs: Vec<f64> = columns.iter().map(|c| dot(c, generation)).collect();
    solve_dense(normal, rhs).unwrap_or_else(|| vec![f64::NAN; r])
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Equilibrium outlet flows at `t` K and `p` Pa.
fn equilibrium_flows(
    constraint: &Constraint,
    chemistry: &Chemistry,
    pkg: &dyn FlowsheetPackage,
    inlet: &[f64],
    t: f64,
    p: f64,
    fugacity: FugacityBasis,
) -> Result<Vec<f64>, String> {
    let s = inlet.len();
    let g = gibbs_over_rt(pkg, t)?;
    let ln_p = (p / P0).ln();
    let (mut c, invariants) = match constraint {
        Constraint::Elements { inert } => {
            let mut rows: Vec<Vec<f64>> = (0..chemistry.elements().len())
                .map(|k| {
                    chemistry
                        .element_column(k)
                        .into_iter()
                        .zip(inert.iter())
                        .map(|(a, &i)| if i { 0.0 } else { a })
                        .collect()
                })
                .collect();
            for (i, &is_inert) in inert.iter().enumerate() {
                if is_inert {
                    let mut row = vec![0.0; s];
                    row[i] = 1.0;
                    rows.push(row);
                }
            }
            (g, rows)
        }
        Constraint::Reactions { nu, constants } => {
            // ln K per reaction, then the least-norm standard potentials that
            // reproduce them: c = -nu (nu^T nu)^-1 ln K.
            let ln_k: Vec<f64> = nu
                .iter()
                .zip(constants.iter())
                .map(|(col, k)| match k {
                    EquilibriumConstant::GibbsEnergy => -dot(col, &g),
                    EquilibriumConstant::Correlation(c) => c.ln_value(t),
                })
                .collect();
            let r = nu.len();
            let normal: Vec<Vec<f64>> = (0..r)
                .map(|a| (0..r).map(|b| dot(&nu[a], &nu[b])).collect())
                .collect();
            let y = solve_dense(normal, ln_k).ok_or("reactions are linearly dependent")?;
            let c: Vec<f64> = (0..s)
                .map(|i| -(0..r).map(|j| nu[j][i] * y[j]).sum::<f64>())
                .collect();
            // Each column of nu is a row of nu^T; its null space holds the
            // quantities the reactions conserve.
            (c, null_space(nu, s))
        }
    };
    for v in &mut c {
        *v += ln_p;
    }
    minimize_gibbs(&c, &invariants, inlet, |x| match fugacity {
        FugacityBasis::IdealGas => Ok(vec![0.0; x.len()]),
        FugacityBasis::Package => pkg
            .ln_fugacity_coefficients(x, t, p, RootSelection::MinimumGibbs)
            .map_err(|e| format!("equilibrium fugacity coefficients: {e}")),
    })
}

/// Reactor context shared by every calculation.
pub(crate) struct Feed<'a> {
    pub inlet: &'a MaterialStream,
    pub state: &'a EquilibriumState,
    pub flows: Vec<f64>,
    pub pressure: f64,
    pub pkg: &'a dyn FlowsheetPackage,
    pub chemistry: Chemistry,
}

impl Feed<'_> {
    fn terms(&self, outlet: &[f64], extents: Vec<f64>) -> ReactionTerms {
        let generation: Vec<f64> = outlet.iter().zip(&self.flows).map(|(o, i)| o - i).collect();
        ReactionTerms {
            standard_reaction_enthalpy_w: self.chemistry.formation_flow_w(&generation),
            generation_mol_s: generation,
            extents_mol_s: extents,
            element_balanced: true,
        }
    }

    /// Energy residual of outlet flows at `t`: `H_out + formation - H_in - Q`.
    fn energy_residual(&self, flows: &[f64], t: f64, duty: f64) -> Result<f64, String> {
        let n: f64 = flows.iter().sum();
        let z: Vec<f64> = flows.iter().map(|v| v / n).collect();
        let out = self
            .pkg
            .tp(&z, t, self.pressure)
            .map_err(|e| format!("reactor outlet: {e}"))?;
        let generation: Vec<f64> = flows.iter().zip(&self.flows).map(|(o, i)| o - i).collect();
        Ok(
            n * out.enthalpy() + self.chemistry.formation_flow_w(&generation)
                - self.inlet.molar_flow_mol_s * self.state.enthalpy()
                - duty,
        )
    }

    /// Outlet for a temperature-dependent outlet composition. A specified
    /// outlet temperature is used directly; adiabatic and duty specifications
    /// solve for the temperature that closes the energy balance.
    fn thermal(
        &self,
        thermal: SeparatorThermal,
        mut flows_at: impl FnMut(f64) -> Result<(Vec<f64>, Vec<f64>), String>,
        mut warnings: Vec<String>,
    ) -> Result<UnitOutput, String> {
        let (flows, extents) = match thermal {
            SeparatorThermal::OutletTemperature(Some(t)) => flows_at(t)?,
            SeparatorThermal::Adiabatic | SeparatorThermal::Duty(Some(_)) => {
                let q = match thermal {
                    SeparatorThermal::Duty(Some(q)) => q,
                    _ => 0.0,
                };
                let [lo, hi] = self.pkg.temperature_interval();
                let t = bracketed_root(
                    |t| {
                        let (f, _) = flows_at(t)?;
                        self.energy_residual(&f, t, q)
                    },
                    self.state.temperature_k,
                    lo,
                    hi,
                    1e-9,
                )
                .map_err(|e| format!("reactor outlet temperature: {e}"))?;
                flows_at(t)?
            }
            _ => return Err("missing reactor thermal specification".into()),
        };
        let terms = self.terms(&flows, extents);
        let out = react(
            self.inlet,
            flows,
            terms,
            std::mem::take(&mut warnings),
            thermal,
            self.pressure,
            self.pkg,
        )?;
        Ok(out)
    }
}

fn single_phase_warning(out: &mut UnitOutput, fugacity: FugacityBasis) {
    let split = out.outlets[0]
        .state
        .as_ref()
        .is_some_and(|s| s.phases.iter().filter(|p| p.phase_fraction > 0.0).count() > 1);
    if split {
        out.warnings.push(format!(
            "The outlet splits into two phases; chemical equilibrium was computed for a single {} fluid phase.",
            if fugacity == FugacityBasis::IdealGas { "ideal-gas" } else { "homogeneous" }
        ));
    }
}

pub(crate) fn equilibrium_reactor(
    feed: &Feed,
    reactions: &[EquilibriumReaction],
    options: EquilibriumOptions,
    thermal: SeparatorThermal,
) -> Result<UnitOutput, String> {
    let nu: Vec<Vec<f64>> = reactions
        .iter()
        .enumerate()
        .map(|(k, r)| feed.chemistry.stoichiometry(k + 1, &r.coefficients))
        .collect::<Result<_, _>>()?;
    if rank(&nu) < nu.len() {
        return Err(
            "equilibrium reactions are linearly dependent; remove the redundant ones".into(),
        );
    }
    let constants: Vec<EquilibriumConstant> = reactions.iter().map(|r| r.constant).collect();
    let constraint = Constraint::Reactions {
        nu: &nu,
        constants: &constants,
    };
    let mut out = feed.thermal(
        thermal,
        |t| {
            let flows = equilibrium_flows(
                &constraint,
                &feed.chemistry,
                feed.pkg,
                &feed.flows,
                t + options.temperature_approach_k,
                feed.pressure,
                options.fugacity,
            )?;
            let generation: Vec<f64> = flows.iter().zip(&feed.flows).map(|(o, i)| o - i).collect();
            Ok((flows, extents_of(&nu, &generation)))
        },
        Vec::new(),
    )?;
    single_phase_warning(&mut out, options.fugacity);
    Ok(out)
}

pub(crate) fn gibbs_reactor(
    feed: &Feed,
    inerts: &[String],
    options: EquilibriumOptions,
    thermal: SeparatorThermal,
) -> Result<UnitOutput, String> {
    let mut inert = vec![false; feed.flows.len()];
    for id in inerts {
        inert[feed.chemistry.index(id)?] = true;
    }
    let constraint = Constraint::Elements { inert: &inert };
    let mut out = feed.thermal(
        thermal,
        |t| {
            let flows = equilibrium_flows(
                &constraint,
                &feed.chemistry,
                feed.pkg,
                &feed.flows,
                t + options.temperature_approach_k,
                feed.pressure,
                options.fugacity,
            )?;
            Ok((flows, Vec::new()))
        },
        Vec::new(),
    )?;
    single_phase_warning(&mut out, options.fugacity);
    Ok(out)
}

/// Dense stoichiometry (columns) of kinetic reactions, element-checked.
fn kinetic_stoichiometry(
    chemistry: &Chemistry,
    kinetics: &Kinetics,
) -> Result<Vec<Vec<f64>>, String> {
    kinetics
        .reactions
        .iter()
        .enumerate()
        .map(|(k, r)| chemistry.stoichiometry(k + 1, &r.coefficients))
        .collect()
}

/// Basis quantity of every component in the reacting phase.
fn basis_values(
    basis: ConcentrationBasis,
    x: &[f64],
    molar_volume: f64,
    ln_phi: &[f64],
    pressure: f64,
    molar_mass: &[f64],
) -> Vec<f64> {
    let mean_mass: f64 = dot(x, molar_mass);
    x.iter()
        .enumerate()
        .map(|(i, &xi)| match basis {
            ConcentrationBasis::MolarConcentration => xi / molar_volume,
            ConcentrationBasis::MassConcentration => xi * molar_mass[i] / molar_volume,
            ConcentrationBasis::MoleFraction => xi,
            ConcentrationBasis::MassFraction => xi * molar_mass[i] / mean_mass,
            ConcentrationBasis::PartialPressure => xi * pressure,
            ConcentrationBasis::Fugacity => xi * ln_phi[i].exp() * pressure,
        })
        .collect()
}

/// `prod_i q_i^a_i` over the listed orders.
fn product(chemistry: &Chemistry, q: &[f64], orders: &[(String, f64)]) -> Result<f64, String> {
    let mut value = 1.0;
    for (id, a) in orders {
        let qi = q[chemistry.index(id)?];
        if *a == 0.0 {
            continue;
        }
        if qi <= 0.0 {
            if *a > 0.0 {
                return Ok(0.0);
            }
            return Err(format!("negative order of {id} with zero {id} present"));
        }
        value *= qi.powf(*a);
    }
    Ok(value)
}

/// Rate of each reaction per m³ of reactor at the given equilibrium state,
/// mol/(m³ s). The reacting phase's share of the volume (no-slip, from the
/// phase molar volumes) scales volume-based rates; catalyst-based rates are
/// multiplied by the catalyst bulk density.
pub(crate) fn rates(
    chemistry: &Chemistry,
    kinetics: &Kinetics,
    state: &EquilibriumState,
) -> Result<Vec<f64>, String> {
    let t = state.temperature_k;
    let total_volume: f64 = state
        .phases
        .iter()
        .map(|p| p.phase_fraction * p.state.molar_volume)
        .sum();
    kinetics
        .reactions
        .iter()
        .map(|r| {
            let phase = state.phases.iter().find(|p| {
                p.phase_fraction > 0.0
                    && matches!(
                        (r.phase, p.phase),
                        (_, Phase::Supercritical)
                            | (ReactingPhase::Vapor, Phase::Vapor)
                            | (ReactingPhase::Liquid, Phase::Liquid)
                    )
            });
            let Some(phase) = phase else {
                return Ok(0.0);
            };
            let q = basis_values(
                r.basis,
                &phase.composition,
                phase.state.molar_volume,
                &phase.state.ln_fugacity_coefficients,
                state.pressure_pa,
                chemistry.molar_masses(),
            );
            let intrinsic = rate_value(chemistry, r, &q, t)?;
            let scale = match r.rate_basis {
                RateBasis::ReactorVolume => {
                    phase.phase_fraction * phase.state.molar_volume / total_volume
                }
                RateBasis::CatalystMass => kinetics
                    .catalyst_bulk_density_kg_m3
                    .ok_or("catalyst bulk density is missing")?,
            };
            Ok(intrinsic * scale)
        })
        .collect()
}

fn rate_value(
    chemistry: &Chemistry,
    r: &KineticReaction,
    q: &[f64],
    t: f64,
) -> Result<f64, String> {
    let k = |a: &crate::Arrhenius| {
        a.value(t)
            .ok_or_else(|| "pre-exponential factor is missing".to_string())
    };
    match &r.rate {
        RateLaw::PowerLaw { forward, reverse } => {
            let mut v = k(&forward.constant)? * product(chemistry, q, &forward.orders)?;
            if let Some(rev) = reverse {
                v -= k(&rev.constant)? * product(chemistry, q, &rev.orders)?;
            }
            Ok(v)
        }
        RateLaw::Lhhw {
            kinetic,
            driving_forward,
            driving_reverse,
            adsorption,
            exponent,
        } => {
            let term = |l: &crate::LhhwTerm| -> Result<f64, String> {
                Ok(l.constant.ln_value(t).exp() * product(chemistry, q, &l.orders)?)
            };
            let mut driving = term(driving_forward)?;
            if let Some(rev) = driving_reverse {
                driving -= term(rev)?;
            }
            let denominator = if adsorption.is_empty() {
                1.0
            } else {
                adsorption.iter().map(term).sum::<Result<f64, String>>()?
            };
            if denominator.is_nan() || denominator <= 0.0 {
                return Err("LHHW adsorption term is not positive".into());
            }
            Ok(k(kinetic)? * driving / denominator.powf(*exponent))
        }
    }
}

fn fractions(flows: &[f64]) -> Result<Vec<f64>, String> {
    let n: f64 = flows.iter().sum();
    if !(n.is_finite() && n > 0.0) {
        return Err("reactor contents have no positive total amount".into());
    }
    Ok(flows.iter().map(|v| v.max(0.0) / n).collect())
}

/// CSTR outlet flows and extents at temperature `t`: `xi = V R(x(xi))`.
fn cstr_at(
    feed: &Feed,
    kinetics: &Kinetics,
    nu: &[Vec<f64>],
    volume: f64,
    t: f64,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let s = feed.flows.len();
    let r = nu.len();
    let scale: f64 = feed.flows.iter().sum();
    let flows_of = |xi: &[f64]| -> Vec<f64> {
        (0..s)
            .map(|i| feed.flows[i] + (0..r).map(|j| nu[j][i] * xi[j]).sum::<f64>())
            .collect()
    };
    let residual = |xi: &[f64], v: f64| -> Result<Vec<f64>, String> {
        let flows = flows_of(xi);
        let state = feed
            .pkg
            .tp(&fractions(&flows)?, t, feed.pressure)
            .map_err(|e| format!("CSTR contents: {e}"))?;
        let rate = rates(&feed.chemistry, kinetics, &state)?;
        Ok((0..r).map(|j| xi[j] - v * rate[j]).collect())
    };
    // Largest step keeping every flow nonnegative.
    let feasible = |xi: &[f64], step: &[f64]| -> f64 {
        let n = flows_of(xi);
        let mut alpha: f64 = 1.0;
        for i in 0..s {
            let dn: f64 = (0..r).map(|j| nu[j][i] * step[j]).sum();
            if dn < 0.0 {
                alpha = alpha.min(0.99 * n[i].max(0.0) / -dn);
            }
        }
        alpha
    };
    let newton = |xi0: Vec<f64>, v: f64| -> Result<Vec<f64>, String> {
        let mut xi = xi0;
        let mut f = residual(&xi, v)?;
        for _ in 0..100 {
            let norm = f.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            if norm <= 1e-12 * scale {
                return Ok(xi);
            }
            let mut jac = vec![vec![0.0; r]; r];
            for k in 0..r {
                let h = 1e-7 * xi[k].abs().max(1e-6 * scale);
                let mut shifted = xi.clone();
                shifted[k] += h;
                if flows_of(&shifted).iter().any(|n| *n < 0.0) {
                    shifted[k] = xi[k] - h;
                }
                let fs = residual(&shifted, v)?;
                for j in 0..r {
                    jac[j][k] = (fs[j] - f[j]) / (shifted[k] - xi[k]);
                }
            }
            let step = solve_dense(jac, f.iter().map(|x| -x).collect())
                .ok_or("singular CSTR Jacobian (dependent reactions?)")?;
            let mut alpha = feasible(&xi, &step);
            let mut improved = false;
            for _ in 0..40 {
                let trial: Vec<f64> = xi.iter().zip(&step).map(|(a, d)| a + alpha * d).collect();
                if let Ok(ft) = residual(&trial, v) {
                    let nt = ft.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
                    if nt < norm {
                        xi = trial;
                        f = ft;
                        improved = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if !improved {
                break;
            }
        }
        Err("CSTR material balance did not converge".into())
    };
    // Direct solve from no reaction; otherwise continue in volume from zero,
    // which follows the low-conversion branch deterministically.
    let xi = match newton(vec![0.0; r], volume) {
        Ok(xi) => xi,
        Err(_) => {
            let mut xi = vec![0.0; r];
            for k in 1..=20 {
                xi = newton(xi, volume * k as f64 / 20.0)?;
            }
            xi
        }
    };
    Ok((flows_of(&xi), xi))
}

pub(crate) fn cstr(
    feed: &Feed,
    kinetics: &Kinetics,
    volume: f64,
    thermal: SeparatorThermal,
) -> Result<UnitOutput, String> {
    let nu = kinetic_stoichiometry(&feed.chemistry, kinetics)?;
    feed.thermal(
        thermal,
        |t| cstr_at(feed, kinetics, &nu, volume, t),
        Vec::new(),
    )
}

/// State of reactor contents from component amounts and, for nonisothermal
/// operation, the formation-basis enthalpy `h_total` (per the same amounts).
fn contents(
    feed: &Feed,
    amounts: &[f64],
    temperature: Option<f64>,
    h_total: f64,
    pressure: f64,
) -> Result<EquilibriumState, String> {
    let n: f64 = amounts.iter().sum();
    let z = fractions(amounts)?;
    match temperature {
        Some(t) => feed.pkg.tp(&z, t, pressure),
        None => {
            // Remove the formation part of the amounts present.
            let h = (h_total - feed.chemistry.formation_flow_w(amounts)) / n;
            feed.pkg.ph(&z, pressure, h)
        }
    }
    .map_err(|e| e.to_string())
}

/// Integration tolerance of the kinetic reactors (relative).
const RTOL: f64 = 1e-9;

pub(crate) fn plug_flow(
    feed: &Feed,
    kinetics: &Kinetics,
    length: f64,
    diameter: f64,
    thermal: PlugFlowThermal,
    pressure_drop: f64,
) -> Result<UnitOutput, String> {
    let nu = kinetic_stoichiometry(&feed.chemistry, kinetics)?;
    let s = feed.flows.len();
    let r = nu.len();
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let p_in = feed.state.pressure_pa;
    let pressure_at = |z: f64| p_in - pressure_drop * z / length;
    let isothermal = match thermal {
        PlugFlowThermal::Isothermal(Some(t)) => Some(t),
        PlugFlowThermal::Isothermal(None) => return Err("missing reactor temperature".into()),
        _ => None,
    };
    let (coolant_t, u) = match thermal {
        PlugFlowThermal::Coolant {
            temperature_k: Some(t),
            u_w_m2_k: Some(u),
        } => (t, u),
        PlugFlowThermal::Coolant { .. } => return Err("missing coolant specification".into()),
        _ => (0.0, 0.0),
    };
    let formation_in = feed.chemistry.formation_flow_w(&feed.flows);
    let h_in = feed.inlet.molar_flow_mol_s * feed.state.enthalpy() + formation_in;
    // State: component flows, extents, total formation-basis enthalpy flow.
    let mut y0 = feed.flows.clone();
    y0.extend(vec![0.0; r]);
    y0.push(h_in);
    let scale: f64 = feed.flows.iter().sum();
    let thermal_scale = feed.inlet.thermal_scale_w().max(1.0);
    let mut atol = vec![1e-12 * scale; s + r];
    atol.push(1e-10 * (thermal_scale + h_in.abs()));
    let rhs = |z: f64, y: &[f64]| -> Result<Vec<f64>, String> {
        let state = contents(feed, &y[..s], isothermal, y[s + r], pressure_at(z))?;
        let rate = rates(&feed.chemistry, kinetics, &state)?;
        let mut d: Vec<f64> = (0..s)
            .map(|i| area * (0..r).map(|j| nu[j][i] * rate[j]).sum::<f64>())
            .collect();
        d.extend(rate.iter().map(|v| area * v));
        d.push(match thermal {
            PlugFlowThermal::Coolant { .. } => {
                u * std::f64::consts::PI * diameter * (coolant_t - state.temperature_k)
            }
            _ => 0.0,
        });
        Ok(d)
    };
    let points = integrate(rhs, 0.0, length, y0, RTOL, &atol, 20_000)
        .map_err(|e| format!("plug-flow integration: {e}"))?;
    let (_, end) = points.last().expect("integration start point");
    let flows: Vec<f64> = end[..s].iter().map(|v| v.max(0.0)).collect();
    let extents = end[s..s + r].to_vec();
    let spec = match thermal {
        PlugFlowThermal::Isothermal(t) => SeparatorThermal::OutletTemperature(t),
        PlugFlowThermal::Adiabatic => SeparatorThermal::Adiabatic,
        PlugFlowThermal::Coolant { .. } => SeparatorThermal::Duty(Some(end[s + r] - h_in)),
    };
    let feed_out = Feed {
        pressure: pressure_at(length),
        ..feed.clone_shallow()
    };
    let terms = feed_out.terms(&flows, extents);
    react(
        feed.inlet,
        flows,
        terms,
        Vec::new(),
        spec,
        feed_out.pressure,
        feed.pkg,
    )
}

pub(crate) fn batch(
    feed: &Feed,
    kinetics: &Kinetics,
    reaction_time: f64,
    thermal: BatchThermal,
) -> Result<UnitOutput, String> {
    let nu = kinetic_stoichiometry(&feed.chemistry, kinetics)?;
    let s = feed.flows.len();
    let r = nu.len();
    let p = feed.state.pressure_pa;
    let isothermal = match thermal {
        BatchThermal::Isothermal(Some(t)) => Some(t),
        BatchThermal::Isothermal(None) => return Err("missing batch temperature".into()),
        BatchThermal::Adiabatic => None,
    };
    // Per mole charged: amounts, extents and formation-basis enthalpy.
    let z = feed.state.mole_fractions.clone();
    let h0 = feed.state.enthalpy() + feed.chemistry.formation_flow_w(&z);
    let mut y0 = z;
    y0.extend(vec![0.0; r]);
    y0.push(h0);
    let mut atol = vec![1e-12; s + r];
    atol.push(1e-10 * (molar_scale(feed.state) + h0.abs()));
    let rhs = |_t: f64, y: &[f64]| -> Result<Vec<f64>, String> {
        let state = contents(feed, &y[..s], isothermal, y[s + r], p)?;
        let rate = rates(&feed.chemistry, kinetics, &state)?;
        let n: f64 = y[..s].iter().sum();
        let volume: f64 = n * state
            .phases
            .iter()
            .map(|ph| ph.phase_fraction * ph.state.molar_volume)
            .sum::<f64>();
        let mut d: Vec<f64> = (0..s)
            .map(|i| volume * (0..r).map(|j| nu[j][i] * rate[j]).sum::<f64>())
            .collect();
        d.extend(rate.iter().map(|v| volume * v));
        d.push(0.0);
        Ok(d)
    };
    let points = integrate(rhs, 0.0, reaction_time, y0, RTOL, &atol, 20_000)
        .map_err(|e| format!("batch integration: {e}"))?;
    let (_, end) = points.last().expect("integration start point");
    let f = feed.inlet.molar_flow_mol_s;
    let flows: Vec<f64> = end[..s].iter().map(|v| f * v.max(0.0)).collect();
    let extents: Vec<f64> = end[s..s + r].iter().map(|v| f * v).collect();
    let terms = feed.terms(&flows, extents);
    let spec = match thermal {
        BatchThermal::Isothermal(t) => SeparatorThermal::OutletTemperature(t),
        BatchThermal::Adiabatic => SeparatorThermal::Adiabatic,
    };
    react(feed.inlet, flows, terms, Vec::new(), spec, p, feed.pkg)
}

fn molar_scale(s: &EquilibriumState) -> f64 {
    let cp: f64 = s.phases.iter().map(|p| p.phase_fraction * p.state.cp).sum();
    (cp * s.temperature_k).abs().max(1.0)
}

impl<'a> Feed<'a> {
    fn clone_shallow(&self) -> Feed<'a> {
        Feed {
            inlet: self.inlet,
            state: self.state,
            flows: self.flows.clone(),
            pressure: self.pressure,
            pkg: self.pkg,
            chemistry: self.chemistry.clone(),
        }
    }
}
