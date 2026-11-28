//! Rate-based (mass-transfer) packed absorber or stripper (#241).
//!
//! The packed height is divided into equal mixed-flow segments, numbered
//! from the top. Each segment has bulk vapor and bulk liquid at their
//! outlet states, an interface at phase equilibrium and component fluxes
//! from vapor to liquid. Two-film model with effective diffusivities:
//!
//! \[ N_i = c_V k_{V,i} A (y_i - y_{I,i}) + \bar y_i N_t, \qquad
//!    N_i = c_L k_{L,i} A (x_{I,i} - x_i) + \bar x_i N_t, \]
//!
//! for every component but one reference per phase (the most abundant in
//! that phase's feed), with \(N_t = \sum N_i\), film-average fractions
//! \(\bar y = (y + y_I)/2\), \(\bar x = (x + x_I)/2\) and interfacial area
//! \(A = f\,a_w A_c \Delta z\). Energy crosses the films as
//! \(E_V = h_V A (T_V - T_I) + \sum N_i \bar H^V_i\) and
//! \(E_L = h_L A (T_I - T_L) + \sum N_i \bar H^L_i\) with partial molar
//! enthalpies at the interface; \(E_V = E_L\) there. Heat-transfer
//! coefficients follow the Chilton–Colburn analogy \(j_H = j_D\) in both
//! phases, \(h = \bar k\,c\,c_p\,\mathrm{Le}^{2/3}\) with
//! \(\mathrm{Le} = \lambda/(c\,c_p \bar D)\) (Chilton and Colburn, Ind.
//! Eng. Chem. 26 (1934) 1183; Geankoplis, *Transport Processes*, 4th ed.,
//! p. 475). Component and
//! energy balances of each phase close every segment. With very large
//! transfer rates every segment becomes an ideal stage, so the seed is the
//! equilibrium-stage absorber with one stage per segment.
//!
//! The unknowns of all segments are solved simultaneously by Newton's method
//! with a finite-difference Jacobian assembled segment by segment.

use property_packages::{EquilibriumPhase, PhaseProperties};
use thermo_models::{DiffusivityModel, RootSelection};
use thermo_types::Phase;

use super::mesh::{self, TEMPERATURE_SCALE};
use super::packing::{MAX_NEGLECTED_TRANSPORT_FRACTION, RateBasedSpec};
use super::packing_data::{self as pd, FilmProperties, Packing, PackingData};
use super::stages::{Stage, StageCascade, StageFeed};
use crate::MaterialStream;
use crate::numerics::solve_dense;
use crate::units::CalcContext;

const ENTHALPY_SCALE: f64 = 10_000.0;
const MAX_ITERATIONS: usize = 60;
const RESIDUAL_TOLERANCE: f64 = 1e-10;
const MAX_PRESSURE_ITERATIONS: usize = 30;
const PRESSURE_TOLERANCE_PA: f64 = 1e-6;
/// Relative mole-number step of the partial molar enthalpies.
const PARTIAL_MOLAR_STEP: f64 = 1e-6;
/// Area factors of the continuation from near-equilibrium rates, as
/// multiples of the specified factor (used only when the direct solve from
/// the equilibrium-stage seed fails).
const CONTINUATION: [f64; 5] = [1e4, 1e3, 1e2, 1e1, 1.0];

/// Profiles and hydraulics of a solved rate-based packed section, segments
/// from the top (liquid inlet) to the bottom (gas inlet).
#[derive(Debug, Clone, PartialEq)]
pub struct RateBasedProfile {
    pub packing: Packing,
    pub diameter_m: f64,
    pub segment_height_m: f64,
    pub vapor_temperatures_k: Vec<f64>,
    pub liquid_temperatures_k: Vec<f64>,
    pub interface_temperatures_k: Vec<f64>,
    /// Vapor leaving each segment upward, mol/s.
    pub vapor_flows_mol_s: Vec<f64>,
    /// Liquid leaving each segment downward, mol/s.
    pub liquid_flows_mol_s: Vec<f64>,
    /// Bulk vapor mole fractions, package order.
    pub vapor_mole_fractions: Vec<Vec<f64>>,
    /// Bulk liquid mole fractions, package order.
    pub liquid_mole_fractions: Vec<Vec<f64>>,
    /// Component transfer from vapor to liquid in each segment, mol/s.
    pub fluxes_mol_s: Vec<Vec<f64>>,
    /// Effective (wetted) interfacial area including the area factor, m²/m³.
    pub interfacial_area_m2_m3: Vec<f64>,
    pub flooding_fractions: Vec<f64>,
    pub pressure_drops_pa_m: Vec<f64>,
    /// Whether the boundary pressures were used in the segment equations.
    pub pressure_drop_coupled: bool,
    /// Estimated absolute pressure at each segment boundary, from the top to
    /// the bottom. The top is anchored at the specified column pressure.
    pub estimated_boundary_pressures_pa: Vec<f64>,
    /// Total bed pressure drop, Pa.
    pub pressure_drop_pa: f64,
    pub max_flooding_fraction: f64,
    /// Largest fraction of a component neglected in a phase's transport
    /// properties, 0 when none is neglected.
    pub max_neglected_transport_fraction: f64,
    pub iterations: usize,
}

/// Inputs shared by every segment.
pub(super) struct RateColumn<'a, 'p> {
    pub ctx: &'a CalcContext<'p>,
    pub spec: &'a RateBasedSpec,
    pub pressure: f64,
    pub temperature_range: [f64; 2],
    pub height_m: f64,
    pub gas: StageFeed,
    pub liquid: StageFeed,
    pub active: Vec<usize>,
    pub total_flow: f64,
    diffusivity: DiffusivityModel,
    gas_neglected: Vec<usize>,
    liquid_neglected: Vec<usize>,
    boundary_pressures_pa: Vec<f64>,
}

/// Inputs of [`RateColumn::new`].
pub(super) struct RateInputs<'a, 'p> {
    pub ctx: &'a CalcContext<'p>,
    pub spec: &'a RateBasedSpec,
    pub pressure: f64,
    pub temperature_range: [f64; 2],
    pub height_m: f64,
    pub gas: StageFeed,
    pub liquid: StageFeed,
    pub active: Vec<usize>,
    pub total_flow: f64,
}

fn indices(ctx: &CalcContext<'_>, ids: &[String]) -> Result<Vec<usize>, String> {
    let names = ctx.component_ids();
    ids.iter()
        .map(|id| {
            names.iter().position(|n| n == id).ok_or_else(|| {
                format!("rate-based transport basis: {id} is not a package component")
            })
        })
        .collect()
}

impl<'a, 'p> RateColumn<'a, 'p> {
    pub(super) fn new(inputs: RateInputs<'a, 'p>) -> Result<Self, String> {
        let ctx = inputs.ctx;
        let diffusivity = DiffusivityModel::new(&ctx.component_ids())
            .map_err(|e| format!("rate-based diffusivities: {e}"))?;
        Ok(Self {
            gas_neglected: indices(ctx, &inputs.spec.gas_transport_neglected)?,
            liquid_neglected: indices(ctx, &inputs.spec.liquid_transport_neglected)?,
            boundary_pressures_pa: vec![inputs.pressure; usize::from(inputs.spec.segments) + 1],
            diffusivity,
            ctx,
            spec: inputs.spec,
            pressure: inputs.pressure,
            temperature_range: inputs.temperature_range,
            height_m: inputs.height_m,
            gas: inputs.gas,
            liquid: inputs.liquid,
            active: inputs.active,
            total_flow: inputs.total_flow,
        })
    }
}

#[derive(Debug, Clone)]
struct Segment {
    v: Vec<f64>,
    l: Vec<f64>,
    t_v: f64,
    t_l: f64,
    t_i: f64,
    x_i: Vec<f64>,
    y_i: Vec<f64>,
    flux: Vec<f64>,
}

impl Segment {
    fn vapor(&self) -> f64 {
        self.v.iter().sum()
    }
    fn liquid(&self) -> f64 {
        self.l.iter().sum()
    }
}

#[derive(Debug, Clone)]
struct Props {
    y: Vec<f64>,
    x: Vec<f64>,
    h_v: f64,
    h_l: f64,
    ln_k: Vec<f64>,
    hbar_v: Vec<f64>,
    hbar_l: Vec<f64>,
    c_v: f64,
    c_l: f64,
    /// k_V,i·A (m³/s) per component; zero for the reference and inactive.
    ka_v: Vec<f64>,
    ka_l: Vec<f64>,
    /// h·A, W/K.
    ha_v: f64,
    ha_l: f64,
    area_m2_m3: f64,
    gas_film: FilmProperties,
    liquid_film: FilmProperties,
    liquid_velocity_m_s: f64,
    surface_tension_n_m: f64,
    neglected: f64,
}

/// Decoded segments, their properties and the scaled residual.
type Evaluation = (Vec<Segment>, Vec<Props>, Vec<f64>);
/// Outlets, profile and warnings of a solved section.
pub(super) type Solved = (Vec<MaterialStream>, RateBasedProfile, Vec<String>);

struct Unknowns {
    c: usize,
}

impl Unknowns {
    fn size(&self) -> usize {
        5 * self.c + 1
    }
}

fn phase_of(
    ctx: &CalcContext<'_>,
    composition: &[f64],
    t: f64,
    p: f64,
    phase: Phase,
) -> Result<EquilibriumPhase, String> {
    let root = match phase {
        Phase::Liquid => RootSelection::SmallestVolume,
        _ => RootSelection::LargestVolume,
    };
    Ok(EquilibriumPhase {
        phase,
        phase_fraction: 1.0,
        composition: composition.to_vec(),
        state: ctx.homogeneous(composition, t, p, root)?,
    })
}

/// `composition` without the neglected components, renormalised, and the
/// largest neglected fraction.
fn transport_basis(
    composition: &[f64],
    neglected: &[usize],
    phase: &str,
    segment: usize,
    ids: &[&str],
) -> Result<(Vec<f64>, f64), String> {
    let mut basis = composition.to_vec();
    let mut largest: f64 = 0.0;
    for &i in neglected {
        if composition[i] > MAX_NEGLECTED_TRANSPORT_FRACTION {
            return Err(format!(
                "rate-based segment {}: {} is {:.4} of the {phase}, more than the {} allowed for a component neglected in {phase} transport properties",
                segment + 1,
                ids[i],
                composition[i],
                MAX_NEGLECTED_TRANSPORT_FRACTION
            ));
        }
        largest = largest.max(composition[i]);
        basis[i] = 0.0;
    }
    let total: f64 = basis.iter().sum();
    if total <= 0.0 {
        return Err(format!(
            "rate-based segment {}: every {phase} component is neglected in transport properties",
            segment + 1
        ));
    }
    basis.iter_mut().for_each(|v| *v /= total);
    Ok((basis, largest))
}

fn partial_molar_enthalpies(
    ctx: &CalcContext<'_>,
    state: &PhaseProperties,
    composition: &[f64],
    active: &[usize],
    root: RootSelection,
) -> Result<Vec<f64>, String> {
    let mut out = vec![0.0; composition.len()];
    for &i in active {
        let mut shifted = composition.to_vec();
        shifted[i] += PARTIAL_MOLAR_STEP;
        let total = 1.0 + PARTIAL_MOLAR_STEP;
        shifted.iter_mut().for_each(|v| *v /= total);
        let h = ctx
            .homogeneous(&shifted, state.temperature_k, state.pressure_pa, root)?
            .enthalpy;
        out[i] = (total * h - state.enthalpy) / PARTIAL_MOLAR_STEP;
    }
    Ok(out)
}

impl RateColumn<'_, '_> {
    fn segments(&self) -> usize {
        usize::from(self.spec.segments)
    }

    fn unknowns(&self) -> Unknowns {
        Unknowns {
            c: self.active.len(),
        }
    }

    fn cross_section_m2(&self) -> f64 {
        let d = self.spec.diameter_m.expect("checked diameter");
        std::f64::consts::PI * d * d / 4.0
    }

    fn packing(&self) -> PackingData {
        self.spec.packing.data()
    }

    fn reference(&self, z: &[f64]) -> usize {
        *self
            .active
            .iter()
            .max_by(|&&a, &&b| z[a].total_cmp(&z[b]))
            .expect("active components")
    }

    fn decode(&self, block: &[f64]) -> Segment {
        let c = self.active.len();
        let n = self.ctx.component_count();
        let f = self.total_flow;
        let mut v = vec![0.0; n];
        let mut l = vec![0.0; n];
        let mut flux = vec![0.0; n];
        for (k, &i) in self.active.iter().enumerate() {
            v[i] = f * block[k].exp();
            l[i] = f * block[c + k].exp();
            flux[i] = f * block[4 * c + 1 + k];
        }
        Segment {
            v,
            l,
            t_v: block[2 * c] * TEMPERATURE_SCALE,
            t_l: block[2 * c + 1] * TEMPERATURE_SCALE,
            t_i: block[2 * c + 2] * TEMPERATURE_SCALE,
            x_i: mesh::fractions(&block[2 * c + 3..3 * c + 2], &self.active, n),
            y_i: mesh::fractions(&block[3 * c + 2..4 * c + 1], &self.active, n),
            flux,
        }
    }

    fn encode(&self, s: &Segment) -> Vec<f64> {
        let f = self.total_flow;
        let mut out: Vec<f64> = self.active.iter().map(|&i| (s.v[i] / f).ln()).collect();
        out.extend(self.active.iter().map(|&i| (s.l[i] / f).ln()));
        out.extend([s.t_v, s.t_l, s.t_i].map(|t| t / TEMPERATURE_SCALE));
        out.extend(mesh::log_ratios(&s.x_i, &self.active));
        out.extend(mesh::log_ratios(&s.y_i, &self.active));
        out.extend(self.active.iter().map(|&i| s.flux[i] / f));
        out
    }

    fn props(&self, j: usize, s: &Segment, factor: f64) -> Result<Props, String> {
        let ctx = self.ctx;
        let p_v = self.boundary_pressures_pa[j];
        let p_l = self.boundary_pressures_pa[j + 1];
        let p_i = p_v.midpoint(p_l);
        let [lo, hi] = self.temperature_range;
        for (name, t) in [("vapor", s.t_v), ("liquid", s.t_l), ("interface", s.t_i)] {
            if !(lo..=hi).contains(&t) {
                return Err(format!(
                    "rate-based segment {} {name} temperature is outside {lo}..={hi} K",
                    j + 1
                ));
            }
        }
        let (vapor, liquid) = (s.vapor(), s.liquid());
        if !(vapor > 0.0 && liquid > 0.0 && vapor.is_finite() && liquid.is_finite()) {
            return Err(format!(
                "rate-based segment {} has no liquid or vapor flow",
                j + 1
            ));
        }
        let y: Vec<f64> = s.v.iter().map(|v| v / vapor).collect();
        let x: Vec<f64> = s.l.iter().map(|l| l / liquid).collect();
        let bulk_v = ctx.homogeneous(&y, s.t_v, p_v, RootSelection::LargestVolume)?;
        let bulk_l = ctx.homogeneous(&x, s.t_l, p_l, RootSelection::SmallestVolume)?;
        let int_v = ctx.homogeneous(&s.y_i, s.t_i, p_i, RootSelection::LargestVolume)?;
        let int_l = ctx.homogeneous(&s.x_i, s.t_i, p_i, RootSelection::SmallestVolume)?;
        if int_v.molar_volume <= 1.05 * int_l.molar_volume
            || bulk_v.molar_volume <= 1.05 * bulk_l.molar_volume
        {
            return Err(format!(
                "rate-based segment {} lacks distinct liquid and vapor roots",
                j + 1
            ));
        }
        let ln_k = int_l
            .ln_fugacity_coefficients
            .iter()
            .zip(&int_v.ln_fugacity_coefficients)
            .map(|(l, v)| l - v)
            .collect();
        let hbar_v = partial_molar_enthalpies(
            ctx,
            &int_v,
            &s.y_i,
            &self.active,
            RootSelection::LargestVolume,
        )?;
        let hbar_l = partial_molar_enthalpies(
            ctx,
            &int_l,
            &s.x_i,
            &self.active,
            RootSelection::SmallestVolume,
        )?;

        // Transport properties on the stated basis.
        let ids = ctx.component_ids();
        let (gas_basis, gas_dropped) = transport_basis(&y, &self.gas_neglected, "gas", j, &ids)?;
        let (liquid_basis, liquid_dropped) =
            transport_basis(&x, &self.liquid_neglected, "liquid", j, &ids)?;
        let context = |e: String| format!("rate-based segment {}: {e}", j + 1);
        let gas_phase = phase_of(ctx, &gas_basis, s.t_v, p_v, Phase::Vapor).map_err(context)?;
        let liquid_phase =
            phase_of(ctx, &liquid_basis, s.t_l, p_l, Phase::Liquid).map_err(context)?;
        let gas_t = ctx.transport(&gas_phase).map_err(context)?;
        let liquid_t = ctx.transport(&liquid_phase).map_err(context)?;
        let unavailable =
            |what: &str, e: thermo_types::CalculationError| context(format!("{what}: {e}"));
        let mu_v = gas_t
            .viscosity
            .map_err(|e| unavailable("gas viscosity", e))?;
        let lambda_v = gas_t
            .thermal_conductivity
            .map_err(|e| unavailable("gas thermal conductivity", e))?;
        let mu_l = liquid_t
            .viscosity
            .map_err(|e| unavailable("liquid viscosity", e))?;
        let lambda_l = liquid_t
            .thermal_conductivity
            .map_err(|e| unavailable("liquid thermal conductivity", e))?;
        let sigma = liquid_t
            .surface_tension
            .map_err(|e| unavailable("liquid surface tension", e))?;
        let density = ctx
            .package()
            .engineering_liquid_density(&x, s.t_l, p_l)
            .map_err(|e| unavailable("engineering liquid density", e))?;

        let d_v = self
            .diffusivity
            .gas_effective(&y, s.t_v, p_v)
            .map_err(|e| context(e.to_string()))?;
        let d_l = self
            .diffusivity
            .liquid_in_solvent(self.reference(&self.liquid.z), s.t_l, mu_l)
            .map_err(|e| context(e.to_string()))?;

        // Hydraulic loads.
        let area = self.cross_section_m2();
        let m_v = bulk_v.mass_density * bulk_v.molar_volume;
        let m_l = density.mass_density * density.molar_volume;
        let gas_film = FilmProperties {
            mass_flux_kg_m2_s: vapor * m_v / area,
            density_kg_m3: bulk_v.mass_density,
            viscosity_pa_s: mu_v,
        };
        let liquid_film = FilmProperties {
            mass_flux_kg_m2_s: liquid * m_l / area,
            density_kg_m3: density.mass_density,
            viscosity_pa_s: mu_l,
        };
        let packing = self.packing();
        let sigma_c = self
            .spec
            .critical_surface_tension_n_m
            .expect("checked critical surface tension");
        let wetted = pd::onda_wetted_area(&packing, &liquid_film, sigma, sigma_c);
        let area_m2_m3 = factor * wetted;
        let interface_m2 = area_m2_m3 * area * self.height_m / self.segments() as f64;
        let ref_v = self.reference(&self.gas.z);
        let ref_l = self.reference(&self.liquid.z);
        let n = ids.len();
        let mut ka_v = vec![0.0; n];
        let mut ka_l = vec![0.0; n];
        let (mut dv_mean, mut dl_mean) = (0.0, 0.0);
        let c = self.active.len() as f64;
        for &i in &self.active {
            if i != ref_v {
                ka_v[i] = pd::onda_gas_coefficient(&packing, &gas_film, d_v[i]) * interface_m2;
                dv_mean += d_v[i] / (c - 1.0);
            }
            if i != ref_l {
                ka_l[i] = pd::onda_liquid_coefficient(&packing, &liquid_film, wetted, d_l[i])
                    * interface_m2;
                dl_mean += d_l[i] / (c - 1.0);
            }
        }
        let c_v = 1.0 / bulk_v.molar_volume;
        let c_l = 1.0 / density.molar_volume;
        let k_v = pd::onda_gas_coefficient(&packing, &gas_film, dv_mean);
        let k_l = pd::onda_liquid_coefficient(&packing, &liquid_film, wetted, dl_mean);
        let le_v = lambda_v / (c_v * bulk_v.cp * dv_mean);
        let le_l = lambda_l / (c_l * bulk_l.cp * dl_mean);
        let ha_v = k_v * c_v * bulk_v.cp * le_v.powf(2.0 / 3.0) * interface_m2;
        let ha_l = k_l * c_l * bulk_l.cp * le_l.powf(2.0 / 3.0) * interface_m2;
        let values = [ha_v, ha_l, area_m2_m3, gas_film.mass_flux_kg_m2_s];
        if values.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(context("nonfinite transfer coefficient".into()));
        }
        Ok(Props {
            y,
            x,
            h_v: bulk_v.enthalpy,
            h_l: bulk_l.enthalpy,
            ln_k,
            hbar_v,
            hbar_l,
            c_v,
            c_l,
            ka_v,
            ka_l,
            ha_v,
            ha_l,
            area_m2_m3,
            gas_film,
            liquid_film,
            liquid_velocity_m_s: liquid * density.molar_volume / area,
            surface_tension_n_m: sigma,
            neglected: gas_dropped.max(liquid_dropped),
        })
    }

    /// Residual block of segment `j` from its own state and its neighbours'.
    fn block(
        &self,
        above: Option<(&Segment, &Props)>,
        (s, pr): (&Segment, &Props),
        below: Option<(&Segment, &Props)>,
    ) -> Vec<f64> {
        let f = self.total_flow;
        let scale = f * ENTHALPY_SCALE;
        let (l_in, hl_in): (Vec<f64>, f64) = match above {
            Some((a, ap)) => (a.l.clone(), a.liquid() * ap.h_l),
            None => (
                self.liquid
                    .z
                    .iter()
                    .map(|z| z * self.liquid.molar_flow_mol_s)
                    .collect(),
                self.liquid.molar_flow_mol_s * self.liquid.molar_enthalpy_j_mol,
            ),
        };
        let (v_in, hv_in): (Vec<f64>, f64) = match below {
            Some((b, bp)) => (b.v.clone(), b.vapor() * bp.h_v),
            None => (
                self.gas
                    .z
                    .iter()
                    .map(|z| z * self.gas.molar_flow_mol_s)
                    .collect(),
                self.gas.molar_flow_mol_s * self.gas.molar_enthalpy_j_mol,
            ),
        };
        let mut out = Vec::with_capacity(self.unknowns().size());
        for &i in &self.active {
            out.push((v_in[i] - s.v[i] - s.flux[i]) / f);
        }
        for &i in &self.active {
            out.push((l_in[i] - s.l[i] + s.flux[i]) / f);
        }
        let carried_v: f64 = self.active.iter().map(|&i| s.flux[i] * pr.hbar_v[i]).sum();
        let carried_l: f64 = self.active.iter().map(|&i| s.flux[i] * pr.hbar_l[i]).sum();
        let e_v = pr.ha_v * (s.t_v - s.t_i) + carried_v;
        let e_l = pr.ha_l * (s.t_i - s.t_l) + carried_l;
        // Film terms are scaled by their own conductance as well, so very
        // large transfer rates (the equilibrium limit) stay well scaled.
        let scale_v = scale + pr.ha_v * TEMPERATURE_SCALE;
        let scale_l = scale + pr.ha_l * TEMPERATURE_SCALE;
        out.push((hv_in - s.vapor() * pr.h_v - e_v) / scale_v);
        out.push((hl_in - s.liquid() * pr.h_l + e_l) / scale_l);
        out.push((e_v - e_l) / (scale_v + scale_l - scale));
        for &i in &self.active {
            out.push(s.y_i[i].ln() - s.x_i[i].ln() - pr.ln_k[i]);
        }
        let total: f64 = s.flux.iter().sum();
        let ref_v = self.reference(&self.gas.z);
        let ref_l = self.reference(&self.liquid.z);
        for &i in self.active.iter().filter(|&&i| i != ref_v) {
            let mean = 0.5 * (pr.y[i] + s.y_i[i]);
            let conductance = pr.c_v * pr.ka_v[i];
            out.push(
                (s.flux[i] - conductance * (pr.y[i] - s.y_i[i]) - mean * total) / (f + conductance),
            );
        }
        for &i in self.active.iter().filter(|&&i| i != ref_l) {
            let mean = 0.5 * (pr.x[i] + s.x_i[i]);
            let conductance = pr.c_l * pr.ka_l[i];
            out.push(
                (s.flux[i] - conductance * (s.x_i[i] - pr.x[i]) - mean * total) / (f + conductance),
            );
        }
        out
    }

    fn evaluate(&self, values: &[f64], factor: f64) -> Result<Evaluation, String> {
        let size = self.unknowns().size();
        let segments: Vec<Segment> = values.chunks_exact(size).map(|b| self.decode(b)).collect();
        let props = segments
            .iter()
            .enumerate()
            .map(|(j, s)| self.props(j, s, factor))
            .collect::<Result<Vec<_>, _>>()?;
        let residual = self.residual(&segments, &props)?;
        Ok((segments, props, residual))
    }

    fn residual(&self, segments: &[Segment], props: &[Props]) -> Result<Vec<f64>, String> {
        let n = segments.len();
        let mut out = Vec::with_capacity(n * self.unknowns().size());
        for j in 0..n {
            out.extend(self.block_at(j, segments, props, None));
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err("rate-based column has a nonfinite residual".into());
        }
        Ok(out)
    }

    /// Block `j` with segment `swap.0` replaced by `swap.1`.
    fn block_at(
        &self,
        j: usize,
        segments: &[Segment],
        props: &[Props],
        swap: Option<(usize, &Segment, &Props)>,
    ) -> Vec<f64> {
        let get = |k: usize| match swap {
            Some((m, s, p)) if m == k => (s, p),
            _ => (&segments[k], &props[k]),
        };
        let above = (j > 0).then(|| get(j - 1));
        let below = (j + 1 < segments.len()).then(|| get(j + 1));
        self.block(above, get(j), below)
    }

    fn solve_at(&self, mut values: Vec<f64>, factor: f64) -> Result<(Vec<f64>, usize), String> {
        let size = self.unknowns().size();
        let n = self.segments();
        let (mut segments, mut props, mut r) = self.evaluate(&values, factor)?;
        for iteration in 0..=MAX_ITERATIONS {
            let error = norm(&r);
            if error <= RESIDUAL_TOLERANCE {
                return Ok((values, iteration));
            }
            if iteration == MAX_ITERATIONS {
                break;
            }
            let total = n * size;
            let mut jacobian = vec![vec![0.0; total]; total];
            for j in 0..n {
                let base = &values[j * size..(j + 1) * size];
                for k in 0..size {
                    let mut trial = base.to_vec();
                    let h = 1e-6 * base[k].abs().max(1.0);
                    let mut evaluate = |delta: f64| {
                        trial[k] = base[k] + delta;
                        let s = self.decode(&trial);
                        self.props(j, &s, factor).map(|p| (s, p))
                    };
                    let (delta, (s, p)) = match evaluate(h) {
                        Ok(v) => (h, v),
                        Err(_) => (
                            -h,
                            evaluate(-h).map_err(|e| format!("rate-based Jacobian: {e}"))?,
                        ),
                    };
                    let col = j * size + k;
                    for row_segment in j.saturating_sub(1)..(j + 2).min(n) {
                        let shifted =
                            self.block_at(row_segment, &segments, &props, Some((j, &s, &p)));
                        for (q, value) in shifted.iter().enumerate() {
                            let row = row_segment * size + q;
                            jacobian[row][col] = (value - r[row]) / delta;
                        }
                    }
                }
            }
            let step = solve_dense(jacobian, r.iter().map(|v| -v).collect())
                .ok_or("rate-based Jacobian is singular")?;
            let mut accepted = None;
            let mut trial = values.clone();
            for backtrack in 0..30 {
                let scale = 0.5_f64.powi(backtrack);
                for ((t, v), d) in trial.iter_mut().zip(&values).zip(&step) {
                    *t = v + scale * d;
                }
                if let Ok(next) = self.evaluate(&trial, factor)
                    && norm(&next.2) < error
                {
                    accepted = Some(next);
                    break;
                }
            }
            let Some(next) = accepted else {
                return Err(format!(
                    "rate-based line search stalled (scaled residual {error:.3e})"
                ));
            };
            values = trial;
            (segments, props, r) = next;
        }
        Err(format!(
            "rate-based column did not converge in {MAX_ITERATIONS} iterations (scaled residual {:.3e})",
            norm(&r)
        ))
    }

    /// Initial segments from an equilibrium cascade with one stage per segment.
    fn seed(&self, cascade: &StageCascade) -> Vec<f64> {
        let n = cascade.len();
        let mut out = Vec::new();
        for (j, stage) in cascade.stages.iter().enumerate() {
            let v: Vec<f64> = stage.y.iter().map(|y| y * stage.vapor_mol_s).collect();
            let l: Vec<f64> = stage.x.iter().map(|x| x * stage.liquid_mol_s).collect();
            let v_in: Vec<f64> = if j + 1 < n {
                let below = &cascade.stages[j + 1];
                below.y.iter().map(|y| y * below.vapor_mol_s).collect()
            } else {
                self.gas
                    .z
                    .iter()
                    .map(|z| z * self.gas.molar_flow_mol_s)
                    .collect()
            };
            let flux = v_in.iter().zip(&v).map(|(a, b)| a - b).collect();
            out.extend(self.encode(&Segment {
                v,
                l,
                t_v: stage.temperature_k,
                t_l: stage.temperature_k,
                t_i: stage.temperature_k,
                x_i: stage.x.clone(),
                y_i: stage.y.clone(),
                flux,
            }));
        }
        out
    }

    pub(super) fn solve(
        mut self,
        seed: &StageCascade,
        gas_inlet_pressure_pa: f64,
    ) -> Result<Solved, String> {
        // Hydraulic correlations outside their stated domain fail before
        // the solve (Robbins: at most 3 bar; Kister-Gill: Fp <= 60 1/ft).
        if self.pressure > pd::ROBBINS_MAX_PRESSURE_PA {
            return Err(format!(
                "rate-based hydraulics: the Robbins pressure-drop correlation is not suitable above 3 bar absolute (column pressure {} Pa)",
                self.pressure
            ));
        }
        pd::flood_pressure_drop_pa_m(&self.packing())
            .map_err(|e| format!("rate-based hydraulics: {e}"))?;
        let target = self.spec.interfacial_area_factor;
        let initial = self.seed(seed);
        let (mut values, mut iterations) = match self.solve_at(initial.clone(), target) {
            Ok(solved) => solved,
            Err(direct) => {
                let mut values = initial;
                let mut total = 0;
                for multiple in CONTINUATION {
                    let (next, used) = self
                        .solve_at(values, target * multiple)
                        .map_err(|e| format!("{e}; direct solve: {direct}"))?;
                    values = next;
                    total += used;
                }
                (values, total)
            }
        };
        let mut pressure_updates = 0;
        loop {
            let (segments, props, _) = self.evaluate(&values, target)?;
            let (flooding, drops) = self.hydraulics(&props)?;
            let mut pressures = Vec::with_capacity(segments.len() + 1);
            pressures.push(self.pressure);
            let dz = self.height_m / segments.len() as f64;
            let mut cumulative_drop = 0.0;
            for drop_per_m in &drops {
                cumulative_drop += drop_per_m * dz;
                pressures.push(self.pressure + cumulative_drop);
            }
            if !self.spec.couple_pressure_drop {
                return self.finish(&segments, &props, flooding, drops, pressures, iterations);
            }
            let bottom = *pressures.last().expect("bed boundary");
            // Keep trial pressures in the correlation's domain. A capped
            // fixed point is admissible only if the uncapped drop agrees.
            for pressure in &mut pressures {
                *pressure = pressure.min(pd::ROBBINS_MAX_PRESSURE_PA);
            }
            let change = pressures
                .iter()
                .zip(&self.boundary_pressures_pa)
                .map(|(next, current)| (next - current).abs())
                .fold(0.0, f64::max);
            if change <= PRESSURE_TOLERANCE_PA {
                if bottom > pd::ROBBINS_MAX_PRESSURE_PA + PRESSURE_TOLERANCE_PA {
                    return Err(format!(
                        "rate-based hydraulics: bed bottom pressure {bottom:.2} Pa exceeds the Robbins limit of 3 bar absolute"
                    ));
                }
                if bottom > gas_inlet_pressure_pa + PRESSURE_TOLERANCE_PA {
                    return Err(format!(
                        "rate-based bed bottom pressure {bottom:.2} Pa exceeds the gas inlet pressure {gas_inlet_pressure_pa:.2} Pa"
                    ));
                }
                return self.finish(
                    &segments,
                    &props,
                    flooding,
                    drops,
                    self.boundary_pressures_pa.clone(),
                    iterations,
                );
            }
            if pressure_updates == MAX_PRESSURE_ITERATIONS {
                return Err(format!(
                    "rate-based pressure profile did not converge in {MAX_PRESSURE_ITERATIONS} updates (maximum pressure change {change:.3e} Pa)"
                ));
            }
            self.boundary_pressures_pa = pressures;
            pressure_updates += 1;
            let (next, used) = self
                .solve_at(values, target)
                .map_err(|e| format!("rate-based pressure update {pressure_updates}: {e}"))?;
            values = next;
            iterations += used;
        }
    }

    fn hydraulics(&self, props: &[Props]) -> Result<(Vec<f64>, Vec<f64>), String> {
        let packing = self.packing();
        let mut flooding = Vec::with_capacity(props.len());
        let mut drops = Vec::with_capacity(props.len());
        for (j, pr) in props.iter().enumerate() {
            let p = self.boundary_pressures_pa[j].midpoint(self.boundary_pressures_pa[j + 1]);
            let fraction = pd::flooding_fraction(&packing, &pr.gas_film, &pr.liquid_film, p)?;
            if !(fraction.is_finite() && fraction < 1.0) {
                return Err(format!(
                    "rate-based segment {}: flooding exceeded (fraction of flood {fraction:.3})",
                    j + 1
                ));
            }
            let load = pr.liquid_velocity_m_s * 3600.0;
            let minimum = pd::minimum_wetting_rate_m_h(
                &packing,
                self.spec.wetting_surface,
                pr.surface_tension_n_m,
            );
            if load < minimum {
                return Err(format!(
                    "rate-based segment {}: liquid load {load:.2} m³/(m² h) is below the minimum wetting rate {minimum:.2} m³/(m² h) for {} ({})",
                    j + 1,
                    packing.name,
                    self.spec.wetting_surface.name()
                ));
            }
            flooding.push(fraction);
            drops.push(pd::robbins_pressure_drop_pa_m(
                &packing,
                &pr.gas_film,
                &pr.liquid_film,
                p,
            )?);
        }
        Ok((flooding, drops))
    }

    fn finish(
        &self,
        segments: &[Segment],
        props: &[Props],
        flooding: Vec<f64>,
        drops: Vec<f64>,
        boundary_pressures: Vec<f64>,
        iterations: usize,
    ) -> Result<Solved, String> {
        let ctx = self.ctx;
        let n = segments.len();
        let dz = self.height_m / n as f64;
        // The interface of every segment must be a stable vapor–liquid split.
        let interfaces = StageCascade {
            stages: segments
                .iter()
                .enumerate()
                .map(|(j, s)| {
                    Stage::equilibrium(
                        s.t_i,
                        self.boundary_pressures_pa[j].midpoint(self.boundary_pressures_pa[j + 1]),
                        s.liquid(),
                        s.vapor(),
                        s.x_i.clone(),
                        s.y_i.clone(),
                    )
                })
                .collect(),
            feeds: Vec::new(),
        };
        let mut warnings = mesh::check_liquid_liquid(ctx, &interfaces)?;
        mesh::verify_phases(ctx, &interfaces).map_err(|e| format!("rate-based interface: {e}"))?;
        let top = &segments[0];
        let bottom = &segments[n - 1];
        let gas_out = Stage::equilibrium(
            top.t_v,
            self.pressure,
            top.liquid(),
            top.vapor(),
            props[0].x.clone(),
            props[0].y.clone(),
        );
        let liquid_out = Stage::equilibrium(
            bottom.t_l,
            self.boundary_pressures_pa[n],
            bottom.liquid(),
            bottom.vapor(),
            props[n - 1].x.clone(),
            props[n - 1].y.clone(),
        );
        let outlets = vec![
            mesh::product(ctx, &gas_out, top.vapor(), Phase::Vapor)?,
            mesh::product(ctx, &liquid_out, bottom.liquid(), Phase::Liquid)?,
        ];
        let neglected = props.iter().map(|p| p.neglected).fold(0.0, f64::max);
        if !self.spec.gas_transport_neglected.is_empty()
            || !self.spec.liquid_transport_neglected.is_empty()
        {
            warnings.push(format!(
                "rate-based transport properties neglect trace components (gas: {:?}, liquid: {:?}; largest neglected fraction {neglected:.2e})",
                self.spec.gas_transport_neglected, self.spec.liquid_transport_neglected
            ));
        }
        let pressure_drop_pa = drops.iter().sum::<f64>() * dz;
        let max_flooding_fraction = flooding.iter().copied().fold(0.0, f64::max);
        let profile = RateBasedProfile {
            packing: self.spec.packing,
            diameter_m: self.spec.diameter_m.expect("checked diameter"),
            segment_height_m: dz,
            vapor_temperatures_k: segments.iter().map(|s| s.t_v).collect(),
            liquid_temperatures_k: segments.iter().map(|s| s.t_l).collect(),
            interface_temperatures_k: segments.iter().map(|s| s.t_i).collect(),
            vapor_flows_mol_s: segments.iter().map(Segment::vapor).collect(),
            liquid_flows_mol_s: segments.iter().map(Segment::liquid).collect(),
            vapor_mole_fractions: props.iter().map(|p| p.y.clone()).collect(),
            liquid_mole_fractions: props.iter().map(|p| p.x.clone()).collect(),
            fluxes_mol_s: segments.iter().map(|s| s.flux.clone()).collect(),
            interfacial_area_m2_m3: props.iter().map(|p| p.area_m2_m3).collect(),
            flooding_fractions: flooding,
            pressure_drops_pa_m: drops,
            pressure_drop_coupled: self.spec.couple_pressure_drop,
            estimated_boundary_pressures_pa: boundary_pressures,
            pressure_drop_pa,
            max_flooding_fraction,
            max_neglected_transport_fraction: neglected,
            iterations,
        };
        Ok((outlets, profile, warnings))
    }
}

fn norm(values: &[f64]) -> f64 {
    values.iter().map(|x| x.abs()).fold(0.0, f64::max)
}
