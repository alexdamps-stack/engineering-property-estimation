//! Native GERG adapter for the shared finite-trial TPD / two-phase TP algorithm.
use crate::equilibrium_model::EquilibriumModel;
use crate::{FlashOptions, StabilityAssessment, StabilityOptions, TpFlash};
use thermo_data::ComponentDatabase;
use thermo_models::{
    FugacityState, GERG_CALORIC_REFERENCE, GERG_DATA_REVISION, GERG_GAS_CONSTANT as R,
    GERG_MODEL_VERSION, GergFugacityState, GergModel, HomogeneousState, RootSelection,
};
use thermo_types::{
    CalculationError as Error, CaloricReference, DataMetadata, Diagnostic, DiagnosticSeverity,
    ModelMetadata, Phase, Pressure, Temperature, ValidityStatus,
};

const COMPONENTS: [&str; 5] = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"];

pub(crate) fn validate(model: &GergModel, t: Temperature, p: Pressure) -> Result<(), Error> {
    if model
        .component_ids()
        .iter()
        .any(|id| !COMPONENTS.contains(id))
    {
        return Err(Error::UnsupportedCapability("GERG equilibrium is qualified only for methane/ethane/propane/nitrogen/carbon-dioxide; no aqueous or heavier-fluid equilibrium".into()));
    }
    if !(230.0..=350.0).contains(&t.value()) || !(1e5..=1e7).contains(&p.value()) {
        return Err(Error::OutsideValidityDomain(
            "GERG equilibrium computational envelope: 230..=350 K, 0.1..=10 MPa".into(),
        ));
    }
    Ok(())
}

/// Locate regular pressure roots on the two exterior branches of a sampled
/// GERG isotherm. Interior roots between sampled negative-stiffness regions are
/// excluded because multiparameter EOS can have nonphysical two-phase loops.
/// Finite sampling does not prove completeness; near-critical roots are unqualified.
pub fn gerg_outer_roots(
    model: &GergModel,
    t: Temperature,
    p: Pressure,
) -> Result<Vec<GergFugacityState>, Error> {
    validate(model, t, p)?;
    let target = p.value();
    let rt = R * t.value();
    let sample = |rho: f64| -> Result<(f64, f64), Error> {
        let a = model.helmholtz(t, rho)?;
        Ok((
            rho * rt * (1.0 + a.ar_d) - target,
            rt * (1.0 + 2.0 * a.ar_d + a.ar_dd),
        ))
    };
    let mut grid = Vec::with_capacity(386);
    for i in 0..=128 {
        grid.push(0.001 * (i as f64 / 128.0 * (30000.0_f64 / 0.001).ln()).exp());
    }
    for i in 0..=256 {
        grid.push(0.001 + i as f64 / 256.0 * (30000.0 - 0.001));
    }
    // Clamp endpoint roundoff so the homogeneous envelope is never exceeded.
    grid.iter_mut()
        .for_each(|rho| *rho = rho.clamp(0.001, 30000.0));
    grid.sort_by(f64::total_cmp);
    grid.dedup();
    let mut roots = Vec::new();
    let mut first_negative = f64::INFINITY;
    let mut last_negative = f64::NEG_INFINITY;
    let (mut previous, mut left) = (grid[0], sample(grid[0])?.0);
    for &current in grid.iter().skip(1) {
        let (right, stiffness) = sample(current)?;
        if stiffness <= 0.0 {
            first_negative = first_negative.min(current);
            last_negative = current;
        }
        if left <= 0.0 && right >= 0.0 {
            let (mut lo, mut hi) = (previous, current);
            let mut best = (f64::INFINITY, lo);
            for _ in 0..100 {
                let rho = lo + 0.5 * (hi - lo);
                let (f, _) = sample(rho)?;
                if f.abs() < best.0 {
                    best = (f.abs(), rho);
                }
                if f.abs() <= (5e-13 * target).max(2e-7) {
                    break;
                }
                if rho == lo || rho == hi {
                    break;
                }
                if f < 0.0 {
                    lo = rho;
                } else {
                    hi = rho;
                }
            }
            if roots
                .last()
                .is_none_or(|last: &(f64, f64)| (best.1 - last.0).abs() > 1e-8 * best.1)
            {
                roots.push((best.1, best.0));
            }
        }
        previous = current;
        left = right;
    }
    let mut states = Vec::new();
    for (rho, residual) in roots {
        if rho > first_negative && rho < last_negative {
            continue;
        }
        if residual > (5e-12 * target).max(2e-7) {
            return Err(Error::NonConvergence(format!(
                "GERG exterior density root failed pressure convergence: T={} K, p={target} Pa, rho={rho} mol/m³, residual={residual} Pa",
                t.value()
            )));
        }
        // On a candidate exterior crossing local admissibility must hold;
        // unresolved thermal/mechanical failures are not silently skipped.
        states.push(model.fugacity(t, rho)?);
    }
    if states.is_empty() {
        return Err(Error::OutsideValidityDomain(
            "no regular exterior GERG pressure root found inside the bounded density search".into(),
        ));
    }
    Ok(states)
}

#[derive(Clone)]
pub(crate) struct GergEquilibrium {
    pub(crate) model: GergModel,
    constants: Vec<(f64, f64, f64)>,
}

impl GergEquilibrium {
    pub(crate) fn new(model: &GergModel) -> Result<Self, Error> {
        let db = ComponentDatabase::embedded()?;
        let constants = model
            .component_ids()
            .iter()
            .map(|id| {
                let c = db.get(id)?;
                Ok((
                    c.critical_temperature(),
                    c.critical_pressure(),
                    c.acentric_factor(),
                ))
            })
            .collect::<Result<_, Error>>()?;
        Ok(Self {
            model: model.clone(),
            constants,
        })
    }
    fn selected(
        &self,
        t: Temperature,
        p: Pressure,
        selection: RootSelection,
    ) -> Result<(GergFugacityState, usize), Error> {
        let states = gerg_outer_roots(&self.model, t, p)?;
        let count = states.len();
        let index = match selection {
            RootSelection::LargestVolume => 0,
            RootSelection::SmallestVolume => count - 1,
            RootSelection::MinimumGibbs => {
                states
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.state.gibbs_energy.total_cmp(&b.state.gibbs_energy))
                    .unwrap()
                    .0
            }
        };
        Ok((states.into_iter().nth(index).unwrap(), count))
    }
}

impl EquilibriumModel for GergEquilibrium {
    fn fractions(&self) -> &[f64] {
        self.model.fractions()
    }
    fn with_fractions(&self, x: &[f64]) -> Result<Self, Error> {
        Ok(Self {
            model: self.model.with_fractions(x)?,
            constants: self.constants.clone(),
        })
    }
    fn ln_fugacity(
        &self,
        t: Temperature,
        p: Pressure,
        selection: RootSelection,
    ) -> Result<FugacityState, Error> {
        let (f, count) = self.selected(t, p, selection)?;
        let reduced_residual_gibbs = self
            .fractions()
            .iter()
            .zip(&f.ln_fugacity_coefficients)
            .map(|(x, phi)| x * phi)
            .sum();
        Ok(FugacityState {
            compressibility: f.state.compressibility,
            molar_volume: 1.0 / f.state.molar_density,
            ln_fugacity_coefficients: f.ln_fugacity_coefficients,
            reduced_residual_gibbs,
            root_count: count,
        })
    }
    fn evaluate(
        &self,
        t: Temperature,
        p: Pressure,
        selection: RootSelection,
    ) -> Result<HomogeneousState, Error> {
        let (f, _) = self.selected(t, p, selection)?;
        let s = f.state;
        let a = s.derivatives;
        Ok(HomogeneousState {temperature_k:s.temperature_k,pressure_pa:s.pressure_pa,compressibility:s.compressibility,
            molar_volume:1.0/s.molar_density,mass_density:s.mass_density,ln_fugacity_coefficients:f.ln_fugacity_coefficients,
            residual_enthalpy:R*t.value()*(a.ar_t+a.ar_d),residual_entropy:R*(a.ar_t-a.ar+s.compressibility.ln()),
            enthalpy:s.enthalpy,entropy:s.entropy,internal_energy:s.internal_energy,cp:s.cp,cv:s.cv,speed_of_sound:s.speed_of_sound,
            dp_dt:s.dp_dt,dp_dv:-s.dp_drho*s.molar_density*s.molar_density,
            model:ModelMetadata {identifier:"GERG-2008".into(),implementation_version:GERG_MODEL_VERSION.into(),
                formulation:"native GERG residual and ideal; R=8.314472; finite exterior-branch pressure search; minimum Gibbs among exterior roots".into()},
            data:DataMetadata {identifier:"NIST AGA8 GERG2008 2.01".into(),version:GERG_DATA_REVISION.into(),source:"data/source-extracts/gerg; native parameters and ideal reference".into()},
            caloric_reference:CaloricReference {identifier:"gerg-native-ideal".into(),definition:GERG_CALORIC_REFERENCE.into()},
            validity:ValidityStatus::Unassessed,selection})
    }
    fn wilson_values(&self, t: Temperature, p: Pressure) -> Vec<f64> {
        self.constants
            .iter()
            .map(|(tc, pc, w)| pc / p.value() * (5.373 * (1.0 + w) * (1.0 - tc / t.value())).exp())
            .collect()
    }
    fn phase_label(&self, s: &HomogeneousState) -> Phase {
        let pc: f64 = self
            .constants
            .iter()
            .zip(self.fractions())
            .map(|((_, pc, _), x)| pc * x)
            .sum();
        if s.temperature_k > self.model.reducing_temperature() && s.pressure_pa > pc {
            Phase::Supercritical
        } else if 1.0 / s.molar_volume < self.model.reducing_density() {
            Phase::Vapor
        } else {
            Phase::Liquid
        }
    }
}

/// Finite-trial stability on GERG exterior density branches in the dry-gas domain.
pub fn gerg_assess_stability(
    model: &GergModel,
    t: Temperature,
    p: Pressure,
    options: StabilityOptions,
) -> Result<StabilityAssessment, Error> {
    validate(model, t, p)?;
    options.validate()?;
    crate::stability::assess_stability_impl(&GergEquilibrium::new(model)?, t, p, options)
}

/// Native dry-gas GERG TP flash; at most two fluid phases, finite-trial stability.
pub fn gerg_tp_flash(
    model: &GergModel,
    t: Temperature,
    p: Pressure,
    options: FlashOptions,
) -> Result<TpFlash, Error> {
    validate(model, t, p)?;
    options.validate()?;
    if model.pure_critical_temperature().is_some() {
        let roots = gerg_outer_roots(model, t, p)?;
        if roots.len() > 1 {
            let active = model.fractions().iter().position(|x| *x > 0.0).unwrap();
            let difference = roots[0].ln_fugacity_coefficients[active]
                - roots.last().unwrap().ln_fugacity_coefficients[active];
            if difference.abs() <= options.tolerance {
                return Err(Error::UnsupportedCapability("pure GERG TP coexistence cannot determine a phase fraction; supply enthalpy through GERG PH flash".into()));
            }
        }
    }
    let mut result = crate::flash::tp_flash_impl(&GergEquilibrium::new(model)?, t, p, options)?;
    result.diagnostics.push(Diagnostic {code:"A6-GERG-EXTERIOR-ROOTS".into(),severity:DiagnosticSeverity::Information,
        message:"GERG dry five-component fluid domain; exterior density branches from finite log/linear isotherm sampling. Interior nonphysical loops excluded; root completeness, solids and global stability are unassessed. Wilson constants initialize trials only; native GERG supplies fugacity and calorics.".into()});
    Ok(result)
}
