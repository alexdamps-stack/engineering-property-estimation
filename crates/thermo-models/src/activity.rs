//! Gamma–phi (modified Raoult) phase model: NRTL liquid with pure-liquid
//! vapor-pressure standard states and an explicitly selected vapor model.
//!
//! Equilibrium: y_i φ_i^V p = x_i γ_i p_i^sat (no Poynting factor, no φ^sat).
//! With the ideal-gas vapor model (φ^V = 1, liquid volume neglected in the Gibbs
//! energy), the consistent pure-liquid Gibbs energy is
//! g_i^L(T) = g_i^ig(T, p0) + RT ln(p_i^sat/p0), giving
//! h_i^L = h_i^ig − R T² dln p^sat/dT (Clausius–Clapeyron),
//! s_i^L = s_i^ig° − R ln(p^sat/p0) − R T dln p^sat/dT, and
//! c_p,i^L = c_p,i^ig − R (2T L' + T² L''). Liquid h, s, g are pressure independent.
//! Saturated-liquid density is a separately labelled correlation, combined by ideal
//! (Amagat) mixing, and is not coupled into fugacities or caloric properties.
use crate::nrtl::{ExcessState, NrtlModel, NrtlParameters};
use thermo_data::{Component, ComponentDatabase, LiquidDatabase, PureLiquid, REFERENCE_PRESSURE};
use thermo_types::{
    COMPOSITION_SUM_TOLERANCE, CalculationError as Error, CaloricReference, ComponentId,
    Composition, CompositionBasis, DataMetadata, InputError, MOLAR_GAS_CONSTANT as R,
    ModelMetadata, Phase, Pressure, Temperature, ValidityStatus,
};

/// Computational envelope of the gamma–phi package (K, Pa). The ideal-gas vapor
/// assumption restricts it to low pressure; it is not an accuracy statement.
pub const ACTIVITY_TEMPERATURE_RANGE_K: [f64; 2] = [273.16, 450.0];
pub const ACTIVITY_PRESSURE_RANGE_PA: [f64; 2] = [100.0, 500_000.0];

/// Vapor-phase model. Selected explicitly: there is no default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VaporModel {
    /// φ_i^V = 1; NASA7 ideal-gas caloric properties.
    IdealGas,
}

/// Pure-liquid quantities at T. Units: Pa, 1/K, J/mol, J/(mol K).
#[derive(Debug, Clone, Copy)]
pub struct PureLiquidState {
    pub vapor_pressure: f64,
    pub d_ln_vapor_pressure_dt: f64,
    pub d2_ln_vapor_pressure_dt2: f64,
    /// Clausius–Clapeyron enthalpy of vaporization implied by the model, R T² L'.
    pub enthalpy_of_vaporization: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub heat_capacity: f64,
}

/// Homogeneous phase state of the gamma–phi package. Vector order matches the model.
#[derive(Debug, Clone)]
pub struct ActivityPhaseState {
    pub phase: Phase,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub composition: Vec<f64>,
    /// ln φ_i with f_i = z_i φ_i p: liquid ln(γ_i p_i^sat/p), vapor 0.
    pub ln_fugacity_coefficients: Vec<f64>,
    pub enthalpy: f64,
    pub entropy: f64,
    pub gibbs_energy: f64,
    pub cp: f64,
    /// Vapor: RT/p. Liquid: Σ x_i/ρ_i^sat(T) from the density correlation (uncoupled).
    pub molar_volume: f64,
    pub mass_density: f64,
    /// Liquid only.
    pub excess: Option<ExcessState>,
    pub model: ModelMetadata,
    pub data: DataMetadata,
    pub caloric_reference: CaloricReference,
    /// Liquid: NRTL regression-interval status. Vapor: unassessed.
    pub validity: ValidityStatus,
}

#[derive(Debug, Clone)]
pub struct ActivityModel {
    vapor: VaporModel,
    components: Vec<Component>,
    liquids: Vec<PureLiquid>,
    nrtl: NrtlModel,
    fractions: Vec<f64>,
    data: DataMetadata,
}

impl ActivityModel {
    /// Parameter component order must equal the composition order exactly.
    pub fn new(
        database: &ComponentDatabase,
        liquids: &LiquidDatabase,
        composition: &Composition,
        vapor: VaporModel,
        parameters: NrtlParameters,
    ) -> Result<Self, Error> {
        if composition.basis() != CompositionBasis::Mole {
            return Err(Error::UnsupportedCapability(
                "convert mass composition explicitly before activity evaluation".into(),
            ));
        }
        let ids: Vec<&ComponentId> = composition.entries().iter().map(|(id, _)| id).collect();
        if parameters.component_ids().len() != ids.len()
            || parameters
                .component_ids()
                .iter()
                .zip(&ids)
                .any(|(a, b)| a != *b)
        {
            return Err(InputError::ComponentMismatch.into());
        }
        let components = ids
            .iter()
            .map(|id| database.get(id.as_str()).cloned())
            .collect::<Result<Vec<_>, _>>()?;
        let pure = ids
            .iter()
            .map(|id| liquids.pure(id.as_str()).cloned())
            .collect::<Result<Vec<_>, _>>()?;
        let (d, l) = (database.metadata(), liquids.metadata());
        Ok(Self {
            vapor,
            components,
            liquids: pure,
            nrtl: NrtlModel::new(parameters),
            fractions: composition.entries().iter().map(|(_, x)| *x).collect(),
            data: DataMetadata {
                identifier: format!("{}+{}", d.identifier, l.identifier),
                version: format!("{}+{}", d.version, l.version),
                source: format!("{}; {}", d.source, l.source),
            },
        })
    }

    /// NRTL interactions for every pair from the embedded dataset; missing pairs fail.
    pub fn from_databases(
        database: &ComponentDatabase,
        liquids: &LiquidDatabase,
        composition: &Composition,
        vapor: VaporModel,
    ) -> Result<Self, Error> {
        let ids: Vec<ComponentId> = composition
            .entries()
            .iter()
            .map(|(id, _)| id.clone())
            .collect();
        let parameters = NrtlParameters::from_database(liquids, &ids)?;
        Self::new(database, liquids, composition, vapor, parameters)
    }

    /// Same model at another mole-fraction vector in model order (finite, ≥ 0, Σ = 1 ± 1e-10).
    pub fn with_fractions(&self, fractions: &[f64]) -> Result<Self, Error> {
        if fractions.len() != self.components.len() {
            return Err(InputError::ComponentMismatch.into());
        }
        if fractions.iter().any(|x| !x.is_finite()) {
            return Err(InputError::NonFinite("mole fraction").into());
        }
        if fractions.iter().any(|x| *x < 0.0) {
            return Err(InputError::Negative("mole fraction").into());
        }
        if (fractions.iter().sum::<f64>() - 1.0).abs() > COMPOSITION_SUM_TOLERANCE {
            return Err(InputError::NotNormalized.into());
        }
        Ok(Self {
            fractions: fractions.to_vec(),
            ..self.clone()
        })
    }

    pub fn fractions(&self) -> &[f64] {
        &self.fractions
    }
    pub fn components(&self) -> &[Component] {
        &self.components
    }
    pub fn component_ids(&self) -> Vec<&str> {
        self.components.iter().map(Component::id).collect()
    }
    pub fn nrtl(&self) -> &NrtlModel {
        &self.nrtl
    }
    pub fn vapor_model(&self) -> VaporModel {
        self.vapor
    }

    pub fn metadata(&self) -> ModelMetadata {
        let nrtl = self.nrtl.metadata();
        ModelMetadata {
            identifier: format!("gamma-phi/{}/{:?}", nrtl.identifier, self.vapor),
            implementation_version: "a7-v1".into(),
            formulation: format!(
                "modified Raoult y_i p = x_i gamma_i psat_i; vapor {:?} (phi=1); no Poynting, no phi_sat; liquid g_i = g_i^ig(T,p0) + RT ln(psat_i/p0); psat: CoolProp pV ancillaries; liquid volume: saturated-liquid ancillary with ideal mixing, uncoupled; {}",
                self.vapor, nrtl.formulation
            ),
        }
    }

    pub fn caloric_reference() -> CaloricReference {
        CaloricReference {
            identifier: "ideal-pure-zero-298.15K-1bar-cp-integrated-v1".into(),
            definition: "NASA7 cp integrals; pure ideal h=s=0 at 298.15 K, 100000 Pa; ideal mixing entropy; liquid pure states from the vapor-pressure correlation via Clausius-Clapeyron with an ideal-gas vapor and negligible liquid volume, plus NRTL excess; J/mol and J/(mol K)".into(),
        }
    }

    /// Package envelope: 273.16..=450 K and 100..=5e5 Pa; every pure correlation must
    /// also cover T. Returns `OutsideValidityDomain` otherwise.
    pub fn check_domain(&self, temperature: Temperature, pressure: Pressure) -> Result<(), Error> {
        self.check_temperature(temperature)?;
        let [lo, hi] = ACTIVITY_PRESSURE_RANGE_PA;
        if !(lo..=hi).contains(&pressure.value()) {
            return Err(Error::OutsideValidityDomain(format!(
                "gamma-phi ideal-gas vapor package: {lo}..={hi} Pa"
            )));
        }
        Ok(())
    }

    pub fn check_temperature(&self, temperature: Temperature) -> Result<(), Error> {
        let t = temperature.value();
        let [lo, hi] = ACTIVITY_TEMPERATURE_RANGE_K;
        if !(lo..=hi).contains(&t) {
            return Err(Error::OutsideValidityDomain(format!(
                "gamma-phi package: {lo}..={hi} K"
            )));
        }
        for l in &self.liquids {
            let [a, b] = l.vapor_pressure_range();
            let [c, d] = l.liquid_density_range();
            if !(a..=b).contains(&t) || !(c..=d).contains(&t) {
                return Err(Error::OutsideValidityDomain(format!(
                    "{} liquid correlations: {}..={} K",
                    l.id(),
                    a.max(c),
                    b.min(d)
                )));
            }
        }
        Ok(())
    }

    /// Pure-liquid standard state of component `i` at T.
    pub fn pure_liquid(
        &self,
        i: usize,
        temperature: Temperature,
    ) -> Result<PureLiquidState, Error> {
        self.check_temperature(temperature)?;
        let t = temperature.value();
        let ideal = self.components[i].ideal_gas(temperature)?;
        let vp = self.liquids[i].vapor_pressure(temperature)?;
        let dh = R * t * t * vp.d_ln_pressure_dt;
        Ok(PureLiquidState {
            vapor_pressure: vp.pressure,
            d_ln_vapor_pressure_dt: vp.d_ln_pressure_dt,
            d2_ln_vapor_pressure_dt2: vp.d2_ln_pressure_dt2,
            enthalpy_of_vaporization: dh,
            enthalpy: ideal.enthalpy - dh,
            entropy: ideal.entropy
                - R * (vp.ln_pressure - REFERENCE_PRESSURE.ln())
                - R * t * vp.d_ln_pressure_dt,
            heat_capacity: ideal.heat_capacity
                - R * (2.0 * t * vp.d_ln_pressure_dt + t * t * vp.d2_ln_pressure_dt2),
        })
    }

    /// Liquid state at the model's fractions.
    pub fn liquid(
        &self,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<ActivityPhaseState, Error> {
        self.check_domain(temperature, pressure)?;
        let (t, p) = (temperature.value(), pressure.value());
        let x = &self.fractions;
        let ex = self.nrtl.excess(t, x)?;
        let (mut h, mut s, mut cp, mut v, mut mass) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let mut ln_phi = Vec::with_capacity(x.len());
        for (i, xi) in x.iter().enumerate() {
            let pure = self.pure_liquid(i, temperature)?;
            ln_phi.push(ex.ln_activity_coefficients[i] + pure.vapor_pressure.ln() - p.ln());
            h += xi * pure.enthalpy;
            s += xi * pure.entropy;
            cp += xi * pure.heat_capacity;
            if *xi > 0.0 {
                s -= R * xi * xi.ln();
            }
            v += xi / self.liquids[i].saturated_liquid_molar_density(temperature)?;
            mass += xi * self.components[i].molar_mass().value();
        }
        h += ex.excess_enthalpy;
        s += ex.excess_entropy;
        cp += ex.excess_heat_capacity;
        let validity = ex.validity;
        self.finish(
            Phase::Liquid,
            t,
            p,
            ln_phi,
            h,
            s,
            cp,
            v,
            mass,
            Some(ex),
            validity,
        )
    }

    /// Vapor state at the model's fractions.
    pub fn vapor(
        &self,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<ActivityPhaseState, Error> {
        self.check_domain(temperature, pressure)?;
        let (t, p) = (temperature.value(), pressure.value());
        let VaporModel::IdealGas = self.vapor;
        let (mut h, mut s, mut cp, mut mass) = (0.0, -R * (p / REFERENCE_PRESSURE).ln(), 0.0, 0.0);
        for (c, y) in self.components.iter().zip(&self.fractions) {
            let ideal = c.ideal_gas(temperature)?;
            h += y * ideal.enthalpy;
            s += y * ideal.entropy;
            cp += y * ideal.heat_capacity;
            if *y > 0.0 {
                s -= R * y * y.ln();
            }
            mass += y * c.molar_mass().value();
        }
        let ln_phi = vec![0.0; self.fractions.len()];
        self.finish(
            Phase::Vapor,
            t,
            p,
            ln_phi,
            h,
            s,
            cp,
            R * t / p,
            mass,
            None,
            ValidityStatus::Unassessed,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        &self,
        phase: Phase,
        t: f64,
        p: f64,
        ln_phi: Vec<f64>,
        h: f64,
        s: f64,
        cp: f64,
        v: f64,
        mass: f64,
        excess: Option<ExcessState>,
        validity: ValidityStatus,
    ) -> Result<ActivityPhaseState, Error> {
        if ![h, s, cp, v, mass / v]
            .iter()
            .chain(&ln_phi)
            .all(|x| x.is_finite())
            || v <= 0.0
        {
            return Err(Error::NumericalFailure(
                "nonfinite gamma-phi phase state".into(),
            ));
        }
        if cp <= 0.0 {
            return Err(Error::NumericalFailure(format!(
                "nonpositive {phase:?} heat capacity {cp} J/(mol K)"
            )));
        }
        Ok(ActivityPhaseState {
            phase,
            temperature_k: t,
            pressure_pa: p,
            composition: self.fractions.clone(),
            ln_fugacity_coefficients: ln_phi,
            enthalpy: h,
            entropy: s,
            gibbs_energy: h - t * s,
            cp,
            molar_volume: v,
            mass_density: mass / v,
            excess,
            model: self.metadata(),
            data: self.data.clone(),
            caloric_reference: Self::caloric_reference(),
            validity,
        })
    }
}
