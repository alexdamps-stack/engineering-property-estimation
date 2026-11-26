use std::sync::Arc;

use thermo_models::EosCg2021;
use thermo_types::{
    CalculationError as Error, CaloricReference, DataMetadata, Diagnostic, DiagnosticSeverity,
    ModelMetadata, Phase, ValidityStatus,
};

use crate::{
    Convergence, CubicPackage, EquilibriumPhase, EquilibriumState, FlowsheetPackage,
    PhaseProperties, RootSelection, SaturationKind, SaturationState,
};

/// A bounded single-fluid-phase CO₂–normal H₂ flowsheet package.
#[derive(Debug, Clone)]
pub struct EosCgPackage {
    identity: Arc<()>,
    components: Vec<String>,
    temperature: [f64; 2],
}

impl EosCgPackage {
    pub const COMPONENTS: [&str; 16] = EosCg2021::COMPONENTS;
    pub const TEMPERATURE_RANGE_K: [f64; 2] = [330.0, 350.0];
    pub const PRESSURE_RANGE_PA: [f64; 2] = [100_000.0, 2_000_000.0];
    pub const VERSION: &str = "multifluid-gas-v2";
    pub const DATA: &str = "EOS-CG-2021 CO2–H2 coefficients";
    pub const DATA_VERSION: &str = "a65b099b891df2c8f315cfbfc57ba820d6d4ca0f";
    pub const REFERENCE: &str = "eos-cg-2021-embedded-pure-fluid-offsets";

    pub fn new(components: Vec<String>, temperature: [f64; 2]) -> Result<Self, Error> {
        let ids: Vec<_> = components.iter().map(String::as_str).collect();
        let mut seed = vec![0.0; ids.len()];
        if let Some(first) = seed.first_mut() {
            *first = 1.0;
        }
        EosCg2021::new(&ids, &seed)?;
        let [lo, hi] = Self::TEMPERATURE_RANGE_K;
        if !temperature.into_iter().all(|t| (lo..=hi).contains(&t))
            || temperature[0] >= temperature[1]
        {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 requires ordered temperature bounds within 330–350 K".into(),
            ));
        }
        Ok(Self {
            identity: Arc::new(()),
            components,
            temperature,
        })
    }

    fn model(&self, fractions: &[f64]) -> Result<EosCg2021, Error> {
        EosCg2021::new(&self.component_ids(), fractions)
    }

    fn state(
        &self,
        model: &EosCg2021,
        fractions: &[f64],
        t: f64,
        p: f64,
    ) -> Result<EquilibriumState, Error> {
        if !(self.temperature[0]..=self.temperature[1]).contains(&t)
            || !(Self::PRESSURE_RANGE_PA[0]..=Self::PRESSURE_RANGE_PA[1]).contains(&p)
        {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 gas calculation requires the package temperature interval and 0.1–2 MPa".into(),
            ));
        }
        // Search from low density and stop at the first pressure bracket on
        // the rising branch. A dense endpoint can be beyond the gas spinodal.
        let [low, mut high] = [
            0.001,
            (p / (thermo_types::MOLAR_GAS_CONSTANT * t) / 4.0).max(0.002),
        ];
        let mut previous = low;
        loop {
            let (actual, slope, _) = model.isotherm(t, high)?;
            if slope > 0.0 && actual >= p {
                break;
            }
            if slope <= 0.0 || high >= 1000.0 {
                return Err(Error::OutsideValidityDomain(
                    "no bracket for a mechanically admissible EOS-CG-2021 gas root".into(),
                ));
            }
            previous = high;
            high = (high * 1.25).min(1000.0);
        }
        let (density, iterations) = bracket([previous, high], p * 1e-11, |density| {
            let (actual, slope, _) = model.isotherm(t, density)?;
            if slope <= 0.0 {
                return Err(Error::OutsideValidityDomain(
                    "EOS-CG-2021 density solver left the gas branch".into(),
                ));
            }
            Ok((actual - p, density))
        })?;
        let calculated = model.evaluate(t, density)?;
        Ok(EquilibriumState {
            temperature_k: t,
            pressure_pa: p,
            mole_fractions: fractions.to_vec(),
            phases: vec![EquilibriumPhase {
                phase: Phase::Vapor,
                phase_fraction: 1.0,
                composition: fractions.to_vec(),
                state: PhaseProperties {
                    temperature_k: t,
                    pressure_pa: p,
                    compressibility: calculated.compressibility,
                    molar_volume: density.recip(),
                    mass_density: calculated.mass_density,
                    ln_fugacity_coefficients: calculated.ln_fugacity_coefficients,
                    residual_enthalpy: calculated.residual_enthalpy,
                    enthalpy: calculated.enthalpy,
                    entropy: calculated.entropy,
                    cp: calculated.cp,
                    cv: Some(calculated.cv),
                    speed_of_sound: Some(calculated.speed_of_sound),
                    dynamic_viscosity: None,
                    thermal_conductivity: None,
                    transport: None,
                    model: ModelMetadata {
                        identifier: "EOS-CG-2021".into(),
                        implementation_version: Self::VERSION.into(),
                        formulation: "CO2–normal H2 multiparameter Helmholtz mixture".into(),
                    },
                    data: DataMetadata {
                        identifier: Self::DATA.into(),
                        version: Self::DATA_VERSION.into(),
                        source: "Neumann et al. (2023), DOI 10.1007/s10765-023-03263-6".into(),
                    },
                    caloric_reference: CaloricReference {
                        identifier: Self::REFERENCE.into(),
                        definition: "Embedded pure-fluid enthalpy and entropy offsets, with ideal entropy of mixing".into(),
                    },
                    validity: ValidityStatus::Unassessed,
                    inferred_regression_interval: false,
                },
            }],
            iterations,
            convergence: Convergence {
                max_ln_fugacity_residual: 0.0,
                max_component_balance_residual: None,
                scaled_energy_residual: None,
                scaled_entropy_residual: None,
            },
            diagnostics: vec![Diagnostic {
                code: "EOS-CG-2021-GAS-BRANCH".into(),
                severity: DiagnosticSeverity::Information,
                message: "Single gas branch, CO₂/normal H₂, 330–350 K and 0.1–2 MPa. Global phase stability is unassessed.".into(),
            }],
        })
    }
}

impl FlowsheetPackage for EosCgPackage {
    fn identity(&self) -> &Arc<()> {
        &self.identity
    }
    fn component_ids(&self) -> Vec<&str> {
        self.components.iter().map(String::as_str).collect()
    }
    fn temperature_interval(&self) -> [f64; 2] {
        self.temperature
    }
    fn validate_composition(&self, fractions: &[f64]) -> Result<(), Error> {
        self.model(fractions).map(|_| ())
    }

    fn tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumState, Error> {
        self.state(
            &self.model(fractions)?,
            fractions,
            temperature_k,
            pressure_pa,
        )
    }

    fn ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<EquilibriumState, Error> {
        if !enthalpy.is_finite() {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 requires a finite molar enthalpy".into(),
            ));
        }
        let model = self.model(fractions)?;
        let scale = enthalpy.abs().max(1.0);
        let (mut state, iterations) = bracket(self.temperature, scale * 1e-10, |t| {
            let state = self.state(&model, fractions, t, pressure_pa)?;
            Ok((state.enthalpy() - enthalpy, state))
        })?;
        state.convergence.scaled_energy_residual = Some((state.enthalpy() - enthalpy) / scale);
        state.iterations += iterations;
        Ok(state)
    }

    fn cubic(&self) -> Result<&CubicPackage, Error> {
        Err(Error::UnsupportedCapability(
            "this unit requires a cubic property package; selected EOS-CG-2021".into(),
        ))
    }

    /// PS on the single gas branch, with the structure of PH: a bracketed
    /// temperature solve over the package interval. Entropy increases with T at
    /// fixed p on the branch (c_p > 0).
    fn ps(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        entropy: f64,
    ) -> Result<EquilibriumState, Error> {
        if !entropy.is_finite() {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 requires a finite molar entropy".into(),
            ));
        }
        let model = self.model(fractions)?;
        let scale = entropy.abs().max(1.0);
        let (mut state, iterations) = bracket(self.temperature, scale * 1e-10, |t| {
            let state = self.state(&model, fractions, t, pressure_pa)?;
            Ok((state.entropy() - entropy, state))
        })?;
        state.convergence.scaled_entropy_residual = Some((state.entropy() - entropy) / scale);
        state.iterations += iterations;
        Ok(state)
    }

    fn homogeneous(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<PhaseProperties, Error> {
        gas_branch_only(root)?;
        let mut state = self.tp(fractions, temperature_k, pressure_pa)?;
        Ok(state.phases.remove(0).state)
    }

    fn ln_fugacity_coefficients(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<Vec<f64>, Error> {
        Ok(self
            .homogeneous(fractions, temperature_k, pressure_pa, root)?
            .ln_fugacity_coefficients)
    }

    fn saturation_temperature_interval(&self) -> Result<[f64; 2], Error> {
        Err(gas_branch_error())
    }

    fn saturation_pressure(
        &self,
        _fractions: &[f64],
        _temperature_k: f64,
        _kind: SaturationKind,
    ) -> Result<SaturationState, Error> {
        Err(gas_branch_error())
    }

    fn saturation_temperature(
        &self,
        _fractions: &[f64],
        _pressure_pa: f64,
        _kind: SaturationKind,
    ) -> Result<SaturationState, Error> {
        Err(gas_branch_error())
    }

    fn formation_components(&self) -> Result<&[thermo_data::Component], Error> {
        Err(Error::UnsupportedCapability(
            "reactors need the ideal-gas formation basis (h = s = 0 at 298.15 K and 1 bar); EOS-CG-2021 uses embedded pure-fluid offsets and supplies no conversion offset".into(),
        ))
    }
}

/// The package solves only the gas branch: a liquid root, a minimum-Gibbs
/// selection (which needs every root) and phase boundaries are unsupported.
fn gas_branch_only(root: RootSelection) -> Result<(), Error> {
    if root == RootSelection::LargestVolume {
        Ok(())
    } else {
        Err(gas_branch_error())
    }
}

fn gas_branch_error() -> Error {
    Error::UnsupportedCapability(
        "EOS-CG-2021 package evaluates only its single gas branch (330–350 K, 0.1–2 MPa); liquid roots, minimum-Gibbs selection, K-values and bubble/dew points are unsupported".into(),
    )
}

fn bracket<T>(
    [mut lo, mut hi]: [f64; 2],
    tolerance: f64,
    evaluate: impl Fn(f64) -> Result<(f64, T), Error>,
) -> Result<(T, usize), Error> {
    let (low, low_state) = evaluate(lo)?;
    let (high, high_state) = evaluate(hi)?;
    if low.abs() <= tolerance {
        return Ok((low_state, 0));
    }
    if high.abs() <= tolerance {
        return Ok((high_state, 0));
    }
    if low >= 0.0 || high <= 0.0 {
        return Err(Error::OutsideValidityDomain(
            "EOS-CG-2021 target falls outside the gas-phase solve bracket".into(),
        ));
    }
    for iteration in 1..=64 {
        let mid = lo.midpoint(hi);
        let (residual, state) = evaluate(mid)?;
        if residual.abs() <= tolerance {
            return Ok((state, iteration));
        }
        if residual < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(
        "EOS-CG-2021 gas-phase solve did not converge".into(),
    ))
}
