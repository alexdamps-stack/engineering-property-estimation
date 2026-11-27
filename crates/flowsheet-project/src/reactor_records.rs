//! Format 4 records of the equilibrium, Gibbs and kinetic reactor
//! specifications (JSON spelling in SI units; `null` means unset).

use flowsheet_core::{
    Arrhenius, BatchThermal, ConcentrationBasis, Correlation, EquilibriumConstant,
    EquilibriumOptions, EquilibriumReaction, FugacityBasis, KineticReaction, Kinetics, LhhwTerm,
    LiquidVolumeBasis, PlugFlowThermal, RateBasis, RateLaw, RateTerm, ReactingPhase,
};
use serde::{Deserialize, Serialize};

use crate::records::CoefficientRecord;

pub(crate) fn coefficients_to(c: Vec<(String, f64)>) -> Vec<CoefficientRecord> {
    c.into_iter()
        .map(|(component, coefficient)| CoefficientRecord {
            component,
            coefficient,
        })
        .collect()
}

pub(crate) fn coefficients_from(c: Vec<CoefficientRecord>) -> Vec<(String, f64)> {
    c.into_iter()
        .map(|r| (r.component, r.coefficient))
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CorrelationRecord {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl From<Correlation> for CorrelationRecord {
    fn from(c: Correlation) -> Self {
        Self {
            a: c.a,
            b: c.b,
            c: c.c,
            d: c.d,
        }
    }
}

impl From<CorrelationRecord> for Correlation {
    fn from(c: CorrelationRecord) -> Self {
        Self {
            a: c.a,
            b: c.b,
            c: c.c,
            d: c.d,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ConstantRecord {
    GibbsEnergy,
    Correlation(CorrelationRecord),
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct EquilibriumReactionRecord {
    pub coefficients: Vec<CoefficientRecord>,
    pub constant: ConstantRecord,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum FugacityRecord {
    Package,
    IdealGas,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct EquilibriumOptionsRecord {
    pub fugacity: FugacityRecord,
    pub temperature_approach_k: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ArrheniusRecord {
    pub pre_exponential: Option<f64>,
    pub activation_energy_j_mol: f64,
    pub temperature_exponent: f64,
    pub reference_temperature_k: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct OrderRecord {
    pub component: String,
    pub order: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RateTermRecord {
    pub constant: ArrheniusRecord,
    pub orders: Vec<OrderRecord>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct LhhwTermRecord {
    pub constant: CorrelationRecord,
    pub orders: Vec<OrderRecord>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum RateLawRecord {
    PowerLaw {
        forward: RateTermRecord,
        reverse: Option<RateTermRecord>,
    },
    Lhhw {
        kinetic: ArrheniusRecord,
        driving_forward: LhhwTermRecord,
        driving_reverse: Option<LhhwTermRecord>,
        adsorption: Vec<LhhwTermRecord>,
        exponent: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum PhaseRecord {
    Vapor,
    Liquid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum BasisRecord {
    MolarConcentrationMolM3,
    MassConcentrationKgM3,
    MoleFraction,
    MassFraction,
    PartialPressurePa,
    FugacityPa,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum RateBasisRecord {
    ReactorVolume,
    CatalystMass,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct KineticReactionRecord {
    pub coefficients: Vec<CoefficientRecord>,
    pub phase: PhaseRecord,
    pub basis: BasisRecord,
    pub rate_basis: RateBasisRecord,
    pub rate: RateLawRecord,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct KineticsRecord {
    pub reactions: Vec<KineticReactionRecord>,
    pub catalyst_bulk_density_kg_m3: Option<f64>,
    /// Liquid volume basis (#228). Omitted when it is the EOS default, so
    /// older projects and their fingerprints are unchanged.
    #[serde(default, skip_serializing_if = "LiquidVolumeRecord::is_default")]
    pub liquid_volume: LiquidVolumeRecord,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum LiquidVolumeRecord {
    #[default]
    EquationOfState,
    EngineeringDensity,
}

impl LiquidVolumeRecord {
    fn is_default(&self) -> bool {
        *self == Self::default()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum PlugFlowThermalRecord {
    IsothermalK(Option<f64>),
    Adiabatic,
    Coolant {
        temperature_k: Option<f64>,
        u_w_m2_k: Option<f64>,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum BatchThermalRecord {
    IsothermalK(Option<f64>),
    Adiabatic,
}

// ----- conversions -----

impl From<EquilibriumReaction> for EquilibriumReactionRecord {
    fn from(r: EquilibriumReaction) -> Self {
        Self {
            coefficients: coefficients_to(r.coefficients),
            constant: match r.constant {
                EquilibriumConstant::GibbsEnergy => ConstantRecord::GibbsEnergy,
                EquilibriumConstant::Correlation(c) => ConstantRecord::Correlation(c.into()),
            },
        }
    }
}

impl From<EquilibriumReactionRecord> for EquilibriumReaction {
    fn from(r: EquilibriumReactionRecord) -> Self {
        Self {
            coefficients: coefficients_from(r.coefficients),
            constant: match r.constant {
                ConstantRecord::GibbsEnergy => EquilibriumConstant::GibbsEnergy,
                ConstantRecord::Correlation(c) => EquilibriumConstant::Correlation(c.into()),
            },
        }
    }
}

impl From<EquilibriumOptions> for EquilibriumOptionsRecord {
    fn from(o: EquilibriumOptions) -> Self {
        Self {
            fugacity: match o.fugacity {
                FugacityBasis::Package => FugacityRecord::Package,
                FugacityBasis::IdealGas => FugacityRecord::IdealGas,
            },
            temperature_approach_k: o.temperature_approach_k,
        }
    }
}

impl From<EquilibriumOptionsRecord> for EquilibriumOptions {
    fn from(o: EquilibriumOptionsRecord) -> Self {
        Self {
            fugacity: match o.fugacity {
                FugacityRecord::Package => FugacityBasis::Package,
                FugacityRecord::IdealGas => FugacityBasis::IdealGas,
            },
            temperature_approach_k: o.temperature_approach_k,
        }
    }
}

impl From<Arrhenius> for ArrheniusRecord {
    fn from(a: Arrhenius) -> Self {
        Self {
            pre_exponential: a.pre_exponential,
            activation_energy_j_mol: a.activation_energy_j_mol,
            temperature_exponent: a.temperature_exponent,
            reference_temperature_k: a.reference_temperature_k,
        }
    }
}

impl From<ArrheniusRecord> for Arrhenius {
    fn from(a: ArrheniusRecord) -> Self {
        Self {
            pre_exponential: a.pre_exponential,
            activation_energy_j_mol: a.activation_energy_j_mol,
            temperature_exponent: a.temperature_exponent,
            reference_temperature_k: a.reference_temperature_k,
        }
    }
}

fn orders_to(o: Vec<(String, f64)>) -> Vec<OrderRecord> {
    o.into_iter()
        .map(|(component, order)| OrderRecord { component, order })
        .collect()
}

fn orders_from(o: Vec<OrderRecord>) -> Vec<(String, f64)> {
    o.into_iter().map(|r| (r.component, r.order)).collect()
}

fn term_to(t: RateTerm) -> RateTermRecord {
    RateTermRecord {
        constant: t.constant.into(),
        orders: orders_to(t.orders),
    }
}

fn term_from(t: RateTermRecord) -> RateTerm {
    RateTerm {
        constant: t.constant.into(),
        orders: orders_from(t.orders),
    }
}

fn lhhw_to(t: LhhwTerm) -> LhhwTermRecord {
    LhhwTermRecord {
        constant: t.constant.into(),
        orders: orders_to(t.orders),
    }
}

fn lhhw_from(t: LhhwTermRecord) -> LhhwTerm {
    LhhwTerm {
        constant: t.constant.into(),
        orders: orders_from(t.orders),
    }
}

impl From<Kinetics> for KineticsRecord {
    fn from(k: Kinetics) -> Self {
        Self {
            catalyst_bulk_density_kg_m3: k.catalyst_bulk_density_kg_m3,
            liquid_volume: match k.liquid_volume {
                LiquidVolumeBasis::EquationOfState => LiquidVolumeRecord::EquationOfState,
                LiquidVolumeBasis::EngineeringDensity => LiquidVolumeRecord::EngineeringDensity,
            },
            reactions: k
                .reactions
                .into_iter()
                .map(|r| KineticReactionRecord {
                    coefficients: coefficients_to(r.coefficients),
                    phase: match r.phase {
                        ReactingPhase::Vapor => PhaseRecord::Vapor,
                        ReactingPhase::Liquid => PhaseRecord::Liquid,
                    },
                    basis: match r.basis {
                        ConcentrationBasis::MolarConcentration => {
                            BasisRecord::MolarConcentrationMolM3
                        }
                        ConcentrationBasis::MassConcentration => BasisRecord::MassConcentrationKgM3,
                        ConcentrationBasis::MoleFraction => BasisRecord::MoleFraction,
                        ConcentrationBasis::MassFraction => BasisRecord::MassFraction,
                        ConcentrationBasis::PartialPressure => BasisRecord::PartialPressurePa,
                        ConcentrationBasis::Fugacity => BasisRecord::FugacityPa,
                    },
                    rate_basis: match r.rate_basis {
                        RateBasis::ReactorVolume => RateBasisRecord::ReactorVolume,
                        RateBasis::CatalystMass => RateBasisRecord::CatalystMass,
                    },
                    rate: match r.rate {
                        RateLaw::PowerLaw { forward, reverse } => RateLawRecord::PowerLaw {
                            forward: term_to(forward),
                            reverse: reverse.map(term_to),
                        },
                        RateLaw::Lhhw {
                            kinetic,
                            driving_forward,
                            driving_reverse,
                            adsorption,
                            exponent,
                        } => RateLawRecord::Lhhw {
                            kinetic: kinetic.into(),
                            driving_forward: lhhw_to(driving_forward),
                            driving_reverse: driving_reverse.map(lhhw_to),
                            adsorption: adsorption.into_iter().map(lhhw_to).collect(),
                            exponent,
                        },
                    },
                })
                .collect(),
        }
    }
}

impl From<KineticsRecord> for Kinetics {
    fn from(k: KineticsRecord) -> Self {
        Self {
            catalyst_bulk_density_kg_m3: k.catalyst_bulk_density_kg_m3,
            liquid_volume: match k.liquid_volume {
                LiquidVolumeRecord::EquationOfState => LiquidVolumeBasis::EquationOfState,
                LiquidVolumeRecord::EngineeringDensity => LiquidVolumeBasis::EngineeringDensity,
            },
            reactions: k
                .reactions
                .into_iter()
                .map(|r| KineticReaction {
                    coefficients: coefficients_from(r.coefficients),
                    phase: match r.phase {
                        PhaseRecord::Vapor => ReactingPhase::Vapor,
                        PhaseRecord::Liquid => ReactingPhase::Liquid,
                    },
                    basis: match r.basis {
                        BasisRecord::MolarConcentrationMolM3 => {
                            ConcentrationBasis::MolarConcentration
                        }
                        BasisRecord::MassConcentrationKgM3 => ConcentrationBasis::MassConcentration,
                        BasisRecord::MoleFraction => ConcentrationBasis::MoleFraction,
                        BasisRecord::MassFraction => ConcentrationBasis::MassFraction,
                        BasisRecord::PartialPressurePa => ConcentrationBasis::PartialPressure,
                        BasisRecord::FugacityPa => ConcentrationBasis::Fugacity,
                    },
                    rate_basis: match r.rate_basis {
                        RateBasisRecord::ReactorVolume => RateBasis::ReactorVolume,
                        RateBasisRecord::CatalystMass => RateBasis::CatalystMass,
                    },
                    rate: match r.rate {
                        RateLawRecord::PowerLaw { forward, reverse } => RateLaw::PowerLaw {
                            forward: term_from(forward),
                            reverse: reverse.map(term_from),
                        },
                        RateLawRecord::Lhhw {
                            kinetic,
                            driving_forward,
                            driving_reverse,
                            adsorption,
                            exponent,
                        } => RateLaw::Lhhw {
                            kinetic: kinetic.into(),
                            driving_forward: lhhw_from(driving_forward),
                            driving_reverse: driving_reverse.map(lhhw_from),
                            adsorption: adsorption.into_iter().map(lhhw_from).collect(),
                            exponent,
                        },
                    },
                })
                .collect(),
        }
    }
}

impl From<PlugFlowThermal> for PlugFlowThermalRecord {
    fn from(t: PlugFlowThermal) -> Self {
        match t {
            PlugFlowThermal::Isothermal(t) => Self::IsothermalK(t),
            PlugFlowThermal::Adiabatic => Self::Adiabatic,
            PlugFlowThermal::Coolant {
                temperature_k,
                u_w_m2_k,
            } => Self::Coolant {
                temperature_k,
                u_w_m2_k,
            },
        }
    }
}

impl From<PlugFlowThermalRecord> for PlugFlowThermal {
    fn from(t: PlugFlowThermalRecord) -> Self {
        match t {
            PlugFlowThermalRecord::IsothermalK(t) => Self::Isothermal(t),
            PlugFlowThermalRecord::Adiabatic => Self::Adiabatic,
            PlugFlowThermalRecord::Coolant {
                temperature_k,
                u_w_m2_k,
            } => Self::Coolant {
                temperature_k,
                u_w_m2_k,
            },
        }
    }
}

impl From<BatchThermal> for BatchThermalRecord {
    fn from(t: BatchThermal) -> Self {
        match t {
            BatchThermal::Isothermal(t) => Self::IsothermalK(t),
            BatchThermal::Adiabatic => Self::Adiabatic,
        }
    }
}

impl From<BatchThermalRecord> for BatchThermal {
    fn from(t: BatchThermalRecord) -> Self {
        match t {
            BatchThermalRecord::IsothermalK(t) => Self::Isothermal(t),
            BatchThermalRecord::Adiabatic => Self::Adiabatic,
        }
    }
}

// ----- numbers for the finite-value check -----

pub(crate) fn kinetics_numbers(k: &KineticsRecord) -> Vec<Option<f64>> {
    let mut v = vec![k.catalyst_bulk_density_kg_m3];
    let arrhenius = |v: &mut Vec<Option<f64>>, a: &ArrheniusRecord| {
        v.extend([
            a.pre_exponential,
            Some(a.activation_energy_j_mol),
            Some(a.temperature_exponent),
            a.reference_temperature_k,
        ])
    };
    let orders =
        |v: &mut Vec<Option<f64>>, o: &[OrderRecord]| v.extend(o.iter().map(|o| Some(o.order)));
    let lhhw = |v: &mut Vec<Option<f64>>, t: &LhhwTermRecord| {
        v.extend([t.constant.a, t.constant.b, t.constant.c, t.constant.d].map(Some));
        orders(v, &t.orders);
    };
    for r in &k.reactions {
        v.extend(r.coefficients.iter().map(|c| Some(c.coefficient)));
        match &r.rate {
            RateLawRecord::PowerLaw { forward, reverse } => {
                for t in std::iter::once(forward).chain(reverse) {
                    arrhenius(&mut v, &t.constant);
                    orders(&mut v, &t.orders);
                }
            }
            RateLawRecord::Lhhw {
                kinetic,
                driving_forward,
                driving_reverse,
                adsorption,
                exponent,
            } => {
                arrhenius(&mut v, kinetic);
                for t in std::iter::once(driving_forward)
                    .chain(driving_reverse)
                    .chain(adsorption)
                {
                    lhhw(&mut v, t);
                }
                v.push(Some(*exponent));
            }
        }
    }
    v
}

pub(crate) fn equilibrium_numbers(
    reactions: &[EquilibriumReactionRecord],
    options: &EquilibriumOptionsRecord,
) -> Vec<Option<f64>> {
    let mut v = vec![Some(options.temperature_approach_k)];
    for r in reactions {
        v.extend(r.coefficients.iter().map(|c| Some(c.coefficient)));
        if let ConstantRecord::Correlation(c) = r.constant {
            v.extend([c.a, c.b, c.c, c.d].map(Some));
        }
    }
    v
}
