//! Project records of the columns group (#185–#187, #190); owned by the
//! columns work package. Add fields with `#[serde(default)]` behind format
//! version 4.

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::packing::{
    FilmModel, PackedCalculation, PackedSectionBasis, PackedSectionSpec, RateBasedSpec,
    SpecifiedPackedHydraulics, StageRounding,
};
use flowsheet_core::units::columns::packing_data::{
    FloodingCorrelation, MassTransferCorrelation, Packing, PressureDropCorrelation, WettingSurface,
};
use flowsheet_core::units::columns::{
    AbsorberSpec, BottomSpecification, ColumnInternals, ColumnPressureDrop, ComponentTarget,
    CurveBasis, CutPoint, DistillationColumnSpec, DistillationCondenser, DrawRate,
    FractionatorBottom, LiquidSideDraw, MAX_COLUMN_PORTS, PetroleumFractionatorSpec, PumpAround,
    PumpAroundCooling, ShortcutColumnSpec, ShortcutCondenser, ShortcutReflux, SideProduct,
    SideStripper, SideStripperBottom, TopSpecification, VaporSideDraw,
};
use serde::{Deserialize, Serialize};

use super::count_in;

fn one() -> u16 {
    1
}

fn unit_factor() -> f64 {
    1.0
}

fn is_unit_factor(v: &f64) -> bool {
    *v == 1.0
}

/// Packing of a rate-based section (#241).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum PackingTypeRecord {
    MetalPallRing25,
    MetalPallRing50,
}

/// Packing surface of the minimum-wetting rule (#241).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum WettingSurfaceRecord {
    UnglazedCeramic,
    OxidizedMetal,
    SurfaceTreatedMetal,
    GlazedCeramic,
    Glass,
    #[default]
    BrightMetal,
    PvcCpvc,
    Polypropylene,
    Fluoropolymer,
}

impl From<WettingSurface> for WettingSurfaceRecord {
    fn from(s: WettingSurface) -> Self {
        match s {
            WettingSurface::UnglazedCeramic => Self::UnglazedCeramic,
            WettingSurface::OxidizedMetal => Self::OxidizedMetal,
            WettingSurface::SurfaceTreatedMetal => Self::SurfaceTreatedMetal,
            WettingSurface::GlazedCeramic => Self::GlazedCeramic,
            WettingSurface::Glass => Self::Glass,
            WettingSurface::BrightMetal => Self::BrightMetal,
            WettingSurface::PvcCpvc => Self::PvcCpvc,
            WettingSurface::Polypropylene => Self::Polypropylene,
            WettingSurface::Fluoropolymer => Self::Fluoropolymer,
        }
    }
}

impl From<WettingSurfaceRecord> for WettingSurface {
    fn from(s: WettingSurfaceRecord) -> Self {
        match s {
            WettingSurfaceRecord::UnglazedCeramic => Self::UnglazedCeramic,
            WettingSurfaceRecord::OxidizedMetal => Self::OxidizedMetal,
            WettingSurfaceRecord::SurfaceTreatedMetal => Self::SurfaceTreatedMetal,
            WettingSurfaceRecord::GlazedCeramic => Self::GlazedCeramic,
            WettingSurfaceRecord::Glass => Self::Glass,
            WettingSurfaceRecord::BrightMetal => Self::BrightMetal,
            WettingSurfaceRecord::PvcCpvc => Self::PvcCpvc,
            WettingSurfaceRecord::Polypropylene => Self::Polypropylene,
            WettingSurfaceRecord::Fluoropolymer => Self::Fluoropolymer,
        }
    }
}

fn is_default_surface(s: &WettingSurfaceRecord) -> bool {
    *s == WettingSurfaceRecord::default()
}

/// Mass-transfer correlation of a rate-based section (#241).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum CorrelationRecord {
    #[default]
    Onda1968,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum PressureDropCorrelationRecord {
    #[default]
    Robbins1991,
    BilletSchultes1999,
    VendorRating,
}

impl From<PressureDropCorrelation> for PressureDropCorrelationRecord {
    fn from(value: PressureDropCorrelation) -> Self {
        match value {
            PressureDropCorrelation::Robbins1991 => Self::Robbins1991,
            PressureDropCorrelation::BilletSchultes1999 => Self::BilletSchultes1999,
            PressureDropCorrelation::VendorRating => Self::VendorRating,
        }
    }
}

impl From<PressureDropCorrelationRecord> for PressureDropCorrelation {
    fn from(value: PressureDropCorrelationRecord) -> Self {
        match value {
            PressureDropCorrelationRecord::Robbins1991 => Self::Robbins1991,
            PressureDropCorrelationRecord::BilletSchultes1999 => Self::BilletSchultes1999,
            PressureDropCorrelationRecord::VendorRating => Self::VendorRating,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum FloodingCorrelationRecord {
    #[default]
    KisterGill,
    VendorRating,
}

impl From<FloodingCorrelation> for FloodingCorrelationRecord {
    fn from(value: FloodingCorrelation) -> Self {
        match value {
            FloodingCorrelation::KisterGill => Self::KisterGill,
            FloodingCorrelation::VendorRating => Self::VendorRating,
        }
    }
}

impl From<FloodingCorrelationRecord> for FloodingCorrelation {
    fn from(value: FloodingCorrelationRecord) -> Self {
        match value {
            FloodingCorrelationRecord::KisterGill => Self::KisterGill,
            FloodingCorrelationRecord::VendorRating => Self::VendorRating,
        }
    }
}

/// Film model of a rate-based section (#241).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum FilmModelRecord {
    #[default]
    EffectiveDiffusivity,
}

/// Rate-based calculation mode of a packed absorber (#241), added behind
/// format version 4. A record without it is the equilibrium-stage mode.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RateBasedRecord {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    diameter_m: Option<f64>,
    segments: u16,
    packing: PackingTypeRecord,
    #[serde(default)]
    correlation: CorrelationRecord,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pressure_drop_correlation: Option<PressureDropCorrelationRecord>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    flooding_correlation: Option<FloodingCorrelationRecord>,
    #[serde(default)]
    film_model: FilmModelRecord,
    #[serde(default, skip_serializing_if = "is_false")]
    couple_pressure_drop: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    specified_hydraulics: Option<SpecifiedPackedHydraulicsRecord>,
    #[serde(default = "unit_factor", skip_serializing_if = "is_unit_factor")]
    interfacial_area_factor: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    critical_surface_tension_n_m: Option<f64>,
    #[serde(default, skip_serializing_if = "is_default_surface")]
    wetting_surface: WettingSurfaceRecord,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    gas_transport_neglected: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    liquid_transport_neglected: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SpecifiedPackedHydraulicsRecord {
    pressure_drop_pa: Option<f64>,
    flood_gas_flow_kg_s: Option<f64>,
}

impl From<&RateBasedSpec> for RateBasedRecord {
    fn from(r: &RateBasedSpec) -> Self {
        Self {
            diameter_m: r.diameter_m,
            segments: r.segments,
            packing: match r.packing {
                Packing::MetalPallRing25 => PackingTypeRecord::MetalPallRing25,
                Packing::MetalPallRing50 => PackingTypeRecord::MetalPallRing50,
            },
            correlation: match r.correlation {
                MassTransferCorrelation::Onda1968 => CorrelationRecord::Onda1968,
            },
            pressure_drop_correlation: Some(r.pressure_drop_correlation.into()),
            flooding_correlation: Some(r.flooding_correlation.into()),
            film_model: match r.film_model {
                FilmModel::EffectiveDiffusivity => FilmModelRecord::EffectiveDiffusivity,
            },
            interfacial_area_factor: r.interfacial_area_factor,
            couple_pressure_drop: r.couple_pressure_drop,
            specified_hydraulics: r.specified_hydraulics.map(|rating| {
                SpecifiedPackedHydraulicsRecord {
                    pressure_drop_pa: rating.pressure_drop_pa,
                    flood_gas_flow_kg_s: rating.flood_gas_flow_kg_s,
                }
            }),
            critical_surface_tension_n_m: r.critical_surface_tension_n_m,
            wetting_surface: r.wetting_surface.into(),
            gas_transport_neglected: r.gas_transport_neglected.clone(),
            liquid_transport_neglected: r.liquid_transport_neglected.clone(),
        }
    }
}

impl From<RateBasedRecord> for RateBasedSpec {
    fn from(r: RateBasedRecord) -> Self {
        Self {
            diameter_m: r.diameter_m,
            segments: r.segments,
            packing: match r.packing {
                PackingTypeRecord::MetalPallRing25 => Packing::MetalPallRing25,
                PackingTypeRecord::MetalPallRing50 => Packing::MetalPallRing50,
            },
            correlation: match r.correlation {
                CorrelationRecord::Onda1968 => MassTransferCorrelation::Onda1968,
            },
            pressure_drop_correlation: r.pressure_drop_correlation.map(Into::into).unwrap_or_else(
                || {
                    if r.specified_hydraulics.is_some() {
                        PressureDropCorrelation::VendorRating
                    } else {
                        PressureDropCorrelation::Robbins1991
                    }
                },
            ),
            flooding_correlation: r.flooding_correlation.map(Into::into).unwrap_or_else(|| {
                if r.specified_hydraulics.is_some() {
                    FloodingCorrelation::VendorRating
                } else {
                    FloodingCorrelation::KisterGill
                }
            }),
            film_model: match r.film_model {
                FilmModelRecord::EffectiveDiffusivity => FilmModel::EffectiveDiffusivity,
            },
            interfacial_area_factor: r.interfacial_area_factor,
            couple_pressure_drop: r.couple_pressure_drop
                || (r.pressure_drop_correlation.is_none() && r.specified_hydraulics.is_some()),
            specified_hydraulics: r
                .specified_hydraulics
                .map(|rating| SpecifiedPackedHydraulics {
                    pressure_drop_pa: rating.pressure_drop_pa,
                    flood_gas_flow_kg_s: rating.flood_gas_flow_kg_s,
                }),
            critical_surface_tension_n_m: r.critical_surface_tension_n_m,
            wetting_surface: r.wetting_surface.into(),
            gas_transport_neglected: r.gas_transport_neglected,
            liquid_transport_neglected: r.liquid_transport_neglected,
        }
    }
}

impl RateBasedRecord {
    fn numbers(&self) -> impl Iterator<Item = f64> {
        self.diameter_m
            .into_iter()
            .chain(self.critical_surface_tension_n_m)
            .chain([self.interfacial_area_factor])
            .chain(self.specified_hydraulics.into_iter().flat_map(|rating| {
                [rating.pressure_drop_pa, rating.flood_gas_flow_kg_s]
                    .into_iter()
                    .flatten()
            }))
    }
}

fn unset_reflux() -> ShortcutRefluxRecord {
    ShortcutRefluxRecord::MinimumMultiple(None)
}

/// Shortcut column reflux or stage-count choice (#185).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ShortcutRefluxRecord {
    Ratio(Option<f64>),
    MinimumMultiple(Option<f64>),
    /// Added behind format version 4.
    TheoreticalStages(Option<f64>),
}

impl From<ShortcutReflux> for ShortcutRefluxRecord {
    fn from(r: ShortcutReflux) -> Self {
        match r {
            ShortcutReflux::Ratio(v) => Self::Ratio(v),
            ShortcutReflux::MinimumMultiple(v) => Self::MinimumMultiple(v),
            ShortcutReflux::TheoreticalStages(v) => Self::TheoreticalStages(v),
        }
    }
}

impl From<ShortcutRefluxRecord> for ShortcutReflux {
    fn from(r: ShortcutRefluxRecord) -> Self {
        match r {
            ShortcutRefluxRecord::Ratio(v) => Self::Ratio(v),
            ShortcutRefluxRecord::MinimumMultiple(v) => Self::MinimumMultiple(v),
            ShortcutRefluxRecord::TheoreticalStages(v) => Self::TheoreticalStages(v),
        }
    }
}

/// Shortcut column condenser (#185), added behind format version 4. Records
/// without it load as a total condenser, and a total condenser is not written.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ShortcutCondenserRecord {
    #[default]
    Total,
    Partial,
}

impl ShortcutCondenserRecord {
    fn is_total(&self) -> bool {
        *self == Self::Total
    }
}

impl From<ShortcutCondenser> for ShortcutCondenserRecord {
    fn from(c: ShortcutCondenser) -> Self {
        match c {
            ShortcutCondenser::Total => Self::Total,
            ShortcutCondenser::Partial => Self::Partial,
        }
    }
}

impl From<ShortcutCondenserRecord> for ShortcutCondenser {
    fn from(c: ShortcutCondenserRecord) -> Self {
        match c {
            ShortcutCondenserRecord::Total => Self::Total,
            ShortcutCondenserRecord::Partial => Self::Partial,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum InternalsRecord {
    Trays(u16),
    Packed {
        hetp_m: Option<f64>,
        basis: PackingRecord,
    },
}

impl InternalsRecord {
    fn numbers(self) -> impl Iterator<Item = f64> {
        let (hetp, height) = match self {
            Self::Trays(_) => (None, None),
            Self::Packed { hetp_m, basis } => (
                hetp_m,
                match basis {
                    PackingRecord::Height { height_m, .. } => height_m,
                    PackingRecord::TheoreticalStages(_) => None,
                },
            ),
        };
        [hetp, height].into_iter().flatten()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum PackingRecord {
    Height {
        height_m: Option<f64>,
        rounding: RoundingRecord,
    },
    TheoreticalStages(u16),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum RoundingRecord {
    Nearest,
    Down,
    Up,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct LiquidDrawRecord {
    stage: u16,
    molar_flow_mol_s: Option<f64>,
}

fn is_zero(v: &f64) -> bool {
    *v == 0.0
}

fn is_zero_count(v: &u16) -> bool {
    *v == 0
}

fn is_one(v: &u16) -> bool {
    *v == 1
}

fn is_false(v: &bool) -> bool {
    !*v
}

/// Flow specification of a petroleum fractionator side product (#190).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum DrawRateRecord {
    MolarFlowMolS(Option<f64>),
    StandardLiquidVolumeFlowM3S(Option<f64>),
    /// TBP or ASTM D86 cut point (#190, added behind format version 4).
    CutPoint(CutPointRecord),
}

impl DrawRateRecord {
    fn numbers(self) -> Vec<f64> {
        match self {
            Self::MolarFlowMolS(v) | Self::StandardLiquidVolumeFlowM3S(v) => {
                v.into_iter().collect()
            }
            Self::CutPoint(c) => [Some(c.volume_percent), c.temperature_k, c.estimate_mol_s]
                .into_iter()
                .flatten()
                .collect(),
        }
    }
}

/// Distillation curve of a cut-point specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum CurveBasisRecord {
    Tbp,
    AstmD86,
}

/// Cut-point specification of a side product (#190).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CutPointRecord {
    basis: CurveBasisRecord,
    volume_percent: f64,
    temperature_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    estimate_mol_s: Option<f64>,
}

impl From<CutPoint> for CutPointRecord {
    fn from(c: CutPoint) -> Self {
        Self {
            basis: match c.basis {
                CurveBasis::Tbp => CurveBasisRecord::Tbp,
                CurveBasis::AstmD86 => CurveBasisRecord::AstmD86,
            },
            volume_percent: c.volume_percent,
            temperature_k: c.temperature_k,
            estimate_mol_s: c.estimate_mol_s,
        }
    }
}

impl From<CutPointRecord> for CutPoint {
    fn from(c: CutPointRecord) -> Self {
        Self {
            basis: match c.basis {
                CurveBasisRecord::Tbp => CurveBasis::Tbp,
                CurveBasisRecord::AstmD86 => CurveBasis::AstmD86,
            },
            volume_percent: c.volume_percent,
            temperature_k: c.temperature_k,
            estimate_mol_s: c.estimate_mol_s,
        }
    }
}

/// Cooling of a petroleum fractionator pump-around (#190).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum PumpAroundCoolingRecord {
    ReturnTemperatureK(Option<f64>),
    DutyW(Option<f64>),
}

/// Pump-around of a petroleum fractionator (#190), behind format version 4.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct PumpAroundRecord {
    draw_stage: u16,
    return_stage: u16,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    molar_flow_mol_s: Option<f64>,
    cooling: PumpAroundCoolingRecord,
}

impl PumpAroundRecord {
    fn numbers(&self) -> impl Iterator<Item = f64> {
        let cooling = match self.cooling {
            PumpAroundCoolingRecord::ReturnTemperatureK(v) | PumpAroundCoolingRecord::DutyW(v) => v,
        };
        self.molar_flow_mol_s.into_iter().chain(cooling)
    }
}

impl From<PumpAround> for PumpAroundRecord {
    fn from(p: PumpAround) -> Self {
        Self {
            draw_stage: p.draw_stage,
            return_stage: p.return_stage,
            molar_flow_mol_s: p.molar_flow_mol_s,
            cooling: match p.cooling {
                PumpAroundCooling::ReturnTemperature(v) => {
                    PumpAroundCoolingRecord::ReturnTemperatureK(v)
                }
                PumpAroundCooling::Duty(v) => PumpAroundCoolingRecord::DutyW(v),
            },
        }
    }
}

impl From<PumpAroundRecord> for PumpAround {
    fn from(p: PumpAroundRecord) -> Self {
        Self {
            draw_stage: p.draw_stage,
            return_stage: p.return_stage,
            molar_flow_mol_s: p.molar_flow_mol_s,
            cooling: match p.cooling {
                PumpAroundCoolingRecord::ReturnTemperatureK(v) => {
                    PumpAroundCooling::ReturnTemperature(v)
                }
                PumpAroundCoolingRecord::DutyW(v) => PumpAroundCooling::Duty(v),
            },
        }
    }
}

/// Stripping agent of a petroleum fractionator side stripper (#190).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum SideStripperBottomRecord {
    Steam,
    ReboilerDutyW(Option<f64>),
    BoilupRatio(Option<f64>),
}

/// Side stripper of a petroleum fractionator (#190), behind format version 4.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SideStripperRecord {
    product: u16,
    stages: u16,
    return_stage: u16,
    bottom: SideStripperBottomRecord,
}

impl SideStripperRecord {
    fn number(&self) -> Option<f64> {
        match self.bottom {
            SideStripperBottomRecord::Steam => None,
            SideStripperBottomRecord::ReboilerDutyW(v)
            | SideStripperBottomRecord::BoilupRatio(v) => v,
        }
    }
}

impl From<SideStripper> for SideStripperRecord {
    fn from(s: SideStripper) -> Self {
        Self {
            product: s.product,
            stages: s.stages,
            return_stage: s.return_stage,
            bottom: match s.bottom {
                SideStripperBottom::Steam => SideStripperBottomRecord::Steam,
                SideStripperBottom::ReboilerDuty(v) => SideStripperBottomRecord::ReboilerDutyW(v),
                SideStripperBottom::BoilupRatio(v) => SideStripperBottomRecord::BoilupRatio(v),
            },
        }
    }
}

impl From<SideStripperRecord> for SideStripper {
    fn from(s: SideStripperRecord) -> Self {
        Self {
            product: s.product,
            stages: s.stages,
            return_stage: s.return_stage,
            bottom: match s.bottom {
                SideStripperBottomRecord::Steam => SideStripperBottom::Steam,
                SideStripperBottomRecord::ReboilerDutyW(v) => SideStripperBottom::ReboilerDuty(v),
                SideStripperBottomRecord::BoilupRatio(v) => SideStripperBottom::BoilupRatio(v),
            },
        }
    }
}

/// Side product of a petroleum fractionator (#190).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SideProductRecord {
    stage: u16,
    rate: DrawRateRecord,
}

impl From<SideProduct> for SideProductRecord {
    fn from(p: SideProduct) -> Self {
        Self {
            stage: p.stage,
            rate: match p.rate {
                DrawRate::Molar(v) => DrawRateRecord::MolarFlowMolS(v),
                DrawRate::StandardLiquidVolume(v) => DrawRateRecord::StandardLiquidVolumeFlowM3S(v),
                DrawRate::CutPoint(c) => DrawRateRecord::CutPoint(c.into()),
            },
        }
    }
}

impl From<SideProductRecord> for SideProduct {
    fn from(p: SideProductRecord) -> Self {
        Self {
            stage: p.stage,
            rate: match p.rate {
                DrawRateRecord::MolarFlowMolS(v) => DrawRate::Molar(v),
                DrawRateRecord::StandardLiquidVolumeFlowM3S(v) => DrawRate::StandardLiquidVolume(v),
                DrawRateRecord::CutPoint(c) => DrawRate::CutPoint(c.into()),
            },
        }
    }
}

/// Rigorous column condenser (#186), added behind format version 4. Records
/// without it load as a total condenser, and a total condenser is not written.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum DistillationCondenserRecord {
    #[default]
    Total,
    PartialVapor,
    PartialVaporLiquid {
        liquid_distillate_mol_s: Option<f64>,
    },
}

impl DistillationCondenserRecord {
    fn is_total(&self) -> bool {
        *self == Self::Total
    }

    fn liquid_distillate(self) -> Option<f64> {
        match self {
            Self::PartialVaporLiquid {
                liquid_distillate_mol_s,
            } => liquid_distillate_mol_s,
            _ => None,
        }
    }
}

impl From<DistillationCondenser> for DistillationCondenserRecord {
    fn from(c: DistillationCondenser) -> Self {
        match c {
            DistillationCondenser::Total => Self::Total,
            DistillationCondenser::PartialVapor => Self::PartialVapor,
            DistillationCondenser::PartialVaporLiquid {
                liquid_distillate_mol_s,
            } => Self::PartialVaporLiquid {
                liquid_distillate_mol_s,
            },
        }
    }
}

impl From<DistillationCondenserRecord> for DistillationCondenser {
    fn from(c: DistillationCondenserRecord) -> Self {
        match c {
            DistillationCondenserRecord::Total => Self::Total,
            DistillationCondenserRecord::PartialVapor => Self::PartialVapor,
            DistillationCondenserRecord::PartialVaporLiquid {
                liquid_distillate_mol_s,
            } => Self::PartialVaporLiquid {
                liquid_distillate_mol_s,
            },
        }
    }
}

/// Stage pressure drop below the first tray (#186), added behind format
/// version 4. Records without it load as zero drop (uniform pressure).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum PressureDropRecord {
    PerStagePa(f64),
    ColumnPa(f64),
}

impl Default for PressureDropRecord {
    fn default() -> Self {
        Self::PerStagePa(0.0)
    }
}

impl PressureDropRecord {
    fn is_zero(&self) -> bool {
        *self == Self::PerStagePa(0.0)
    }

    fn value(self) -> f64 {
        let (Self::PerStagePa(drop) | Self::ColumnPa(drop)) = self;
        drop
    }
}

impl From<ColumnPressureDrop> for PressureDropRecord {
    fn from(d: ColumnPressureDrop) -> Self {
        match d {
            ColumnPressureDrop::PerStage(v) => Self::PerStagePa(v),
            ColumnPressureDrop::Column(v) => Self::ColumnPa(v),
        }
    }
}

impl From<PressureDropRecord> for ColumnPressureDrop {
    fn from(d: PressureDropRecord) -> Self {
        match d {
            PressureDropRecord::PerStagePa(v) => Self::PerStage(v),
            PressureDropRecord::ColumnPa(v) => Self::Column(v),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ComponentTargetRecord {
    component: Option<String>,
    value: Option<f64>,
}

impl From<&ComponentTarget> for ComponentTargetRecord {
    fn from(t: &ComponentTarget) -> Self {
        Self {
            component: t.component.clone(),
            value: t.value,
        }
    }
}

impl From<ComponentTargetRecord> for ComponentTarget {
    fn from(t: ComponentTargetRecord) -> Self {
        Self {
            component: t.component,
            value: t.value,
        }
    }
}

/// First column specification other than the reflux ratio, which keeps the
/// original `reflux_ratio` field (#186, added behind format version 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum TopSpecificationRecord {
    RefluxRateMolS(Option<f64>),
    DistillateRateMolS(Option<f64>),
    CondenserDutyW(Option<f64>),
    DistillatePurity(ComponentTargetRecord),
    DistillateRecovery(ComponentTargetRecord),
}

/// Second column specification other than the distillate fraction, which
/// keeps the original `distillate_fraction` field (#186, format version 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum BottomSpecificationRecord {
    BottomsRateMolS(Option<f64>),
    ReboilerDutyW(Option<f64>),
    BoilupRatio(Option<f64>),
    BottomsPurity(ComponentTargetRecord),
    BottomsRecovery(ComponentTargetRecord),
}

impl TopSpecificationRecord {
    fn number(&self) -> Option<f64> {
        match self {
            Self::RefluxRateMolS(v) | Self::DistillateRateMolS(v) | Self::CondenserDutyW(v) => *v,
            Self::DistillatePurity(t) | Self::DistillateRecovery(t) => t.value,
        }
    }
}

impl BottomSpecificationRecord {
    fn number(&self) -> Option<f64> {
        match self {
            Self::BottomsRateMolS(v) | Self::ReboilerDutyW(v) | Self::BoilupRatio(v) => *v,
            Self::BottomsPurity(t) | Self::BottomsRecovery(t) => t.value,
        }
    }
}

/// Splits a top specification into the original reflux-ratio field and the
/// specification record.
fn top_record(top: &TopSpecification) -> (Option<f64>, Option<TopSpecificationRecord>) {
    use TopSpecificationRecord as R;
    match top {
        TopSpecification::RefluxRatio(v) => (*v, None),
        TopSpecification::RefluxRate(v) => (None, Some(R::RefluxRateMolS(*v))),
        TopSpecification::DistillateRate(v) => (None, Some(R::DistillateRateMolS(*v))),
        TopSpecification::CondenserDuty(v) => (None, Some(R::CondenserDutyW(*v))),
        TopSpecification::DistillatePurity(t) => (None, Some(R::DistillatePurity(t.into()))),
        TopSpecification::DistillateRecovery(t) => (None, Some(R::DistillateRecovery(t.into()))),
    }
}

fn bottom_record(bottom: &BottomSpecification) -> (Option<f64>, Option<BottomSpecificationRecord>) {
    use BottomSpecificationRecord as R;
    match bottom {
        BottomSpecification::DistillateFraction(v) => (*v, None),
        BottomSpecification::BottomsRate(v) => (None, Some(R::BottomsRateMolS(*v))),
        BottomSpecification::ReboilerDuty(v) => (None, Some(R::ReboilerDutyW(*v))),
        BottomSpecification::BoilupRatio(v) => (None, Some(R::BoilupRatio(*v))),
        BottomSpecification::BottomsPurity(t) => (None, Some(R::BottomsPurity(t.into()))),
        BottomSpecification::BottomsRecovery(t) => (None, Some(R::BottomsRecovery(t.into()))),
    }
}

fn top_spec(
    reflux_ratio: Option<f64>,
    record: Option<TopSpecificationRecord>,
) -> Result<TopSpecification, String> {
    use TopSpecificationRecord as R;
    Ok(match record {
        None => TopSpecification::RefluxRatio(reflux_ratio),
        Some(_) if reflux_ratio.is_some() => {
            return Err("column has both a reflux ratio and another top specification".into());
        }
        Some(R::RefluxRateMolS(v)) => TopSpecification::RefluxRate(v),
        Some(R::DistillateRateMolS(v)) => TopSpecification::DistillateRate(v),
        Some(R::CondenserDutyW(v)) => TopSpecification::CondenserDuty(v),
        Some(R::DistillatePurity(t)) => TopSpecification::DistillatePurity(t.into()),
        Some(R::DistillateRecovery(t)) => TopSpecification::DistillateRecovery(t.into()),
    })
}

fn bottom_spec(
    distillate_fraction: Option<f64>,
    record: Option<BottomSpecificationRecord>,
) -> Result<BottomSpecification, String> {
    use BottomSpecificationRecord as R;
    Ok(match record {
        None => BottomSpecification::DistillateFraction(distillate_fraction),
        Some(_) if distillate_fraction.is_some() => {
            return Err(
                "column has both a distillate fraction and another bottom specification".into(),
            );
        }
        Some(R::BottomsRateMolS(v)) => BottomSpecification::BottomsRate(v),
        Some(R::ReboilerDutyW(v)) => BottomSpecification::ReboilerDuty(v),
        Some(R::BoilupRatio(v)) => BottomSpecification::BoilupRatio(v),
        Some(R::BottomsPurity(t)) => BottomSpecification::BottomsPurity(t.into()),
        Some(R::BottomsRecovery(t)) => BottomSpecification::BottomsRecovery(t.into()),
    })
}

fn draw_records(draws: impl Iterator<Item = (u16, Option<f64>)>) -> Vec<LiquidDrawRecord> {
    draws
        .map(|(stage, molar_flow_mol_s)| LiquidDrawRecord {
            stage,
            molar_flow_mol_s,
        })
        .collect()
}

impl From<ColumnInternals> for InternalsRecord {
    fn from(value: ColumnInternals) -> Self {
        match value {
            ColumnInternals::Trays(count) => Self::Trays(count),
            ColumnInternals::Packed(spec) => Self::Packed {
                hetp_m: spec.hetp_m,
                basis: match spec.basis {
                    PackedSectionBasis::TheoreticalStages(count) => {
                        PackingRecord::TheoreticalStages(count)
                    }
                    PackedSectionBasis::PackedHeight { height_m, rounding } => {
                        PackingRecord::Height {
                            height_m,
                            rounding: match rounding {
                                StageRounding::Nearest => RoundingRecord::Nearest,
                                StageRounding::Down => RoundingRecord::Down,
                                StageRounding::Up => RoundingRecord::Up,
                            },
                        }
                    }
                },
            },
        }
    }
}

impl From<InternalsRecord> for ColumnInternals {
    fn from(value: InternalsRecord) -> Self {
        match value {
            InternalsRecord::Trays(count) => Self::Trays(count),
            InternalsRecord::Packed { hetp_m, basis } => Self::Packed(PackedSectionSpec {
                hetp_m,
                basis: match basis {
                    PackingRecord::TheoreticalStages(count) => {
                        PackedSectionBasis::TheoreticalStages(count)
                    }
                    PackingRecord::Height { height_m, rounding } => {
                        PackedSectionBasis::PackedHeight {
                            height_m,
                            rounding: match rounding {
                                RoundingRecord::Nearest => StageRounding::Nearest,
                                RoundingRecord::Down => StageRounding::Down,
                                RoundingRecord::Up => StageRounding::Up,
                            },
                        }
                    }
                },
            }),
        }
    }
}

/// Column specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ColumnRecord {
    /// Shortcut fields were added behind format version 4; the earlier empty
    /// record loads with every value unset.
    ShortcutColumn {
        #[serde(default)]
        light_key: Option<String>,
        #[serde(default)]
        heavy_key: Option<String>,
        #[serde(default)]
        light_key_recovery: Option<f64>,
        #[serde(default)]
        heavy_key_recovery: Option<f64>,
        #[serde(default)]
        pressure_pa: Option<f64>,
        #[serde(default = "unset_reflux")]
        reflux: ShortcutRefluxRecord,
        #[serde(default, skip_serializing_if = "ShortcutCondenserRecord::is_total")]
        condenser: ShortcutCondenserRecord,
    },
    DistillationColumn {
        #[serde(default = "one")]
        feeds: u16,
        #[serde(default)]
        side_draws: u16,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        feed_stages: Vec<u16>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        liquid_draws: Vec<LiquidDrawRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        internals: Option<InternalsRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pressure_pa: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reflux_ratio: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        distillate_fraction: Option<f64>,
        // The fields below were added behind format version 4 (#186
        // specification set); records without them load as the original
        // total-condenser, uniform-pressure column.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        top_specification: Option<TopSpecificationRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        bottom_specification: Option<BottomSpecificationRecord>,
        #[serde(default, skip_serializing_if = "DistillationCondenserRecord::is_total")]
        condenser: DistillationCondenserRecord,
        #[serde(default, skip_serializing_if = "is_zero")]
        condenser_pressure_drop_pa: f64,
        #[serde(default, skip_serializing_if = "PressureDropRecord::is_zero")]
        pressure_drop: PressureDropRecord,
        #[serde(default, skip_serializing_if = "is_zero_count")]
        vapor_side_draws: u16,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        vapor_draws: Vec<LiquidDrawRecord>,
    },
    Absorber {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pressure_pa: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        internals: Option<InternalsRecord>,
        /// Rate-based mode (#241), behind format version 4.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        rate_based: Option<RateBasedRecord>,
    },
    /// The drawn-only fractionator stored only `side_draws`; the other
    /// fields were added behind format version 4 (#190) and load as one
    /// feed, ten trays, a total condenser and a reboiler with unset values.
    PetroleumFractionator {
        #[serde(default)]
        side_draws: u16,
        #[serde(default = "one", skip_serializing_if = "is_one")]
        feeds: u16,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        feed_stages: Vec<u16>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        products: Vec<SideProductRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        internals: Option<InternalsRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pressure_pa: Option<f64>,
        #[serde(default, skip_serializing_if = "is_zero")]
        condenser_pressure_drop_pa: f64,
        #[serde(default, skip_serializing_if = "PressureDropRecord::is_zero")]
        pressure_drop: PressureDropRecord,
        #[serde(default, skip_serializing_if = "DistillationCondenserRecord::is_total")]
        condenser: DistillationCondenserRecord,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reflux_ratio: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        top_specification: Option<TopSpecificationRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        distillate_fraction: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        bottom_specification: Option<BottomSpecificationRecord>,
        /// Stripping steam instead of a reboiler; no bottom specification.
        #[serde(default, skip_serializing_if = "is_false")]
        stripping_steam: bool,
        /// Pump-arounds; added behind format version 4 after the first
        /// fractionator records, absent means none.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        pump_arounds: Vec<PumpAroundRecord>,
        /// Side strippers; added like `pump_arounds`.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        side_strippers: Vec<SideStripperRecord>,
    },
}

impl ColumnRecord {
    /// Every number set in the record (counts are integers and always finite).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        match self {
            Self::ShortcutColumn {
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux,
                ..
            } => {
                let reflux = match reflux {
                    ShortcutRefluxRecord::Ratio(v)
                    | ShortcutRefluxRecord::MinimumMultiple(v)
                    | ShortcutRefluxRecord::TheoreticalStages(v) => *v,
                };
                [
                    *light_key_recovery,
                    *heavy_key_recovery,
                    *pressure_pa,
                    reflux,
                ]
                .into_iter()
                .flatten()
                .collect()
            }
            Self::DistillationColumn {
                liquid_draws,
                internals,
                pressure_pa,
                reflux_ratio,
                distillate_fraction,
                top_specification,
                bottom_specification,
                condenser,
                condenser_pressure_drop_pa,
                pressure_drop,
                vapor_draws,
                ..
            } => {
                let liquid_distillate = match condenser {
                    DistillationCondenserRecord::PartialVaporLiquid {
                        liquid_distillate_mol_s,
                    } => *liquid_distillate_mol_s,
                    _ => None,
                };
                let (PressureDropRecord::PerStagePa(drop) | PressureDropRecord::ColumnPa(drop)) =
                    *pressure_drop;
                let mut values: Vec<f64> = [
                    *pressure_pa,
                    *reflux_ratio,
                    *distillate_fraction,
                    top_specification.as_ref().and_then(|t| t.number()),
                    bottom_specification.as_ref().and_then(|b| b.number()),
                    liquid_distillate,
                    Some(*condenser_pressure_drop_pa),
                    Some(drop),
                ]
                .into_iter()
                .flatten()
                .collect();
                values.extend(liquid_draws.iter().filter_map(|draw| draw.molar_flow_mol_s));
                values.extend(vapor_draws.iter().filter_map(|draw| draw.molar_flow_mol_s));
                values.extend(internals.iter().flat_map(|i| i.numbers()));
                values
            }
            Self::Absorber {
                pressure_pa,
                internals,
                rate_based,
            } => pressure_pa
                .iter()
                .copied()
                .chain(internals.iter().flat_map(|i| i.numbers()))
                .chain(rate_based.iter().flat_map(|r| r.numbers()))
                .collect(),
            Self::PetroleumFractionator {
                products,
                internals,
                pressure_pa,
                condenser_pressure_drop_pa,
                pressure_drop,
                condenser,
                reflux_ratio,
                top_specification,
                distillate_fraction,
                bottom_specification,
                pump_arounds,
                side_strippers,
                ..
            } => {
                let mut values: Vec<f64> = [
                    *pressure_pa,
                    *reflux_ratio,
                    *distillate_fraction,
                    top_specification.as_ref().and_then(|t| t.number()),
                    bottom_specification.as_ref().and_then(|b| b.number()),
                    condenser.liquid_distillate(),
                    Some(*condenser_pressure_drop_pa),
                    Some(pressure_drop.value()),
                ]
                .into_iter()
                .flatten()
                .collect();
                values.extend(products.iter().flat_map(|p| p.rate.numbers()));
                values.extend(internals.iter().flat_map(|i| i.numbers()));
                values.extend(pump_arounds.iter().flat_map(PumpAroundRecord::numbers));
                values.extend(side_strippers.iter().filter_map(SideStripperRecord::number));
                values
            }
        }
    }
}

impl From<&ColumnSpec> for ColumnRecord {
    fn from(s: &ColumnSpec) -> Self {
        match s {
            ColumnSpec::Shortcut(x) => Self::ShortcutColumn {
                light_key: x.light_key.clone(),
                heavy_key: x.heavy_key.clone(),
                light_key_recovery: x.light_key_recovery,
                heavy_key_recovery: x.heavy_key_recovery,
                pressure_pa: x.pressure_pa,
                reflux: x.reflux.into(),
                condenser: x.condenser.into(),
            },
            ColumnSpec::Distillation(x) => {
                let (reflux_ratio, top_specification) = top_record(&x.top);
                let (distillate_fraction, bottom_specification) = bottom_record(&x.bottom);
                Self::DistillationColumn {
                    feeds: x.feeds,
                    side_draws: x.side_draws,
                    feed_stages: x.feed_stages.clone(),
                    liquid_draws: draw_records(
                        x.liquid_draws.iter().map(|d| (d.stage, d.molar_flow_mol_s)),
                    ),
                    internals: (x.internals != DistillationColumnSpec::default().internals)
                        .then(|| x.internals.into()),
                    pressure_pa: x.pressure_pa,
                    reflux_ratio,
                    distillate_fraction,
                    top_specification,
                    bottom_specification,
                    condenser: match x.condenser {
                        DistillationCondenser::Total => DistillationCondenserRecord::Total,
                        DistillationCondenser::PartialVapor => {
                            DistillationCondenserRecord::PartialVapor
                        }
                        DistillationCondenser::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        } => DistillationCondenserRecord::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        },
                    },
                    condenser_pressure_drop_pa: x.condenser_pressure_drop_pa,
                    pressure_drop: match x.pressure_drop {
                        ColumnPressureDrop::PerStage(v) => PressureDropRecord::PerStagePa(v),
                        ColumnPressureDrop::Column(v) => PressureDropRecord::ColumnPa(v),
                    },
                    vapor_side_draws: x.vapor_side_draws,
                    vapor_draws: draw_records(
                        x.vapor_draws.iter().map(|d| (d.stage, d.molar_flow_mol_s)),
                    ),
                }
            }
            ColumnSpec::Absorber(x) => Self::Absorber {
                pressure_pa: x.pressure_pa,
                internals: (x.internals != AbsorberSpec::default().internals)
                    .then(|| x.internals.into()),
                rate_based: match &x.calculation {
                    PackedCalculation::Equilibrium => None,
                    PackedCalculation::RateBased(r) => Some(r.into()),
                },
            },
            ColumnSpec::PetroleumFractionator(x) => {
                let (reflux_ratio, top_specification) = top_record(&x.top);
                let (distillate_fraction, bottom_specification) = match &x.bottom {
                    FractionatorBottom::Reboiler(b) => bottom_record(b),
                    FractionatorBottom::StrippingSteam => (None, None),
                };
                Self::PetroleumFractionator {
                    side_draws: x.side_draws,
                    feeds: x.feeds,
                    feed_stages: x.feed_stages.clone(),
                    products: x.products.iter().map(|&p| p.into()).collect(),
                    internals: (x.internals != PetroleumFractionatorSpec::default().internals)
                        .then(|| x.internals.into()),
                    pressure_pa: x.pressure_pa,
                    condenser_pressure_drop_pa: x.condenser_pressure_drop_pa,
                    pressure_drop: x.pressure_drop.into(),
                    condenser: x.condenser.into(),
                    reflux_ratio,
                    top_specification,
                    distillate_fraction,
                    bottom_specification,
                    stripping_steam: x.steam(),
                    pump_arounds: x.pump_arounds.iter().map(|&p| p.into()).collect(),
                    side_strippers: x.side_strippers.iter().map(|&s| s.into()).collect(),
                }
            }
        }
    }
}

impl TryFrom<ColumnRecord> for ColumnSpec {
    type Error = String;

    fn try_from(r: ColumnRecord) -> Result<Self, String> {
        Ok(match r {
            ColumnRecord::ShortcutColumn {
                light_key,
                heavy_key,
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux,
                condenser,
            } => Self::Shortcut(ShortcutColumnSpec {
                light_key,
                heavy_key,
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux: reflux.into(),
                condenser: condenser.into(),
            }),
            ColumnRecord::DistillationColumn {
                feeds,
                side_draws,
                feed_stages,
                liquid_draws,
                internals,
                pressure_pa,
                reflux_ratio,
                distillate_fraction,
                top_specification,
                bottom_specification,
                condenser,
                condenser_pressure_drop_pa,
                pressure_drop,
                vapor_side_draws,
                vapor_draws,
            } => {
                if feed_stages.len() > usize::from(feeds)
                    || liquid_draws.len() > usize::from(side_draws)
                    || vapor_draws.len() > usize::from(vapor_side_draws)
                {
                    return Err("column has more feed or draw specifications than ports".into());
                }
                Self::Distillation(DistillationColumnSpec {
                    feeds: count_in(feeds, "feeds", 1..=MAX_COLUMN_PORTS)?,
                    side_draws: count_in(side_draws, "side_draws", 0..=MAX_COLUMN_PORTS)?,
                    feed_stages,
                    liquid_draws: liquid_draws
                        .into_iter()
                        .map(|d| LiquidSideDraw {
                            stage: d.stage,
                            molar_flow_mol_s: d.molar_flow_mol_s,
                        })
                        .collect(),
                    internals: internals
                        .map_or_else(|| DistillationColumnSpec::default().internals, Into::into),
                    vapor_side_draws: count_in(
                        vapor_side_draws,
                        "vapor_side_draws",
                        0..=MAX_COLUMN_PORTS,
                    )?,
                    vapor_draws: vapor_draws
                        .into_iter()
                        .map(|d| VaporSideDraw {
                            stage: d.stage,
                            molar_flow_mol_s: d.molar_flow_mol_s,
                        })
                        .collect(),
                    pressure_pa,
                    condenser_pressure_drop_pa,
                    pressure_drop: match pressure_drop {
                        PressureDropRecord::PerStagePa(v) => ColumnPressureDrop::PerStage(v),
                        PressureDropRecord::ColumnPa(v) => ColumnPressureDrop::Column(v),
                    },
                    condenser: match condenser {
                        DistillationCondenserRecord::Total => DistillationCondenser::Total,
                        DistillationCondenserRecord::PartialVapor => {
                            DistillationCondenser::PartialVapor
                        }
                        DistillationCondenserRecord::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        } => DistillationCondenser::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        },
                    },
                    top: top_spec(reflux_ratio, top_specification)?,
                    bottom: bottom_spec(distillate_fraction, bottom_specification)?,
                })
            }
            ColumnRecord::Absorber {
                pressure_pa,
                internals,
                rate_based,
            } => Self::Absorber(AbsorberSpec {
                pressure_pa,
                internals: internals.map_or_else(|| AbsorberSpec::default().internals, Into::into),
                calculation: rate_based.map_or(PackedCalculation::Equilibrium, |r| {
                    PackedCalculation::RateBased(r.into())
                }),
            }),
            ColumnRecord::PetroleumFractionator {
                side_draws,
                feeds,
                feed_stages,
                products,
                internals,
                pressure_pa,
                condenser_pressure_drop_pa,
                pressure_drop,
                condenser,
                reflux_ratio,
                top_specification,
                distillate_fraction,
                bottom_specification,
                stripping_steam,
                pump_arounds,
                side_strippers,
            } => {
                if feed_stages.len() > usize::from(feeds)
                    || products.len() > usize::from(side_draws)
                    || pump_arounds.len() > usize::from(MAX_COLUMN_PORTS)
                    || side_strippers.len() > usize::from(side_draws)
                {
                    return Err(
                        "petroleum fractionator has more feed or product specifications than ports"
                            .into(),
                    );
                }
                let bottom = if stripping_steam {
                    if distillate_fraction.is_some() || bottom_specification.is_some() {
                        return Err(
                            "petroleum fractionator has both stripping steam and a reboiler specification"
                                .into(),
                        );
                    }
                    FractionatorBottom::StrippingSteam
                } else {
                    FractionatorBottom::Reboiler(bottom_spec(
                        distillate_fraction,
                        bottom_specification,
                    )?)
                };
                Self::PetroleumFractionator(PetroleumFractionatorSpec {
                    side_draws: count_in(side_draws, "side_draws", 0..=MAX_COLUMN_PORTS)?,
                    feeds: count_in(feeds, "feeds", 1..=MAX_COLUMN_PORTS)?,
                    feed_stages,
                    products: products.into_iter().map(Into::into).collect(),
                    internals: internals.map_or_else(
                        || PetroleumFractionatorSpec::default().internals,
                        Into::into,
                    ),
                    pressure_pa,
                    condenser_pressure_drop_pa,
                    pressure_drop: pressure_drop.into(),
                    condenser: condenser.into(),
                    top: top_spec(reflux_ratio, top_specification)?,
                    bottom,
                    pump_arounds: pump_arounds.into_iter().map(Into::into).collect(),
                    side_strippers: side_strippers.into_iter().map(Into::into).collect(),
                })
            }
        })
    }
}
