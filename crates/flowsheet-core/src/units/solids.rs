//! Solids group (#198–#211): cyclone, hydrocyclone, baghouse, rotary
//! filter, simple solids separator, crystallizer, centrifuge, electrostatic
//! precipitator, venturi scrubber, crusher, screen, dryer, granulator and
//! solids washer.
//!
//! Owned by the solids work package (`docs/a17-architecture.md`). Every kind
//! needs a solids substream / particle-size-distribution model, which this
//! group designs; `MaterialStream` carries no solids yet.

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info};
use crate::unit::{CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, UnitKind};

const fn split(a: &'static str, b: &'static str) -> [PortSpec; 3] {
    [
        PortSpec::inlet("in"),
        PortSpec::outlet(a),
        PortSpec::outlet(b),
    ]
}

const fn contact(
    in_a: &'static str,
    in_b: &'static str,
    out_a: &'static str,
    out_b: &'static str,
) -> [PortSpec; 4] {
    [
        PortSpec::inlet(in_a),
        PortSpec::inlet(in_b),
        PortSpec::outlet(out_a),
        PortSpec::outlet(out_b),
    ]
}

/// Metadata of a solids kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    const GAS_SOLIDS: [PortSpec; 3] = split("gas", "solids");
    const HYDROCYCLONE: [PortSpec; 3] = split("overflow", "underflow");
    const FILTER: [PortSpec; 3] = split("filtrate", "cake");
    const SEPARATOR: [PortSpec; 3] = split("fluid", "solids");
    const CRYSTALLIZER: [PortSpec; 3] = split("vapor", "slurry");
    const CENTRIFUGE: [PortSpec; 3] = split("liquid", "solids");
    const SCREEN: [PortSpec; 3] = split("oversize", "undersize");
    const SINGLE: [PortSpec; 2] = [PortSpec::inlet("in"), PortSpec::outlet("out")];
    const VENTURI: [PortSpec; 4] = contact("gas in", "liquid in", "gas out", "slurry");
    const DRYER: [PortSpec; 4] = contact("wet solids", "gas in", "dry solids", "gas out");
    const WASHER: [PortSpec; 4] = contact("solids in", "wash in", "solids out", "wash out");
    const DUTY: [PortSpec; 1] = [PortSpec::energy_inlet("duty")];
    const POWER: [PortSpec; 1] = [PortSpec::energy_inlet("power")];
    let solids = |label, tag, ports, energy, issue| {
        kind_info(label, tag, PaletteCategory::Solids, ports, energy, issue)
    };
    match kind {
        UnitKind::Cyclone => solids("Cyclone", "CY", &GAS_SOLIDS, &[], 198),
        UnitKind::Hydrocyclone => solids("Hydrocyclone", "HC", &HYDROCYCLONE, &[], 199),
        UnitKind::Baghouse => solids("Baghouse filter", "BH", &GAS_SOLIDS, &[], 200),
        UnitKind::RotaryFilter => solids("Rotary vacuum filter", "F", &FILTER, &[], 201),
        UnitKind::SolidsSeparator => solids("Solids separator", "SS", &SEPARATOR, &[], 202),
        UnitKind::Crystallizer => solids("Crystallizer", "CR", &CRYSTALLIZER, &DUTY, 203),
        UnitKind::Centrifuge => solids("Centrifuge", "CF", &CENTRIFUGE, &POWER, 204),
        UnitKind::ElectrostaticPrecipitator => solids(
            "Electrostatic precipitator",
            "ESP",
            &GAS_SOLIDS,
            &POWER,
            205,
        ),
        UnitKind::VenturiScrubber => solids("Venturi scrubber", "VS", &VENTURI, &[], 206),
        UnitKind::Crusher => solids("Crusher", "CRU", &SINGLE, &POWER, 207),
        UnitKind::Screen => solids("Screen", "SCR", &SCREEN, &[], 208),
        UnitKind::Dryer => solids("Dryer", "DR", &DRYER, &DUTY, 209),
        UnitKind::Granulator => solids("Granulator", "GR", &SINGLE, &[], 210),
        UnitKind::SolidsWasher => solids("Solids washer", "SW", &WASHER, &[], 211),
        _ => unreachable!("{kind:?} is not a solids kind"),
    }
}

/// Calculation status of a solids kind.
pub(crate) fn status(_kind: UnitKind) -> CalculationStatus {
    DRAWN_ONLY
}

/// Specifications of the solids kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum SolidsSpec {
    /// Cyclone (#198).
    Cyclone(CycloneSpec),
    /// Hydrocyclone (#199).
    Hydrocyclone(HydrocycloneSpec),
    /// Baghouse / fabric filter (#200).
    Baghouse(BaghouseSpec),
    /// Rotary vacuum / cake filter (#201).
    RotaryFilter(RotaryFilterSpec),
    /// Simple solids separator (#202).
    SolidsSeparator(SolidsSeparatorSpec),
    /// Crystallizer (#203).
    Crystallizer(CrystallizerSpec),
    /// Centrifuge (#204).
    Centrifuge(CentrifugeSpec),
    /// Electrostatic precipitator (#205).
    ElectrostaticPrecipitator(ElectrostaticPrecipitatorSpec),
    /// Venturi scrubber (#206).
    VenturiScrubber(VenturiScrubberSpec),
    /// Crusher (#207).
    Crusher(CrusherSpec),
    /// Screen / classifier (#208).
    Screen(ScreenSpec),
    /// Dryer (#209).
    Dryer(DryerSpec),
    /// Granulator / fluidised bed (#210).
    Granulator(GranulatorSpec),
    /// Solids washer / CCD (#211).
    SolidsWasher(SolidsWasherSpec),
}

macro_rules! empty_specs {
    ($($name:ident: $doc:literal),* $(,)?) => {
        $(
            #[doc = $doc]
            #[derive(Debug, Clone, Default, PartialEq)]
            pub struct $name {}
        )*
    };
}

empty_specs! {
    CycloneSpec: "Cyclone specification (#198).",
    HydrocycloneSpec: "Hydrocyclone specification (#199).",
    BaghouseSpec: "Baghouse specification (#200).",
    RotaryFilterSpec: "Rotary filter specification (#201).",
    SolidsSeparatorSpec: "Simple solids separator specification (#202).",
    CrystallizerSpec: "Crystallizer specification (#203).",
    CentrifugeSpec: "Centrifuge specification (#204).",
    ElectrostaticPrecipitatorSpec: "Electrostatic precipitator specification (#205).",
    VenturiScrubberSpec: "Venturi scrubber specification (#206).",
    CrusherSpec: "Crusher specification (#207).",
    ScreenSpec: "Screen / classifier specification (#208).",
    DryerSpec: "Dryer specification (#209).",
    GranulatorSpec: "Granulator specification (#210).",
    SolidsWasherSpec: "Solids washer specification (#211).",
}

/// Default specification of a solids kind.
pub(crate) fn default_spec(kind: UnitKind) -> SolidsSpec {
    use SolidsSpec as S;
    match kind {
        UnitKind::Cyclone => S::Cyclone(Default::default()),
        UnitKind::Hydrocyclone => S::Hydrocyclone(Default::default()),
        UnitKind::Baghouse => S::Baghouse(Default::default()),
        UnitKind::RotaryFilter => S::RotaryFilter(Default::default()),
        UnitKind::SolidsSeparator => S::SolidsSeparator(Default::default()),
        UnitKind::Crystallizer => S::Crystallizer(Default::default()),
        UnitKind::Centrifuge => S::Centrifuge(Default::default()),
        UnitKind::ElectrostaticPrecipitator => S::ElectrostaticPrecipitator(Default::default()),
        UnitKind::VenturiScrubber => S::VenturiScrubber(Default::default()),
        UnitKind::Crusher => S::Crusher(Default::default()),
        UnitKind::Screen => S::Screen(Default::default()),
        UnitKind::Dryer => S::Dryer(Default::default()),
        UnitKind::Granulator => S::Granulator(Default::default()),
        UnitKind::SolidsWasher => S::SolidsWasher(Default::default()),
        _ => unreachable!("{kind:?} has no SolidsSpec"),
    }
}

impl GroupSpec for SolidsSpec {
    fn kind(&self) -> UnitKind {
        use SolidsSpec as S;
        match self {
            S::Cyclone(_) => UnitKind::Cyclone,
            S::Hydrocyclone(_) => UnitKind::Hydrocyclone,
            S::Baghouse(_) => UnitKind::Baghouse,
            S::RotaryFilter(_) => UnitKind::RotaryFilter,
            S::SolidsSeparator(_) => UnitKind::SolidsSeparator,
            S::Crystallizer(_) => UnitKind::Crystallizer,
            S::Centrifuge(_) => UnitKind::Centrifuge,
            S::ElectrostaticPrecipitator(_) => UnitKind::ElectrostaticPrecipitator,
            S::VenturiScrubber(_) => UnitKind::VenturiScrubber,
            S::Crusher(_) => UnitKind::Crusher,
            S::Screen(_) => UnitKind::Screen,
            S::Dryer(_) => UnitKind::Dryer,
            S::Granulator(_) => UnitKind::Granulator,
            S::SolidsWasher(_) => UnitKind::SolidsWasher,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        Vec::new()
    }
}

/// Calculates a solids unit.
pub(crate) fn calculate(
    kind: UnitKind,
    _spec: &SolidsSpec,
    _inputs: &UnitInputs,
    _ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    Err(super::not_implemented(kind))
}

#[cfg(test)]
mod tests {
    use crate::UnitGroup;

    #[test]
    fn solids_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::Solids);
    }
}
