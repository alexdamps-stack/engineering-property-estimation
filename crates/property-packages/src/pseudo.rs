//! Cubic packages with petroleum pseudo-components (#229).
//!
//! Pseudo-components are characterised in `thermo-data` from stored inputs
//! (Tb, SG, optional M or specified constants) and appended to the embedded
//! component database. The package is an ordinary [`CubicPackage`]: the
//! interaction matrix must cover every pair, including pseudo-component
//! pairs, with an explicit source or approximation note. Nothing is filled
//! in by default.
use thermo_data::{
    ComponentDatabase, PSEUDO_COMPONENT_PREFIX, PSEUDO_IDEAL_GAS_RANGE_K, PseudoComponentInput,
};
use thermo_models::{CUBIC_TEMPERATURE_RANGE_K, CubicModel, CubicVariant};
use thermo_solvers::PhFlashOptions;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, InputError, NormalizationPolicy,
};

use crate::CubicPackage;

/// Temperature interval in which a cubic package with pseudo-components can
/// evaluate every component: the cubic envelope intersected with the
/// declared Kesler–Lee ideal-gas domain, K.
pub const PSEUDO_CUBIC_TEMPERATURE_RANGE_K: [f64; 2] = [
    if PSEUDO_IDEAL_GAS_RANGE_K[0] > CUBIC_TEMPERATURE_RANGE_K[0] {
        PSEUDO_IDEAL_GAS_RANGE_K[0]
    } else {
        CUBIC_TEMPERATURE_RANGE_K[0]
    },
    if PSEUDO_IDEAL_GAS_RANGE_K[1] < CUBIC_TEMPERATURE_RANGE_K[1] {
        PSEUDO_IDEAL_GAS_RANGE_K[1]
    } else {
        CUBIC_TEMPERATURE_RANGE_K[1]
    },
];

/// The embedded database extended by the characterised pseudo-components,
/// after checking that `component_ids` (package order) and the definitions
/// match one to one: every `pseudo:` identifier must be defined exactly
/// once, and every definition must be selected.
pub fn pseudo_component_database(
    component_ids: &[String],
    pseudo_components: &[PseudoComponentInput],
) -> Result<ComponentDatabase, CalculationError> {
    for id in component_ids {
        if id.starts_with(PSEUDO_COMPONENT_PREFIX) && !pseudo_components.iter().any(|p| &p.id == id)
        {
            return Err(CalculationError::MissingParameter(format!(
                "pseudo-component {id} is selected but has no characterisation input"
            )));
        }
    }
    for (i, p) in pseudo_components.iter().enumerate() {
        if pseudo_components[..i].iter().any(|q| q.id == p.id) {
            return Err(InputError::DuplicateComponent.into());
        }
        if !component_ids.contains(&p.id) {
            return Err(CalculationError::UnsupportedCapability(format!(
                "pseudo-component {} is defined but not selected; remove it or add it to the components",
                p.id
            )));
        }
    }
    ComponentDatabase::embedded()?.with_pseudo_components(pseudo_components)
}

/// A cubic model over database components and pseudo-components in
/// `component_ids` order, with the explicit interaction matrix and source.
pub fn cubic_model_with_pseudo_components(
    variant: CubicVariant,
    component_ids: &[String],
    pseudo_components: &[PseudoComponentInput],
    kij: Vec<Vec<f64>>,
    interaction_source: impl Into<String>,
) -> Result<CubicModel, CalculationError> {
    let db = pseudo_component_database(component_ids, pseudo_components)?;
    let entries = component_ids
        .iter()
        .map(|id| Ok((ComponentId::new(id.clone())?, 1.0)))
        .collect::<Result<Vec<_>, CalculationError>>()?;
    // The fractions only seed the model; each calculation supplies its own.
    let composition = Composition::new(
        CompositionBasis::Mole,
        entries,
        NormalizationPolicy::Normalize,
    )?;
    CubicModel::new(&db, &composition, variant, kij, interaction_source)
}

/// A cubic package with pseudo-components. The PH/PS search interval in
/// `options` must lie inside [`PSEUDO_CUBIC_TEMPERATURE_RANGE_K`], so that no
/// search evaluates a pseudo-component outside its ideal-gas domain.
pub fn cubic_package_with_pseudo_components(
    variant: CubicVariant,
    component_ids: &[String],
    pseudo_components: &[PseudoComponentInput],
    kij: Vec<Vec<f64>>,
    interaction_source: impl Into<String>,
    options: PhFlashOptions,
) -> Result<CubicPackage, CalculationError> {
    let [lo, hi] = PSEUDO_CUBIC_TEMPERATURE_RANGE_K;
    let (t_lo, t_hi) = (options.minimum_temperature_k, options.maximum_temperature_k);
    if pseudo_components.is_empty() {
        return Err(CalculationError::MissingParameter(
            "no pseudo-components; build an ordinary cubic package instead".into(),
        ));
    }
    if !(t_lo.is_finite() && t_hi.is_finite()) || t_lo < lo || t_hi > hi || t_lo >= t_hi {
        return Err(CalculationError::OutsideValidityDomain(format!(
            "a package with pseudo-components needs an ordered search interval within {lo}..={hi} K (Kesler-Lee ideal-gas domain)"
        )));
    }
    let model = cubic_model_with_pseudo_components(
        variant,
        component_ids,
        pseudo_components,
        kij,
        interaction_source,
    )?;
    Ok(CubicPackage::new(model, options))
}
