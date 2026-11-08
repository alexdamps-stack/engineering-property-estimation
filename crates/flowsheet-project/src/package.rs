//! Exact cubic package selection stored in a project.

use property_packages::CubicPackage;
use serde::{Deserialize, Serialize};
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::PhFlashOptions;
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

use crate::ProjectError;

/// Cubic equation of state variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CubicVariantChoice {
    /// Peng–Robinson (1976).
    #[serde(rename = "PR")]
    PengRobinson,
    /// Peng–Robinson with the Boston–Mathias alpha extrapolation.
    #[serde(rename = "PR-BM")]
    PengRobinsonBostonMathias,
}

impl CubicVariantChoice {
    /// Both variants, in menu order.
    pub const ALL: [CubicVariantChoice; 2] = [
        CubicVariantChoice::PengRobinsonBostonMathias,
        CubicVariantChoice::PengRobinson,
    ];

    /// Short label, `PR` or `PR-BM`.
    pub fn label(self) -> &'static str {
        match self {
            CubicVariantChoice::PengRobinson => "PR",
            CubicVariantChoice::PengRobinsonBostonMathias => "PR-BM",
        }
    }

    fn model(self) -> CubicVariant {
        match self {
            CubicVariantChoice::PengRobinson => CubicVariant::PengRobinson,
            CubicVariantChoice::PengRobinsonBostonMathias => {
                CubicVariant::PengRobinsonBostonMathias
            }
        }
    }
}

/// One immutable cubic package choice for the whole flowsheet.
///
/// The interaction matrix and its source are always explicit; no missing
/// pair is filled in. The PH interval bounds the temperature search of
/// PH flashes and must lie within the 200..600 K cubic domain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PackageSelection {
    /// PR or PR-BM.
    pub variant: CubicVariantChoice,
    /// Ordered component identifiers from the embedded database.
    pub components: Vec<String>,
    /// Full symmetric binary interaction matrix in component order.
    pub kij: Vec<Vec<f64>>,
    /// Source of the interaction values, or an explicit approximation note.
    pub interaction_source: String,
    /// Lower PH search temperature, K.
    pub ph_minimum_temperature_k: f64,
    /// Upper PH search temperature, K.
    pub ph_maximum_temperature_k: f64,
}

/// Model, data and engine identity recorded with results. Results computed
/// under a different provenance are never shown as current.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    /// [`crate::ENGINE_REVISION`] of the producing program.
    pub engine: String,
    /// Model identifier, for example `PR-exact-BM-Clapeyron-0.6.30`.
    pub model: String,
    /// Model implementation version.
    pub model_version: String,
    /// Component dataset identifier.
    pub data: String,
    /// Component dataset version.
    pub data_version: String,
}

impl PackageSelection {
    /// A selection with an explicitly chosen all-zero interaction matrix.
    /// `source` must say that this is an approximation and why it was chosen.
    pub fn zero_interactions(
        variant: CubicVariantChoice,
        components: Vec<String>,
        source: impl Into<String>,
    ) -> Self {
        let n = components.len();
        Self {
            variant,
            components,
            kij: vec![vec![0.0; n]; n],
            interaction_source: source.into(),
            ph_minimum_temperature_k: 200.0,
            ph_maximum_temperature_k: 600.0,
        }
    }

    /// Identifiers of all components in the embedded database.
    pub fn available_components() -> Result<Vec<String>, ProjectError> {
        let db = ComponentDatabase::embedded().map_err(|e| ProjectError::Package(e.to_string()))?;
        Ok(db.components().iter().map(|c| c.id().to_string()).collect())
    }

    /// Builds the immutable package. Each call creates a new package
    /// identity; result currency in this crate is decided by fingerprints.
    pub fn build(&self) -> Result<CubicPackage, ProjectError> {
        Ok(self.build_with_model()?.0)
    }

    /// Model/data/engine identity of results computed with this selection.
    pub fn provenance(&self) -> Result<Provenance, ProjectError> {
        let model = self.build_with_model()?.1;
        let m = model.metadata();
        let d = ComponentDatabase::embedded()
            .map_err(|e| ProjectError::Package(e.to_string()))?
            .metadata();
        Ok(Provenance {
            engine: crate::ENGINE_REVISION.into(),
            model: m.identifier,
            model_version: m.implementation_version,
            data: d.identifier,
            data_version: d.version,
        })
    }

    fn build_with_model(&self) -> Result<(CubicPackage, CubicModel), ProjectError> {
        let package = |e: &dyn std::fmt::Display| ProjectError::Package(e.to_string());
        if self.components.is_empty() {
            return Err(ProjectError::Package(
                "select at least one component".into(),
            ));
        }
        let (lo, hi) = (self.ph_minimum_temperature_k, self.ph_maximum_temperature_k);
        if !(lo.is_finite() && hi.is_finite() && (200.0..=600.0).contains(&lo))
            || !(200.0..=600.0).contains(&hi)
            || lo >= hi
        {
            return Err(ProjectError::Package(
                "PH search interval must be ordered and within 200..600 K".into(),
            ));
        }
        let db = ComponentDatabase::embedded().map_err(|e| package(&e))?;
        let entries = self
            .components
            .iter()
            .map(|id| Ok((ComponentId::new(id.clone()).map_err(|e| package(&e))?, 1.0)))
            .collect::<Result<Vec<_>, ProjectError>>()?;
        // The fractions only seed the model; each calculation supplies its own.
        let composition = Composition::new(
            CompositionBasis::Mole,
            entries,
            NormalizationPolicy::Normalize,
        )
        .map_err(|e| package(&e))?;
        let model = CubicModel::new(
            &db,
            &composition,
            self.variant.model(),
            self.kij.clone(),
            self.interaction_source.clone(),
        )
        .map_err(|e| package(&e))?;
        let options = PhFlashOptions {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            ..Default::default()
        };
        Ok((CubicPackage::new(model.clone(), options), model))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_selections_are_rejected_without_defaults() {
        let ok = PackageSelection::zero_interactions(
            CubicVariantChoice::PengRobinson,
            vec!["methane".into(), "ethane".into()],
            "test approximation",
        );
        assert!(ok.build().is_ok());

        let mut empty = ok.clone();
        empty.components.clear();
        empty.kij.clear();
        assert!(empty.build().is_err());

        let mut missing_pair = ok.clone();
        missing_pair.kij = vec![vec![0.0; 2]];
        assert!(missing_pair.build().is_err());

        let mut no_source = ok.clone();
        no_source.interaction_source = " ".into();
        assert!(no_source.build().is_err());

        let mut unknown = ok.clone();
        unknown.components[1] = "unobtainium".into();
        assert!(unknown.build().is_err());

        let mut bounds = ok.clone();
        bounds.ph_maximum_temperature_k = 700.0;
        assert!(bounds.build().is_err());
        bounds.ph_maximum_temperature_k = 150.0;
        assert!(bounds.build().is_err());
    }

    #[test]
    fn provenance_names_model_and_data() {
        let sel = PackageSelection::zero_interactions(
            CubicVariantChoice::PengRobinsonBostonMathias,
            vec!["methane".into()],
            "test approximation",
        );
        let p = sel.provenance().unwrap();
        assert_eq!(p.engine, crate::ENGINE_REVISION);
        assert!(p.model.contains("BM"));
        assert!(!p.data_version.is_empty());
        assert!(
            PackageSelection::available_components()
                .unwrap()
                .contains(&"methane".to_string())
        );
    }
}
