//! Property package selection stored in a project.

use std::sync::Arc;

use property_packages::{CubicPackage, EosCgPackage, FlowsheetPackage, GergPackage, NrtlPackage};
use serde::{Deserialize, Serialize};
use thermo_data::{COMPONENT_DATASET_ID, ComponentDatabase, LiquidDatabase};
use thermo_models::{
    ActivityModel, CUBIC_TEMPERATURE_RANGE_K, CubicModel, CubicVariant, GERG_CALORIC_REFERENCE_ID,
    GERG_DATA_REVISION, GERG_MODEL_VERSION, GergModel, VaporModel,
};
use thermo_solvers::PhFlashOptions;
use thermo_types::{ComponentId, Composition, CompositionBasis, DataMetadata, NormalizationPolicy};

use crate::{Fingerprint, ProjectError};

const LEGACY_PR_BM_ID: &str = "PR-exact-BM-Clapeyron-0.6.30";
const LEGACY_COMPONENT_DATASET_ID: &str = "coolprop-7.2-cubic-cantera-3.2-nasa7";

/// Property model supported by the flowsheet.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PackageChoice {
    /// Peng–Robinson (1976).
    #[serde(rename = "PR")]
    PengRobinson,
    /// Peng–Robinson with the Boston–Mathias alpha extrapolation.
    #[serde(rename = "PR-BM")]
    PengRobinsonBostonMathias,
    #[serde(rename = "SRK")]
    SoaveRedlichKwong,
    #[serde(rename = "GERG-2008")]
    Gerg2008,
    #[serde(rename = "EOS-CG-2021")]
    EosCg2021,
    #[serde(rename = "NRTL-IG")]
    NrtlIdealGas,
}

/// Existing name for the PR and PR-BM selection API.
pub type CubicVariantChoice = PackageChoice;

impl PackageChoice {
    pub const ALL: [Self; 6] = [
        Self::PengRobinsonBostonMathias,
        Self::PengRobinson,
        Self::SoaveRedlichKwong,
        Self::Gerg2008,
        Self::EosCg2021,
        Self::NrtlIdealGas,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::PengRobinson => "PR",
            Self::PengRobinsonBostonMathias => "PR-BM",
            Self::SoaveRedlichKwong => "SRK",
            Self::Gerg2008 => "GERG-2008",
            Self::EosCg2021 => "EOS-CG-2021",
            Self::NrtlIdealGas => "NRTL-IG",
        }
    }

    pub fn temperature_range(self) -> [f64; 2] {
        match self {
            Self::Gerg2008 => GergPackage::TEMPERATURE_RANGE_K,
            Self::EosCg2021 => EosCgPackage::TEMPERATURE_RANGE_K,
            Self::NrtlIdealGas => NrtlPackage::TEMPERATURE_RANGE_K,
            _ => CUBIC_TEMPERATURE_RANGE_K,
        }
    }

    fn cubic_model(self) -> Option<CubicVariant> {
        match self {
            Self::PengRobinson => Some(CubicVariant::PengRobinson),
            Self::PengRobinsonBostonMathias => Some(CubicVariant::PengRobinsonBostonMathias),
            Self::SoaveRedlichKwong => Some(CubicVariant::SoaveRedlichKwong),
            Self::Gerg2008 | Self::EosCg2021 | Self::NrtlIdealGas => None,
        }
    }

    pub fn uses_kij(self) -> bool {
        self.cubic_model().is_some()
    }
}

/// One immutable package choice for the whole flowsheet.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PackageSelection {
    pub variant: PackageChoice,
    /// Ordered component identifiers from the embedded database.
    pub components: Vec<String>,
    /// Full symmetric cubic interaction matrix in component order. Empty for embedded packages.
    #[serde(default)]
    pub kij: Vec<Vec<f64>>,
    /// Source of the interaction values, or an explicit approximation note.
    #[serde(default)]
    pub interaction_source: String,
    /// Lower PH search temperature, K.
    pub ph_minimum_temperature_k: f64,
    /// Upper PH search temperature, K.
    pub ph_maximum_temperature_k: f64,
}

/// Model, data and engine identity recorded with results. The fingerprint
/// treats aliases of the same model and dataset as the same identity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    /// [`crate::ENGINE_REVISION`] of the producing program.
    pub engine: String,
    /// Model identifier, for example `PR-exact-BM`.
    pub model: String,
    /// Model implementation version.
    pub model_version: String,
    /// Component dataset identifier.
    pub data: String,
    /// Component dataset version.
    pub data_version: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub caloric_reference: Option<String>,
}

impl Provenance {
    pub(crate) fn update_identifiers(&mut self) {
        if self.model == LEGACY_PR_BM_ID {
            self.model = CubicVariant::PengRobinsonBostonMathias.identifier().into();
        }
        if self.data == LEGACY_COMPONENT_DATASET_ID {
            self.data = COMPONENT_DATASET_ID.into();
        }
    }

    /// Use the original identifiers in the serialized hash input so a label
    /// change alone does not make existing results stale.
    pub(crate) fn fingerprint_identity(&self) -> Self {
        let mut identity = self.clone();
        if self.model == CubicVariant::PengRobinsonBostonMathias.identifier() {
            identity.model = LEGACY_PR_BM_ID.into();
        }
        if self.data == COMPONENT_DATASET_ID {
            identity.data = LEGACY_COMPONENT_DATASET_ID.into();
        }
        identity
    }
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
            ph_minimum_temperature_k: CUBIC_TEMPERATURE_RANGE_K[0],
            ph_maximum_temperature_k: CUBIC_TEMPERATURE_RANGE_K[1],
        }
    }

    /// Identifiers of all components in the embedded database.
    pub fn available_components() -> Result<Vec<String>, ProjectError> {
        let db = ComponentDatabase::embedded().map_err(|e| ProjectError::Package(e.to_string()))?;
        Ok(db.components().iter().map(|c| c.id().to_string()).collect())
    }

    pub fn available_components_for(variant: PackageChoice) -> Result<Vec<String>, ProjectError> {
        match variant {
            PackageChoice::Gerg2008 => Ok(GergPackage::COMPONENTS.map(String::from).to_vec()),
            PackageChoice::EosCg2021 => Ok(EosCgPackage::COMPONENTS.map(String::from).to_vec()),
            PackageChoice::NrtlIdealGas => Ok(NrtlPackage::COMPONENTS.map(String::from).to_vec()),
            _ => Self::available_components(),
        }
    }

    pub fn with_model(&self, variant: PackageChoice, zero_kij_source: &str) -> Self {
        use PackageChoice::{PengRobinson, PengRobinsonBostonMathias};
        if variant == self.variant
            || matches!(
                (self.variant, variant),
                (PengRobinson, PengRobinsonBostonMathias)
                    | (PengRobinsonBostonMathias, PengRobinson)
            )
        {
            return Self {
                variant,
                ..self.clone()
            };
        }
        Self::default_for(variant, self.components.clone(), zero_kij_source)
    }

    pub fn default_for(
        variant: PackageChoice,
        components: Vec<String>,
        zero_kij_source: &str,
    ) -> Self {
        if variant.uses_kij() {
            return Self::zero_interactions(variant, components, zero_kij_source);
        }
        let [lo, hi] = variant.temperature_range();
        Self {
            variant,
            components,
            kij: Vec::new(),
            interaction_source: String::new(),
            ph_minimum_temperature_k: lo,
            ph_maximum_temperature_k: hi,
        }
    }

    /// Builds the immutable package. Each call creates a new package
    /// identity; result currency in this crate is decided by fingerprints.
    pub fn build(&self) -> Result<Arc<dyn FlowsheetPackage>, ProjectError> {
        Ok(self.build_with_model()?.0)
    }

    /// Builds the package once together with its provenance and package
    /// fingerprint, for [`crate::solve_project`]. The prepared package
    /// remembers the selection it was built from, so a solve refuses a
    /// package that does not belong to the project being solved.
    pub fn prepare(&self) -> Result<PreparedPackage, ProjectError> {
        let (package, provenance) = self.build_with_model()?;
        let fingerprint = crate::fingerprint::package_fingerprint(self, &provenance);
        Ok(PreparedPackage {
            selection: self.clone(),
            package,
            provenance,
            fingerprint,
        })
    }

    /// Model/data/engine identity of results computed with this selection.
    pub fn provenance(&self) -> Result<Provenance, ProjectError> {
        Ok(self.build_with_model()?.1)
    }

    fn build_with_model(&self) -> Result<(Arc<dyn FlowsheetPackage>, Provenance), ProjectError> {
        let package = |e: &dyn std::fmt::Display| ProjectError::Package(e.to_string());
        if self.components.is_empty() {
            return Err(ProjectError::Package(
                "select at least one component".into(),
            ));
        }
        let (lo, hi) = (self.ph_minimum_temperature_k, self.ph_maximum_temperature_k);
        let [t_lo, t_hi] = self.variant.temperature_range();
        if !(lo.is_finite() && hi.is_finite() && (t_lo..=t_hi).contains(&lo))
            || !(t_lo..=t_hi).contains(&hi)
            || lo >= hi
        {
            return Err(ProjectError::Package(format!(
                "PH search interval must be ordered and within {t_lo}..{t_hi} K"
            )));
        }
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
        let Some(cubic_variant) = self.variant.cubic_model() else {
            if !self.kij.is_empty() || !self.interaction_source.is_empty() {
                return Err(ProjectError::Package(format!(
                    "{} uses embedded interactions. Clear the cubic matrix and source.",
                    self.variant.label()
                )));
            }
            if self.variant == PackageChoice::EosCg2021 {
                let model = EosCgPackage::new(self.components.clone(), [lo, hi])
                    .map_err(|e| package(&e))?;
                let provenance = Provenance {
                    engine: crate::ENGINE_REVISION.into(),
                    model: self.variant.label().into(),
                    model_version: EosCgPackage::VERSION.into(),
                    data: EosCgPackage::DATA.into(),
                    data_version: EosCgPackage::DATA_VERSION.into(),
                    caloric_reference: Some(EosCgPackage::REFERENCE.into()),
                };
                return Ok((Arc::new(model), provenance));
            }
            if self.variant == PackageChoice::NrtlIdealGas {
                let components = ComponentDatabase::embedded().map_err(|e| package(&e))?;
                let liquids = LiquidDatabase::embedded().map_err(|e| package(&e))?;
                let model = ActivityModel::from_databases(
                    &components,
                    &liquids,
                    &composition,
                    VaporModel::IdealGas,
                )
                .map_err(|e| package(&e))?;
                let metadata = model.metadata();
                let component_data =
                    components.metadata_for(self.components.iter().map(String::as_str));
                let liquid_data = liquids.metadata();
                let provenance = Provenance {
                    engine: crate::ENGINE_REVISION.into(),
                    model: self.variant.label().into(),
                    model_version: metadata.implementation_version,
                    data: format!("{} + {}", component_data.identifier, liquid_data.identifier),
                    data_version: format!("{} + {}", component_data.version, liquid_data.version),
                    caloric_reference: Some(ActivityModel::caloric_reference().identifier),
                };
                let nrtl = NrtlPackage::new(model, [lo, hi]).map_err(|e| package(&e))?;
                return Ok((Arc::new(nrtl), provenance));
            }
            let model = GergModel::new(&composition).map_err(|e| package(&e))?;
            let gerg = GergPackage::new(model, [lo, hi]).map_err(|e| package(&e))?;
            let provenance = Provenance {
                engine: crate::ENGINE_REVISION.into(),
                model: self.variant.label().into(),
                model_version: GERG_MODEL_VERSION.into(),
                data: "NIST AGA8 GERG2008 2.01".into(),
                data_version: GERG_DATA_REVISION.into(),
                caloric_reference: Some(GERG_CALORIC_REFERENCE_ID.into()),
            };
            return Ok((Arc::new(gerg), provenance));
        };
        let db = ComponentDatabase::embedded().map_err(|e| package(&e))?;
        let model = CubicModel::new(
            &db,
            &composition,
            cubic_variant,
            self.kij.clone(),
            self.interaction_source.clone(),
        )
        .map_err(|e| package(&e))?;
        let options = PhFlashOptions {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            ..Default::default()
        };
        let data = db.metadata_for(self.components.iter().map(String::as_str));
        let provenance = provenance_of(&model, data);
        Ok((Arc::new(CubicPackage::new(model, options)), provenance))
    }
}

fn provenance_of(model: &CubicModel, data: DataMetadata) -> Provenance {
    let m = model.metadata();
    Provenance {
        engine: crate::ENGINE_REVISION.into(),
        model: m.identifier,
        model_version: m.implementation_version,
        data: data.identifier,
        data_version: data.version,
        caloric_reference: None,
    }
}

/// A package built from one [`PackageSelection`], with the provenance and
/// package fingerprint of the results it produces. Prepare it once per
/// selection with [`PackageSelection::prepare`] and reuse it for every solve
/// of that selection; [`crate::solve_project`] builds nothing itself.
#[derive(Debug, Clone)]
pub struct PreparedPackage {
    selection: PackageSelection,
    package: Arc<dyn FlowsheetPackage>,
    provenance: Provenance,
    fingerprint: Fingerprint,
}

impl PreparedPackage {
    /// The selection this package was built from.
    pub fn selection(&self) -> &PackageSelection {
        &self.selection
    }

    /// The built property package.
    pub fn package(&self) -> &dyn FlowsheetPackage {
        self.package.as_ref()
    }

    /// Provenance recorded with results of this package.
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }

    /// Package fingerprint; equal to [`crate::Fingerprints::package`] of
    /// the selection.
    pub fn fingerprint(&self) -> Fingerprint {
        self.fingerprint
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fingerprint::package_fingerprint;
    use crate::{Currency, Fingerprints, Project, demonstration, solve_project};
    use std::ops::ControlFlow;

    #[test]
    fn previous_engine_results_load_as_stale() {
        let mut project = demonstration(CubicVariantChoice::PengRobinson);
        let package = project.package.prepare().unwrap();
        let mut results =
            solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
        results.provenance.engine = "flowsheet-core a10-tears-v1".into();
        let old_package = package_fingerprint(&project.package, &results.provenance);
        let old = Fingerprints::compute_with(&project.flowsheet, &project.feeds, Some(old_package));
        for (id, result) in &mut results.units {
            result.fingerprint = old.units[id];
        }
        for (id, result) in &mut results.streams {
            result.fingerprint = old.stream(&project.flowsheet, *id).unwrap();
        }
        project.results = Some(results);

        let loaded = Project::from_json(&project.to_json().unwrap()).unwrap();
        let results = loaded.results.as_ref().unwrap();
        let current = loaded.fingerprints();
        for unit in loaded.flowsheet.units() {
            assert_eq!(results.unit_currency(unit.id, &old), Currency::Current);
            assert_eq!(results.unit_currency(unit.id, &current), Currency::Stale);
        }
        for stream in loaded.flowsheet.streams() {
            assert_eq!(
                results.stream_currency(&loaded.flowsheet, stream.id, &old),
                Currency::Current
            );
            assert_eq!(
                results.stream_currency(&loaded.flowsheet, stream.id, &current),
                Currency::Stale
            );
        }
    }

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
