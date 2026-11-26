//! Pinned component constants and ideal-gas NASA7 properties.
//! Source coefficients and notices are retained in `data/source-extracts`.
mod liquid;
pub use liquid::*;
mod transport;
pub use transport::*;
mod assay;
pub use assay::*;
mod pseudo;
pub use pseudo::*;

use serde::Deserialize;
use thermo_types::{
    CalculationError, ComponentId, DataMetadata, InputError, MOLAR_GAS_CONSTANT as R, MolarMass,
    Temperature,
};

pub const DATA_VERSION: &str = "2026-11-07.1";
/// Identifier of the embedded component constants and ideal-gas NASA7 data.
pub const COMPONENT_DATASET_ID: &str = "cubic-constants-nasa7";
/// Version of the A17 reacting-component extension
/// (`components-reacting.provisional.json`: hydrogen, carbon monoxide).
pub const REACTING_DATA_VERSION: &str = "2026-11-23.1";
const COMPONENT_SCHEMA: &str = "a2-source-extract-provisional-v1";
pub const REFERENCE_TEMPERATURE: f64 = 298.15;
pub const REFERENCE_PRESSURE: f64 = 100_000.0;

#[derive(Debug, Clone, Deserialize)]
struct Constants {
    #[serde(rename = "critical_temperature_K")]
    tc: f64,
    #[serde(rename = "critical_pressure_Pa")]
    pc: f64,
    acentric_factor: f64,
    molar_mass_kg_per_mol: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct Nasa7 {
    #[serde(rename = "temperature_breaks_K")]
    breaks: [f64; 3],
    coefficients: [[f64; 7]; 2],
}

#[derive(Debug, Clone, Deserialize)]
struct Record {
    id: String,
    cas: String,
    formula: String,
    constants: Constants,
    ideal_gas: Nasa7,
    /// From the reacting-component extension rather than the base extract.
    #[serde(skip)]
    extension: bool,
    /// Characterisation of a petroleum pseudo-component (#229); `None` for
    /// embedded database records.
    #[serde(skip)]
    pseudo: Option<std::sync::Arc<PseudoCharacterisation>>,
}

/// Immutable component from the checked embedded dataset.
#[derive(Debug, Clone)]
pub struct Component(Record);

/// Molar ideal-gas quantities, h and s anchored to zero at 298.15 K.
/// Pure entropy excludes pressure and mixture terms. Units: J/mol, J/(mol K).
#[derive(Debug, Clone, Copy)]
pub struct IdealGasProperties {
    pub enthalpy: f64,
    pub entropy: f64,
    pub heat_capacity: f64,
}

impl Component {
    pub fn id(&self) -> &str {
        &self.0.id
    }
    pub fn cas(&self) -> &str {
        &self.0.cas
    }
    /// Molecular formula as recorded in the source extract, for example `CH3OH`.
    pub fn formula(&self) -> &str {
        &self.0.formula
    }
    /// Atoms per molecule by element symbol, in first-appearance order of the
    /// formula. The embedded formulas are checked to parse when loading.
    /// Pseudo-components have no formula and return an empty list; reactors
    /// reject them before using element balances.
    pub fn elements(&self) -> Vec<(String, u32)> {
        if self.is_pseudo() {
            return Vec::new();
        }
        parse_formula(&self.0.formula).expect("validated embedded formula")
    }
    /// True for a characterised petroleum pseudo-component (#229).
    pub fn is_pseudo(&self) -> bool {
        self.0.pseudo.is_some()
    }
    /// Inputs, derived constants and methods of a pseudo-component.
    pub fn pseudo_characterisation(&self) -> Option<&PseudoCharacterisation> {
        self.0.pseudo.as_deref()
    }
    /// Ideal-gas standard enthalpy of formation at 298.15 K, J/mol, from the
    /// retained NASA7 integration constant a6 of the lower interval:
    /// `R (a1 T + a2 T^2/2 + ... + a5 T^5/5 + a6)` at T = 298.15 K. The
    /// engine caloric scale ([`Self::ideal_gas`]) is zero at that state; this
    /// value restores the source formation basis for reacting balances.
    /// Pseudo-components have no formation basis and return NaN.
    pub fn standard_formation_enthalpy(&self) -> f64 {
        if self.is_pseudo() {
            return f64::NAN;
        }
        let low = &self.0.ideal_gas.coefficients[0];
        R * (primitive(low, REFERENCE_TEMPERATURE, false) + low[5])
    }
    /// Ideal-gas absolute (third-law) entropy at 298.15 K and 1 bar,
    /// J/(mol K), from the retained NASA7 integration constant a7 of the
    /// lower interval: `R (a1 ln T + a2 T + ... + a5 T^4/4 + a7)`.
    /// Pseudo-components have no third-law basis and return NaN.
    pub fn standard_entropy(&self) -> f64 {
        if self.is_pseudo() {
            return f64::NAN;
        }
        let low = &self.0.ideal_gas.coefficients[0];
        R * (primitive(low, REFERENCE_TEMPERATURE, true) + low[6])
    }
    /// Ideal-gas standard Gibbs energy at T and 1 bar on the formation basis
    /// (elements in their reference states at 298.15 K), J/mol:
    /// `dHf + h(T) - T (S0 + s(T))`, with h and s the cp-integrated values
    /// of [`Self::ideal_gas`]. Differences over a balanced reaction give
    /// `-R T ln K`.
    pub fn standard_gibbs_energy(&self, temperature: Temperature) -> Result<f64, CalculationError> {
        if self.is_pseudo() {
            return Err(CalculationError::MissingParameter(format!(
                "pseudo-component {} has no formation properties",
                self.id()
            )));
        }
        let ig = self.ideal_gas(temperature)?;
        let t = temperature.value();
        Ok(self.standard_formation_enthalpy() + ig.enthalpy
            - t * (self.standard_entropy() + ig.entropy))
    }
    pub fn critical_temperature(&self) -> f64 {
        self.0.constants.tc
    }
    pub fn critical_pressure(&self) -> f64 {
        self.0.constants.pc
    }
    pub fn acentric_factor(&self) -> f64 {
        self.0.constants.acentric_factor
    }
    pub fn molar_mass(&self) -> MolarMass {
        MolarMass::new(self.0.constants.molar_mass_kg_per_mol).expect("validated embedded mass")
    }

    pub fn ideal_gas(
        &self,
        temperature: Temperature,
    ) -> Result<IdealGasProperties, CalculationError> {
        let t = temperature.value();
        let [lo, join, hi] = self.0.ideal_gas.breaks;
        if !(lo..=hi).contains(&t) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "{} {}: {lo}..={hi} K",
                self.id(),
                if self.is_pseudo() {
                    "Kesler-Lee ideal-gas cp"
                } else {
                    "NASA7"
                }
            )));
        }
        // Upper interval owns the join. Integrating cp across it avoids a jump
        // from the independently fitted source integration constants a6/a7.
        let low = &self.0.ideal_gas.coefficients[0];
        let high = &self.0.ideal_gas.coefficients[1];
        let integral = |entropy| {
            if t < join {
                primitive(low, t, entropy) - primitive(low, REFERENCE_TEMPERATURE, entropy)
            } else {
                primitive(low, join, entropy) - primitive(low, REFERENCE_TEMPERATURE, entropy)
                    + primitive(high, t, entropy)
                    - primitive(high, join, entropy)
            }
        };
        let a = if t < join { low } else { high };
        Ok(IdealGasProperties {
            enthalpy: R * integral(false),
            entropy: R * integral(true),
            heat_capacity: R * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4])))),
        })
    }
}

/// Parses a plain molecular formula (`CH3OH`, `C2H6`, `N2`): element symbols
/// of one capital and optional lowercase letters, each followed by an optional
/// positive count. Repeated symbols are summed. Charges, parentheses and
/// hydrates are not supported and are rejected.
pub fn parse_formula(formula: &str) -> Result<Vec<(String, u32)>, CalculationError> {
    let bad = || CalculationError::MissingParameter(format!("unsupported formula {formula:?}"));
    let mut out: Vec<(String, u32)> = Vec::new();
    let mut chars = formula.chars().peekable();
    while let Some(c) = chars.next() {
        if !c.is_ascii_uppercase() {
            return Err(bad());
        }
        let mut symbol = c.to_string();
        while let Some(&l) = chars.peek().filter(|l| l.is_ascii_lowercase()) {
            symbol.push(l);
            chars.next();
        }
        let mut digits = String::new();
        while let Some(&d) = chars.peek().filter(|d| d.is_ascii_digit()) {
            digits.push(d);
            chars.next();
        }
        let count = if digits.is_empty() {
            1
        } else {
            digits.parse::<u32>().map_err(|_| bad())?
        };
        if count == 0 {
            return Err(bad());
        }
        match out.iter_mut().find(|(s, _)| *s == symbol) {
            Some((_, n)) => *n = n.checked_add(count).ok_or_else(bad)?,
            None => out.push((symbol, count)),
        }
    }
    if out.is_empty() {
        return Err(bad());
    }
    Ok(out)
}

fn primitive(a: &[f64; 7], t: f64, entropy: bool) -> f64 {
    if entropy {
        a[0] * t.ln()
            + a[1] * t
            + a[2] * t * t / 2.0
            + a[3] * t.powi(3) / 3.0
            + a[4] * t.powi(4) / 4.0
    } else {
        a[0] * t
            + a[1] * t * t / 2.0
            + a[2] * t.powi(3) / 3.0
            + a[3] * t.powi(4) / 4.0
            + a[4] * t.powi(5) / 5.0
    }
}

/// Immutable, checked embedded dataset. This does not qualify experimental accuracy.
#[derive(Debug, Clone)]
pub struct ComponentDatabase {
    components: Vec<Component>,
}

impl ComponentDatabase {
    /// The embedded dataset. It is parsed and checked once per process; later calls
    /// return a clone of that result.
    pub fn embedded() -> Result<Self, CalculationError> {
        static EMBEDDED: std::sync::OnceLock<Result<ComponentDatabase, CalculationError>> =
            std::sync::OnceLock::new();
        EMBEDDED
            .get_or_init(|| {
                let mut components = Self::from_extract(
                    include_str!("../../../data/source-extracts/components.provisional.json"),
                    DATA_VERSION,
                    false,
                )?;
                components.extend(Self::from_extract(
                    include_str!(
                        "../../../data/source-extracts/components-reacting.provisional.json"
                    ),
                    REACTING_DATA_VERSION,
                    true,
                )?);
                for (i, c) in components.iter().enumerate() {
                    if components[..i].iter().any(|o| o.id() == c.id()) {
                        return Err(CalculationError::MissingParameter(format!(
                            "component {} is defined in both extracts",
                            c.id()
                        )));
                    }
                }
                Ok(Self { components })
            })
            .clone()
    }

    /// The header must match the expected version (`DATA_VERSION` or
    /// `REACTING_DATA_VERSION`), which the metadata reports, so a regenerated
    /// extract cannot silently carry a stale version label.
    fn from_extract(
        text: &str,
        version: &str,
        extension: bool,
    ) -> Result<Vec<Component>, CalculationError> {
        #[derive(Deserialize)]
        struct Extract {
            schema: String,
            version: String,
            components: Vec<Record>,
        }
        let extract: Extract = serde_json::from_str(text)
            .map_err(|e| CalculationError::MissingParameter(e.to_string()))?;
        if extract.schema != COMPONENT_SCHEMA || extract.version != version {
            return Err(CalculationError::MissingParameter(format!(
                "invalid embedded component extract header: schema {} version {} (expected {COMPONENT_SCHEMA} {version})",
                extract.schema, extract.version
            )));
        }
        for (i, c) in extract.components.iter().enumerate() {
            ComponentId::new(&c.id)?;
            if c.id.starts_with(PSEUDO_COMPONENT_PREFIX) {
                return Err(CalculationError::MissingParameter(format!(
                    "embedded record {} uses the reserved pseudo-component prefix",
                    c.id
                )));
            }
            if extract.components[..i].iter().any(|other| other.id == c.id)
                || !c.constants.tc.is_finite()
                || c.constants.tc <= 0.0
                || !c.constants.pc.is_finite()
                || c.constants.pc <= 0.0
                || !c.constants.acentric_factor.is_finite()
                || c.ideal_gas.breaks != [200.0, 1000.0, 6000.0]
                || c.ideal_gas
                    .coefficients
                    .iter()
                    .flatten()
                    .any(|v| !v.is_finite())
            {
                return Err(CalculationError::MissingParameter(format!(
                    "invalid embedded record {}",
                    c.id
                )));
            }
            MolarMass::new(c.constants.molar_mass_kg_per_mol)?;
            parse_formula(&c.formula).map_err(|_| {
                CalculationError::MissingParameter(format!(
                    "invalid embedded formula {} for {}",
                    c.formula, c.id
                ))
            })?;
        }
        Ok(extract
            .components
            .into_iter()
            .map(|mut record| {
                record.extension = extension;
                Component(record)
            })
            .collect())
    }

    pub fn get(&self, id: &str) -> Result<&Component, CalculationError> {
        self.components
            .iter()
            .find(|c| c.id() == id)
            .ok_or_else(|| CalculationError::MissingParameter(format!("component {id}")))
    }
    pub fn components(&self) -> &[Component] {
        &self.components
    }
    /// Metadata of the base extract. Use [`Self::metadata_for`] for a
    /// selection that may include extension components.
    pub fn metadata(&self) -> DataMetadata {
        DataMetadata { identifier: COMPONENT_DATASET_ID.into(), version: DATA_VERSION.into(),
            source: "data/source-extracts/components.provisional.json; upstream URLs, hashes and notices retained".into() }
    }

    /// Metadata of the data a component selection uses: the base metadata
    /// unless a selected component comes from the reacting extension, in
    /// which case version and source name both extracts. Selections of base
    /// components keep their provenance unchanged.
    pub fn metadata_for<'a>(&self, ids: impl IntoIterator<Item = &'a str>) -> DataMetadata {
        let selected: Vec<&Component> =
            ids.into_iter().filter_map(|id| self.get(id).ok()).collect();
        let extended = selected.iter().any(|c| c.0.extension);
        let pseudo = selected.iter().any(|c| c.is_pseudo());
        let mut data = if extended {
            DataMetadata {
                identifier: COMPONENT_DATASET_ID.into(),
                version: format!("{DATA_VERSION}+reacting-{REACTING_DATA_VERSION}"),
                source: "data/source-extracts/components.provisional.json and components-reacting.provisional.json; upstream URLs, hashes and notices retained".into(),
            }
        } else {
            self.metadata()
        };
        if pseudo {
            // The characterisation inputs themselves are stored with the
            // selection (and its fingerprint); the version names the method set.
            data.version = format!("{}+pseudo-{PSEUDO_CHARACTERISATION_VERSION}", data.version);
            data.source = format!(
                "{}; pseudo-components characterised by {PSEUDO_CHARACTERISATION_ID} (Fahim et al. 2010, ch. 3 transcription)",
                data.source
            );
        }
        data
    }

    /// The database extended by characterised pseudo-components, appended in
    /// input order. Identifiers must carry [`PSEUDO_COMPONENT_PREFIX`] and be
    /// unique; each input is characterised and domain-checked here.
    pub fn with_pseudo_components(
        &self,
        inputs: &[PseudoComponentInput],
    ) -> Result<Self, CalculationError> {
        let mut components = self.components.clone();
        for input in inputs {
            let component = input.component()?;
            if components.iter().any(|c| c.id() == component.id()) {
                return Err(InputError::DuplicateComponent.into());
            }
            components.push(component);
        }
        Ok(Self { components })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EMBEDDED: &str =
        include_str!("../../../data/source-extracts/components.provisional.json");

    #[test]
    fn embedded_header_matches_reported_version() {
        let db = ComponentDatabase::embedded().unwrap();
        assert_eq!(db.metadata().version, DATA_VERSION);
        assert!(!db.components().is_empty());
    }

    #[test]
    fn formulas_parse_to_element_counts() {
        let e = |f| parse_formula(f).unwrap();
        assert_eq!(
            e("CH3OH"),
            [("C".into(), 1), ("H".into(), 4), ("O".into(), 1)]
        );
        assert_eq!(e("C3H8"), [("C".into(), 3), ("H".into(), 8)]);
        assert_eq!(e("N2"), [("N".into(), 2)]);
        assert_eq!(e("NaCl"), [("Na".into(), 1), ("Cl".into(), 1)]);
        for bad in ["", "ch4", "C0", "H2O+", "Ca(OH)2", "2H"] {
            assert!(parse_formula(bad).is_err(), "{bad}");
        }
        let db = ComponentDatabase::embedded().unwrap();
        assert_eq!(
            db.get("methanol").unwrap().elements(),
            e("CH3OH"),
            "embedded formula"
        );
    }

    /// The pinned Cantera NASA7 constants reproduce tabulated ideal-gas
    /// formation enthalpies (CODATA/JANAF/ATcT values as compiled in the NIST
    /// Chemistry WebBook) to well within their source uncertainty. Independent
    /// Cantera values for the exact reacting balances are in the
    /// `flowsheet-core` reaction reference test.
    #[test]
    fn formation_enthalpies_match_tabulated_values() {
        let db = ComponentDatabase::embedded().unwrap();
        for (id, kj_per_mol) in [
            ("methane", -74.87),
            ("ethane", -84.0),
            ("propane", -104.7),
            ("nitrogen", 0.0),
            ("oxygen", 0.0),
            ("carbon-dioxide", -393.51),
            ("water", -241.826),
            ("methanol", -201.0),
            ("hydrogen", 0.0),
            ("carbon-monoxide", -110.53),
        ] {
            let h = db.get(id).unwrap().standard_formation_enthalpy() / 1000.0;
            assert!((h - kj_per_mol).abs() < 0.6, "{id}: {h} kJ/mol");
        }
    }

    #[test]
    fn reacting_extension_is_loaded_and_named_in_provenance() {
        let db = ComponentDatabase::embedded().unwrap();
        for id in ["hydrogen", "carbon-monoxide"] {
            let c = db.get(id).unwrap();
            assert!(c.0.extension, "{id}");
            assert!(c.ideal_gas(Temperature::new(500.0).unwrap()).is_ok());
        }
        assert_eq!(db.components().len(), 10);
        assert_eq!(db.metadata_for(["methane", "water"]), db.metadata());
        let mixed = db.metadata_for(["methane", "hydrogen"]);
        assert_eq!(
            mixed.version,
            format!("{DATA_VERSION}+reacting-{REACTING_DATA_VERSION}")
        );
        assert_eq!(mixed.identifier, db.metadata().identifier);
    }

    #[test]
    fn stale_version_or_schema_is_rejected() {
        let stale = EMBEDDED.replacen(DATA_VERSION, "1999-01-01.1", 1);
        assert_ne!(stale, EMBEDDED);
        assert!(matches!(
            ComponentDatabase::from_extract(&stale, DATA_VERSION, false),
            Err(CalculationError::MissingParameter(_))
        ));
        let schema = EMBEDDED.replacen(COMPONENT_SCHEMA, "a2-source-extract-provisional-v0", 1);
        assert_ne!(schema, EMBEDDED);
        assert!(ComponentDatabase::from_extract(&schema, DATA_VERSION, false).is_err());
    }
}
