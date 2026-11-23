//! Reactor chemistry: element and formation data in package order, conversion
//! and yield calculations, and the reacting terms of the unit balances.
//!
//! The package caloric scale is zero for every pure ideal gas at 298.15 K and
//! 1 bar (`thermo_models::IDEAL_GAS_CALORIC_REFERENCE_ID`), so stream
//! enthalpies carry no formation contribution. A reactor therefore adds the
//! ideal-gas standard formation enthalpies of what it generates and consumes:
//!
//! `n_out h_out - n_in h_in + sum_i g_i dHf_i = Q`,
//!
//! where `g_i` is the net generation of component `i` (mol/s) and `dHf_i` is
//! `thermo_data::Component::standard_formation_enthalpy`. Non-reacting units
//! have `g = 0` and keep their usual balance.

use crate::{Reaction, ReactionExtent};
use property_packages::CubicPackage;

/// Relative tolerance on the element balance of one reaction,
/// `|sum_i nu_i a_ik| <= tol * max(1, sum_i |nu_i a_ik|)`.
pub(crate) const ELEMENT_BALANCE_TOLERANCE: f64 = 1e-9;
/// Relative tolerance on the closure of total mass across a yield reactor.
/// Conversion reactors conserve elements instead: the tabulated molar masses
/// are not atom-consistent beyond about 1e-6 relative.
pub(crate) const MASS_BALANCE_TOLERANCE: f64 = 1e-8;
/// Scaled element imbalance above which a yield reactor reports a warning.
pub(crate) const YIELD_ELEMENT_WARNING: f64 = 1e-6;

/// Element, formation and mass data of the package components, in package order.
#[derive(Debug, Clone)]
pub(crate) struct Chemistry {
    ids: Vec<String>,
    formation_j_mol: Vec<f64>,
    molar_mass_kg_mol: Vec<f64>,
    elements: Vec<String>,
    /// `atoms[i][k]`: atoms of element `k` in one molecule of component `i`.
    atoms: Vec<Vec<f64>>,
}

impl Chemistry {
    pub(crate) fn of(package: &CubicPackage) -> Self {
        let components = package.model().components();
        let mut elements: Vec<String> = Vec::new();
        let parsed: Vec<_> = components.iter().map(|c| c.elements()).collect();
        for (symbol, _) in parsed.iter().flatten() {
            if !elements.contains(symbol) {
                elements.push(symbol.clone());
            }
        }
        let atoms = parsed
            .iter()
            .map(|counts| {
                elements
                    .iter()
                    .map(|e| {
                        counts
                            .iter()
                            .find(|(s, _)| s == e)
                            .map_or(0.0, |(_, n)| f64::from(*n))
                    })
                    .collect()
            })
            .collect();
        Self {
            ids: components.iter().map(|c| c.id().to_owned()).collect(),
            formation_j_mol: components
                .iter()
                .map(|c| c.standard_formation_enthalpy())
                .collect(),
            molar_mass_kg_mol: components.iter().map(|c| c.molar_mass().value()).collect(),
            elements,
            atoms,
        }
    }

    pub(crate) fn molar_masses(&self) -> &[f64] {
        &self.molar_mass_kg_mol
    }

    /// Atoms of element `k` per molecule of each component.
    pub(crate) fn element_column(&self, k: usize) -> Vec<f64> {
        self.atoms.iter().map(|a| a[k]).collect()
    }

    pub(crate) fn elements(&self) -> &[String] {
        &self.elements
    }

    pub(crate) fn index(&self, id: &str) -> Result<usize, String> {
        self.ids
            .iter()
            .position(|c| c == id)
            .ok_or_else(|| format!("component {id:?} is not in the property package"))
    }

    /// Dense coefficients in package order; rejects unknown components and
    /// reactions that do not conserve every element.
    pub(crate) fn stoichiometry(
        &self,
        number: usize,
        coefficients: &[(String, f64)],
    ) -> Result<Vec<f64>, String> {
        let mut nu = vec![0.0; self.ids.len()];
        for (id, v) in coefficients {
            nu[self
                .index(id)
                .map_err(|e| format!("reaction {number}: {e}"))?] += v;
        }
        for (k, element) in self.elements.iter().enumerate() {
            let (net, scale) = nu
                .iter()
                .zip(&self.atoms)
                .fold((0.0, 0.0), |(net, scale), (v, a)| {
                    (net + v * a[k], scale + (v * a[k]).abs())
                });
            if net.abs() > ELEMENT_BALANCE_TOLERANCE * scale.max(1.0) {
                return Err(format!(
                    "reaction {number} does not conserve {element}: {net:+.6} atoms per unit extent"
                ));
            }
        }
        Ok(nu)
    }

    /// Standard formation enthalpy flow of a generation vector, W.
    pub(crate) fn formation_flow_w(&self, generation: &[f64]) -> f64 {
        generation
            .iter()
            .zip(&self.formation_j_mol)
            .map(|(g, h)| g * h)
            .sum()
    }

    /// Mass flow of component flows, kg/s.
    pub(crate) fn mass_flow_kg_s(&self, flows: &[f64]) -> f64 {
        flows
            .iter()
            .zip(&self.molar_mass_kg_mol)
            .map(|(n, m)| n * m)
            .sum()
    }

    /// Atom flow per element, mol/s.
    pub(crate) fn element_flows(&self, flows: &[f64]) -> Vec<f64> {
        (0..self.elements.len())
            .map(|k| flows.iter().zip(&self.atoms).map(|(n, a)| n * a[k]).sum())
            .collect()
    }
}

/// Reacting terms of one reactor calculation, retained with its result so
/// that recycle closure can re-evaluate the balance.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionTerms {
    /// Net generation per component in package order, mol/s
    /// (negative when consumed).
    pub generation_mol_s: Vec<f64>,
    /// Extent of each conversion-reactor reaction, mol/s; empty for a yield
    /// reactor, which has no stoichiometry.
    pub extents_mol_s: Vec<f64>,
    /// `sum_i g_i dHf_i` at 298.15 K, W: the standard heat of reaction of the
    /// generation, negative when exothermic.
    pub standard_reaction_enthalpy_w: f64,
    /// Whether the reactor conserves elements (stoichiometric, equilibrium,
    /// Gibbs and kinetic reactors) or only mass (yield reactor).
    pub element_balanced: bool,
}

/// Reacting part of a reactor [`crate::Balance`].
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionBalance {
    /// The terms the component and energy residuals account for.
    pub terms: ReactionTerms,
    /// Outlet minus inlet mass flow, kg/s.
    pub mass_residual_kg_s: f64,
    /// Element symbols of the package components, in first-appearance order.
    pub elements: Vec<String>,
    /// Outlet minus inlet atom flow per element, mol/s.
    pub element_residuals_mol_s: Vec<f64>,
}

/// Applies the reactions in order to the inlet component flows.
pub(crate) fn convert(
    chemistry: &Chemistry,
    reactions: &[Reaction],
    inlet: &[f64],
) -> Result<(Vec<f64>, ReactionTerms), String> {
    let total: f64 = inlet.iter().sum();
    let mut flows = inlet.to_vec();
    let mut extents = Vec::with_capacity(reactions.len());
    for (k, r) in reactions.iter().enumerate() {
        let number = k + 1;
        let nu = chemistry.stoichiometry(number, &r.coefficients)?;
        let extent = match &r.extent {
            ReactionExtent::Conversion {
                key: Some(key),
                fraction: Some(x),
            } => {
                let j = chemistry.index(key)?;
                if nu[j] >= 0.0 {
                    return Err(format!("reaction {number}: key {key} is not a reactant"));
                }
                x * flows[j] / -nu[j]
            }
            ReactionExtent::MolarExtent(Some(e)) => *e,
            _ => {
                return Err(format!(
                    "reaction {number}: incomplete extent specification"
                ));
            }
        };
        for (n, v) in flows.iter_mut().zip(&nu) {
            *n += v * extent;
        }
        // Complete conversion of a key can leave a rounding residue of either
        // sign; anything beyond it is a genuine shortage.
        for (i, n) in flows.iter_mut().enumerate() {
            if *n < 0.0 {
                if *n >= -1e-12 * total {
                    *n = 0.0;
                } else {
                    return Err(format!(
                        "reaction {number} consumes more {} than is available ({:.6e} mol/s short)",
                        chemistry.ids[i], -*n
                    ));
                }
            }
        }
        extents.push(extent);
    }
    let generation: Vec<f64> = flows.iter().zip(inlet).map(|(o, i)| o - i).collect();
    let terms = ReactionTerms {
        standard_reaction_enthalpy_w: chemistry.formation_flow_w(&generation),
        generation_mol_s: generation,
        extents_mol_s: extents,
        element_balanced: true,
    };
    Ok((flows, terms))
}

/// Distributes the mass of the non-inert inlet components over the yield
/// components; inerts pass through. Returns outlet flows, terms and warnings.
pub(crate) fn yields(
    chemistry: &Chemistry,
    mass_yields: &[(String, f64)],
    inerts: &[String],
    inlet: &[f64],
) -> Result<(Vec<f64>, ReactionTerms, Vec<String>), String> {
    let mut inert = vec![false; inlet.len()];
    for id in inerts {
        inert[chemistry.index(id)?] = true;
    }
    let mut weight = vec![0.0; inlet.len()];
    for (id, y) in mass_yields {
        weight[chemistry.index(id)?] += y;
    }
    let weight_sum: f64 = weight.iter().sum();
    if !(weight_sum.is_finite() && weight_sum > 0.0) {
        return Err("mass yields must have a positive finite sum".into());
    }
    let reacting_mass: f64 = inlet
        .iter()
        .zip(&chemistry.molar_mass_kg_mol)
        .zip(&inert)
        .filter(|(_, inert)| !**inert)
        .map(|((n, m), _)| n * m)
        .sum();
    let mut warnings = Vec::new();
    let flows: Vec<f64> = (0..inlet.len())
        .map(|i| {
            if inert[i] {
                inlet[i]
            } else {
                reacting_mass * weight[i] / weight_sum / chemistry.molar_mass_kg_mol[i]
            }
        })
        .collect();
    if reacting_mass == 0.0 {
        warnings.push(
            "No reacting mass: the feed contains only inert components and passes through unchanged."
                .into(),
        );
    }
    let generation: Vec<f64> = flows.iter().zip(inlet).map(|(o, i)| o - i).collect();
    let net = chemistry.element_flows(&generation);
    let entering = chemistry.element_flows(inlet);
    let leaving = chemistry.element_flows(&flows);
    let unbalanced: Vec<String> = chemistry
        .elements
        .iter()
        .enumerate()
        .filter_map(|(k, e)| {
            let scale = entering[k].max(leaving[k]);
            (scale > 0.0 && net[k].abs() > YIELD_ELEMENT_WARNING * scale)
                .then(|| format!("{e} {:+.3} %", 100.0 * net[k] / scale))
        })
        .collect();
    if !unbalanced.is_empty() {
        warnings.push(format!(
            "Mass yields conserve mass but not elements (outlet minus inlet atoms: {}); adjust the yields if the products must be atom-balanced.",
            unbalanced.join(", ")
        ));
    }
    let terms = ReactionTerms {
        standard_reaction_enthalpy_w: chemistry.formation_flow_w(&generation),
        generation_mol_s: generation,
        extents_mol_s: Vec::new(),
        element_balanced: false,
    };
    Ok((flows, terms, warnings))
}
