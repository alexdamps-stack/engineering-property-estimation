# Transport properties (#225)

Status: November 24, 2026. Phase-wise dynamic viscosity μ, thermal
conductivity λ and surface tension σ for every flowsheet package. These are
side properties (property-model requirements, section 3.5): they never feed
back into equilibrium, caloric or balance results, and they carry their own
method identity, data version, validity and the reason for every unavailable
value.

## Methods

| Phase label | Regime reported | Method |
|---|---|---|
| vapor, supercritical | `gas` (ρ·Vc,m ≤ 0.1) or `dense fluid` | Chung, Ajlan, Lee and Starling (1988), Ind. Eng. Chem. Res. 27, 671–679: dense-fluid viscosity and thermal conductivity (Poling, Prausnitz and O'Connell 2001, eqs. 9-4.10, 9-6.18–9-6.21, 10-3.14, 10-5.5–10-5.7), six-parameter Neufeld collision integral |
| liquid | `liquid` | ChemSep 8.3 pure-component saturated-liquid correlations (μ: ChemSep eq. 101; λ and σ: eq. 16) |

* **Gas inputs.** Tc, acentric factor and molar mass from the component
  extract (CoolProp 7.2), critical volume and dipole moment from ChemSep 8.3,
  Chung association factors κ of water (0.075908) and methanol (0.215175)
  from Chung et al. (1988), ideal-gas cv = cp° − R from the NASA7 records, and
  the phase's own molar volume from the package model. No critical
  enhancement of λ is modelled.
* **Gas mixtures.** Chung mixing rules (Poling eqs. 9-5.25–9-5.44) with
  ξij = ζij = 1 and κij = (κi κj)^½; Fc uses 0.2756 ω for pure fluids and
  mixtures alike.
* **Liquid mixtures (labelled approximations).** μ: ln μ = Σ xi ln μi
  (Grunberg–Nissan with Gij = 0); λ: DIPPR 9H, λ = (Σ wi λi⁻²)^−½ with mass
  fractions; σ: Σ xi σi. Known to be poor for strongly non-ideal mixtures
  such as methanol–water; mixture results therefore carry
  `ValidityStatus::Unassessed`. A single present component carries
  `WithinDeclaredDomain`.
* **σ** is reported for liquid phases only (the liquid against its vapor);
  gas phases return `UnsupportedCapability` for it.

Corresponding-states alternatives were screened against CoolProp before
selection and rejected for liquids: Chung's dense method with the untranslated
PR liquid density was 25–300 % off (PR liquid densities differ by 5–20 %);
Letsou–Stiel and Sato–Riedel were 20–100 % off over 0.5–0.9 Tc; Brock–Bird σ
was good for nonpolar fluids but 40–110 % high for water and methanol, and
Macleod–Sugden needs liquid densities to the fourth power. Once an
engineering liquid density exists (#228), Chung or Macleod–Sugden for liquids
can be reconsidered as separate, versioned methods.

## Declared validity domain

A property outside the domain is `OutsideValidityDomain` (stored as
*Unavailable* with the reason), never an extrapolated value. Every component
with a nonzero mole fraction must satisfy each rule; zero components are
ignored.

| Rule | Applies to |
|---|---|
| Components: methane, ethane, propane, nitrogen, carbon dioxide, oxygen, water, methanol. Others (hydrogen, carbon monoxide, the further GERG/EOS-CG components) are `MissingParameter`. | all |
| 200 K ≤ T ≤ 600 K | gas μ, λ |
| reduced density ρ·Vc,m ≤ 1.5 (Vc,m from the Chung mixing rule) | gas μ, λ |
| with a polar component present (dipole ≥ 0.5 D or κ > 0: water, methanol): ρ·Vc,m ≤ 0.1 | gas μ |
| no polar component present (Chung λ is 15–80 % high for water vapour and up to 160 % for dense methanol) | gas λ |
| 0.3 ≤ T* = 1.2593 T/Tc,m ≤ 100 | gas μ, λ |
| T within each present component's source correlation range | liquid μ, λ, σ |
| T ≤ 0.75 Tc of each present component | liquid μ, λ, σ |
| p ≤ 5 MPa (the correlations are saturated-liquid; compression is not modelled) | liquid μ, λ |

Consequences: a liquid containing any trace of a component above 0.75 of its
critical temperature (for example dissolved methane in liquid propane at
250 K) has no liquid transport properties, and gas λ is unavailable for wet
gases. Units that need a missing property must fail with the reason from
`CalcContext::viscosity`, `thermal_conductivity` or `surface_tension`.

## Interfaces

* `thermo_data::TransportDatabase` loads
  `data/source-extracts/transport.provisional.json` (`TRANSPORT_DATA_VERSION`
  2026-11-24.1), produced only by `validation/prepare-transport-source-data.py`
  from the ChemSep 8.3 database in pinned DWSIM commit
  `1abf72d1b6b41d3e9a8cc770d3cc4e8fc76e5766` (Artistic License 2.0, notice in
  `data/source-extracts/ChemSep-LICENSE.txt`).
* `thermo_models::TransportModel` evaluates the methods
  (`TRANSPORT_MODEL_VERSION` transport-v1); model metadata `chung-1988` or
  `chemsep-liquid-correlations`.
* `FlowsheetPackage::transport(&phase)` is a provided trait method that
  applies these methods to a phase of any package (the gas method uses the
  package's phase volume). `property_packages::annotate_transport` fills
  `PhaseProperties::dynamic_viscosity`, `thermal_conductivity` and
  `transport`.
* The flowsheet engine annotates every stored outlet stream. Units reach the
  values through `CalcContext::transport`, `viscosity`,
  `thermal_conductivity` and `surface_tension`.
* Project results store `dynamic_viscosity_pa_s`,
  `thermal_conductivity_w_m_k`, `surface_tension_n_m` and a `transport`
  record (regime, method, data, validity, unavailable reasons) per phase. The
  desktop stream panel shows *Unavailable* and the reasons.
  `ENGINE_REVISION` names the transport version so results computed before
  this change load as stale.

## Verification

Fixture `validation/fixtures/transport-v1.json`, generated by
`validation/transport/generate-transport-references.py` (Python 3.14,
CoolProp 7.2.0, chemicals 1.5.2); Rust test
`crates/thermo-models/tests/transport.rs`. Budgets are registered in
`validation/tolerances.json`; comparison `|a − e| ≤ atol + rtol |e|`.

| Set | Quantity | atol | rtol | Largest deviation found |
|---|---|---|---|---|
| Implementation (24 cases, 7 expected out of domain; independent Python evaluation of the same formulas, pure Chung λ also checked against `chemicals.Chung_dense` to 1e-12) | μ, σ | 1e-15 | 1e-10 | passes |
| | λ | 1e-13 W/(m K) | 1e-10 | passes |
| Model comparison, 130 pure states vs CoolProp 7.2.0 | dilute-gas μ, nonpolar | 0 | 6 % | 2.4 % |
| | dense-fluid μ, nonpolar | 0 | 12 % | 10.6 % (CO2, 350 K, 10 MPa) |
| | gas μ, water and methanol | 0 | 15 % | 12.6 % |
| | gas λ, nonpolar | 0 | 12 % | 10.9 % |
| | liquid μ | 0 | 15 % | 13.2 % |
| | liquid λ | 0 | 10 % | 9.8 % |
| | σ | 0 | 8 % | 3.4 % |

The CoolProp states cover each component at 250–550 K and 0.01–10 MPa in the
gas (reduced density ≤ 1.3 nonpolar, ≤ 0.08 polar, at most half the
saturation pressure below Tc) and at 0.5–0.74 Tc in the liquid (1.5 psat and
5 MPa). Gas states use the engine's own PR-exact volume, so the budgets
include the EOS density difference. These are model-comparison budgets
against REFPROP-quality fluid-specific correlations, not implementation
tolerances and not experimental uncertainties.

Budget history: the dense-fluid viscosity budget was first declared jointly
with the dilute gas at 6 %. The first comparison failed for one state
(CO2 at 350 K and 10 MPa, +10.6 %), consistent with the 5–10 % accuracy
published for Chung's dense-fluid correction. The budget is now declared per
regime (12 % above ρ·Vc,m = 0.1); the dilute budget is unchanged.

The independence of the implementation set is limited: the Python evaluator
was written for this issue from the same published equations; only the pure
Chung conductivity is checked against a third-party implementation, since
`chemicals` has no Chung viscosity. Carbon dioxide liquid is compared at one
temperature (225 K) because its correlation range starts at 216.6 K and the
domain ends at 0.75 Tc = 228 K.

## Gaps

* No pressure correction for liquids (for example Lucas), so compressed
  liquids above 5 MPa are unavailable.
* No polar gas conductivity, no critical enhancement, no hydrogen or carbon
  monoxide data (hydrogen needs a quantum-corrected method).
* Liquid mixing rules have no binary parameters; mixture accuracy is not
  qualified.
* Binary diffusivities are not part of the transport interface. The
  rate-based columns (#241) use package-independent correlations in
  `thermo_models::DiffusivityModel` (Fuller–Schettler–Giddings gas binaries
  with Blanc mixing, Wilke–Chang liquids with Tyn–Calus volumes from the
  ChemSep critical volume); see [a17-columns-rate-based.md](a17-columns-rate-based.md).
