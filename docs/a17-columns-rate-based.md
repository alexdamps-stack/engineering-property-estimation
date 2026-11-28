# A17 columns: rate-based packed sections (#241)

Part of the columns group; see [a17-implementation-brief.md](a17-implementation-brief.md)
§ "Columns and reactors as packed beds" and [a17-columns-packing.md](a17-columns-packing.md)
for the equilibrium (HETP) path. Code:
`crates/flowsheet-core/src/units/columns/{packing,packing_data,rate_based,absorber}.rs`
and `crates/thermo-models/src/diffusivity.rs`.

Status: the **adiabatic absorber / stripper** has a rate-based mode. The
rigorous distillation column's packed bed is not yet rate-based (next step,
below).

## Calculation mode

`AbsorberSpec::calculation` is a `PackedCalculation`:

* `Equilibrium` (default): ideal trays or an HETP bed, unchanged. Records
  without the new field load as this mode.
* `RateBased(RateBasedSpec)`: needs packed internals with a **height** basis
  (the HETP field is not used; trays or the stage basis report
  "packed height is missing").

| `RateBasedSpec` field | Meaning | Validation |
|---|---|---|
| `diameter_m` | column inside diameter, m | required, > 0 |
| `segments` | mixed-flow segments over the height | 1–24 (`MAX_RATE_SEGMENTS`) |
| `packing` | library entry | see below |
| `correlation` | mass-transfer correlation | `Onda1968` |
| `film_model` | film model | `EffectiveDiffusivity` |
| `interfacial_area_factor` | multiplier of the correlated area | 0 < f ≤ 1e6, default 1 |
| `gas_transport_neglected`, `liquid_transport_neglected` | components left out of that phase's transport properties | package ids, no duplicates; each ≤ 1 mol% in the phase |

Project records: `rate_based` on the absorber record (snake-case packing,
correlation and film-model names), behind format version 4. Results:
`AbsorberProfile::rate_based` (`RateBasedSegmentsProfile`), absent for the
equilibrium mode and for older results.

## Segment model

The bed of height \(Z\) is divided into \(n\) equal segments of height
\(\Delta z\), numbered from the top. Each segment is a mixed-flow cell (the
"mixed" flow model of RateSep-type nonequilibrium models; Krishnamurthy and
Taylor, AIChE J. 31 (1985) 449): bulk vapor and bulk liquid are at the
segment's outlet states, and an interface at phase equilibrium separates two
films. Unknowns per segment: vapor and liquid component flows, \(T_V\),
\(T_L\), \(T_I\), interface \(x_I\) and \(y_I\), and component fluxes
\(N_i\) (vapor → liquid).

* Component balances of each phase:
  \(v_{i,j+1} - v_{i,j} - N_{i,j} = 0\), \(l_{i,j-1} - l_{i,j} + N_{i,j} = 0\).
* Energy balances of each phase with the film energy fluxes
  \(E_V = h_V A (T_V - T_I) + \sum_i N_i \bar H_i^V(T_I, y_I)\),
  \(E_L = h_L A (T_I - T_L) + \sum_i N_i \bar H_i^L(T_I, x_I)\), and
  \(E_V = E_L\) at the interface. Partial molar enthalpies are forward
  differences of the package phase enthalpy in mole numbers.
* Interface equilibrium \(\ln y_{I,i} - \ln x_{I,i} = \ln K_i(T_I, x_I, y_I)\)
  with the package fugacity coefficients.
* **Film model (effective diffusivity)**, for every component except one
  reference per phase (the most abundant in that phase's feed):
  \(N_i = c_V k_{V,i} A (y_i - y_{I,i}) + \bar y_i N_t\),
  \(N_i = c_L k_{L,i} A (x_{I,i} - x_i) + \bar x_i N_t\), with
  \(N_t = \sum N_i\), \(\bar y = (y + y_I)/2\) and \(A = f a_w A_c \Delta z\).
  The bulk-flow term carries the reference components (for example the
  carrier gas). This is a documented first step: Maxwell–Stefan
  multicomponent films with matrix coefficients (Taylor and Krishna,
  *Multicomponent Mass Transfer*, Wiley 1993) are not implemented, and there
  is no high-flux correction.
* **Heat transfer.** Chilton–Colburn analogy in the vapor,
  \(h_V = \bar k_V c_V c_{p,V}\,\mathrm{Le}_V^{2/3}\), and the penetration
  analogy in the liquid, \(h_L = \bar k_L c_L c_{p,L}\,\mathrm{Le}_L^{1/2}\)
  (Taylor and Krishna 1993), with \(\mathrm{Le} = \lambda/(c\,c_p \bar D)\)
  and \(\bar k\) the correlation's coefficient at the mean effective
  diffusivity \(\bar D\) of the non-reference components.

Concentrations: vapor \(1/v\) from the package EOS; liquid from the
package's **engineering liquid density** (#228), which fails explicitly when
unavailable. Mass fluxes use the phase molar masses.

### Solver

All segment unknowns are solved simultaneously by Newton's method with a
finite-difference Jacobian assembled segment by segment (a perturbation of
segment \(j\) changes only the residuals of segments \(j-1\) to \(j+1\)),
dense LU and a max-norm line search, to a scaled residual of 1e-10. Seed:
the equilibrium-stage absorber with one ideal stage per segment (the
infinite-rate limit of the segment equations). If the package's PH flash of
the feed blend fails, the seed cascade starts from a TP flash of the blend at
the flow-weighted inlet temperature instead (a starting point only; the
equilibrium mode itself is unchanged). If the direct solve from the seed
fails, a continuation in the area factor (10⁴, 10³, 10², 10, 1 times the
specified value) is used. Nonconvergence is an explicit error.

After convergence every interface is checked: the tie-line midpoint must
flash to the solved vapor and liquid (package TP), and the package's VLLE
flash must not find a second liquid (warning when unsupported), as for the
equilibrium stages.

## Correlations and data

| Item | Source |
|---|---|
| Gas binary diffusivity | Fuller, Schettler and Giddings (1966), diffusion volumes of Fuller, Ensley and Giddings (1969), as Poling, Prausnitz and O'Connell (2001) eq. 11-4.4 / Table 11-1 |
| Gas mixture diffusivity | Blanc's law, Poling et al. (2001) |
| Liquid diffusivity | Wilke and Chang (1955), φ = 2.6 (water), 1.9 (methanol), 1.0 otherwise; mixed solvent φM = Σ x_j φ_j M_j / Σ x_j (Poling et al. 2001 §11-10); liquid-mixture viscosity |
| Solute molar volume at the normal boiling point | Tyn and Calus (1975), \(V_b = 0.285 V_c^{1.048}\), Vc from the ChemSep transport extract |
| Wetted area, k_L, k_G | Onda, Takeuchi and Okumoto (1968), J. Chem. Eng. Japan 1, 56; σc of steel 0.075 N/m; gas constant 5.23 above 15 mm |
| Pressure drop | Robbins (1991), Chem. Eng. Prog. 87(5) 87 (analytic GPDC; C3 = 7.4e-8, C4 = 2.7e-5, US units) |
| Flooding | Kister and Gill (1991) flood pressure drop 0.115 Fp^0.7 inH2O/ft; fraction = G / G_flood at the same L |
| Minimum wetting | Morris and Jackson (1953), 0.08 m³/(m h) for rings ≤ 75 mm (Sinnott, Coulson & Richardson Vol. 6, §11.14) |
| Viscosity, conductivity, surface tension | the package transport interface (#225) |

Diffusivity domain (`thermo_models::DiffusivityModel`): the eight transport
components (others `MissingParameter`), gas 200–600 K and at most 2 MPa (the
Fuller method is a low-pressure correlation; no dense-gas correction),
liquid 250–500 K; outside it `OutsideValidityDomain`.

### Packing library

| Packing | Size | a, m²/m³ | ε | Fp, 1/m | Source |
|---|---|---|---|---|---|
| Pall rings, metal | 25 mm | 215 | 0.95 | 157 | Raschig GmbH, *RASCHIG PALL-Ring* brochure 350 (v2025010), operational data and tower packing factor Table 1 |
| Pall rings, metal | 50 mm | 105 | 0.96 | 66 | same |

Not included: ceramic Raschig rings and structured packings (Mellapak-type)
with Bravo–Rocha–Fair, and Billet–Schultes (1999). Their packing-specific
constants (Billet–Schultes C_h, C_P, C_L, C_V, C_Fl; corrugation geometry S,
B, h for Bravo–Rocha–Fair) were not taken from a verified source in this
batch; they are not invented here.

### Hydraulics

Per segment: wetted area, fraction of flood and Robbins pressure drop at the
segment's loads. The column is solved at **one pressure**: the total bed
pressure drop is reported only and does not change the segment pressures.
The Robbins correlation is used with the library's tower packing factor in
place of Robbins's dry packing factor (the bulletin gives one factor);
this is an approximation of the hydraulics, not of the mass transfer.

## Explicit failures

* flooding: a segment at or above the Kister–Gill flood pressure drop;
* wetting below the Morris–Jackson minimum;
* transport properties unavailable through the package (#225 domain, for
  example dissolved light gases in a liquid above 0.75 Tc or a wet gas's
  conductivity) — unless the component is listed as neglected and stays at
  or below 1 mol%; a listed component above that fails;
* engineering liquid density or diffusivity unavailable;
* temperatures outside the package bubble/dew interval, missing distinct
  roots, interface not reproduced by the package flash or a second liquid;
* nonconvergence of the seed or the segment equations.

**Trace transport basis.** The #225 transport domain rejects any liquid
containing a trace of a supercritical component and the conductivity of any
gas containing water or methanol, so no practical absorber has all its
transport properties. The neglected-component lists are an explicit,
opt-in, reported approximation: the listed components are removed from that
phase's transport-property evaluation (the rest renormalised; the phase's
own EOS state at that composition), each must stay at or below 1 mol% in
the phase, and the result reports the largest neglected fraction and a
warning. Diffusivities and densities always use the actual phase.

## Validation

* `validation/columns/rate-based-v1.json` from
  `validation/columns/generate-rate-based-references.py` (Python standard
  library + NumPy/SciPy with the teqp/Cantera, transport and liquid-density
  reference adapters; never calls Rust). Two nitrogen/methanol/water
  scrubbing columns (PR, 50 mm Pall rings, 2 m, 6 segments; PR-BM, 25 mm
  Pall rings, 1.5 m, 5 segments). The generator uses a different
  formulation (bulk mole fractions and total flows with summations, partial
  molar enthalpies from \(\partial \ln\varphi/\partial T\) and NASA7) and
  checks both phases of every interface with the independent multi-start
  TPD stability search. Budgets are declared
  in the generator before any comparison.
  `flowsheet-project/tests/rate_based_absorber.rs` compares every segment
  temperature, flow, composition, flux, area, flooding fraction and pressure
  drop, the outlets, balance closure and the project round trip.
* Limiting behaviour: with an area factor of 1e5 the nitrogen/CO2/water wash
  reproduces the equilibrium-stage absorber with one tray per segment
  (temperatures within 2 mK, CO2 slip within 1e-4 relative); with factor 1
  less CO2 is absorbed.
* Correlation checks: Fuller CO2–N2 at 298 K against the measured 0.167
  cm²/s (Cussler, *Diffusion*, Table 5.1-1) within 5 %; Wilke–Chang CO2 in
  water against about 1.9e-9 m²/s within 20 %; Robbins dry-bed form and the
  flooding inversion; Onda magnitudes. No published worked example of the
  complete Onda/Robbins chain was reproduced (not available in this
  environment); the reference fixture checks implementation agreement only.
* Failure tests: flooding, minimum wetting, strict transport domain, a
  neglected component above 1 mol%, an unknown neglected component, spec
  issues; legacy records load as the equilibrium mode.

## Next steps

* Rate-based packed bed of the rigorous distillation column (#186): the same
  segment equations with condenser/reboiler end stages and feeds by
  elevation (`PackedStages::stage_at_depth` analogue for segments).
* Billet–Schultes and Bravo–Rocha–Fair with verified packing constants,
  structured packings, Maxwell–Stefan films, high-flux corrections, pressure
  profile from the computed pressure drop, HETP from the rate model.
