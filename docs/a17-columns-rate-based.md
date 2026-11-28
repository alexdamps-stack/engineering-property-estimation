# A17 columns: rate-based packed sections (#241)

Part of the columns group; see [a17-implementation-brief.md](a17-implementation-brief.md)
§ "Columns and reactors as packed beds" and [a17-columns-packing.md](a17-columns-packing.md)
for the equilibrium (HETP) path. Code:
`crates/flowsheet-core/src/units/columns/{packing,packing_data,rate_based,absorber}.rs`
and `crates/thermo-models/src/diffusivity.rs`.

Status: the **adiabatic absorber / stripper** has a rate-based mode, limited
to **at most 3 bar absolute** by the stated domain of the Robbins
pressure-drop correlation. The rigorous distillation column's packed bed is
not yet rate-based (next step, below).

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
| `couple_pressure_drop` | solve the absolute pressure profile from the specified bed top pressure and the Robbins pressure drop | opt-in, default false |
| `interfacial_area_factor` | multiplier of the correlated area | 0 < f ≤ 1e6, default 1 |
| `critical_surface_tension_n_m` | σc of the packing material (Onda wetted area), N/m | **required user input**, > 0 |
| `wetting_surface` | packing surface class of the minimum-wetting rule | Perry's Table 14-16 classes, default bright metal |
| `gas_transport_neglected`, `liquid_transport_neglected` | components left out of that phase's transport properties (approximation, opt-in) | package ids, no duplicates; each ≤ 1 mol% in the phase |

Project records: `rate_based` on the absorber record (snake-case packing,
correlation, film-model and surface names), behind format version 4.
Results: `AbsorberProfile::rate_based` (`RateBasedSegmentsProfile`), absent
for the equilibrium mode and for older results. Rate-based results include
absolute pressures at segment boundaries. With fixed pressure, the profile is
an estimate obtained by adding the correlated pressure drops downward from
the specified pressure. All thermodynamic and transport calculations and both
outlet states use the specified pressure in this mode.

With coupled pressure, each segment's bulk vapor uses the upper boundary
pressure, its bulk liquid uses the lower boundary pressure, and interface
equilibrium and the Robbins pressure correction use the mean of those
pressures. Both outlet streams use their corresponding bed boundary pressure.
Inlet enthalpy flows enter the existing adiabatic energy balances at their
stream values. The gas feed must supply at least the computed bottom pressure,
the liquid feed must supply at least the specified top pressure, and the
whole bed must stay within the 3 bar absolute Robbins limit.

The solver starts at fixed pressure. It then integrates the segment pressure
drops from the top down, resolves the transfer equations at those boundary
pressures, and repeats until each boundary pressure changes by at most
1e-6 Pa. It reports failure if the profile fails to converge within 30 updates.

## Segment model

The bed of height \(Z\) is divided into \(n\) equal segments of height
\(\Delta z\), numbered from the top. Each segment is a mixed-flow cell (the
"mixed" flow model of nonequilibrium-stage models after Krishnamurthy and
Taylor, 1985): bulk vapor and bulk liquid are at the
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
  This is a documented first step: no Maxwell–Stefan multicomponent films,
  no high-flux correction.
* **Heat transfer.** Chilton–Colburn analogy \(j_H = j_D\) in both phases,
  \(h = \bar k\,c\,c_p\,\mathrm{Le}^{2/3}\), \(\mathrm{Le} = \lambda/(c\,c_p \bar D)\),
  with \(\bar k\) the correlation's coefficient at the mean effective
  diffusivity \(\bar D\) of the non-reference components. (The first version
  used a penetration analogy with \(\mathrm{Le}^{1/2}\) in the liquid; it was
  replaced because its source could not be checked.)

Concentrations: vapor \(1/v\) from the package EOS; liquid from the
package's **engineering liquid density** (#228), which fails explicitly when
unavailable. Diffusivities: gas Fuller binaries combined by Wilke's
stagnant-mixture rule; liquid Wilke–Chang with every component a dilute
solute in the liquid's reference component (the solvent) and the liquid's
viscosity as the solvent viscosity.

### Solver

All segment unknowns are solved simultaneously by Newton's method with a
finite-difference Jacobian assembled segment by segment, dense LU and a
max-norm line search, to a scaled residual of 1e-10 (film rate and energy
residuals are scaled by their own conductances so that the equilibrium limit
stays well scaled). Seed: the equilibrium-stage absorber with one ideal
stage per segment (the infinite-rate limit of the segment equations). If the
package's PH flash of the feed blend fails, the seed cascade starts from a
TP flash of the blend at the flow-weighted inlet temperature instead (a
starting point only). The nitrogen/methanol/water case in #243 now solves
in equilibrium mode using the corrected cubic PH flash. If the direct solve
fails, a continuation in the area factor (10⁴ … 1 times the specified value)
is used. Nonconvergence is an explicit error. After convergence every
interface is checked with the package TP flash of its tie-line midpoint and
the VLLE flash, as for equilibrium stages.

## Verification of constants

Every constant in the code was checked against the source listed, digit by
digit, from a copy that could be read in this work (2026-11-28). "Primary"
means the original publication; "secondary" means a reproduction that cites
the original.

| Constant / form | Value in code | Source read | Status |
|---|---|---|---|
| Onda wetted area, eq. (1) | \(a_w/a_t = 1-\exp[-1.45(\sigma_c/\sigma)^{0.75}Re_L^{0.1}Fr_L^{-0.05}We_L^{0.2}]\), \(Re_L = L/(a_t\mu_L)\), \(Fr_L = L^2a_t/(\rho_L^2 g)\), \(We_L = L^2/(\rho_L\sigma a_t)\) | Onda, Takeuchi, Okumoto, J. Chem. Eng. Japan 1(1) (1968) 56, eq. (1) p. 56 and nomenclature p. 61 (primary, <https://www.jstage.jst.go.jp/article/jcej1968/1/1/1_1_56/_article>) | verified |
| Onda k_L, eq. (2) | \(k_L(\rho_L/\mu_L g)^{1/3} = 0.0051(L/a_w\mu_L)^{2/3}Sc_L^{-1/2}(a_tD_p)^{0.4}\) | same, eq. (2) p. 57 | verified |
| Onda k_G, eq. (3) | \(k_GRT/(a_tD_G) = 5.23(G/a_t\mu_G)^{0.7}Sc_G^{1/3}(a_tD_p)^{-2.0}\); 2.00 for packings smaller than 15 mm | same, eq. (3) and text p. 58 | verified |
| Onda σc of steel (0.075 N/m) | **removed**: user input | the σc table is in Onda, Takeuchi, Koyama, Kagaku Kōgaku 31 (1967) 126, which could not be obtained | **not verified → user input** |
| Metal Pall rings 25 mm: a, ε, F_p, F_pd | 205 m²/m³, 0.94, 183 1/m, 174 1/m | Perry's Chemical Engineers' Handbook, 8th ed. (2008), Table 14-13, p. 14-60 (secondary; F_p from Kister 1992, F_pd from Robbins 1991) | verified |
| Metal Pall rings 50 mm: a, ε, F_p, F_pd | 105 m²/m³, 0.96, 89 1/m, 79 1/m | same | verified |
| Previous library values (Raschig brochure 350: 215/0.95/157 and 105/0.96/66) | **replaced** | Raschig "PALL-Ring 350" p. 3 and Table 1 p. 10: values confirmed, but they are wet "tower packing factors" of a different correlation basis, and no dry factor F_pd is given | replaced by Perry's Table 14-13 |
| Robbins ΔP, eqs. 14-145 to 14-151 | \(\Delta P = C_3G_f^210^{C_4L_f} + 0.4(L_f/20000)^{0.1}[C_3G_f^210^{C_4L_f}]^4\), C3 = 7.4e-8, C4 = 2.7e-5; \(G_f = 986F_s(F_{pd}/20)^{0.5}\), × \(10^{0.3\rho_G}\) above atmospheric pressure; \(L_f = L(62.4/\rho_L)(F_{pd}/20)^{0.5}\mu_L^{0.1}\) (\(\mu_L^{0.2}\) for F_pd > 200, \((20/F_{pd})^{0.5}\) for F_pd < 15) | Perry's 8th ed. p. 14-59 (secondary, citing Robbins, Chem. Eng. Prog. 87(5) 87, 1991) | verified; the pressure factor and F_pd regimes were **added** (missing before) and the dry factor F_pd is now used (the wet factor was used before) |
| Robbins domain | ≤ 3 bar absolute, **enforced** | Perry's 8th ed. p. 14-61 | verified; the previous 8–10 bar cases were outside it |
| Kister–Gill flood ΔP, eq. 14-142 | 0.12 F_p^0.7 inH2O/ft (F_p in 1/ft); not for F_p > 60 1/ft (**enforced**) | Perry's 8th ed. p. 14-58 (secondary, citing Kister and Gill, Chem. Eng. Prog. 87(2) 32, 1991) | verified; **changed from 0.115** |
| Flooding fraction with Robbins + Kister–Gill | G / G at flood ΔP | Perry's 8th ed. p. 14-61 ("For flood … Robbins recommends his pressure drop method together with Eq. (14-142)") | verified |
| Minimum wetting | Q_MW = Q_table (200/a_p)^0.5, Table 14-16 (bright metal 3.0, oxidized metal 0.7, …, m³/(m² h)), doubled for liquids above 25 mN/m | Perry's 8th ed. eq. 14-166 and Table 14-16, p. 14-79 (secondary, Glitsch Bulletin 345 and Kister 1992) | verified; **replaces** the Morris–Jackson 0.08 m³/(m h) rule, which could not be checked |
| Fuller–Schettler–Giddings | \(D = 10^{-3}T^{1.75}(1/M_A+1/M_B)^{1/2}/(P[\Sigma_{vA}^{1/3}+\Sigma_{vB}^{1/3}]^2)\), cm²/s, P in atm | Cussler, *Diffusion*, 3rd ed. (2009) Eq. 5.1-9 / Example 5.1-2, reproduced with citation of Fuller et al., Ind. Eng. Chem. 58(5) 18 (1966) in University of Jordan "Transport Phenomena II, 03 Diffusion coefficients" slides 19–21 (<https://www.madar-ju.com/storage/images/files/file_1731705111VrPpO.pdf>) | verified; **changed** from Poling's 0.00143/M_AB form with the 1969 volumes to the 1966 form and volumes |
| Fuller diffusion volumes | N2 17.9, O2 16.6, CO2 26.9, H2O 12.7; C 16.5, H 1.98, O 5.48 | same slides, slide 20 (FSG 1966 table) | verified; **changed** (were the 1969 set, with CO2 26.7 wrong even for that set, 26.9 per COMSOL's reproduction) |
| Wilke (1950) mixture rule | \(D_{i,m} = (1-y_i)/\sum_{j\ne i}y_j/D_{ij}\) | Getzinger, "Steady-state diffusion in ternary gas mixtures", LBNL 1962, eq. 13 p. 5 (<https://escholarship.org/uc/item/0w40t47j>), citing Wilke (1950) | verified (secondary) (previously labelled "Blanc's law") |
| Wilke–Chang | \(D = 7.4\times10^{-8}(\phi_BM_B)^{1/2}T/(\mu_BV_A^{0.6})\), cm²/s; φ water 2.6, methanol 1.9, others 1.0 | U. Jordan slides 34 and 37, citing Wilke and Chang, AIChE J. 1 (1955) 264 | verified; the unverified mixed-solvent φM rule was **removed** (single solvent = liquid reference component) |
| Tyn–Calus V_b | 0.285 V_c^1.048 (cm³/mol) | Thesis, Al-Nahrain University (<https://nahrainuniv.edu.iq/sites/default/files/Thesis_68.pdf>), eq. 2-36 §2.10.2, citing Tyn and Calus (1975) | verified (secondary) |
| Chilton–Colburn | \(j_H = (h/c_pG)Pr^{2/3} = j_D = (k_c/\bar v)Sc^{2/3}\) | Wikipedia "Chilton and Colburn J-factor analogy", citing Geankoplis 4th ed. p. 475 (secondary) | verified; liquid penetration analogy **removed** |
| Fuller domain ≤ 2 MPa, liquid 250–500 K | declared limits | — | declared (not literature constants) |

Published worked examples reproduced (in `thermo-models` and
`packing_data` unit tests and recorded in the fixture's
`published_examples`):

* Fuller: Cussler Example 5.1-2, H2–N2 at 294 K and 2 atm, 0.37 cm²/s
  (computed 0.3697); n-butanol–air at 273 K, 1 atm, 7.73e-6 m²/s (computed
  7.70e-6; the printed value carries intermediate rounding).
* Wilke–Chang: ethylbenzene in water at 293 K, 0.77e-5 cm²/s (computed 0.766e-5).
* Robbins: Perry's 8th ed. Example 13 (2-in metal Pall rings, air/water):
  G_f = 1653 (computed 1654.4; the example rounds F_s to 1.53), L_f = 9859
  (9859.0), ΔP 0.40 inH2O/ft read from the chart (equations 0.381).

No worked example of the Onda correlations with all inputs stated was
available; Onda is checked for implementation agreement only.

### Packing library

| Packing | Size | a, m²/m³ | ε | F_p, 1/m | F_pd, 1/m |
|---|---|---|---|---|---|
| Pall rings, metal | 25 mm | 205 | 0.94 | 183 | 174 |
| Pall rings, metal | 50 mm | 105 | 0.96 | 89 | 79 |

Source: Perry's 8th ed. Table 14-13. Not included: ceramic Raschig rings
(F_p 587 1/m exceeds the Kister–Gill limit of 60 1/ft), structured packings
(Bravo–Rocha–Fair geometry not verified) and Billet–Schultes (constants not
verified).

### Hydraulics

Per segment: wetted area, fraction of flood and Robbins pressure drop at the
segment's loads. Fixed pressure is the default, and coupled pressure applies
the drop to the segment equations. Specified pressures above 3 bar absolute and
packings above F_p = 60 1/ft fail before the segment solve.

## Explicit failures

* column pressure above 3 bar absolute (Robbins domain); F_p > 60 1/ft;
* in coupled mode, bed bottom pressure above 3 bar absolute or the gas inlet
  pressure, or pressure profile nonconvergence;
* flooding: a segment at or above the Kister–Gill flood pressure drop;
* liquid load below the minimum wetting rate;
* transport properties unavailable through the package (#225 domain) —
  unless the component is listed as neglected and stays at or below 1 mol%;
  a listed component above that fails;
* engineering liquid density or diffusivity unavailable;
* temperatures outside the package bubble/dew interval, missing distinct
  roots, interface not reproduced by the package flash or a second liquid;
* nonconvergence of the seed or the segment equations.

**Trace transport basis (approximation, opt-in).** The #225 transport
domain rejects any liquid containing a trace of a supercritical component
and the conductivity of any gas containing water or methanol, so no
practical absorber has all its transport properties. The neglected-component
lists are empty by default. When used, the listed components are removed
from that phase's viscosity, conductivity and surface tension (the rest
renormalised; the phase's own EOS state at that composition), each must
stay at or below 1 mol% in the phase, the unit publishes a warning naming
the components and the largest neglected fraction, and the result stores
that fraction (shown in the desktop). Diffusivities and densities always use
the actual phase.

## Validation

* `validation/columns/rate-based-v1.json` from
  `validation/columns/generate-rate-based-references.py` (NumPy/SciPy with
  the teqp/Cantera, transport and liquid-density reference adapters; never
  calls Rust). Two nitrogen/methanol/water scrubbing columns (PR, 50 mm Pall
  rings, 2.5 bar, 2 m, 6 segments; PR-BM, 25 mm Pall rings, 2 bar, 1.5 m,
  5 segments), σc = 0.075 N/m as a case input. The generator uses a
  different formulation (bulk mole fractions and total flows with
  summations, partial molar enthalpies from \(\partial \ln\varphi/\partial T\)
  and NASA7), checks both phases of every interface with the multi-start
  TPD stability search, and first reproduces the published worked examples
  above. Budgets are declared in the generator before any comparison.
  `flowsheet-project/tests/rate_based_absorber.rs` compares every segment
  temperature, flow, composition, flux, area, flooding fraction and
  pressure drop, the outlets, balance closure and the project round trip.
* Limiting behaviour: with an area factor of 1e5 the nitrogen/CO2/water wash
  at 3 bar reproduces the equilibrium-stage absorber with one tray per
  segment; with factor 1 less CO2 is absorbed.
* Failure tests: pressure above 3 bar, flooding, minimum wetting, strict
  transport domain, a neglected component above 1 mol%, an unknown
  neglected component, spec issues (including the missing σc); legacy
  records load as the equilibrium mode.

## Next steps

* Rate-based packed bed of the rigorous distillation column (#186).
* A verified σc table (Onda et al. 1967) as library data; Billet–Schultes
  and Bravo–Rocha–Fair with verified packing constants; structured packings;
  Maxwell–Stefan films;
  hydraulics above 3 bar (a correlation with a stated high-pressure domain).
