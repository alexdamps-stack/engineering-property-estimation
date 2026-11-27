# Petroleum pseudo-components — status (#229)

Status: November 24, 2026. This record follows the
[property-model requirements](property-model-requirements.md) for a
pseudo-component capability of the cubic packages (PR, PR-BM, SRK).

## What is implemented

| Layer | Content |
|---|---|
| `thermo-data` (`pseudo.rs`, `assay.rs`) | Stored input `PseudoComponentInput` (id, Tb, SG 60 °F/60 °F, optional measured M, optional specified critical constants with a source); the correlations; characterisation into a `Component`; `ComponentDatabase::with_pseudo_components`; ASTM D86 → TBP, VABP/MeABP, Watson K, TBP cutting into equal-temperature cuts at constant Watson K, and cut mole fractions. |
| `property-packages` (`pseudo.rs`) | `cubic_package_with_pseudo_components`: database components and pseudo-components in one `CubicPackage`, with a complete explicit interaction matrix. |
| `flowsheet-core` | Reactors reject pseudo-components in the preflight (no formula, no formation properties). |
| `flowsheet-project` | `PackageSelection::pseudo_components` (serde default, omitted when empty, unknown fields refused) and `PackageSelection::with_pseudo_components`. The stored inputs are part of the package fingerprint, so editing a cut makes results stale. |

## Correlations and sources

Coefficients are transcribed from M. A. Fahim, T. A. Al-Sahhaf, A. S. Elkilani,
*Fundamentals of Petroleum Refining*, Elsevier 2010 (ISBN 978-0-444-52785-1),
chapter 3.

| Quantity | Correlation | Book equation |
|---|---|---|
| D86 → TBP | API method of Riazi and Daubert (1980) | (3.3), Table 3.2 |
| VABP, MeABP | API method | (3.6), (3.7) |
| Watson K | K = (1.8 Tb)^(1/3)/SG | (3.5) |
| M (unless measured) | Riazi and Daubert (1987) | (3.11) |
| Tc, Pc | Kesler and Lee (1976), SI form | (3.25), (3.26) |
| ω | Lee and Kesler (1975) for Tb/Tc < 0.8; Kesler and Lee (1976) above | (3.27), (3.28) |
| Ideal-gas cp | Kesler and Lee (1976), quadratic in T | (3.38)–(3.41) |

Two print issues are recorded rather than copied:

- Eq. (3.25) prints the Tb coefficient as 0.422. The exact SI conversion of
  the original °R coefficients gives 0.4244, and only 0.4244 reproduces the
  book's own E3.11 value Tc = 765.5 K. Example E3.9 prints 764.17 K for the
  same oil, which neither coefficient reproduces; it is excluded with this reason.
- E3.11 prints an ideal-gas enthalpy change of 312,088.8 kJ/kmol between 298
  and 723.15 K. That value does not follow from the book's own A, B and C (they
  integrate to about 250,200 kJ/kmol), so it is excluded. The engine's integral
  is checked exactly against the closed form instead.

The cp correlation is implemented in its SI form (kJ/(kg K), T in K). The
book's per-kmol form folds conversion factors 4.1843 and 1.8 into the
coefficients; every coefficient agrees with the SI form to 2e-4 relative.

## Declared domain

| Quantity | Domain | Failure outside |
|---|---|---|
| Tb | 300–850 K | `OutsideValidityDomain` |
| SG | 0.60–1.00 | `OutsideValidityDomain` |
| M (correlated or measured) | 70–700 g/mol | `OutsideValidityDomain` |
| Ideal-gas cp, h, s | 250–900 K | `OutsideValidityDomain` |
| Package search interval | 250–600 K (cubic envelope ∩ cp domain) | the package is not built |

These are engine-declared limits that cover the worked examples and
kerosene/diesel cuts. They were not taken from the original papers, which
were not re-audited here, so they make no claim to be the published ranges.
Derived constants that are not physical (Tc ≤ Tb, Pc ≤ 1 atm, ω ≤ 0, or cp ≤ 0
anywhere in the cp domain) are `NumericalFailure`. Nothing is clipped or
extrapolated.

Pseudo-component identifiers start with `pseudo:`. Embedded records may not
use that prefix, so the two cannot collide. Pseudo-components have no
formation basis: `standard_gibbs_energy` returns `MissingParameter`, and
reactors refuse them before solving.

## Interactions, caloric reference and provenance

- **Interactions:** the cubic matrix must cover every pair, including
  pseudo-component pairs, with a non-empty source or approximation note.
  Missing pairs are an error. The tests use an explicitly labelled kij = 0
  approximation for hydrocarbon pairs (Fahim et al. 2010, sec. 3.6.1). No
  values are filled in by default.
- **Caloric reference:** unchanged (`ideal-pure-zero-298.15K-1bar-cp-integrated-v1`).
  The Kesler–Lee quadratic is stored as a single-interval NASA7-shaped
  polynomial (a3 = a4 = 0) and integrated from 298.15 K like the database
  components.
- **Provenance:** `DataMetadata.version` gains `+pseudo-2026-11-24.1`, and the
  source names the method set `petroleum-pseudo-kesler-lee-1976-lee-kesler-1975-riazi-daubert-1987`.
  The inputs themselves are stored with the selection.

## Verification

| Evidence | Kind | Where | Budget |
|---|---|---|---|
| Published worked examples E3.1, E3.2, E3.6, E3.7, E3.9, E3.11 | Implementation verification | `validation/pseudo-components/fahim-2010-worked-examples.json`, `crates/thermo-data/tests/pseudo_worked_examples.rs` | Half a unit of the last printed digit; 1e-3 relative for the E3.11 cp constants (unstated intermediate rounding) |
| Independent Python evaluation of every correlation at the surrogate species' Tb and SG | Implementation verification | `validation/fixtures/pseudo-components-v1.json` (generated by `validation/generate-pseudo-component-references.py`) | 1e-12 relative |
| NASA7-shaped integration of the quadratic cp | Consistency | `thermo-data` unit tests | 1e-9 relative |
| PH and PS round trips at a two-phase TP state; cp = ∂h/∂T at liquid, vapor and low-pressure states | Consistency | `crates/property-packages/tests/pseudo_components.rs` | 1e-6 K; 1e-6 relative |
| Kerosene surrogate against its pseudo-component characterisation | Model comparison | same test and fixture | Declared before the first comparison (below) |
| Heater and flash separator with five TBP cuts and methane, save/reload, stale-on-edit | Flowsheet evidence | `crates/flowsheet-project/tests/pseudo_components.rs` | Component residual < 1e-8, scaled energy residual < 2e-9 |
| Failure paths: undefined, unselected or duplicate pseudo-components, partial matrix, empty source, interval outside the cp domain, non-cubic packages, saturation outside the requested bounds, reactors, out-of-domain inputs, unknown stored fields | Failure behaviour | the three test files | — |

### Surrogate versus pseudo-component case

The surrogate is m-xylene, n-decane and n-dodecane at 0.2/0.4/0.4 mol, a
kerosene-range test surrogate defined for this comparison (it is not a
published surrogate formulation). Two PR representations of the same cut are
compared at 1 atm with kij = 0:

- **Surrogate:** defined compounds with CoolProp 7.2.0 HEOS critical constants, acentric factors and molar masses (entered as specified constants);
- **Pseudo-components:** each species characterised from its own Tb and SG alone.

Tb and SG come from the same CoolProp fluids. Both representations use the
Kesler–Lee ideal-gas cp, because the embedded NASA7 data have no C8+ records.
The cp itself is compared separately with CoolProp's ideal-gas cp0.

| Quantity | Pseudo-components | Surrogate | Difference | Budget |
|---|---|---|---|---|
| Bubble temperature at 1 atm | 445.04 K | 445.85 K | 0.82 K | 5 K |
| Dew temperature at 1 atm | 466.36 K | 466.28 K | 0.08 K | 5 K |
| Liquid mass density at 300 K (untranslated PR) | 683.7 kg/m³ | 684.1 kg/m³ | 0.06 % | 6 % |
| h(560 K vapor) − h(300 K liquid) | 887.9 kJ/kg | 896.0 kJ/kg | 0.9 % | 5 % |

Per species, against the HEOS constants:

| Species | Tc | Pc | ω | M | Kesler–Lee cp0, 300–600 K |
|---|---|---|---|---|---|
| n-decane | 0.27 % | 2.7 % | 0.6 % | 1.5 % | 0.05–0.9 % |
| n-dodecane | 0.51 % | 3.9 % | 2.5 % | 2.0 % | 0.3–2.3 % |
| m-xylene | 0.67 % | 4.5 % | 4.3 % | **5.1 %** | **7.7–9.6 %** |

The budgets (Tc 2 %, Pc 8 %, ω 10 %, M 5 %, cp0 5 %) were declared before the
first comparison. The two m-xylene exceedances in bold were not absorbed by
widening a budget. They are recorded with their observed values in the
fixture (`recorded_exceedances`), and the test requires them to reproduce
within 0.001 and to still exceed the budget. Both are expected: Riazi–Daubert
and Kesler–Lee are fitted to petroleum fractions, and m-xylene is a pure
aromatic with Watson K 10.4. With the HEOS ω and M, the m-xylene Kesler–Lee
cp is within 2.3–3.9 %. The mixture properties stay well inside their budgets.

This is a model comparison between two representations. It is not
experimental validation of either one.

## Limits and gaps

- **Bubble/dew solvers:** since #227 the A13 solvers cover pseudo-components.
  Independent Clapeyron references for the kerosene cuts (alone and with methane)
  are in `validation/fixtures/a13-extended-v1.json`, built from the correlated
  constants of this record's Python fixture; the surrogate test also checks that
  the solvers reproduce its TP-flash bisection within 2e-4 K.
- **No independent EOS fixture:** there is no Clapeyron/teqp fixture for PR
  flashes that include pseudo-components. The cubic formulation is shared and
  already qualified; the new paths are the constants and the cp polynomial,
  which are covered above.
- **Characterisation scope:** no Twu, Cavett or Riazi–Daubert Tc/Pc
  alternatives, no D1160/D2887 conversions, no Daubert (1994) D86 method, and
  no extrapolation of TBP curves (cuts span only the supplied points).
  Lumping, binary-interaction correlations (for example Chueh–Prausnitz) and
  hydrogen-sulfide/water pairs are not provided; such pairs need explicit values.
- **Interface:** the desktop package dialog cannot yet edit pseudo-components.
  Projects that carry them load and solve, but changing the package in the
  dialog rebuilds the selection from database components only. There is no
  `thermo-cli` example yet.
- **Reactors:** pseudo-components are refused (no formula or formation properties).
- **Installed prototype:** this change has not refreshed the installed desktop prototype.
