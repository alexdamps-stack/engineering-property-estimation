# A17 multiphase: liquid–liquid extractor (#189)

Aspen Plus *Extract* / HYSYS *Liquid-Liquid Extractor* equivalent. Code:
`crates/flowsheet-core/src/units/multiphase/extractor.rs`; records
`crates/flowsheet-project/src/records/multiphase.rs`; form
`crates/thermo-desktop/src/panels/spec/multiphase.rs`.

## Ports and specification

| Port | Role |
|---|---|
| `feed` (inlet 0) | Raffinate-phase feed, enters stage 1 |
| `solvent` (inlet 1) | Extract-phase solvent, enters stage N |
| `extract` (outlet 2) | Extract phase leaving stage 1 |
| `raffinate` (outlet 3) | Raffinate phase leaving stage N |

No energy ports. An isothermal column reports its duty in the unit balance
and result; an adiabatic column has zero duty.

| Field | Choice / unit | Validation |
|---|---|---|
| Internals | Ideal stages (1–50) **or** one packed section: HETP plus packed height (rounding nearest/down/up) or theoretical stages ([packing rules](a17-columns-packing.md)) | count out of range, missing or non-positive HETP/height |
| Column pressure | Pa, absolute | required, > 0, at most both inlet pressures |
| Thermal mode | Isothermal (temperature, K) **or** adiabatic | temperature required and inside the package interval |
| Feed phase | Detect / heavy / light | a specified phase must agree with the stage densities |

Stages are numbered from the feed end. A heavy feed enters at the top, a light
feed at the bottom. Packed stages are ideal (efficiency 1, HETP includes
non-ideal contacting); no Murphree efficiency.

## Model

Per stage j (R_0 = feed, E_{N+1} = solvent): component balances, iso-fugacity
of the two liquids `ln x^E + ln φ^E = ln x^R + ln φ^R` with the package's liquid
fugacity coefficients (`CalcContext::homogeneous`, smallest-volume root: NRTL
activity coefficients with pure-liquid standard states for NRTL-IG, the cubic
liquid root for PR/PR-BM), and, if adiabatic, the stage energy balance with the
package liquid enthalpies. Inlets keep their enthalpy through any pressure
let-down. No package-interface extension was needed; no `cubic()` call.

Solution: two-liquid VLLE flash of the combined inlets (at the specified
temperature, or at the adiabatic mixing temperature from a VLLE PH flash);
the liquid carrying more solvent-origin material starts the extract profile.
Isothermal sum-rates passes (fixed K, tridiagonal component balances) refine
the profile; a damped Newton iteration on ln component flows (and stage
temperatures) with a block-tridiagonal Jacobian closes the scaled residuals to
1e-12.

Verification before results are published, each an explicit failure:

| Situation | Message contains |
|---|---|
| Inlet missing, zero or non-finite flow | `feed` / `solvent … positive finite flow` |
| Column pressure above an inlet pressure | `pressure exceeds` |
| Package without VLLE/LLE (GERG-2008, EOS-CG-2021) | `needs the VLLE/LLE flash … not supported` |
| Mixed inlets or a stage give one liquid (miscible region, beyond the plait point) | `no liquid-liquid split` |
| Vapor in the mixed inlets or at a stage | `vapor forms` |
| Stage liquids closer than 1e-3 in mole fraction (plait point, trivial solution) | `coincide` |
| Solved tie line differs from the stage VLLE flash by more than 1e-6 | `not the stable liquid-liquid split` |
| Density order of the phases inverts between stages | `invert their density order` |
| Specified feed phase contradicts the densities | `specified heavy` / `specified light` |
| Newton stall or 60 iterations | `stalled` / `did not converge` |

Stage liquids within 0.02 give the warning `LLX-NEAR-PLAIT-POINT`.

## Result

`LiquidExtractorResult` (saved as `results.units[*].liquid_extractor`, added
behind format version 4): stage temperatures, extract and raffinate flows and
compositions per stage, whether the feed is the heavy phase, isothermal duty,
packed height and Newton iterations. A version-4 record written while the kind
was drawn only (`{"kind": "liquid_extractor"}`) opens with the default
specification (5 ideal stages, isothermal, unset pressure and temperature).

## Validation

`validation/multiphase/generate-extractor-references.py` (Python standard
library, no Rust) solves the isothermal countercurrent cascade in composition
and flow coordinates with the NRTL activity coefficients (Renon and
Prausnitz 1968) read directly from the source extracts (ChemSep/DECHEMA
water/1-butanol, methanol/1-butanol and methanol/water pairs, as distributed
with DWSIM). Each stage tie line passes a 1/400 composition-grid tangent-plane
scan. Fixture: `validation/multiphase/extractor-v1.json`.

| Case | Stages | T, p | Feed | Solvent |
|---|---|---|---|---|
| `methanol-from-water-3-stage-298K` | 3 | 298.15 K, 1 atm | 1 mol/s water 0.92 / methanol 0.08 | 0.15 mol/s 1-butanol |
| `wet-solvent-8-stage-310K` | 8 | 310 K, 2 bar | 1 mol/s water 0.95 / methanol 0.05 | 0.2 mol/s water 0.3 / 1-butanol 0.7 |
| `single-stage-298K` | 1 | 298.15 K, 1 atm | 1 mol/s water 0.9 / methanol 0.1 | 0.15 mol/s 1-butanol |

Declared tolerances (before comparison): composition 1e-8 absolute; phase flow
1e-10 mol/s + 1e-8 relative. Tests: `crates/flowsheet-core/tests/multiphase_extractor.rs`
(references, single stage against `vlle_tp`, adiabatic single stage against
`vlle_ph`, adiabatic 4-stage energy balance, packed section, PR liquid-root
cascade with an explicit water–propane kij 0.5 approximation, failure paths,
unsupported package) and `crates/flowsheet-project/tests/multiphase_extractor.rs`
(project solve, save/reload with results, legacy record).

Regenerate (WSL): `python3 validation/multiphase/generate-extractor-references.py --output /tmp/extractor-candidate.json`
(about 15 s). Review differences before replacing the fixture.

## Limits

- The NRTL pairs are 1 atm VLE regressions: the model's water/1-butanol mutual
  solubilities at 298 K are off ([VLLE status](vlle-status.md)), and the
  methanol-containing two-liquid region closes near 12–15 % methanol. The
  references verify the implementation, not experimental extraction data.
- NRTL-IG is limited to the components with complete NRTL pairs; the desktop
  component list for NRTL-IG still offers methanol and water only, so a
  water/methanol/1-butanol project must be created programmatically or by file.
- One pressure, no pressure profile; ideal stages only (no Murphree or overall
  efficiency, no rate-based packing, no hydraulics or flooding).
- No side feeds or draws, no stage heaters; adiabatic or isothermal only.
- The adiabatic reference is the VLLE PH flash for one stage and balance
  closure for several; no independent multi-stage adiabatic reference.
- Up to two liquids per stage; a third liquid or vapor is an explicit failure.
