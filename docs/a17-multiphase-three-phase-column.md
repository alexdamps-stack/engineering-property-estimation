# A17 multiphase: three-phase distillation (#188)

Aspen Plus *RadFrac* with a decanter (three-phase) / HYSYS *Three Phase
Distillation* equivalent, first scope. Code:
`crates/flowsheet-core/src/units/multiphase/three_phase_column.rs`; records
`crates/flowsheet-project/src/records/multiphase.rs`; results
`ThreePhaseColumnProfile` in `crates/flowsheet-project/src/results.rs`; form
`crates/thermo-desktop/src/panels/spec/multiphase.rs`.

## Supported scope

- One column pressure. Stage 1 is a **total condenser with a decanter** at a
  specified temperature; stages 2..N−1 are ideal trays (or the ideal stages of
  one packed bed, HETP rules as in [packing](a17-columns-packing.md)); stage N
  is a **partial reboiler**.
- The vapor from stage 2 condenses completely and splits in the decanter into a
  light and a heavy liquid (liquid–liquid equilibrium at the decanter
  temperature and column pressure). The light liquid is the less dense one, as
  in the VLLE roles ([VLLE status](vlle-status.md)).
- **Two liquids are allowed only in the decanter.** Trays and reboiler carry
  one liquid in vapor–liquid equilibrium. After convergence every stage
  tie-line midpoint is flashed with the package's three-phase operation; a
  second liquid on any stage is an explicit failure. A full VLLE stage model
  (two liquids on any stage) is not implemented: phase appearance on stages
  needs a switching or Michelsen-type stage formulation that this first
  version does not attempt.
- No side feeds or draws, no Murphree efficiency, no pressure profile, no
  stage heaters, no subcooled reflux return to a stage other than stage 2.

## Ports and specification

| Port | Role |
|---|---|
| `feed` (inlet 0) | Enters the specified feed stage as a whole |
| `distillate` (outlet 1) | Non-refluxed part of the reflux liquid (zero for total reflux of that liquid); with *Fractions*, the non-refluxed light liquid |
| `decant` (outlet 2) | The other decanter liquid; with *Fractions*, the non-refluxed heavy liquid |
| `bottoms` (outlet 3) | Reboiler liquid |
| `condenser duty` (energy outlet) | Heat removed, W (positive) |
| `reboiler duty` (energy inlet) | Heat supplied, W (positive) |

Distillate and decant leave at the decanter temperature; bottoms at the
reboiler temperature. The unit duty is `Q_reboiler − Q_condenser`.

| Field | Choice / unit | Validation |
|---|---|---|
| Internals | 1–24 ideal trays **or** one packed bed (HETP plus height or theoretical stages) | count out of range, missing HETP/height |
| Feed stage | 2 (top tray) to N (reboiler) | out of range |
| Column pressure | Pa, absolute, at most the feed pressure | required, > 0 |
| Decanter temperature | K, below the condensate bubble point, inside the package interval | required, > 0 |
| Reflux | *Light liquid* or *Heavy liquid*, each with *whole liquid* (no distillate) or a *reflux ratio* L/D > 0 of that liquid; **or** *Fractions* (0–1) of both liquids, not both 0 and not both 1 | missing / not a fraction / out of range |
| Feed to bottoms | B/F strictly between 0 and 1 | required |

Reflux ratio R maps to the refluxed fraction R/(1+R) of that liquid. The
specification set has two degrees of freedom (reflux, bottoms) as for the
rigorous distillation column; the decanter temperature replaces the
condenser specification.

## Model and solution

Equations (Seader, Henley and Roper, *Separation Process Principles*, 3rd ed.,
2011, §10.1; heterogeneous azeotropic columns: Doherty and Malone, *Conceptual
Design of Distillation Systems*, 2001, ch. 8), in − out:

- decanter: component split of the condensate into a and b, iso-fugacity of
  the two liquids with the package liquid fugacity coefficients
  (`CalcContext::homogeneous`, smallest-volume root: NRTL activity
  coefficients for NRTL-IG, the cubic liquid root for PR/PR-BM), energy
  balance giving the condenser duty;
- stages: component balances with the reflux φ_L·a + φ_H·b entering stage 2,
  ln y − ln x − ln K = 0 with K = φ^L/φ^V of the package's liquid and vapor at
  the stage state, energy balances on stages 2..N−1; the reboiler replaces its
  energy balance by L_N = B and its duty follows from its energy balance.

Unknowns per block: ln component flows of both liquids (decanter) or of liquid
and vapor (stages), plus the scaled condenser duty or the stage temperature.
Start: a single-liquid bubble point of the feed, constant-molar-overflow
passes (tridiagonal component balances, one bubble-point Newton step per
stage) with a lagged decanter VLLE flash of the top vapor. Then a damped
Newton iteration on the block-tridiagonal stage Jacobian (three-color forward
differences, block Thomas solver shared with the extractor) to a scaled
residual of 1e-11. No package-interface extension; no `cubic()` call.

## Explicit failures

| Situation | Message contains |
|---|---|
| Package without VLLE/LLE (GERG-2008, EOS-CG-2021) | `needs the VLLE/LLE flash` |
| Condensate stays one liquid at the decanter temperature | `no liquid-liquid split in the decanter … use the distillation column` |
| Vapor in the decanter (temperature above the condensate bubble point) | `vapor forms in the decanter` |
| Decanter liquids within 1e-3 in mole fraction | `decanter liquids coincide` |
| Solved decanter split differs from the VLLE flash by more than 1e-6 | `not the stable liquid-liquid split` |
| Decanter liquids swap their density order | `exchanged their density order` |
| A tray or the reboiler liquid splits into two liquids | `stage N liquid splits into two liquids` |
| A stage does not reproduce its VLE split in the check flash | `did not reproduce a stable vapor-liquid split` |
| Condenser does not remove or reboiler does not supply heat | `heat input at the reboiler … infeasible` |
| Newton stall or 80 iterations (also an infeasible bottoms/reflux specification) | `stalled` / `did not converge` |
| Column pressure above the feed pressure | `exceeds the feed pressure` |

Decanter liquids within 0.02 give the warning `TPC-NEAR-PLAIT-POINT`.

## Validation

`validation/multiphase/generate-three-phase-column-references.py` (Python
standard library, no Rust) solves the same equilibrium-stage model in
different variables (explicit temperatures, mole fractions and total flows
with summation equations) with Newton's method. Gamma–phi: NRTL (Renon and
Prausnitz 1968) from the ChemSep/DECHEMA water/1-butanol pair, water
Wagner–Pruss vapor-pressure ancillary, ChemSep 1-butanol vapor pressure,
NASA7 ideal-gas enthalpies, Clausius–Clapeyron pure-liquid enthalpies and
the NRTL excess enthalpy, temperature derivatives by complex step. Each stage
liquid passes a 1/4000 composition-grid tangent-plane scan, the decanter a
tie-line scan and a bubble-pressure check. Fixture:
`validation/multiphase/three-phase-column-v1.json`.

| Case | Trays | Feed (330 K, 2 bar, 1 mol/s) | Reflux | B/F | T stage 2 / N, K | Q_C / Q_R, kW |
|---|---|---|---|---|---|---|
| `butanol-side-light-total-6-trays` | 6, feed stage 2 | water 0.3 | light, whole | 0.71 | 367.50 / 388.37 | 41.0 / 49.0 |
| `butanol-side-light-ratio-8-trays` | 8, feed stage 4 | water 0.3 | light, R = 4 | 0.65 | 368.32 / 391.30 | 39.5 / 47.1 |
| `butanol-side-fractions-6-trays` | 6, feed stage 2 | water 0.3 | 0.9 light, 0.2 heavy | 0.68 | 367.29 / 389.63 | 44.1 / 51.8 |
| `water-side-heavy-total-6-trays` | 6, feed stage 2 | water 0.996 | heavy, whole | 0.992 | 366.42 / 372.57 | 0.70 / 3.56 |

All at 101.325 kPa with the decanter at 320 K. Declared tolerances (before
comparison): temperature 1e-6 K; composition 1e-8; flow 1e-9 mol/s + 1e-7
relative; duty 1e-3 W + 1e-7 relative.

Tests: `crates/flowsheet-core/tests/multiphase_three_phase_column.rs`
(references including energy-port values and products; decanter against
`vlle_tp` of the condensate; heavy-liquid reflux ratio; packed bed equal to the
equivalent trays; failure paths: decanter vapor, stage two-liquid split,
negative reboiler duty, infeasible bottoms specification, feed pressure,
unset value, miscible methanol/water condensate, GERG-2008 package) and
`crates/flowsheet-project/tests/multiphase_three_phase_column.rs` (project
solve, save/reload with results, drawn-only record with its connections).

Regenerate (WSL, about 3 s):
`python3 validation/multiphase/generate-three-phase-column-references.py --output /tmp/three-phase-column-candidate.json`.
Review differences before replacing the fixture.

## Project format

The record `{"kind": "three_phase_column", …}` gained `internals`,
`feed_stage`, `pressure_pa`, `decanter_temperature_k`, `reflux` (`mode`:
`light_total`, `heavy_total`, `light_ratio`, `heavy_ratio`, `fractions`) and
`bottoms_fraction`, all with defaults, behind format version 4. A record
written while the kind was drawn only (`{"kind": "three_phase_column"}`) opens
with five ideal trays, feed at stage 2, light-liquid total reflux and unset
pressure, decanter temperature and bottoms fraction. Results are saved as
`results.units[*].three_phase_column`.

## Limits

- Decanter-only liquid–liquid split (see scope). Columns whose trays enter the
  two-liquid region (for example refluxing most of both liquids) fail.
- The water/1-butanol NRTL pair is a 1 atm VLE regression; its ambient mutual
  solubilities are off ([VLLE status](vlle-status.md)). The references verify
  the implementation, not experimental column data. No PR water–hydrocarbon
  column reference yet (PR works through the same package operations).
- The initial profile assumes a liquid feed (constant molar overflow);
  vapor feeds converge in the tested cases but rely on the Newton stage.
- Specifying B/F beyond what the decanter liquids allow (for example more
  aqueous decant than the feed water) has no solution and ends as a Newton
  failure rather than a dedicated feasibility message.
- The ISO 10628 symbol is the existing traced valve-tray column
  (`s1-x8012-column-with-valve-trays`, attribution in the design reference);
  the decanter is not drawn separately.
