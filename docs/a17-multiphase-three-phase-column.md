# A17 multiphase: three-phase distillation (#188)

Aspen Plus *RadFrac* with a decanter and three-phase stages / HYSYS *Three
Phase Distillation* equivalent. Code:
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
- **Every tray and the reboiler carry one or two liquids** (full VLLE stage
  model). A single-liquid stage is in vapor–liquid equilibrium; a three-phase
  stage holds a vapor and two liquids in vapor–liquid–liquid equilibrium, and
  its total liquid flows to the stage below. Which stages carry two liquids is
  found automatically (see *Phase-set rounds*). A three-phase reboiler
  delivers both liquids as the bottoms stream (two liquid phases, light
  first). A column in which no stage splits gives the decanter-only result.
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
Design of Distillation Systems*, 2001, ch. 8; three-phase stages: Cairns and
Furzer, Ind. Eng. Chem. Res. 29 (1990) 1349, here in a Naphtali–Sandholm
simultaneous-correction form), in − out:

- decanter: component split of the condensate into a and b, iso-fugacity of
  the two liquids with the package liquid fugacity coefficients
  (`CalcContext::homogeneous`, smallest-volume root: NRTL activity
  coefficients for NRTL-IG, the cubic liquid root for PR/PR-BM), energy
  balance giving the condenser duty;
- stages: component balances with the reflux φ_L·a + φ_H·b entering stage 2
  and the total liquid of the stage above entering every other stage,
  ln y − ln x − ln K = 0 with K = φ^L/φ^V of the package's liquid and vapor at
  the stage state, energy balances on stages 2..N−1; the reboiler replaces its
  energy balance by L_N + L′_N = B and its duty follows from its energy
  balance;
- three-phase stages add a second liquid l′ with ln x + ln φ^L(x) − ln x′ −
  ln φ^L(x′) = 0 (both with the smallest-volume root); it enters the
  component and energy balances and the liquid flowing down.

Unknowns per block: ln component flows of both liquids (decanter) or of
liquid, vapor and second liquid (stages; the second-liquid variables of a
single-liquid stage are inert and decoupled, so single-liquid columns give the
decanter-only answer), plus the scaled condenser duty or the stage
temperature. Start: a single-liquid bubble point of the feed,
constant-molar-overflow passes (tridiagonal component balances, one
bubble-point Newton step per stage) with a lagged decanter VLLE flash of the
top vapor. Then a damped Newton iteration on the block-tridiagonal stage
Jacobian (three-color forward differences, block Thomas solver shared with the
extractor) to a scaled residual of 1e-11. No package-interface extension; no
`cubic()` call.

### Phase-set rounds

1. Solve with one liquid on every stage. If that fails, the stages whose
   initial-profile liquid splits receive a second liquid and the solve is
   repeated once from the initial profile.
2. After every converged pass each single-liquid stage is flashed at its
   temperature with `vlle_tp` (tie-line midpoint; if that flash fails, the
   stage liquid alone, then a Michelsen liquid tangent-plane test with the
   package liquid fugacity coefficients). A second liquid is added on the tie
   line of that flash (lever rule of the stage liquid, 5..95 %).
3. Each three-phase stage is checked with `vlle_ph` of its contents (same
   temperature within 1e-5 K, compositions within 1e-5, fractions within
   1e-4). If the package finds fewer phases, the second liquid is removed;
   if its PH flash fails (seen for strongly immiscible PR binaries at their
   three-phase temperature), both solved liquids must pass the liquid
   tangent-plane test instead.
4. During Newton a second liquid that falls below 1e-9 of the stage liquid or
   comes within 1e-4 in mole fraction of the first is merged into it and the
   round restarts.
5. Rounds repeat until the phase set no longer changes (at most 12). A phase
   set that recurs is a phase-count oscillation failure.

In a binary at fixed pressure a three-phase stage is invariant (phase rule):
all such stages share the heteroazeotropic temperature and phase
compositions; only the phase amounts differ.

## Explicit failures

| Situation | Message contains |
|---|---|
| Package without VLLE/LLE (GERG-2008, EOS-CG-2021) | `needs the VLLE/LLE flash` |
| Condensate stays one liquid at the decanter temperature | `no liquid-liquid split in the decanter … use the distillation column` |
| Vapor in the decanter (temperature above the condensate bubble point) | `vapor forms in the decanter` |
| Decanter liquids within 1e-3 in mole fraction | `decanter liquids coincide` |
| Solved decanter split differs from the VLLE flash by more than 1e-6 | `not the stable liquid-liquid split` |
| Decanter liquids swap their density order | `exchanged their density order` |
| A single-liquid stage does not reproduce its VLE split in the check flash | `did not reproduce a stable vapor-liquid split` |
| A three-phase stage differs from the VLLE PH flash of its contents | `did not reproduce a stable vapor-liquid-liquid split` |
| The stage phase set recurs | `stage phase counts oscillate` |
| More than 12 phase-set rounds | `did not settle` |
| Condenser does not remove or reboiler does not supply heat | `heat input at the reboiler … infeasible` |
| Newton stall or 80 iterations (also an infeasible bottoms/reflux specification); the message lists the three-phase stages of that round | `stalled` / `did not converge` / `two liquids on stages` |
| Column pressure above the feed pressure | `exceeds the feed pressure` |

Decanter or stage liquids within 0.02 give the warning `TPC-NEAR-PLAIT-POINT`.

## Validation

### Decanter-only columns

`validation/multiphase/generate-three-phase-column-references.py` (Python
standard library, no Rust) solves the same equilibrium-stage model in
different variables (explicit temperatures, mole fractions and total flows
with summation equations) with Newton's method. Gamma–phi: NRTL (Renon and
Prausnitz 1968) from the ChemSep/DECHEMA water/1-butanol pair, water
Wagner–Pruss vapor-pressure ancillary, ChemSep 1-butanol vapor pressure,
NASA7 ideal-gas enthalpies, Clausius–Clapeyron pure-liquid enthalpies and
the NRTL excess enthalpy, temperature derivatives by complex step. Each stage
liquid passes a 1/4000 composition-grid tangent-plane scan, the decanter a
tie-line scan and a bubble-pressure check. Fixture (unchanged by the VLLE stage
model): `validation/multiphase/three-phase-column-v1.json`.

| Case | Trays | Feed (330 K, 2 bar, 1 mol/s) | Reflux | B/F | T stage 2 / N, K | Q_C / Q_R, kW |
|---|---|---|---|---|---|---|
| `butanol-side-light-total-6-trays` | 6, feed stage 2 | water 0.3 | light, whole | 0.71 | 367.50 / 388.37 | 41.0 / 49.0 |
| `butanol-side-light-ratio-8-trays` | 8, feed stage 4 | water 0.3 | light, R = 4 | 0.65 | 368.32 / 391.30 | 39.5 / 47.1 |
| `butanol-side-fractions-6-trays` | 6, feed stage 2 | water 0.3 | 0.9 light, 0.2 heavy | 0.68 | 367.29 / 389.63 | 44.1 / 51.8 |
| `water-side-heavy-total-6-trays` | 6, feed stage 2 | water 0.996 | heavy, whole | 0.992 | 366.42 / 372.57 | 0.70 / 3.56 |

All at 101.325 kPa with the decanter at 320 K.

### Columns with three-phase stages

`validation/multiphase/generate-three-phase-column-vlle-references.py` (Python
standard library, no Rust) solves the VLLE stage model in explicit stage
temperatures, ln mole fractions with summation equations and signed phase
amounts (negative-flash form). Its phase set comes from its own active-set
loop: one liquid everywhere first; a second liquid where a composition-grid
tangent-plane scan of the stage liquid is negative (or, when one liquid
everywhere does not converge, on every stage, started at the binary
three-phase point); removed where a converged liquid amount is not positive.
NRTL-IG exactly as above. PR (Peng and Robinson 1976) for water/propane with
the exact Ω_a, Ω_b, the 1976 m(ω), van der Waals one-fluid mixing, no volume
translation and closed-form departure functions, NASA7 ideal-gas enthalpies.
Water–propane **kij = 0.5** is the explicit water–hydrocarbon approximation of
the #224 VLLE references ([VLLE status](vlle-status.md)), not a regressed
value; it gives water in liquid propane 4.2e-4 at 300 K and propane in water
3e-12. Accepted solutions pass the grid scan on every stage liquid (both
liquids on three-phase stages), the decanter scans, and component and energy
balance closure. Fixture: `validation/multiphase/three-phase-column-vlle-v1.json`.

| Case | Package | Feed | Reflux | B/F | Phases on stages 2..N | T stage 2 / N, K | Q_C / Q_R, kW |
|---|---|---|---|---|---|---|---|
| `nrtl-butanol-fractions-0.95-0.9-split-trays` | NRTL-IG | 330 K, 2 bar, water 0.3 | 0.95 light, 0.9 heavy | 0.71 | 3 3 3 3 3 2 2 | 366.21 / 380.87 | 190.8 / 197.4 |
| `nrtl-butanol-fractions-1-0.8-split-upper-trays` | NRTL-IG | 330 K, 2 bar, water 0.3 | 1.0 light, 0.8 heavy | 0.71 | 3 3 3 2 2 2 2 | 366.21 / 388.37 | 156.1 / 164.0 |
| `pr-propane-drying-free-water-upper-trays` | PR | 400 K, 25 bar, water 0.02 | light, whole | 0.9802 | 3 3 3 3 3 2 2 | 329.80 / 330.17 | 47.5 / 28.6 |
| `pr-propane-drying-dry-reboiler` | PR | 400 K, 25 bar, water 0.05 | light, whole | 0.951 | 3 3 3 3 3 3 2 | 329.80 / 329.97 | 117.5 / 97.7 |

NRTL cases: 6 trays, feed stage 2, 101.325 kPa, decanter 320 K (the
heteroazeotropic temperature of the model is 366.21 K). PR cases (wet propane
drying with free-water decant): 6 trays, feed stage 3, 20 bar, decanter 300 K
(water/propane three-phase temperature 329.80 K). Every case has
single-liquid and three-phase stages. Declared tolerances (before comparison,
both fixtures): temperature 1e-6 K; composition 1e-8; flow 1e-9 mol/s + 1e-7
relative; duty 1e-3 W + 1e-7 relative. Phase counts must agree exactly.

Tests: `crates/flowsheet-core/tests/multiphase_three_phase_column.rs`
(both reference sets including energy-port values, products, phase counts and
both liquids; decanter against `vlle_tp` of the condensate; heavy-liquid reflux
ratio; packed bed equal to the equivalent trays; a decanter-only column reports
one liquid below the decanter; a three-phase stage satisfies both equilibria and
matches the package `vlle_ph` flash; a three-phase reboiler delivers both
liquids as bottoms; failure paths: decanter vapor, negative reboiler duty,
infeasible bottoms specification, nonconvergence naming the three-phase
stages, feed pressure, unset value, miscible methanol/water condensate,
GERG-2008 package; phase-count oscillation as a unit test) and
`crates/flowsheet-project/tests/multiphase_three_phase_column.rs` (project
solve, save/reload with results including three-phase stages, older results
without the new fields, drawn-only record with its connections).

Regenerate (WSL, about 3 s and 20 s):
`python3 validation/multiphase/generate-three-phase-column-references.py --output /tmp/three-phase-column-candidate.json`
and
`python3 validation/multiphase/generate-three-phase-column-vlle-references.py --output /tmp/three-phase-column-vlle-candidate.json`.
Review differences before replacing a fixture.

## Project format

The record `{"kind": "three_phase_column", …}` gained `internals`,
`feed_stage`, `pressure_pa`, `decanter_temperature_k`, `reflux` (`mode`:
`light_total`, `heavy_total`, `light_ratio`, `heavy_ratio`, `fractions`) and
`bottoms_fraction`, all with defaults, behind format version 4. A record
written while the kind was drawn only (`{"kind": "three_phase_column"}`) opens
with five ideal trays, feed at stage 2, light-liquid total reflux and unset
pressure, decanter temperature and bottoms fraction. Results are saved as
`results.units[*].three_phase_column`; `stage_phase_counts` (2 or 3 per stage,
decanter 2) and `stage_liquid_splits` (both liquids per stage, `null` on
single-liquid stages, stage 1 the decanter) were added behind format version 4
with defaults, so earlier results open with empty lists. On a three-phase stage
`liquid_flows_mol_s` and `liquid_compositions` hold the total liquid. The
desktop results list the phases per stage and both liquids (flows and mole
fractions) of every three-phase stage.

## Limits

- One column pressure, no side draws (a free-water side draw from a
  three-phase tray is not offered).
- The water/1-butanol NRTL pair is a 1 atm VLE regression; its ambient mutual
  solubilities are off ([VLLE status](vlle-status.md)). PR water–propane uses
  the kij = 0.5 approximation. The references verify the implementation, not
  experimental column data.
- Robustness rests on the single-liquid start and the phase-set rounds; some
  specifications with free water on most trays do not converge (for example
  the wet-propane column with a light-liquid reflux ratio of 3 at B/F 0.9) and
  end as explicit Newton failures that name the three-phase stages.
- At most two liquids per stage (the package VLLE operations support three
  phases in total).
- The initial profile assumes a liquid feed (constant molar overflow);
  vapor feeds converge in the tested cases but rely on the Newton stage.
- Specifying B/F beyond what the decanter liquids allow (for example more
  aqueous decant than the feed water) has no solution and ends as a Newton
  failure rather than a dedicated feasibility message.
- The ISO 10628 symbol is the existing traced valve-tray column
  (`s1-x8012-column-with-valve-trays`, attribution in the design reference);
  the decanter is not drawn separately.
