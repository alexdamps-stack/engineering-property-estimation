# A17 shortcut column references (#185)

`shortcut-v1.json` holds independent reference values for the
Fenske–Underwood–Gilliland(Molokanov)–Kirkbride shortcut distillation column
(A17 #185 method contract: alpha by successive substitution on the geometric
mean of the top dew-point and bottom bubble-point relative volatilities,
Hengstebeck–Geddes Fenske distribution, Underwood root between the key
volatilities, Molokanov fit of Gilliland, Kirkbride feed location, total or
partial condenser, partial reboiler counted in N, reboiler duty by overall enthalpy
balance). The generator, `generate-shortcut-references.py`, cites the five
correlation sources in its docstring.

Cases (100 mol/s feed, explicit zero kij, components in the A13 set, all
temperatures 200–350 K):

| id | variant | feed (T, p) | column p | keys / recoveries | reflux |
|---|---|---|---|---|---|
| `c2c3-binary-PR-BM` | PR-BM | 295 K, 2.0 MPa (subcooled, q≈1.03) | 2.0 MPa | C2 0.98 / C3 0.97 | 1.3 R_min |
| `c2c3-binary-PR` | PR | same | 2.0 MPa | same | 1.3 R_min |
| `c1c2c3-ternary-PR-BM` | PR-BM | 260 K, 2.5 MPa (two-phase) | 2.5 MPa | C2 0.95 / C3 0.95 | R = 1.1 (≈1.54 R_min) |
| `c1c2c3-letdown-PR-BM` | PR-BM | 270 K, 4.0 MPa (liquid) | 2.5 MPa | C2 0.97 / C3 0.96 | 1.4 R_min |
| `c2c3-partial-PR-BM` | PR-BM | as `c2c3-binary-PR-BM` | 2.0 MPa | C2 0.98 / C3 0.97 | 1.3 R_min, partial condenser |
| `c1c2c3-stages-PR-BM` | PR-BM | as `c1c2c3-letdown-PR-BM` | 2.5 MPa | C2 0.97 / C3 0.96 | N = 25 given |
| `c1c2c3-partial-stages-PR` | PR | 260 K, 2.5 MPa (two-phase) | 2.5 MPa | C2 0.95 / C3 0.95 | N = 18 given, partial condenser |

The last three cases were added for the partial condenser and the stage-count
specification; the first four were regenerated unchanged (identical values).
With `"condenser": "partial"`, the distillate is vapour at the dew point of
x_D (T_D = T_condenser), the reflux is that dew point's incipient liquid x_R
(h_R), and the vapour entering the condenser has y_1 = (x_D + R x_R)/(R + 1)
at its dew point (T_top, h_V,top); Q_C = D[(R + 1) h_V,top − h_D − R h_R].
With reflux kind `stages`, X solves the Molokanov fit for
Y = (N − N_min)/(N + 1) by SciPy brentq and R = (X + R_min)/(1 − X).

Each case records its inputs, the feed's own equilibrium state, q, initial and
converged alpha (and the top/bottom alphas), iteration count, N_min, theta,
R_min, R, X, Y, N, Kirkbride ratio, N_R/N_S and feed stage, D, B, component
product flows, x_D, x_B, T_top (dew of the vapour entering the condenser),
T_D (distillate temperature), T_bot (bubble of x_B), h_F, h_D, h_B, h_V,top,
Q_C, Q_R and the net duty. Budgets
(`tolerances`, `[absolute, relative]`) were fixed with their propagation
argument in the generator before any Rust comparison; do not widen them to
resolve a failure.

## Regenerate

In the pinned reference environment (`validation/reference-requirements.txt`;
on Windows run it in WSL2 Ubuntu):

```sh
python validation/columns/generate-shortcut-references.py --output /tmp/shortcut-candidate.json
```

Review the difference before replacing `shortcut-v1.json`. The run takes
about 15 s and is deterministic.

## Independence

Independent of the Rust engine: no Rust code is called. Fugacity coefficients
and residual properties come from teqp 0.23.2 automatic differentiation;
ideal-gas enthalpies from Cantera NASA7 via the reviewed A3 adapter
(`validation/generate-references.py`, same absolute caloric basis as the other
fixtures); the feed state from the A4 multi-start TPD flash adapter
(`generate-flash-references.py`) through the A14 phase-weighted TP state
(`rotating/generate-ps-references.py`). Bubble and dew temperatures are solved
here by inner ln K successive substitution inside an outer SciPy brentq on
temperature (a 200–350 K scan brackets the root); liquid and vapor use the
smallest and largest mechanically stable roots, as the contract specifies.
As a spot check this solver reproduced the Clapeyron-generated A13 bubble/dew
temperatures (`validation/fixtures/a13-v1.json`) to within 1e-8 K for all
bracketable mixture cases.

Not independent: the component constants and NASA7 coefficients
(`data/source-extracts/components.provisional.json`), the PR/PR-BM
formulation, the zero interaction matrix and the shortcut method contract
itself are shared with the engine. The fixture therefore checks implementation
agreement, not the accuracy of shortcut methods against rigorous columns or
experiment.

## Rigorous distillation references (#186)

`distillation-v1.json` covers a PR ideal-tray column with one feed and
a PR-BM HETP column with two feeds and a liquid side draw. Each case
includes the temperature of every stage, both utility duties and the
flow and composition of each liquid product.

SciPy fits temperatures, binary phase compositions and vapor flows with
a bounded least-squares solve. It uses the reviewed teqp and Cantera
reference adapters for fugacity and enthalpy. The Rust engine uses
different variables and a damped Newton solve. The generator records
the input and adapter hashes and has no Rust dependency. It checks
implementation agreement for the specified equilibrium-stage model.

Generate a candidate in the pinned reference environment and compare it
with the committed reference before any replacement.

```sh
python validation/columns/generate-distillation-references.py --output /tmp/distillation-candidate.json
```

## Adiabatic absorber and stripper references (#187)

`absorber-v1.json` covers propane uptake on PR ideal trays, ethane stripping
in a PR-BM HETP bed and a single PR contact. It records stage temperatures
and both outlet flows and compositions. The cases include pressure let-down
on either inlet. The same component data, cubic formulations and reference
adapters apply as for the rigorous distillation references.

SciPy solves the countercurrent component, enthalpy and equilibrium
equations in direct temperature, composition and vapor-flow coordinates.
teqp checks both inlet states and independently flashes each fitted stage
tie line. The generator runs without the Rust engine and records its input
and adapter hashes.

```sh
python validation/columns/generate-absorber-references.py --output /tmp/absorber-candidate.json
```

## Extended-domain column references (#186, #187)

`extended-v1.json` covers columns outside the original dry ethane/propane,
200–350 K domain, which the columns accept since they take their domain from
the property package:

| id | kind | package | conditions |
|---|---|---|---|
| `methanol-water-trays-PR` | distillation, 6 trays | PR, methanol/water kij = −0.075 | 101.325 kPa, 338–368 K |
| `propane-methanol-trays-PR` | distillation, 4 trays | PR, kij = 0 | 2 MPa, reboiler ≈ 430 K |
| `methanol-solvent-absorber-PR` | absorber, 5 trays | PR, kij = 0 | 3 MPa, methane/CO2 gas, methanol solvent at 250 K |

The methanol/water kij is an explicit approximation recorded in each case
(`kij_note`): with kij = 0 PR splits the bottom liquid into two liquids, which
the column rejects (tested in `flowsheet-project/tests/columns_domain.rs`).
The agreement is implementation verification, not experimental accuracy.

The generator solves the full multicomponent component, energy, logarithmic
equilibrium and explicit summation equations with a damped finite-difference
Newton iteration (SciPy Levenberg–Marquardt as a fallback) in temperature,
log mole-fraction and log vapor-flow coordinates, using the same teqp and
Cantera reference adapters as the other column references (the adapter's
bubble/dew scan is widened from 200–350 K to the cubic package interval
200–600 K and, for the initial profile, starts above the condenser bubble
point). The Rust engine uses different variables (normalized log ratios,
eliminated summations) and its own solver. The independent TPD flash adapter
classifies every inlet as single-phase. At every fitted stage both phases
must pass the adapter's finite-start TPD search, which excludes a second
liquid at the reference solution, and the fitted liquid and vapor must be the
distinct minimum-Gibbs roots of their compositions. The budgets equal those
of `distillation-v1.json` and `absorber-v1.json` and were fixed before any
Rust comparison. Generation takes a few minutes.

```sh
python validation/columns/generate-extended-column-references.py --output /tmp/extended-candidate.json
```

## Column specification set references (#186)

`specs-v1.json` covers the rigorous column features added after the first
references: a PR partial-vapor condenser; a PR-BM partial vapor-and-liquid
condenser (20 mol/s liquid distillate) with a vapor side draw; PR
methanol/water (kij = −0.075) at 101.325 kPa with a 5 kPa condenser drop and
0.7 kPa per stage; and the specification pairs reflux ratio + reboiler duty,
distillate mole fraction + boilup ratio (PR-BM) and condenser duty + bottoms
recovery. Each case records stage temperatures and pressures, both duties,
the distillate, reflux, bottoms and boilup flows, and the flow, phase and
composition of every product in port order.

`generate-column-spec-references.py` imports the stage-property, feed,
tie-line-check and Newton/Levenberg–Marquardt pieces of
`generate-extended-column-references.py` (and through it the reviewed teqp,
Cantera and TPD adapters). It solves temperatures, log mole fractions with
explicit summations, log vapor flows and, for specifications that do not
fix them, the log distillate and reflux, with one residual per such
specification. Starting distillate and reflux for those cases are stated by
hand in the case list. Every fitted stage passes the TPD and root checks of
the extended references at its own pressure. No Rust code or output is used.
The tolerances are the budgets of `extended-v1.json`, fixed before any Rust
comparison. Generation takes about two minutes.

```sh
python validation/columns/generate-column-spec-references.py --output /tmp/specs-candidate.json
```

## Petroleum fractionator references (#190)

`petroleum-v1.json` covers the petroleum fractionator with an 8-cut
crude-like feed (TBP 340–550 K, Watson K 11.8, equal-temperature cuts):
two molar side products with a total condenser and reboiler; two standard
liquid volume side products plus a molar one with a partial condenser and a
pressure profile; and a steam-stripped column (5 mol/s steam at 520 K to the
adiabatic bottom stage, partial condenser, water/cut kij = 0.5 as an
explicit approximation). Each case records the cut inputs and correlated
constants, the feed states, stage temperatures and pressures, duties,
internal flows and every product (flow, phase, composition and, without
water, standard liquid volume flow).

`generate-petroleum-references.py` uses the independent Python correlations
of `validation/generate-pseudo-component-references.py` (CoolProp is stubbed
when absent; only its plain-Python functions are used), teqp 0.23.2 PR for
fugacity coefficients and residual enthalpy with the A4 adapter's root
selection, the Kesler–Lee cp and the water NASA7 lower interval for the
ideal gas. It solves the MESH equations in temperatures, log mole fractions
with explicit summations, log flows, duties and the free distillate and
volume-specified draw rates, started by a Wang–Henke bubble-point iteration.
Every fitted stage phase and both feed phases pass the A4 adapter's
multi-start TPD search. The budgets are those of `extended-v1.json`, plus
1e-8 m³/s on standard liquid volume flows, fixed before any Rust comparison.
No Rust code or output is used.

```sh
python validation/columns/generate-petroleum-references.py --output /tmp/petroleum-candidate.json
```

### Product curves and cut points (#190)

`petroleum-curves-v1.json` comes from
`generate-petroleum-curve-references.py`, Python standard library only and
independent of the Rust engine. It writes its own TBP construction (component
boiling intervals between neighbouring Tb midpoints, linear in cumulative
standard liquid volume) and the inverse API D86 conversion with the Fahim et
al. Table 3.2 constants typed independently, and evaluates them on four
synthetic component sets, on every product of `petroleum-v1.json` at its
reference composition, and on example E3.1 of
`validation/pseudo-components/fahim-2010-worked-examples.json` read backwards
(printed TBP to printed D86 within 0.105 °C). The cut-point case replaces the
side-product flows of `crude-two-draws-reboiled-PR` by the TBP 95 % point of
side draw 1 and the D86 50 % point of side draw 2 read from the reference
products, so the expected column is the reference column. Budgets, fixed
before any Rust comparison: 1e-9 K for the curve functions at given
compositions, and for curves of Rust-solved products the composition budget
2e-7 propagated through Σ|∂T/∂xᵢ| by central differences (recorded per point).

```sh
python validation/columns/generate-petroleum-curve-references.py --output /tmp/petroleum-curves-candidate.json
```

## Petroleum fractionator pump-around and side-stripper references (#190)

`petroleum-v2.json` adds three cases on the same 8-cut crude-like feed
(`petroleum-v1.json` is unchanged):

| id | Column | Pump-arounds | Side strippers |
|---|---|---|---|
| `crude-pump-arounds-reboiled-PR` | total / reboiler (R = 2, D/F = 0.3), feed 480 K | stage 7 → 5 at 30 mol/s returned at 400 K; stage 4 → 2 at 20 mol/s removing 300 kW | — |
| `crude-reboiled-side-strippers-volume-PR` | total / reboiler (R = 3, D/F = 0.3), feed 470 K | — | product 1 (stage 4, 15 mol/s): 2 stages, reboiler duty 150 kW, vapor to stage 3; product 2 (stage 7, 3.0 L/s standard liquid volume): 3 stages, boilup ratio 0.5, vapor to stage 6 |
| `crude-steam-side-stripper-pump-around-PR` | partial vapor / 5 mol/s steam to the bottom (water/cut kij = 0.5) | stage 7 → 6 at 25 mol/s removing 200 kW | product 1 (stage 4, 15 mol/s): 3 stages, 2 mol/s steam at 520 K and 3 bar, vapor to stage 3 |

Each case records everything `petroleum-v1.json` records, plus the
pump-around duties and, per side stripper, its stage temperatures, draw from
the main column, vapor return, bottoms and reboiler duty; a stripped product
is that stripper's bottoms.

`generate-petroleum-side-section-references.py` imports the feed,
characterisation, property, flash, stability, Wang–Henke and Newton pieces
of `generate-petroleum-references.py` unchanged and adds its own network
MESH formulation: every stage of the main column and of each stripper has
its temperature, log mole fractions with explicit summations and log liquid
and vapor flows; the stripper draws, the free distillate, the boilup-ratio
reboiler duty and the volume-specified product rates are further unknowns
with one residual per specification (duties carried in units of 100 kW). The
start is the Wang–Henke iteration of the main column with each stripper
replaced by its product drawn at the draw stage and its steam fed to the
return stage; pump-arounds enter by continuation in their circulation
(0, 25, 50, 75, 100 %). Every fitted phase of every cascade and each
pump-around liquid at its return temperature pass the A4 TPD search. The
budgets are those of `petroleum-v1.json`, fixed before any Rust comparison.
No Rust code or output is used. Each case takes a few minutes; the cases can
be generated separately with `--case` and merged.

```sh
python validation/columns/generate-petroleum-side-section-references.py --output /tmp/petroleum-v2-candidate.json
```

## Rate-based packed absorber references (#241)

`rate-based-v1.json` holds two nitrogen/methanol/water scrubbing columns in
rate-based mode (PR, 50 mm metal Pall rings, 2 m, 0.4 m diameter, 6
segments; PR-BM, 25 mm metal Pall rings, 1.5 m, 0.3 m, 5 segments; methanol
and water neglected in the gas transport properties, nitrogen in the
liquid). Each case records every segment's vapor, liquid and interface
temperatures, flows, bulk compositions, fluxes, wetted area, fraction of
flood and pressure drop, the bed pressure drop and both outlets. The model
contract is [docs/a17-columns-rate-based.md](../../docs/a17-columns-rate-based.md).

The generator solves the segment equations in a different formulation from
the engine (bulk mole fractions and total flows closed by summations;
partial molar enthalpies from the temperature derivative of teqp ln φ plus
Cantera NASA7 ideal-gas enthalpies) by Newton's method with continuation in
the interfacial-area factor, and checks both phases of each interface with
the A4 multi-start TPD stability search. Transport properties come from the #225 reference evaluator
(`validation/transport/generate-transport-references.py`) and the liquid
density from the #228 evaluator (`validation/generate-liquid-density-references.py`,
chemicals 1.5.2); Fuller/Blanc/Wilke–Chang/Tyn–Calus, Onda, Robbins and
Kister–Gill are transcribed separately in the generator. It needs CoolProp
7.2.0 and chemicals 1.5.2 in addition to `validation/reference-requirements.txt`.
Budgets are declared in the generator. The fixture checks implementation
agreement, not the accuracy of the correlations.

```sh
python -I validation/columns/generate-rate-based-references.py --output /tmp/rate-based-candidate.json
```
