# A17 columns: petroleum fractionator (#190)

The petroleum fractionator (`UnitKind::PetroleumFractionator`, Aspen PetroFrac
/ HYSYS crude column equivalent) is calculated since #190: one main column for
feeds characterised as petroleum pseudo-components
([pseudo-component status](pseudo-components-status.md), #229), with
pump-arounds and side strippers solved simultaneously with it. Code:
`crates/flowsheet-core/src/units/columns/petroleum.rs` and
`side_sections.rs`.

## Model

The fractionator is composed on the rigorous MESH solve of the
[distillation column](a17-columns-stages.md#rigorous-distillation) (#186). Its
specification maps onto a `DistillationColumnSpec` with the same stage
numbering, feeds, liquid side draws, pressure profile, condenser and top
specification. Results are reported as the distillation column result (stage
temperatures and pressures, duties, distillate, reflux, bottoms, boilup), so
the project results and the desktop panel show them unchanged.

- **Stages.** Stage 1 is the condenser, stages 2..N+1 the ideal trays (or one
  HETP bed) and the last stage the reboiler or, with stripping steam, the
  bottom tray. 1–24 internal stages.
- **Feeds.** The fixed `feed` port is feed 1; `feed 2..n` are numbered ports.
  Each feed enters its stage as a whole with its own enthalpy and must be at or
  above the stage pressure. The combined feed (without steam) must carry at
  least one characterised pseudo-component.
- **Condenser.** Total (liquid `overhead`), partial vapor (vapor `overhead`),
  or partial vapor and liquid (vapor `overhead` plus `liquid distillate` at a
  molar rate). A three-phase condenser with a free-water decant is not
  modelled: a condenser liquid that splits fails as
  "stage 1 liquid splits into two liquids".
- **Top specification.** The distillation set: reflux ratio, reflux rate,
  distillate rate, condenser duty, distillate mole fraction or recovery.
- **Bottom.** *Reboiler*: a partial reboiler with one of the distillation
  bottom specifications (distillate fraction, bottoms rate, reboiler duty,
  boilup ratio, bottoms mole fraction or recovery) and the `reboiler duty`
  energy inlet. *Stripping steam*: a `stripping steam` inlet enters the bottom
  stage, which is adiabatic and has no reboiler. The zero bottom duty replaces
  the second specification, so the distillate becomes a Newton unknown
  (`BottomEnd::Adiabatic` in `distillation.rs`). The steam can be any stream;
  water needs explicit water/pseudo-component interaction parameters.
- **Side products.** `side draw k` leaves an internal stage at its liquid
  state. The rate is either a molar flow or a **standard liquid volume flow**
  at 60 °F. The standard liquid volume of a pseudo-component is
  M / (SG ρw), with SG its 60 °F/60 °F specific gravity and
  ρw = 999.0170824 kg/m³ (water at 60 °F and 1 atm, IAPWS-95 via CoolProp
  7.2.0); mixtures mix ideally. The same definition gives the assay cut
  mole fractions (nᵢ ∝ VᵢSGᵢ/Mᵢ). A volume-specified product is converted to
  its molar flow by successive substitution on its own composition, each pass
  a converged column, until the molar flow changes by less than 1e-11
  relative (at most 30 passes). Every component that can reach a draw must
  then be a pseudo-component; database components (water, light ends) have
  no specific gravity and fail explicitly.

## Pump-arounds and side strippers

Both are solved **simultaneously** with the main column in one Newton system
(`crates/flowsheet-core/src/units/columns/side_sections.rs`, hooked into
`distillation::calculate_sections`; the plain distillation column passes no
side sections and is unchanged). The unknowns of each side stripper follow
the main-column profile, and its MESH residuals follow the main-column
residuals; the coupling streams are feeds and draws whose flows,
compositions and enthalpies are functions of the unknowns, so the Jacobian
is block-structured (main block, one diagonal block per stripper, coupling
blocks through the draw and return stages). The dense finite-difference
Newton of the column solves it as a whole.

- **Pump-around** (`PumpAround`): a molar flow P of liquid leaves internal
  stage i at the stage liquid state and returns to stage j, 2 ≤ j ≤ i, as a
  feed of the same composition at the return-stage pressure. The cooling is
  either a **return temperature** (enthalpy of the liquid root at that
  temperature) or a **duty** Q removed (h_return = h_L,i − Q/P). The duty
  leaves on the energy outlet `pump-around duty k`; with a return
  temperature it is P (h_L,i − h_return) of the converged column. After
  convergence a return temperature must lie below the draw-stage temperature
  and the returned liquid must be a single liquid at the return pressure (TP
  flash); a duty must be positive.
- **Side stripper** (`SideStripper`, one per side product): 1–24 equilibrium
  stages at the pressure of the product's draw stage. The draw S from main
  stage i enters its top stage with the main stage's liquid enthalpy, its top
  vapor returns to main stage r, 2 ≤ r ≤ i, with its own enthalpy, and its
  bottom liquid is the side product on that product's `side draw` port at
  the product's flow specification (molar or standard liquid volume, the
  latter through the same outer substitution as unstripped products). It is
  stripped by **steam** (port `stripper steam k`, to its adiabatic bottom
  stage; at or above the stripper pressure, positive flow) or by a partial
  **reboiler** (energy inlet `stripper reboiler duty k`) with a **duty** or a
  **boilup ratio** V_N/B. The stripper's total balance gives
  S = B + V_1 − F_steam, so the stripper adds only its stage unknowns
  (temperature, log ratios, log vapor flow per stage) and, with a reboiler,
  its specification residual. Every stripper stage passes the same
  liquid–liquid and phase-split checks as the main column.
- **Start.** A column with side strippers first solves a precursor: the main
  column with each stripper replaced by its bottoms drawn at the draw stage
  and its steam fed to the return stage (the same external balance), with
  the pump-arounds. The coupled system starts from that profile, each
  stripper stage at the draw-stage state, and goes straight to the
  simultaneous Newton solve (free distillate and reflux included); the outer
  flow iteration is only the fallback. Infeasible specifications show up as
  a failed coupled solve, typically a main stage running dry below a draw
  (singular Jacobian) or a stripper whose vapor vanishes.
- **Balances.** The unit balance includes every stripper steam inlet and
  every duty: `duty_w` = reboiler − condenser − Σ pump-around duties +
  Σ stripper reboiler duties. The result reports the pump-around duties and
  each stripper's stage temperatures, draw, vapor return and reboiler duty
  (`DistillationColumnResult::pump_around_duties_w`, `side_strippers`).

**General coupled columns.** The coupling is general in the solver: any
cascade whose top vapor returns to, and whose liquid feed is drawn from,
main-column stages fits the same block structure. It is exposed only as
side strippers of the fractionator. A general MultiFrac-style interface
(several columns with arbitrary interconnecting liquid and vapor streams,
a column's condenser or reboiler feeding another, crude/vacuum or preflash
trains) is not implemented: it needs a column-network specification and
cascades with their own condensers and reboilers.

## Ports

| Index | Port | Present |
|---|---|---|
| 0 | `feed` (inlet) | always |
| 1 | `overhead` (outlet) | always |
| 2 | `bottoms` (outlet) | always |
| 3 | `condenser duty` (energy outlet) | always |
| 4.. | `side draw 1..k` (outlets) | k side products |
| | `feed 2..n` (inlets) | n feeds |
| | `liquid distillate` (outlet) | partial vapor-and-liquid condenser |
| | `stripping steam` (inlet) | steam bottom |
| | `reboiler duty` (energy inlet) | reboiled bottom |
| | `pump-around duty 1..p` (energy outlets) | p pump-arounds |
| | `stripper steam k` (inlet) or `stripper reboiler duty k` (energy inlet) | side stripper k (list order), by stripping agent |

Side draws come first among the variable ports, so a drawn-only fractionator
saved before #190 keeps its side-draw connections at the same indices.

## Failures

Explicit failures, without fallback: a package without characterised
pseudo-components ("the selected package has none") or without the component
characterisation (GERG-2008, EOS-CG-2021: "needs a cubic package with
characterised pseudo-components"); a feed without pseudo-components; side
products that leave no bottoms or distillate; a volume flow the column cannot
deliver; a volume-specified product with a non-pseudo component; a draw on
the condenser or bottom stage; every failure of the distillation solve
(infeasible specifications, loss of two-phase behaviour, a stage liquid that
splits into two liquids, nonconvergence). A selected pseudo-component without
its characterisation input fails when the package is prepared.

Pump-arounds and side strippers add: a pump-around drawn from the condenser
or bottom stage or returned below its draw, a missing flow, return
temperature or duty, or a non-positive duty (specification issues); a return
temperature not below the converged draw-stage temperature, or a returned
liquid that is not a single liquid at the return pressure; a side stripper on
a product that does not exist or a second stripper on one product, a vapor
return below the draw stage, 0 or more than 24 stages, a missing reboiler
duty or boilup ratio (specification issues); a stripper steam inlet without
flow or below the stripper pressure; a stripper reboiler that would remove
heat; a stripper whose top vapor and bottoms are below its steam (no liquid
drawn); and a coupled solve that does not converge (for example a boilup too
small to carry vapor back up the stripper, which runs the main column dry
below the draw).

## Validation

`validation/columns/generate-petroleum-references.py` (teqp 0.23.2 PR,
NumPy, SciPy; no Rust code) characterises an 8-cut crude-like TBP curve
(340–550 K, Watson K 11.8) with the independent Python correlations of
`validation/generate-pseudo-component-references.py` and solves the MESH
equations in its own variables (temperatures, log mole fractions with
explicit summations, log flows, duties, distillate and volume-specified draw
rates as unknowns), started by a Wang–Henke bubble-point iteration. Each
fitted stage phase passes the A4 adapter's multi-start TPD search.
`petroleum-v1.json` holds three cases:

| id | Condenser / bottom | Side products | Notes |
|---|---|---|---|
| `crude-two-draws-reboiled-PR` | total / reboiler (R = 2, D/F = 0.3) | stages 4 and 7, 15 mol/s each | 10 trays, 1.5 bar |
| `crude-volume-draws-partial-pressure-drop-PR` | partial vapor / reboiler (R = 3, D/F = 0.25) | stage 4 at 2.5 L/s and stage 10 at 2.0 L/s standard liquid volume, stage 6 at 10 mol/s | 10 kPa condenser drop, 0.7 kPa per stage |
| `crude-steam-stripped-PR` | partial vapor / 5 mol/s steam at 520 K, 3 bar | stages 4 and 7, 15 mol/s each | water/cut kij = 0.5 (explicit approximation) |

`crates/flowsheet-project/tests/petroleum_fractionator.rs` checks the Rust
assay cuts and standard volumes against the fixture (1e-12), every case
within the declared budgets (temperature 2e-5 K, composition 2e-7, flow
2e-5 mol/s, duty 0.02 W + 2e-6 relative), component and energy balance
closure below 1e-9, the volume specifications to 1e-10, the energy streams,
the failures above, project round trips, and the drawn-only legacy record.
This is implementation agreement for the stated equilibrium-stage model, not
accuracy against plant data or a commercial simulator.

Pump-arounds and side strippers are checked against `petroleum-v2.json`
from `validation/columns/generate-petroleum-side-section-references.py`,
which reuses the property, feed and start pieces of the first generator and
solves its own network MESH formulation (every stage of every cascade with
explicit summations and both flows as unknowns, the stripper draws and
volume-specified product rates as unknowns with residual specifications,
pump-arounds by continuation in their circulation):

| id | Column | Side sections |
|---|---|---|
| `crude-pump-arounds-reboiled-PR` | total / reboiler (R = 2, D/F = 0.3) | pump-arounds 7 → 5 (30 mol/s, return 400 K) and 4 → 2 (20 mol/s, 300 kW) |
| `crude-reboiled-side-strippers-volume-PR` | total / reboiler (R = 3, D/F = 0.3) | reboiled strippers on both products: 2 stages at 150 kW, and 3 stages at boilup ratio 0.5 on a 3.0 L/s volume-specified product |
| `crude-steam-side-stripper-pump-around-PR` | partial vapor / steam bottom | a 3-stage stripper with 2 mol/s steam and a 200 kW pump-around 7 → 6 |

`crates/flowsheet-project/tests/petroleum_side_sections.rs` checks main and
stripper stage temperatures, duties (condenser, reboiler, each pump-around
and stripper reboiler, also on their energy streams), internal and stripper
flows and every product within the same budgets, unit component and energy
balance closure below 1e-9 including all steam inlets, the volume
specification to 1e-10, the failures above, ports, record and result round
trips, and that fractionator records without side sections load unchanged.

## Persistence

`records/columns.rs` adds the fractionator fields behind format version 4
without a version bump (`feeds`, `feed_stages`, `products`, `internals`,
`pressure_pa`, `condenser_pressure_drop_pa`, `pressure_drop`, `condenser`,
`reflux_ratio`/`top_specification`, `distillate_fraction`/
`bottom_specification`, `stripping_steam`). A drawn-only record
`{"kind": "petroleum_fractionator", "side_draws": k}` loads as one feed, ten
trays, a total condenser and a reboiler with unset pressure and
specifications, and is written back unchanged. Stripping steam together with
a reboiler specification is refused.

Pump-arounds and side strippers are two further optional lists behind format
version 4 (no bump), omitted when empty, so every earlier fractionator record
loads without them and is written back unchanged:

```json
"pump_arounds": [{"draw_stage": 7, "return_stage": 6, "molar_flow_mol_s": 25.0,
                  "cooling": {"duty_w": 200000.0}}],
"side_strippers": [{"product": 1, "stages": 3, "return_stage": 3, "bottom": "steam"}]
```

`cooling` is `{"return_temperature_k": T}` or `{"duty_w": Q}`; `bottom` is
`"steam"`, `{"reboiler_duty_w": Q}` or `{"boilup_ratio": r}`. More than 20
pump-arounds or more side strippers than side products are refused on load.
The saved results add `pump_around_duties_w` and `side_strippers` (stage
temperatures, draw, vapor return, reboiler duty) to the column estimate, also
omitted when empty.

## Not yet modelled

- **General multi-column coupling.** Side strippers are solved with the main
  column; other coupled columns (MultiFrac/SCFrac networks, crude/vacuum or
  preflash trains, a column's condenser or reboiler feeding another) are not
  (see [General coupled columns](#pump-arounds-and-side-strippers)). Side
  strippers have no pressure of their own (they run at the draw-stage
  pressure), no condenser and no efficiencies; pump-arounds return liquid
  only (no partial vaporisation on return is checked for a duty
  specification) and have no temperature-drop or return-flow-ratio
  specifications.
- **Cut-point specifications** (ASTM D86 or TBP 5/95 % points, gaps and
  overlaps). The repository converts D86 to TBP for assay input only; it has
  no product distillation-curve calculation from a composition. Products are
  specified by molar or standard liquid volume flow only.
- **Three-phase condenser / free-water decant**, needed for realistic steam
  stripping with a cold reflux drum.
- Tray efficiencies, overflash and heater (furnace)
  specifications, refinery correlations beyond the pseudo-component
  characterisation (for example cold properties, RVP, flash point).
- **Package domain.** Pseudo-components are available in the cubic packages
  only, within 250–600 K; heavier cuts (vacuum residue) fall outside the
  Kesler–Lee ideal-gas and package search domains.
