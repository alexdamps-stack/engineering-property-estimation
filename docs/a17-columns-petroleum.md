# A17 columns: petroleum fractionator (#190)

The petroleum fractionator (`UnitKind::PetroleumFractionator`, Aspen PetroFrac
/ HYSYS crude column equivalent) is calculated since #190. This first scope is
one main column for feeds characterised as petroleum pseudo-components
([pseudo-component status](pseudo-components-status.md), #229). Code:
`crates/flowsheet-core/src/units/columns/petroleum.rs`.

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
  no specific gravity and fail explicitly. A side product may instead be
  specified by a **cut point** (next section).

## Distillation curves and cut-point specifications

Code: `crates/thermo-data/src/distillation_curves.rs` (pure functions) and
`crates/flowsheet-core/src/units/columns/petroleum_curves.rs` (results and
the cut-point iteration).

**TBP curve of a product.** The curve is built from the product composition
and the package's pseudo-components (normal boiling point Tbᵢ, standard
liquid molar volume vᵢ = Mᵢ/(SGᵢ ρw)):

- the components are ordered by Tb; component i boils over the interval
  between the midpoints to its neighbours, [(Tbᵢ₋₁ + Tbᵢ)/2, (Tbᵢ + Tbᵢ₊₁)/2],
  and the end intervals extend half the neighbouring spacing beyond the end
  components (a single component boils at its Tb);
- its share of the curve is its standard liquid volume share xᵢvᵢ/Σxⱼvⱼ, and
  the temperature rises linearly with cumulative volume within the interval.

The intervals are fixed by every pseudo-component of the package, including
those absent from the product, so the curve changes continuously with the
composition. 0 % and 100 % are the lower end of the lightest and the upper
end of the heaviest component present. This is the inverse of the assay
cutting: a product made of equal-temperature cuts in their assay proportions
returns the assay curve at every cut boundary (unit test). Curves are on a
**pseudo-component basis**: database components (water from steam, light
ends) have no boiling-point or standard-density data in the package, are left
out, and their mole fraction is reported with the curve.

**ASTM D86.** D86 = (TBP/a)^(1/b) (°R), the exact inverse of the API
(Riazi–Daubert) correlation TBP = a·D86^b already used for assay input
(M. R. Riazi, T. E. Daubert, Hydrocarbon Processing 59(3), 1980, and Ind.
Eng. Chem. Process Des. Dev. 25, 1986; constants of Fahim, Al-Sahhaf and
Elkilani 2010, Eq. (3.3) and Table 3.2), applied at the tabulated 0, 10, 30,
50, 70, 90 and 95 vol% points of the TBP curve. Between those points D86 is
linear in volume percent; above 95 vol% it is not defined. The correlation was
developed for whole fractions; for narrow cuts the converted D86 points need
not increase with volume and are reported as calculated.

**Results.** Every product with pseudo-components (overhead, bottoms, side
draws, liquid distillate) reports TBP and D86 temperatures at 5, 10, 30, 50,
70, 90 and 95 vol% in `DistillationColumnResult::product_curves`, saved as
`distillation_column.product_curves` behind format version 4 (results
saved before it load with no curves) and shown in the desktop results.

**Cut-point side products.** `DrawRate::CutPoint` specifies the TBP (0–100
vol%) or D86 (0–95 vol%) temperature of a side product at a volume percent,
with an optional starting molar flow (default: hydrocarbon feed flow divided
by twice the number of side products plus one). The molar flows of all
cut-point products are the unknowns of an outer Broyden iteration on their
logarithms, started from a forward-difference Jacobian; each evaluation is a
converged column, including any standard-volume products. Steps are limited
to a factor e^0.5 and halved (at most four times) until the residual norm
decreases. Converged when every cut point is met within 1e-6 K
(`CUT_POINT_TOLERANCE_K`), at most 25 iterations. Explicit failures: a TBP
target outside the pseudo-components' TBP span; a column solve that fails on
two consecutive steps ("not attainable", with the closest temperature and
flow); a residual that no longer decreases after a Jacobian refresh;
interdependent cut points (singular Jacobian); no convergence. Cut points
on the distillate or bottoms, gap/overlap specifications and light ends in
the curves are not implemented.

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

`validation/columns/generate-petroleum-curve-references.py` (Python standard
library, no Rust code) gives `petroleum-curves-v1.json`: its own TBP and
inverse-D86 implementation on four synthetic component sets, on every
product of `petroleum-v1.json` at the reference composition, and on example
E3.1 of Fahim et al. read backwards (printed TBP gives the printed D86 within
0.105 °C, the propagated print rounding). The cut-point case replaces the two
molar side products of `crude-two-draws-reboiled-PR` by cut points read from
its reference products (side draw 1 TBP 95 %, side draw 2 D86 50 %), so the
expected answer is the reference column. `thermo-data/tests/distillation_curves.rs`
checks the curve functions against the fixture (1e-9 K) and E3.1;
`flowsheet-project/tests/petroleum_curves.rs` checks the reported curves of
the solved two-draw and steam columns within the composition budget
propagated per point (Σ|∂T/∂xᵢ|·2e-7, recorded in the fixture), the cut-point
column against the petroleum-v1 budgets (the recovered draw flows agree to
7e-7 mol/s), the cut points met to 1e-6 K, unattainable and invalid cut
points, the stored record, a save/reload round trip and results without
curves.

This is implementation agreement for the stated equilibrium-stage model, not
accuracy against plant data or a commercial simulator.

## Persistence

`records/columns.rs` adds the fractionator fields behind format version 4
without a version bump (`feeds`, `feed_stages`, `products`, `internals`,
`pressure_pa`, `condenser_pressure_drop_pa`, `pressure_drop`, `condenser`,
`reflux_ratio`/`top_specification`, `distillate_fraction`/
`bottom_specification`, `stripping_steam`). A drawn-only record
`{"kind": "petroleum_fractionator", "side_draws": k}` loads as one feed, ten
trays, a total condenser and a reboiler with unset pressure and
specifications, and is written back unchanged. Stripping steam together with
a reboiler specification is refused. A cut-point side product is stored as
`{"rate": {"cut_point": {"basis": "tbp" | "astm_d86", "volume_percent": …,
"temperature_k": …, "estimate_mol_s": …}}}` (estimate omitted when unset),
also behind format version 4.

## Not yet modelled

- **Side strippers** (small columns on side draws returning vapor to the main
  column) and **pump-arounds** (side cooling loops). Both need coupled
  cascades or side heat removal (stage energy ports) in the MESH solve; the
  distillation solver has neither. A side stripper can be drawn today as a
  separate absorber/stripper unit with a recycle, but that is a sequential
  approximation, not a simultaneous solution.
- **Further curve specifications.** Cut points on the distillate or bottoms
  (they would replace a top or bottom specification), 5–95 % gap and overlap
  specifications between products, and light ends in the curves (database
  components lack boiling-point and standard-density data). Side-product cut
  points are implemented (above).
- **Three-phase condenser / free-water decant**, needed for realistic steam
  stripping with a cold reflux drum.
- **Multi-column systems** (MultiFrac, SCFrac, coupled crude/vacuum or
  preflash columns), tray efficiencies, overflash and heater (furnace)
  specifications, refinery correlations beyond the pseudo-component
  characterisation (for example cold properties, RVP, flash point).
- **Package domain.** Pseudo-components are available in the cubic packages
  only, within 250–600 K; heavier cuts (vacuum residue) fall outside the
  Kesler–Lee ideal-gas and package search domains.
