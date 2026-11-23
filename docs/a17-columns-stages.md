# A17 columns: equilibrium-stage foundation

The stage cascade in `crates/flowsheet-core/src/units/columns/stages.rs`
supports rigorous distillation (#186) and adiabatic absorption and stripping
(#187) through the flowsheet engine, project records and desktop.

## Contents

| Item | Purpose |
|---|---|
| `Stage` | T, p, `L_j`, `V_j`, side draws `U_j`/`W_j`, `x_j`, `y_j`, duty `Q_j`, Murphree vapor efficiency |
| `StageFeed` (`from_stream`) | A feed entering one stage as a whole, built from a flowsheet `MaterialStream` |
| `StageProperties`, `CalcContext::stage_properties` | `K_ij` and liquid/vapor molar enthalpies from the cubic package at the stage state |
| `StageCascade::residuals` | MESH residuals per stage (M, E with Murphree, S, H) |
| `StageCascade::total_flow_residuals`, `balance` | Stage total-flow residuals and overall component/energy balance |
| `StageCascade::solve_liquid_fractions` | Thomas-algorithm component balance at fixed flows and K (Wang–Henke inner step) |
| `ColumnEnds`, `Condenser` | Condenser (stage 1) / reboiler (stage N) checks and energy-port values |
| `EndApproach` | Fixed condenser/reboiler terminal approach as an end-stage temperature |

## Conventions

Stages are numbered 1..N from the top. `L_j` leaves downward, `V_j` upward;
`V_1` is the top vapor product and `L_N` the bottom liquid. Side draws leave
at stage composition. Duty `Q_j` is positive into the stage. Residuals are in
minus out (Seader, Henley and Roper, *Separation Process Principles*, 3rd ed.,
Eqs. 10-1 to 10-13, duty sign flipped to the project convention). The bottom
stage is always an equilibrium stage.

## Condenser / reboiler contract

- Condenser is stage 1, reboiler stage N (HYSYS / Aspen Plus numbering).
  Total condenser: `V_1 = 0`, distillate = `U_1`, reflux = `L_1`. Partial
  condenser: vapor distillate `V_1`, optional liquid `U_1`.
- Only end stages with a condenser/reboiler may carry duty; side heaters and
  coolers need new energy ports first.
- `ColumnEnds::port_duties` returns `UnitOutput::energy` in the existing port
  order of `columns::info` (`condenser duty` outlet, `reboiler duty` inlet):
  `−Q_1` and `Q_N`, positive in the port direction. `duty_w = Q_1 + Q_N`.
- A fixed terminal approach follows the exchangers' hot-minus-cold
  definition: condenser `T_1 − T_coolant,in = ΔT`, reboiler
  `T_heating,in − T_N = ΔT`. It becomes an end-stage temperature
  specification replacing a duty or temperature specification.

## Validation

`validation/columns/generate-stage-references.py` (standard library only)
builds a 5-stage, 3-component distillation cascade and a 4-stage,
2-component absorber with a vapor side draw. It solves the component
balances exactly with rational Gauss–Jordan elimination and closes stage
duties. `crates/flowsheet-core/tests/columns_stages.rs` checks the Thomas
solution against it (relative 1e-12), zero MESH residuals, matching
summation residuals, residual sums against the overall balance, end-port
duties, and a single stage built from a cubic TP flash (ethane/propane, PR).

## Rigorous distillation

- Specify one absolute column pressure, a positive molar reflux ratio
  (L/D), and the fraction of total feed withdrawn as distillate. The
  distillate fraction must be strictly between 0 and 1. Each entering feed
  retains its enthalpy through any pressure let-down.
- Choose 1–24 ideal internal trays or one [ideal packed bed](a17-columns-packing.md)
  spanning the internal stages. The total condenser is stage 1 and the
  partial reboiler is the last stage. Numbered feeds can enter from stage 2
  through the reboiler. Liquid side draws specify an internal stage and a
  positive molar flow. The distillate and draws must leave positive bottoms
  flow.
- The vapor flowing to the total condenser is fully condensed to the bubble
  point. Distillate and reflux share that liquid state. Bottoms and liquid
  side draws use the equilibrium liquid of their stage. The condenser removes
  heat and the reboiler adds heat. Their duties are calculated and report on
  the existing optional energy ports.
- The solve eliminates the summation and total-flow equations by expressing
  compositions as normalized log ratios and obtaining each internal liquid
  flow from the corresponding vapor flow and the cumulative material balance.
  A damped Newton solve closes the remaining component, log-fugacity and
  internal energy equations. The column then TP-flashes the midpoint of every
  solved tie line to verify its two phases before publishing results.
- The solve uses PR or PR-BM and the bounded A13 bubble-point domain
  (active methane, ethane, propane, nitrogen and carbon dioxide,
  200–350 K and 10 kPa–10 MPa). It also honors the selected package
  temperature limits. It reports an explicit failure for infeasible duties,
  loss of two-phase behavior or nonconvergence.
- Existing version-4 projects keep their numbered feed, draw and energy
  connections. Legacy columns open with five ideal internal trays, feeds at
  stage 2, and unset pressure, reflux, distillate split and draw flows.
  The desktop shows the stage temperatures, duties and applicable bed height.

The independent [distillation references](../validation/columns/README.md)
cover a single-feed PR tray column and a two-feed PR-BM packed column with a
liquid side draw. Both run through full project solving and save/reload.

## Absorber and stripper

- Connect a single-phase liquid to the top inlet and a single-phase gas to
  the bottom inlet. Specify an absolute pressure at or below both inlet
  pressures. Each inlet retains its enthalpy through any pressure let-down.
- Choose 1–24 ideal trays or one [ideal packed bed](a17-columns-packing.md)
  covering all contacting stages. The liquid enters stage 1 and the gas
  enters the last stage. Gas leaves at the top and liquid leaves at the
  bottom. Both feeds enter the same stage when there is one contact.
- Every stage is adiabatic. The damped Newton solve closes the component,
  equilibrium and energy equations and calculates both outlet flows. It
  starts from a two-phase PH flash of the combined feed. The solver checks
  each stage against a separate TP flash of its tie-line midpoint.
- The column uses the same PR and PR-BM pressure, temperature and active
  component limits as rigorous distillation. It requires two positive
  internal phase flows at each contact and reports initialization or
  convergence failures in the flowsheet.
- Existing version-4 projects keep their four material connections. A legacy
  absorber opens with five ideal trays and unset column pressure. The
  desktop reports stage temperatures, bed height when applicable, and
  the gas and liquid outlet states.

Independent [absorber references](../validation/columns/README.md) cover
propane uptake in a PR tray column, ethane stripping in a PR-BM packed
column, and one ideal contact.

## Further column models

The current calculation uses one pressure and ideal equilibrium stages.
Pressure profiles, tray efficiency, partial condensers, purity or recovery
specifications, and vapor draws are additional solver work. Fixed-approach
ends require a utility inlet temperature. Stacked sections, HETP
correlations, hydraulic rating and rate-based packed segments remain
separate work.

## Packed sections (joint contract with `packing.rs`)

`PlacedPackedSection { top_stage, packed }` in `stages.rs` places the ideal
stages of a resolved `PackedStages` in the cascade, so the two modules agree:

- **Numbering.** `PackedStages::stage_at_depth` counts from 0 within the bed;
  cascade stages are numbered from 1 at the column top. The bed occupies column
  stages `top_stage..=top_stage + stages − 1`, and `column_stage_at_depth`
  returns `top_stage + stage_at_depth(depth)`.
- **Feed placement.** A feed or draw at depth `d` enters that one stage as
  a whole (`StageFeed`). The top face belongs to the first bed stage and the
  bottom face to the last. A point on an internal boundary belongs to the
  lower stage. In an absorber bed spanning the column, lean liquid at the top
  face is stage 1 and rich gas at the bottom face is stage N, matching the
  absorber ports.
- **Ends.** The condenser (stage 1) and reboiler (stage N) are separate
  equilibrium stages and never part of a bed; `check` rejects overlap.
- **Efficiency.** Packed stages carry Murphree efficiency
  `PackedStages::stage_efficiency()` = 1; `check` rejects any other value.
  Trayed stages keep their own Murphree efficiency in `Stage::efficiency`.
  An *overall* (O'Connell-type) efficiency is not a stage property. It
  would convert actual trays to theoretical stages before the cascade is
  built, and it is not implemented.

The distillation solver places one section over all internal stages.
