# A17 columns: equilibrium-stage foundation

Status: foundation only (no unit is enabled by it). Module
`crates/flowsheet-core/src/units/columns/stages.rs`, declared by one line in
`units/columns.rs`. It serves rigorous distillation (#186) and the absorber /
stripper (#187). Shortcut distillation (#185), the whole-column nonlinear
solver, packing / rate-based correlations and project or desktop enablement
are separate work.

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

## Remaining prerequisites

1. Column solver (#186/#187): unknown selection, specifications (reflux
   ratio, distillate rate, purity/recovery), initialization, Newton or
   inside-out iteration, convergence diagnostics, phase-existence handling
   on stages (single-phase stages, supercritical components).
2. Column spec fields and records (stage count, feed and draw stages,
   condenser type, pressure profile, efficiencies) behind format version 4,
   forms and `status()` switch.
3. Fixed-approach ends need a utility inlet temperature: energy streams carry
   duty only, so either a temperature-level field on energy streams or
   utility material ports (shared change through the integrator).
4. Absorber has no energy ports; any stage duty there needs ports added in
   `columns::info` (shared columns file).
5. Packed sections: theoretical stages from height/HETP map onto this cascade;
   rate-based segments need their own residuals.
