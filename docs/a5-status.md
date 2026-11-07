# A5 implementation and verification: PH flash and energy foundations

This change builds on A4 at `ccb3cdf`. The public `thermo_solvers::ph_flash`
accepts a cubic model, absolute pressure, target molar enthalpy in J/mol and
explicit `PhFlashOptions`. `TpFlash::enthalpy()` returns the phase-weighted
molar enthalpy required for heaters and ideal isenthalpic valves.

## Algorithm and energy contract

The solver brackets temperature around equilibrium TP evaluations. Every inner
TP/stability failure propagates. It accepts a state only when
`abs(h - h_target) / max(abs(h_target), 1 J/mol)` meets the requested tolerance
(default `1e-10`; permitted positive values up to `1e-6`). A collapsed
temperature interval alone never establishes convergence. Default bounds are
200 and 600 K; both bounds and every intermediate evaluation must be supported.
Callers can narrow the interval explicitly. Targets outside its endpoint
enthalpy range return `OutsideValidityDomain`.

Pure coexistence is handled before the ordinary temperature iteration: at
fixed pressure, bracket the outer-root log-fugacity difference in temperature,
require its absolute residual to be at most `1e-12`, and calculate the vapor
fraction from `beta = (h_target - h_liquid)/(h_vapor - h_liquid)`.
Subcritical single-root states bracket the liquid/vapor sides by EOS volume
relative to the exact PR critical volume. This avoids evaluating saturation
pressures outside the EOS envelope, such as water at the 200 K lower bound.
The pure path supports a single active species with retained zero components.
Bubble/dew endpoints retain both saturation states with phase fractions 0/1;
the zero-fraction phase is an incipient phase, not a finite material amount.

All phase enthalpies use the existing, unchanged NASA7 ideal contribution plus
same-state EOS residual enthalpy. Pure ideal h/s vanish at 298.15 K and 100000 Pa;
real-fluid total enthalpy need not vanish there. Reference identifiers and
model/data provenance accompany each homogeneous phase state. The input target
must use this same reference and molar basis. No new caloric approximation,
density correction or fitted interaction data is introduced.

Default limits are 100 saturation-temperature iterations and 100 outer energy
iterations, each TP solve retaining A4's own limits. `PhFlash` reports outer
iterations, signed scaled energy residual, phase states, and either the TP
stability/split evidence or pure saturation evidence with its iteration count.
Diagnostics identify the selected A5 path. No TP evidence is fabricated for
pure coexistence.

## Independent evidence

`validation/fixtures/ph-v1.json` stores 28 states (14 per variant): three binary
single/two-phase targets, one natural-gas split, propane and water saturation
at vapor fractions 0, 0.25, 0.75 and 1, one narrow-boiling ethane/propane split,
and one natural-gas isenthalpic valve outlet. The latter uses a 250 K,
3 MPa inlet and a 1.5 MPa outlet; SciPy independently solves its outlet energy.

`validation/generate-ph-references.py` reuses the independent A3 teqp AD/Cantera
NASA7 caloric adapter and A4 SciPy equilibrium/teqp pure-VLE adapter, without
Rust calls. It records pinned package versions, source/adapter hashes and
quantity-specific tolerances selected before Rust comparisons. The existing
reference fixtures and production model/data parameters are unchanged.

Reference budgets are temperature `2e-6 K`, enthalpy `2e-5 J/mol + 2e-9 |h|`,
phase fraction `2e-8`, and composition `2e-8 + 1e-7 |x|`. These budgets test
implementation agreement, not experimental accuracy.

An additional mandatory [Clapeyron standard check](../validation/clapeyron-checks.md)
now compares these 28 PH inputs with Clapeyron.jl 0.6.30's native PH solver,
alongside 60 homogeneous A3 states. Julia uses the same explicit parameters,
NASA7 heat-capacity polynomials and pure ideal reference anchor. The standard
check combines offline Rust regression tests with fresh Julia comparisons in
push/PR and scheduled numerical CI. It supplements the teqp/Cantera evidence.

Additional tests cover uniquely recoverable TP-to-PH states (single vapor,
liquid, supercritical and two-phase), zero target enthalpy and its absolute
scaling floor, pure endpoints at both temperature bounds, inactive components,
component permutation, a narrow coexistence interval, isenthalpic expansion,
nonfinite targets, invalid options, out-of-bracket targets, pressure limits
and forced iteration exhaustion. The CLI accepts `ph-flash` JSON cases and
rejects inconsistent temperature, enthalpy and root specifications.

## Supported domain and remaining gates

Computational limits remain 200..=600 K and 0.01..=1e8 Pa for PR/PR-BM. Pure
coexistence is limited to T/Tc <= 0.95. Verification applies to the fixture
states and tested regular branches, not every state inside that rectangular
envelope. Unsupported inner TP states, unresolved coexistence and failed energy
iterations return explicit errors; the solver does not skip such states.

Finite-start stability evidence is inherited from A4 and does not prove global
stability. Near-critical accuracy, general multiple-liquid phases/VLLE,
aqueous hydrocarbon equilibria and experimental VLE accuracy remain unqualified.
The standalone pure-water coexistence test establishes numerical agreement
for the cubic model, not a qualified water property package. General mixture
bubble/dew curves and phase envelopes remain A13. The narrow-boiling test
verifies a finite coexistence interval; arbitrarily near-pure traces are not
guaranteed to resolve below A4's stability thresholds.

See [executed checks](../validation/a5-checks.md). Hosted Windows/Linux CI and
independent code review remain formal acceptance gates. A8 owns the connected
flowsheet energy balances and the complete mix/cool/expand/separate process.
