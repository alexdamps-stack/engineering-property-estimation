# A8 steady-state flowsheet engine

Implemented November 8, 2026 on `alexdamps-stack/A8`, based on freshly fetched
`origin/main` commit `00cc8fe7b6998cef57306d9465de242dfdb35b07`. No uncommitted
changes were present when work began. The repository-pinned Rust is 1.99.1.

| Issue | Delivered evidence |
|---|---|
| #47: resolve decisions | `flowsheet-and-ui.md` section 5 records the A8 choices, JSON project-format decision and the separate A9 GUI framework gate |
| #48: streams/connections/order | `Flowsheet::calculation_order` and `solve`; explicit feed identity/basis checks; phase-resolved material results; deterministic ordering independent of unit creation order; branched/rejoined graphs; recycle rejection |
| #49: unit operations | Native TP/PH feeds, mixers, heaters/coolers, isenthalpic valves and equilibrium separators; component and energy balance checks; two-outlet splitters; two-stream exchanger specifications already exposed by the prototype |
| #50: headless demonstration | Executable example plus connected PR and PR-BM comparisons against all 14 independent A2 process states and applicable flow scaling |
| #46: parent acceptance | All of the above, with the verification record below |

## Calculation contract

`Flowsheet::solve(&CubicPackage, &FeedCompositions)` returns an atomic
`Calculation` or a `SolveError` naming the failing unit. Structural diagnostics
are checked before any property call; unsupported rotating equipment and missing,
misordered or incompatible feed compositions also fail preflight. Enum-based
specifications prevent conflicting temperature/duty or pressure/drop inputs.
Inlet-dependent constraints, such as a mixer pressure above the lowest flowing
inlet or a valve pressure rise, fail at the affected unit. Inner solver and
domain failures are propagated, never clamped or replaced by previous values.

The immutable package fixes component order, exact PR or PR-BM, explicit full
interaction matrix and its source, data, caloric reference and PH options.
It reuses the existing qualified Rust TP/PH path without modifying its equations,
fixtures or comparison tolerances. The example explicitly selects synthetic
zero interactions; the API never invents missing pairs. NRTL/GERG connected
package adapters and conversion between packages are not implemented in A8.

Material streams contain mol/s and an optional equilibrium state, including
phase fractions/compositions, phase properties, total molar h/s, mass flow,
model/data/reference metadata, convergence residuals and original diagnostics.
Units report signed duty (W, positive into the unit), warnings and balances.
No global stability or experimental-model accuracy is implied.

An unchanged adiabatic separator copies the actual inlet equilibrium phases,
preserving saturation endpoint roots. With pressure drop or thermal specification
it performs the appropriate PH or TP solve. An absent phase outlet has exactly
zero flow and no thermodynamic state, plus a unit warning. A single supercritical
fluid goes to the vapor port by an explicit routing convention, retains its
supercritical label, and emits a warning. Splitter endpoints are supported;
mixers ignore zero-flow branches. All-zero mixers and thermal operations without
positive material flow fail explicitly instead of dividing by zero.

For internal units, component residuals are outlet minus inlet mol/s, scaled
by total inlet mol/s (minimum positive floating-point floor). The acceptance
limit is `1e-8`. Energy residual is outlet enthalpy flow minus inlet enthalpy
flow minus duty, scaled by `max(sum(abs(inlet/outlet enthalpy flows)) + abs(duty),
1 W)`; the limit is `2e-9`. Nonfinite flows, powers and balances fail. Boundary
feeds/products have no internal balance. The passive two-stream exchanger uses
positive hot-to-cold transfer, conserves total energy, and rejects terminal
temperature crossover; no UA/pinch design or sizing model is claimed.

`Calculation::is_current` compares the full topology/specification snapshot,
including its revision, feed compositions and shared immutable package identity.
Cloned packages retain identity; newly selected packages conservatively require
recalculation. A successful report never acquires new results after an input edit.
Consumers must check currency; A9 owns the desktop background worker, visual
invalidation, save/load, cancellation, and connection to the Solve button.

## Reproduce and supported scope

```sh
cargo run -p flowsheet-core --example demonstration --locked -- PR-BM
cargo run -p flowsheet-core --example demonstration --locked -- PR 100
cargo test -p flowsheet-core --all-targets --locked
```

The executable accepts optional PR/PR-BM selection and a positive feed-flow
multiplier. It emits JSON containing eight unit reports and seven material
streams. Canonical inputs are in `crates/flowsheet-core/examples/support/mod.rs`;
no expected outputs are read by the executable. Test expectations come only
from the unchanged `validation/fixtures/process-v1.json`, independently generated
with teqp/Cantera/SciPy. The A2 feed compositions, 330/350 K feeds at 5 MPa,
240 K cooler, and 1.5 MPa valve are retained.

The independent test runs both variants at flow multipliers 0.01, 1 and 100,
checks all seven stream temperatures, absolute pressures, compositions, phase
counts/fractions/compositions, enthalpies, cooler duty and valve enthalpy, and
requires component/energy closure. Duty-specified cooling also recovers the
independent process. Other tests cover thermal modes, pressure drops, branch
rejoining, splitter endpoints, empty outlets, input invalidation, unsupported
equipment, missing/invalid inputs, recycle rejection and PH nonconvergence.

The example's PH interval is 200..400 K. The independently qualified states
are the fixed A2 points and their flow-scaled counterparts; this is not a
continuous phase-map qualification. The package permits an explicit PH search
interval within the existing 200..600 K cubic domain. Existing model/solver
limitations and the A2 exclusions still apply. No recycles, reactions, shaft
work, kinetic/potential energy, solids, VLLE or new critical-region claim is added.

## Verification

The live Julia 1.11.7 standard check passed for PR/PR-BM and NRTL, with all
committed fixtures unchanged. Windows Smart App Control blocked a freshly
compiled dependency build script (OS error 4551), so executable verification
uses WSL Ubuntu against this exact A8 source directory, with the pinned Rust
1.99.1 and build output on the Linux filesystem. Native Windows execution is
not claimed.

Final checks on this source, November 8, 2026:

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | Passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | Passed |
| `cargo test --workspace --locked` | 164 passed, 0 failed, 0 ignored; doc-test targets also passed |
| `cargo test -p flowsheet-core -p thermo-desktop --all-targets --locked` | 62 passed, including all 13 new engine integration tests and 13 desktop interaction/snapshot tests |
| `cargo doc --workspace --no-deps --locked` with `RUSTDOCFLAGS=-D warnings` | Passed |
| `cargo build --workspace --locked` | Passed |
| `python3 validation/consolidate-tolerances.py` | Passed; approved fixture budgets unchanged |
| Julia 1.11.7 `--project=validation/julia validation/julia/check.jl` | Live PR/PR-BM and NRTL comparisons passed; committed fixtures unchanged |
| Headless example executable | PR-BM at 1x and PR at 100x passed; seven streams and eight unit reports |
| `git diff --check` | Passed |

The selected-cooler desktop snapshot was intentionally regenerated and visually
reviewed because the old message incorrectly said A8/A5 did not exist. Its
replacement reports pending desktop integration. All other UI baselines and
all independent thermodynamic references remain unchanged. A second fetch
confirmed that `origin/main` still matched the recorded base commit.
