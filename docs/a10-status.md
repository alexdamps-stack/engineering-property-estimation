# A10 — Recycles and larger networks

Implementation and verification record, November 8, 2026.

A10 implements bounded sequential-modular recycle solving in `flowsheet-core`
and uses it for project and desktop calculations. The defined recycle case,
cancellation, nonconvergence, persistence and result invalidation are executable
tests. The existing A8 acyclic API and unit-operation balance budgets remain.

## Run the demonstration

```sh
cargo run -p thermo-desktop --locked -- --recycle-demo
cargo run -p thermo-desktop --locked -- examples/recycle.flowsheet.json
cargo run -p flowsheet-project --example recycle --locked
```

In the editor, choose **File > Load recycle demonstration**, then **Solve (F5)**.
The project is editable and supports the same background solve, Cancel, save,
reload and auto-solve workflow as A9. The toolbar reports recycle iteration;
Messages reports closure and selected stream diagnostics show the tear residuals.
An optional path passed to the headless example writes an input-only project;
the example prints the calculated result/provenance/closure report to stdout.
The committed JSON contains inputs, not reference answers.

## Calculation contract

- `solve` / `solve_with_progress` and `calculation_order` retain their acyclic
  contract. Their cycle diagnostic points callers to the recycle APIs.
- `solve_network` / `solve_network_with_progress` discover cycles deterministically
  and remove one edge from each remaining cycle, preferring mixer inlets. This
  does not minimize the number of tears. Acyclic graphs use the original path.
- Automatic starts originally used zero recycle flow with the first feed's
  flashed state. Since the #87 fix, each tear uses a feed reachable upstream
  along material connections, preserving exchanger-side separation;
  mixer tears retain zero initial flow, while other tears start with that feed's
  flow. A direct source outlet-temperature specification refines the state.
  `solve_recycles` accepts explicit tear-stream initial states and
  `RecycleOptions`. Unknown tears, incomplete cycle cuts, invalid flows,
  missing inputs, invalid options and unsupported equipment fail explicitly.
  Initial states are re-flashed with the calculation package using their
  composition, pressure and enthalpy; cached phase properties are not reused.
- Direct substitution reads **previous-pass** tear values at every torn inlet,
  independently of whether its source has already run in the current pass.
  Every pass produces a fresh set of unit and stream results.
- Defaults are **500 passes** and **1e-10 maximum scaled closure residual**.
  Closure covers total and component molar flows, enthalpy flow, T, absolute p,
  molar h and mole fractions. A dimensional residual divides by
  `max(abs(old), abs(new), 1 SI unit)`; mole-fraction residuals are absolute.
  Thermodynamic state mismatch is ignored only when both streams have zero
  flow. A small state change alone cannot establish convergence.
- Before accepting a calculation, every internal unit's balance is checked
  again against its **final connected** inlet/outlet stream results. Existing
  budgets remain `1e-8` relative component balance and `2e-9` scaled energy
  balance. Tear closure and the connected-stream balances must both pass.
- Progress includes the 1-based recycle pass (zero for initialization/acyclic
  execution). Cancellation is checked before tear initialization flashes and
  before every unit on every pass. It cannot interrupt a flash already running.
  Cancellation, numerical failure and exhaustion return no partial result.
  Exhaustion reports the limit, last residual, target, tear identifiers and any
  final connected-stream balance failure.

## Defined case and independent conservation solution

The natural-gas composition is `[0.88, 0.06, 0.03, 0.02, 0.01]` in the A8
component order methane/ethane/propane/nitrogen/carbon dioxide. One fresh feed
at 330 K, 5 MPa and 100 mol/s mixes adiabatically with its recycle, heats to
350 K without pressure loss, then splits equally between recycle and product.
The package uses the same explicit synthetic zero-kij approximation and
200..400 K PH interval as A8; this is not measured natural-gas accuracy.

Conservation gives the expected values without running the recycle algorithm:

```text
R = f * (F + R) = f * F / (1 - f)
circulation = F / (1 - f)
product = F
product and recycle composition = fresh composition
heater duty = F * (h(350 K, p, z) - h(330 K, p, z))
```

For `f = 0.5`, recycle and product are 100 mol/s and circulation is 200 mol/s.
The PR-BM example takes **34 passes**, with maximum scaled tear residual
**5.82077e-11**. Reported recycle/product flow is 99.99999999417923 mol/s.
The energy check uses package enthalpies in an independent conservation
expression; it is not an independent thermodynamic model comparison. The
existing Julia/Clapeyron comparisons qualify the thermodynamic implementation.

The core tests additionally cover PR and PR-BM, flow factors 0.01/1/100,
different initial flows/compositions/temperatures, zero recycle, 90% recycle,
reversed unit identifiers, two interacting loops, and two loops with a 100-unit
downstream chain (109 total units). A 100% recycle with positive fresh feed
has no finite steady-state solution and is tested as nonconvergence.
An additional 1e-14 mol/s fresh-feed case verifies that small absolute tear
changes cannot bypass the final relative component-balance budget.

## Projects and interface

Input fingerprints now traverse the complete upstream closure iteratively and
hash each upstream unit and its inlet connections once in identifier order.
Editing any loop member invalidates that loop and its downstream results;
unaffected feeds and independent branches retain their currency. There is no
recursive cycle hashing. Feed composition and connection changes are covered.

Successful results store iteration count and per-tear closure evidence in an
optional `recycle` field. At the A10 milestone, version-1 project files remained readable;
older applications could not read this new result field. Current writers use format
version 2 and readers accept versions 1 and 2 (see `flowsheet-and-ui.md`). The engine provenance was
`flowsheet-core a10-tears-v1`, so results from the previous engine load stale.
Failed or cancelled solves do not replace stored results; input fingerprints
keep obsolete values visibly stale.

No new equipment graphics are introduced. Existing traced ISO source geometry,
attribution/license and separately defined ports are reused.

## Verification

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
cargo doc --workspace --no-deps --locked
julia --project=validation/julia validation/julia/check.jl
```

Clippy and documentation run on Windows/MSVC. The full test suite runs in WSL
Ubuntu with `WGPU_BACKEND=vulkan` and a separate Linux target directory; UI tests
run headlessly with Mesa. UI acceptance covers solving the recycle, editing it,
checking stale statuses, recalculation, save/reload and a visible unit failure.
The solved recycle screenshot is committed. Three A9 baselines were visually
reviewed and refreshed solely for the new engine provenance text.

The live standard check uses **Julia 1.11.7 / Clapeyron.jl 0.6.30** and passes
PR/PR-BM homogeneous and PH comparisons plus NRTL/gamma-phi comparisons.
Committed numerical fixtures and tolerances are unchanged.

All commands above passed. The final focused project/desktop run also passes
22 UI tests (one optional performance measurement ignored), including all
seven reviewed image snapshots. The small-flow balance guard passes separately
after the full workspace run.

## Limits

Direct substitution need not converge for arbitrary loops. There is no
Wegstein/Newton recycle acceleration, optimal tear selection, or desktop tear
configuration form. The explicit API supplies process-specific starting states
and iteration limits when the automatic starting state is unsuitable.
Automatic initialization rejects loops without an upstream reachable feed.
Pressure-loss loops still need a supported pressure-raising operation; drawing
a pump or compressor does not implement it. No new flash domain, equilibrium
stability claim, or equipment capability is introduced by A10.
