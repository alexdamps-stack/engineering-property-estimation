# A13 saturation and phase-envelope acceptance

> **#227 update.** The cubic domain now covers every component of the cubic
> package (water, methanol, 1-butanol, oxygen, hydrogen, carbon monoxide and
> pseudo-components, besides the dry set) over the package temperature interval
> (at most 200–600 K), with multiphase evidence, three-phase bubble points, and
> explicit near-critical and freezing failures. See
> [the extended domain](#227-extended-domain) below; the original dry record
> that follows is unchanged.

A13 implements PR and PR-BM bubble/dew pressure and temperature in
`thermo-solvers`, delegates them through `CubicPackage`, and exposes reproducible
JSON CLI examples. It also samples either envelope boundary over a bounded
temperature interval, returns the actual verified temperature/pressure range,
and stops before returning a failed or near-critical sample. Issues #66–#69
are the acceptance scope; flowsheet persistence and UI specifications are separate.

## Formulation and supported scope

The specified composition is liquid for a bubble point and vapor for a dew
point. The unknowns are ln(p) or ln(T), plus ln(K_i) for active components.
The equations are ln(K_i) + ln(phi_i^V) - ln(phi_i^L) = 0 and
ln(sum(z_i K_i)) = 0 (bubble) or ln(sum(z_i/K_i)) = 0 (dew).
Zero fractions remain zero in both phases. For a single active component,
A13 canonicalizes the accepted near-unit fraction to exactly one before either
pressure or temperature solving, retaining all zero-fraction slots. This is
local to saturation and leaves `CubicModel::with_fractions` unchanged. Thirteen Wilson-initialized seeds
span the requested log-coordinate interval; finite-difference Newton steps use
pivoted elimination and backtracking within the requested bounds. Discovered
multiple physical roots produce an explicit unsupported-branch error rather
than implicit retrograde selection. This finite search does not prove that all
roots have been discovered.

A returned point must pass fugacity equality to 1e-9 in ln(f), mechanically
stable homogeneous evaluation, and the existing finite Wilson/near-pure TPD
trials for both phases. The returned phase states retain model/data/caloric
metadata and the result retains both TPD assessments. Finite trials provide
evidence, not proof of global stability. No solids, water-bearing hydrocarbon
VLLE, or third liquid phase is implemented.

The computational limits remain 200–600 K and 0.01–1e8 Pa absolute. Defaults
are 200–350 K and 1e4–1e7 Pa; callers must explicitly supply other bounds.
Package calls also enforce the package's existing temperature interval.
Active species are restricted to methane, ethane, propane, nitrogen and carbon
dioxide. Independent qualification covers the specific dry mixtures and
states in `validation/a13-inputs.json`, with exact PR factors, no volume
translation, and explicit synthetic zero kij. Broader bounds are search
limits, not experimental accuracy or qualification claims. GERG and NRTL
saturation are not extended by A13.

Near-critical guards require vapor/liquid molar-volume ratio >= 1.5 and,
for mixtures, maximum absolute composition separation >= 0.03. Pure propane
reuses the existing spinodal-bracketed fugacity solver and its T/Tc <= 0.95
limit; temperature is found by a bounded bisection of that pressure function.
These conservative guards do not locate or estimate the critical point.
Trivial single-phase roots and singular/coalescing boundary equations are
rejected. Failure to find a boundary is never converted into a phase result.

The trace samples one requested branch, never interpolates across a failed
point, and never extrapolates through the critical region. It requires a
positive finite temperature step, at most 1000 points, a sufficient point
budget for the requested interval, and 1–2000 iterations per seed. Its stop
is `RequestedRangeComplete`, `NearCritical`, or `UnverifiedBoundary`, carrying
the failed temperature and diagnostic where applicable. An empty trace has
no supported range. Supported ranges describe discrete verified samples;
they do not claim a continuous global envelope or retrograde continuation.

## Independent acceptance evidence

The new generator `validation/julia/generate_a13.jl` uses **Julia 1.11.7** and
locked **Clapeyron.jl 0.6.30** native `bubble_pressure`, `dew_pressure`,
`bubble_temperature`, and `dew_temperature`. It reuses only the previously
qualified independent model-construction adapter from unchanged `generate.jl`.
It never calls Rust, consumes Rust output, or supplies reference solutions as
initial guesses. Input, generator, model-adapter and manifest hashes are
recorded in the new fixture. Existing committed fixtures/generators and their
tolerances remain unchanged.

`validation/fixtures/a13-v1.json` contains 62 independent points: 32 pressure/
temperature cases (both PR and PR-BM for equimolar methane/propane, the original
natural-gas feed, the mixed process feed, and pure propane) plus 30 PR-BM
pressure points along four bounded envelope segments. Offline tests compare
T, p, both compositions, both molar volumes, and fugacity residuals. Trace tests
compare the actual returned samples with those independent pressure points.
The 200-225 K natural-gas dew segment explicitly limits pressure to 4 MPa:
at 225 K the broad default interval contains multiple physical dew roots,
which the solver rejects; the qualified lower-pressure branch is selected
by its stated bounds, not an implicit root choice.

Numerical budgets were fixed before comparing Rust with the new references:
T: 2e-6 K absolute; p: 0.02 Pa + 2e-8 relative; composition: 2e-8 absolute +
1e-7 relative; molar volume: 1e-12 m³/mol + 2e-8 relative. These are matched
implementation budgets, not physical-model accuracy estimates.

The original natural gas is z = [0.88, 0.06, 0.03, 0.02, 0.01] in order
[methane, ethane, propane, nitrogen, carbon-dioxide]. Its independently
calculated hydrocarbon dew temperature at 1.5 MPa absolute is **214.601242213 K
(PR-BM)** and **214.617431732 K (PR)**. The mixed demonstration process feed
[0.85, 0.075, 0.04666666666666667, 0.018333333333333333, 0.01] has dew
temperatures **223.062355807 K (PR-BM)** and **223.086292610 K (PR)** at that
pressure. These are water-free hydrocarbon dew points under the pinned
synthetic interaction approximation, not water-dew-point predictions.

Pure-propane tracing over requested 300–360 K at 5 K steps returns 300–350 K
and stops with a near-critical diagnostic at 355 K. The natural-gas bubble
trace also exercises a mixture near-critical guard. Failure tests cover
invalid and extreme step sizes, inadequate budgets, invalid/nonfinite bounds,
iteration exhaustion, unsupported wet components, component reordering and
retained zero fractions. Package tests check composition and package bounds;
CLI tests reject conflicting specifications.

## Reproduce

```sh
cargo run -p thermo-cli -- examples/a13-bubble-pressure-binary.json
cargo run -p thermo-cli -- examples/a13-dew-temperature-natural-gas.json
cargo run -p thermo-cli -- examples/a13-bubble-temperature-propane.json
cargo run -p thermo-cli -- examples/a13-dew-envelope-natural-gas.json
cargo run -p thermo-cli -- examples/a13-bubble-envelope-propane.json
cargo test -p thermo-solvers --test a13_clapeyron --locked
julia --project=validation/julia validation/julia/check.jl
```

CLI calculations are `bubble-pressure`, `dew-pressure`, `bubble-temperature`,
`dew-temperature`, `bubble-envelope`, and `dew-envelope`. Optional
`saturation_temperature_bounds_K` and `saturation_pressure_bounds_Pa` are
ordered two-element arrays. Envelope inputs also accept `envelope_step_K` and
`envelope_maximum_points`; omit fixed T/p. A failed boundary returns the CLI's
existing nonzero error outcome; a trace returns its verified prefix and
explicit stop, so callers must inspect the stop before using its range.

The standard Julia check now validates A13 provenance, generates a temporary
A13 candidate and runs the same Rust tests using `THERMO_A13_REFERENCE`, in
addition to every existing PR/PR-BM and NRTL check. Reference regeneration is
an explicit review operation using `generate_a13.jl --output candidate.json`;
CI never overwrites committed expected values.

## Executed checks

Executed on 2026-11-09 with Rust 1.99.1, combined A11/A12 main baseline `4887212`,
Windows native scoped tests and WSL Linux full-workspace/reference checks.

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | PASS |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | PASS |
| `cargo test --workspace --locked` (Linux) | PASS: 231 passed, 0 failed, 1 existing ignored performance measurement |
| `cargo test -p thermo-solvers -p property-packages -p thermo-cli --locked` (Windows, pure-fraction correction before A11 rebase) | PASS: 70 passed, 0 failed |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | PASS |
| `python validation/consolidate-tolerances.py` | PASS for all registered fixtures, including unchanged A13 budgets now reconciled in the shared registry/table |
| Julia 1.11.7: `julia --project=validation/julia validation/julia/check.jl` | PASS, exit 0: existing PR/PR-BM and NRTL comparisons plus live A13 candidate (62 points, six tests) |
| All five A13 CLI examples | PASS |
| Git whitespace and A13 provenance against normalized Git-index bytes | PASS |

The six A13 solver tests include all 62 reference points, sampled envelope
comparisons, explicit ambiguity rejection, mixture/pure near-critical stops,
component invariance and failure paths. The existing ignored Linux test is
`canvas_frame_time_with_100_units`, a performance measurement rather than a
correctness gate; no A13 test is ignored. The full live check retains the old
60 homogeneous, 28 PH, and all six NRTL test groups unchanged. Hosted CI and
review remain visible on the PR; shared plans and reference-budget tables are reconciled in this PR.

### Review correction: accepted near-unit pure compositions

The original pure-temperature path searched for an exact fraction of one and
could panic for a pure fraction accepted by `with_fractions` within its 1e-10
unit-sum budget; the pressure path rejected that same accepted input instead.
Both paths now canonicalize a single active fraction to exactly one locally,
including returned phase compositions, while retaining inactive component slots
and leaving the model API unchanged. Solver/package regressions cover fractions
1 +/- 5e-11, PR and PR-BM, one-component and zero-containing vectors, both bubble
and dew calculations, invalid package inputs, and Result-based iteration-limit
failures. All checks above were rerun after the correction; independent inputs,
generators, fixtures and reference tolerances are unchanged.

### Shared-record reconciliation

The delegated documentation follow-up reconciles A11 PR #85 (pinned status at
`c54d3d2`, now merged after final CI/review), A12 merged PR #84, and A13 PR #86 in the shared
plans and reference/tolerance records. The registry now includes A13's four
quantity budgets copied unchanged from `a13-v1.json`; only critical/retrograde
extensions remain deferred there. Existing registry rows, ideal-gas budgets,
independent inputs, generators, fixtures and the Julia manifest are unchanged.
The registry check and deliberate drift rejection for existing cubic and new
A13 rows passed, as did local Markdown-link, provenance and whitespace checks.
After A11 merged, A13 was rebased on main `4887212`, retaining the
`property-tables` and `thermo-bridge` lock entries and the `property-packages`
test dependency. Combined workspace formatting, Clippy, tests and the full Julia
standard check are rerun on this tree; final-head hosted CI is observed before
delivery.

## #227 extended domain

Status: November 27, 2026. Issue #227 removes the dry-component restriction of
the cubic solver. GERG-2008 saturation keeps its own dry 230–350 K, 0.1–10 MPa
domain, and NRTL-IG keeps its modified-Raoult solver.

### Domain

| Item | Before | Now |
|---|---|---|
| Components | methane, ethane, propane, N2, CO2 | every cubic package component: also O2, water, methanol, 1-butanol, H2, CO and pseudo-components |
| Package temperature | 200–350 K, clipped to the package | the package's whole interval (at most 200–600 K) |
| Package pressure | 10 kPa–10 MPa | 1 kPa–10 MPa (`CUBIC_PACKAGE_SATURATION_PRESSURE_RANGE_PA`), so ambient water vapor pressures are reachable |
| Solver defaults (`SaturationOptions::default`, CLI) | 200–350 K, 10 kPa–10 MPa | unchanged; callers pass wider bounds explicitly |
| Three-phase bubble points | not available | returned with the second liquid |

### Method additions (extended sets only)

A set is *extended* when an active component lies outside the dry five.
Dry sets follow exactly the original path, so the `a13-v1.json` results are
unchanged.

- **Seeds.** A water-bearing dew point gets two extra incipient-liquid starts at
  each of the 13 coordinates: aqueous (0.999 water) and hydrocarbon (the
  water-free vapor with 0.001 water), mirroring the #224 trials.
- **Multiphase evidence.** Both phases of a mixture boundary must also pass the
  #224 multiphase TPD trials (Wilson, aqueous, hydrocarbon and near-pure). A root
  whose material phase is unstable, for example a water dew point below a
  hydrocarbon dew point, is rejected. The surviving root is the first liquid to
  appear. `SaturationPoint::multiphase_stability` carries the evidence.
- **Three-phase bubble point.** When no two-phase bubble point is verified, the
  VLLE flash scans 25 coordinates from the liquid side. If a two-liquid state is
  followed by a vapor-bearing one, the bracket is bisected to 1e-6 in ln T or ln p.
  Newton's method then solves the 2m+2 equations (m active components): equal
  ln fugacities of both liquids and the vapor, the liquid–liquid Rachford–Rice
  equation and Σy = 1, for ln T or ln p, ln K₂, the heavy-liquid amount and ln Kᵥ.
  The point passes the fugacity, volume-ratio, root and TPD checks for all three
  phases, plus the multiphase trials of the vapor. `liquid` is then the lighter
  liquid and `second_liquid` the heavier one, with its amount per mole of the
  specified liquid. `SaturationState::second_liquid` and the warning
  `A13-THREE-PHASE-BUBBLE` carry it to package callers. Column models
  (`CalcContext::saturation_temperature`) reject such a state explicitly, because
  their stages hold one liquid.
- **Freezing guard.** Solids are not modelled. A boundary below the water
  (273.16 K) or CO2 (216.592 K) triple point is `OutsideValidityDomain` when a
  liquid phase is at least 0.9 of that component. Hydrates are not detected.
- **Near-critical handling.** The guards are unchanged for two phases: vapor/liquid
  volume ratio ≥ 1.5, mixture max |ln K| ≥ 0.03, pure T/Tc ≤ 0.95. A three-phase
  point adds liquids at least 0.03 apart in mole fraction (plait point) and a
  singular three-phase Jacobian. All of these return `OutsideValidityDomain` with a
  near-critical diagnostic, and traces stop with `NearCritical`. The pure
  temperature bracket is now capped at the largest temperature with T/Tc ≤ 0.95.
  Before, 0.95 Tc could round above that limit and reject a valid bracket
  (found by the ethane references). No critical point is located or estimated.
- **Non-convergence.** Iteration limits (Newton, bracket bisection) return
  `NonConvergence`. Generic VLLE flash failures in the three-phase attempt do not
  replace a more specific two-phase rejection.

### Independent references

`validation/julia/generate_a13_extended.jl` (Julia 1.11.7, Clapeyron 0.6.30,
[method](../validation/clapeyron-checks.md#227-extended-bubbledew-check)) writes
`validation/fixtures/a13-extended-v1.json` from `validation/a13-extended-inputs.json`.
It has 70 cases, all within 200–600 K:

| Set | Cases |
|---|---|
| Pure water (PR bubble pressure 300–575 K every 25 K, sampled as an envelope; PR-BM 600 K to 12.5 MPa; PR-BM dew temperature 10 kPa–5 MPa), methanol, ethane, CO2, 1-butanol | 29 |
| Methanol/water 0.4/0.6 (bubble envelope 320–460 K, dew pressure, bubble and dew temperature at 1 atm and 1 MPa, PR and PR-BM), dilute methanol, 1-butanol/methanol | 17 |
| Water-bearing gases: wet methane (1 and 5 MPa), flue gas N2/O2/CO2/H2O, methane/propane/water with an aqueous first liquid (three water contents, including one where a hydrocarbon dew point lies 36 K lower) and with a hydrocarbon first liquid, H2/CO/methanol syngas, methane/methanol/water | 13 |
| Kerosene pseudo-components (bubble and dew temperature at 1 atm and 1 MPa, 445–582 K) and methane with them (bubble pressure at 450 K) | 5 |
| Three-phase bubble pressure and temperature: propane/water (PR and PR-BM) and methane/propane/water | 6 |

Interactions are explicit. They are synthetic zero, except for two values.
Water/nonpolar pairs use kij = 0.5, the #224 approximation. Methanol/water uses a
synthetic kij = −0.08: the generator's tangent-plane test showed that PR with
kij = 0 predicts a spurious methanol/water liquid split below about 330 K. None
of these values is fitted, and the references verify the implementation, not
experimental accuracy.

The budgets were declared in the generator before the first comparison: the
unchanged A13 budgets (T 2e-6 K; p 0.02 Pa + 2e-8 relative; compositions 2e-8 +
1e-7 relative; volumes 1e-12 m³/mol + 2e-8 relative). The second liquid uses the
same composition and volume budgets, and its amount 2e-8 absolute. The fixture
was committed before the engine was compared with it. In the first comparison
all two-phase cases and both three-phase pressure cases passed. One three-phase
temperature case failed with an `OutsideValidityDomain` error from the scan's
upper end: exp(ln 600 K) rounded above the 600 K EOS envelope. Moving the scan and
Newton limits 1e-12 inside the bounds fixed the engine. All 70 cases now pass.
No fixture value or budget was changed.

`crates/property-packages/tests/a13_extended_clapeyron.rs` compares T, p, both
compositions and volumes, and for three-phase points the second liquid's
composition, volume and amount. It traces the water and methanol/water
envelopes against the reference pressures. It also tests these failures:

- near-critical methanol, ethane and 1-butanol, and a methanol trace that stops
  at 490 K;
- an ethane/propane bubble trace that ends `NearCritical`;
- freezing water and CO2, and a trace-water gas below 273.16 K;
- iteration limits for a wet dew point, a methanol/water bubble point and a
  three-phase bubble point.

The pseudo-component test now checks that the solver reproduces its TP-flash
bisection within 2e-4 K. A package test covers the 200–600 K interval, water at
300 K and a propane/water three-phase bubble point. The CLI example
`examples/a13-bubble-pressure-propane-water.json` reports the heavy liquid.

### Limits

- The methanol/water (0.4/0.6) bubble search stops without a verified point at
  550 K. The reference model gives 9.9 MPa there, with a volume ratio of 4.2, so
  this is not a near-critical state. The trace ends with `UnverifiedBoundary`, and
  no value is returned. The finite 13-seed search over a 1 kPa–100 MPa interval is
  the likely cause; narrower pressure bounds are the workaround.
- The domain is implementation-verified at the listed states only. Other states
  inside the bounds are found by the same bounded search and verified by the same
  checks, but they have no independent reference.
- No retrograde branch selection, critical-point location or envelope continuation.
  Multiple verified roots in the bounds remain `UnsupportedCapability`.
- No three-phase dew point exists (an incipient liquid is one phase), and no
  four-phase or solid boundary is computed. Hydrates are not detected.
- The rigorous column stages still hold one liquid and do not check for a second
  one. The absorber still accepts only the dry set.
