# A13 saturation and phase-envelope acceptance

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
Zero fractions remain zero in both phases. Thirteen Wilson-initialized seeds
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
The 200?225 K natural-gas dew segment explicitly limits pressure to 4 MPa:
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
`saturation_temperature_bounds_k` and `saturation_pressure_bounds_pa` are
ordered two-element arrays. Envelope inputs also accept `envelope_step_k` and
`envelope_maximum_points`; omit fixed T/p. A failed boundary returns the CLI's
existing nonzero error outcome; a trace returns its verified prefix and
explicit stop, so callers must inspect the stop before using its range.

The standard Julia check now validates A13 provenance, generates a temporary
A13 candidate and runs the same Rust tests using `THERMO_A13_REFERENCE`, in
addition to every existing PR/PR-BM and NRTL check. Reference regeneration is
an explicit review operation using `generate_a13.jl --output candidate.json`;
CI never overwrites committed expected values.

## Executed checks

Executed on 2026-11-09 with Rust 1.99.1, latest main baseline `728a3bc`,
Windows native scoped tests and WSL Linux full-workspace/reference checks.

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | PASS |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | PASS |
| `cargo test --workspace --locked` (Linux) | PASS: 219 passed, 0 failed, 1 existing ignored performance measurement |
| `cargo test -p thermo-solvers -p property-packages -p thermo-cli --locked` (Windows) | PASS: 68 passed, 0 failed |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | PASS |
| `python validation/consolidate-tolerances.py` | PASS for all registered prior fixtures; shared A13 table rows remain coordinator-owned, with the new A13 budgets recorded above and in the fixture |
| Julia 1.11.7: `julia --project=validation/julia validation/julia/check.jl` | PASS, exit 0: existing PR/PR-BM and NRTL comparisons plus live A13 candidate (62 points, five tests) |
| All five A13 CLI examples | PASS |
| Git whitespace and A13 provenance against normalized Git-index bytes | PASS |

The five A13 solver tests include all 62 reference points, sampled envelope
comparisons, explicit ambiguity rejection, mixture/pure near-critical stops,
component invariance and failure paths. The existing ignored Linux test is
`canvas_frame_time_with_100_units`, a performance measurement rather than a
correctness gate; no A13 test is ignored. The full live check retains the old
60 homogeneous, 28 PH, and all six NRTL test groups unchanged. Hosted CI and
review remain visible on the PR; shared-plan reconciliation is coordinator-owned.
