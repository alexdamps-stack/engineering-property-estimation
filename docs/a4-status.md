# A4 implementation and verification record — stability and TP flash

Date: 2026-11-07. Starting revision: `538d9f0` (A1–A3 merged). Scope: PR and
PR-BM ([A1–A3](a1-a3-status.md)), at most two phases, inside the existing
200–600 K, 0.01 Pa–100 MPa computational envelope. Check outcomes are recorded in
[the A4 verification report](../validation/a4-checks.md).

Update 2026-11-08: issue #30 adds safeguarded Newton acceleration as described
below. The [A3/A4 acceptance review](a3-a4-acceptance.md) records the milestone
scope and current evidence, superseding the historical pending review below.

## Implemented

All calculation code is native Rust in `thermo-solvers`; there is no external
calculation API, subprocess or Python dependency.

- **Model access (`thermo-models`).** `CubicModel::with_fractions` evaluates the
  same model, data and complete `kij` matrix at another composition. Zero
  entries remain present at infinite dilution. `CubicModel::ln_fugacity`
  returns Z and ln φ of the requested root without caloric evaluation. It
  rejects a selected root with nonnegative (∂p/∂v)_T.
- **Rachford–Rice.** `rachford_rice` classifies `LiquidOnly` (g(0) ≤ 0) and
  `VaporOnly` (g(1) ≥ 0), otherwise it solves 0 < β < 1.
  `rachford_rice_negative` solves between the poles 1/(1−Kmax) < β < 1/(1−Kmin)
  with safeguarded Newton/bisection. Inactive (z = 0) entries stay zero, and
  K on one side of one is an explicit `UnsupportedCapability`. The default
  iteration limit is 200.
- **Stability.** `assess_stability` implements Michelsen modified TPD,
  tm = 1 + Σ W_i (ln W_i + ln φ_i(w) − ln z_i − ln φ_i(z) − 1).
  - Trials: Wilson vapor-like (z·K) and liquid-like (z/K), plus one
    near-pure trial per active component (0.999 of that component), which can
    be disabled.
  - Root selection: successive substitution uses the minimum-Gibbs root at each
    trial composition.
  - Classification: the feed is unstable if a nontrivial trial reaches
    tm < −1e-8.
  - Iteration limit: 2000 per trial. A trial that hits the limit after showing
    tm < −1e-8 still proves instability. Otherwise it is a `NonConvergence`
    error (stability undetermined).
  - A single active component has only the trivial stationary point.
  - **Passing the finite trials is evidence, not a proof of global stability.**
- **TP flash.** `tp_flash` runs these steps:
  1. Runs the feed stability test. A stable feed returns one phase from the
     minimum-Gibbs root.
  2. Takes the initial K from the most negative trial, oriented by its volume
     relative to the feed.
  3. Runs successive substitution on ln K with the negative-flash
     Rachford–Rice inner solve. From iteration 3, acceleration tries Newton on
     the active log-fugacity residual. Central differences of the full residual
     include the implicit phase-fraction dependence; no analytical composition
     derivative API is assumed. A small dense partial-pivoting solve avoids
     forming an inverse. Each ln K step is limited to one, and up to 12
     backtracks require sufficient decrease in the maximum fugacity residual.
     Newton is attempted only for a physical phase fraction, and its trial
     fractions must remain physical. Trials must not increase split Gibbs
     energy beyond 1e-14 in dimensionless per-mole units; this prevents a
     decreasing fugacity residual from steering toward the trivial solution.
     Singular Jacobians, failed model/RR probes, trivial candidates and failed
     line searches fall back to SS. Every fifth fallback step can apply the
     existing dominant-eigenvalue extrapolation (0 < λ < 0.95).
  4. Converges on max_i |ln(y_i φ_i^V) − ln(x_i φ_i^L)| ≤ 1e-10, with a limit
     of 2000 iterations.
  5. Rejects: approach to the trivial solution (all |ln K| < 1e-6), a
     converged β outside (0, 1), or a split that does not lower Gibbs energy.
  6. Runs post-split stability on both phases. An indicated third phase
     returns `UnsupportedCapability`.

  Phases carry full `HomogeneousState` properties (shared caloric reference)
  and a molar phase fraction. The result also reports iterations, fugacity
  residual, component-balance residual, the Δg/RT of the split and every trial.
  `A4-NEWTON-STEPS` reports the accepted Newton step count for a split.
  The CLI's `split_iterations` includes both Newton and SS; the legacy
  `successive_substitution_iterations` key is retained as an alias.
  `FlashOptions::newton=false` retains eigenvalue-only acceleration;
  `acceleration=false` disables both accelerators. Defaults enable Newton.
- **Labels.** These are labels, not phase-boundary results:
  - A split labels the larger molar volume `Vapor` and the other `Liquid`. A
    warning is raised when the molar volumes differ by < 0.1 %.
  - A single phase is `Supercritical` above Kay's pseudo-critical T and p.
    Otherwise it is `Vapor` when v/b ≥ Zc/Bc of the exact-factor PR
    (0.3074013087/OMEGA_B ≈ 3.9514), else `Liquid`.
- **Pure saturation.** `pure_saturation_pressure` uses a bracketed Newton
  solve in ln p between the spinodal pressures (scan plus bisection of
  (∂p/∂v)_T = 0) on ln φ^L − ln φ^V, which is strictly decreasing in p.
  Declared domain: T/Tc ≤ 0.95, 200–600 K, p_sat ≥ 0.01 Pa. Outside it, the
  function returns `OutsideValidityDomain`. The iteration limit is 200.
- **CLI.** `"calculation": "tp-flash"` runs the flash; the `root` field is then
  rejected. `"homogeneous"` (default) keeps the A3 behavior. The output covers
  phases, labels, fractions, compositions, all properties, convergence
  evidence, completed TPD trials and the diagnostics. Since the #88 fix,
  the trial search can stop once a nontrivial negative TPD proves instability;
  inconclusive trials cannot establish stability. See
  `examples/tp-flash-*.json`.

The A1 residual targets are met with tighter internal values:

| Quantity | A1 target | Used here |
|---|---|---|
| Active log-fugacity residual | 1e-7 | 1e-10 |
| Normalized component balance | 1e-8 | Observed ≤ 1.2e-16 absolute |

Component balance is reported as max_i |β y_i + (1−β) x_i − z_i| in
mole-fraction units.

## Independent references

`validation/generate-flash-references.py` is a new adapter (no Rust calls). It
uses teqp 0.23.2 AD fugacities and NumPy roots with SciPy pressure refinement,
the same formulation as the A3 generator. It solves equilibrium with
**different algorithms** from the Rust code:

- **Stability:** global multi-start L-BFGS-B minimization of tm over ln W.
  Starts: 2 Wilson starts plus 48 seeded (0xA4A4) Dirichlet starts.
- **Split:** MINPACK hybrid Powell on ln K with a brentq Rachford–Rice inner
  solve. It asserts fugacity residual < 1e-11, β in (0, 1), Δg < 0 and
  stability of both phases.
- **Pure saturation:** teqp's cubic superancillary, refined with teqp's
  `pure_VLE_T`.

The fixture `validation/fixtures/flash-v1.json` holds 22 TP states and 30
saturation states:

- **TP flash, 11 states per variant:**
  - Methane/propane two-phase, single vapor and single liquid.
  - A zero-fraction nitrogen trace case.
  - Methane/CO2 with synthetic kij = 0.0919.
  - Natural-gas feed A (single phase), the mixed demonstration feed at 240 K /
    5 MPa and 220 K / 1.5 MPa, and feed B at two states.
- **Pure saturation:** ethane, propane, CO2, water and methanol at Tr 0.7/0.8/0.9.
  BM equals PR below Tc, but both are run.

Methane, nitrogen and oxygen saturation lies below 200 K and is tested as an
explicit domain error.

Budgets were fixed in the generator before Rust comparison:

| Quantity | Absolute | Relative |
|---|---|---|
| β | 2e-8 | — |
| Composition | 2e-8 | 1e-7 |
| Phase Z | 2e-9 | 2e-8 |
| Saturation p | 1e-6 Pa | 2e-9 |
| Saturation Z | 2e-10 | 2e-9 |

Observed maximum deviations were:

| Quantity | Maximum deviation |
|---|---|
| β | 2.5e-11 |
| Compositions | 3.5e-11 |
| Phase Z | 3.3e-12 |
| Saturation p | 2.0e-14 relative |

The reference adapter ran in a local venv: Python 3.14.4, teqp 0.23.2,
NumPy 2.5.4, SciPy 1.18.1. These versions are recorded in the fixture and
differ from the A3 pins for NumPy and SciPy. Agreement verifies the
implementation against the shared parameters, not measurements.

Additional Rust tests cover:

- The Rachford–Rice closed forms and all-K>1, all-K<1 and K≈1 cases, negative
  flash, monotone residual, invalid input and forced exhaustion.
- Fugacity equality re-derived from the reported phase properties.
- Permutation and inactive-component invariance.
- A 41-point isotherm (0.2–8.2 MPa), whose phase count runs 1 → 2 → 1 with a
  monotone vapor fraction.
- Single-phase labels and pure feeds.
- Forced iteration limits for the flash, stability and saturation.
- Out-of-envelope, near-critical and non-pure saturation errors.
- The CLI schema and its failures.

A manual probe of methane/propane 75/25 at 250 K converged up to β = 0.007 at
9.025 MPa (~100 iterations) and returned one phase from 9.1 MPa. It produced
no NaN or silent failure. This is not a near-critical accuracy claim.

## Reproduce

```sh
cargo test -p thermo-solvers --locked
cargo run -p thermo-cli --locked -- examples/tp-flash-binary.json
cargo run -p thermo-cli --locked -- examples/tp-flash-natural-gas.json
/path/to/reference-venv/bin/python -I validation/generate-flash-references.py --output /tmp/flash-candidate.json
```

Do not regenerate fixtures to make a failing test pass. Replacing expected
values requires a reviewed change.

## Acceptance boundaries and open items

Delivered A4 evidence:

- Phase classification and stability assessment, tested separately from the
  split.
- Component balances and fugacity equality.
- Boundary cases (appearance and disappearance, trace and zero components).
- Explicit nonconvergence for each iterative solver.

Not claimed:

- Proof of global stability.
- More than two phases, including VLLE (water/methanol hydrocarbon systems
  remain excluded).
- Near-critical accuracy.
- Bubble/dew points and phase envelopes (A13).
- PH flash (A5).
- Fitted kij or experimental accuracy.

Newton convergence is verified on eighteen difficult methane/propane states at
250 K and 8.5–9.3 MPa: 6–20 split iterations versus 63–599 for the previous
eigenvalue-only algorithm, with matching phase fractions and compositions.
This measures iteration reduction, not wall-clock speed or accuracy at the
critical point. Finite-difference Newton requires extra EOS evaluations.

The [A2 process qualification](a2-process-qualification.md) and A5 PH checks
now cover the cooler and valve outlet reference states. A13 phase boundaries
and envelopes and A8 connected flowsheet execution retain their own gates.
