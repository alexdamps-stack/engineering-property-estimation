# Testing Strategy

This is the testing guide for the thermodynamic engine and the connected flowsheet simulator. PR/PR-BM models and TP/PH solvers now have executable numerical comparisons; requirements for future capabilities remain acceptance targets. Evidence and supported domains are recorded in the milestone status documents.

### Mandatory Clapeyron standard check

For the currently supported PR/PR-BM paths, the independent Julia comparison
is a standard correctness check alongside formatting, Clippy and workspace
tests. `cargo test --workspace --locked` always runs the committed offline
Clapeyron references. The live check independently evaluates the same inputs
using locked Julia 1.11.7 / Clapeyron.jl 0.6.30 and compares them with the
current Rust engine:

```sh
julia --project=validation/julia validation/julia/check.jl
```

The live check runs on every push and pull request, plus the scheduled/manual
numerical workflow. Any generation failure, changed fixture provenance,
nonfinite result or out-of-budget comparison fails the job; cases are not
silently skipped. The check generates temporary candidates and never updates
committed expected values. Julia is a validation dependency only.

Changes to cubic models, caloric data, equilibrium solvers or their numerical
contracts must preserve these checks. Adding a supported formulation/domain
requires matched independent cases and explicit budgets before acceptance.
Package/default-data differences must be qualified separately; the standard
check supplies matched parameters, ideal contributions and caloric references.
See [Clapeyron verification](../validation/clapeyron-checks.md) for scope,
reproduction, tolerances and local execution evidence.

## 1. What each kind of evidence establishes

- **Implementation verification:** the code implements the selected formulation correctly.
- **Model comparison:** different formulations predict similar or different results under matched conditions.
- **Experimental validation:** a model predicts measured behavior with sufficient accuracy for its intended use.

Keep these categories separate in tests, reports, and release claims. Passing compilation, formatting, linting, or an empty test harness establishes no numerical verification. Agreement between two implementations does not establish experimental accuracy.

## 2. Test layers and placement

| Layer | Required evidence | Placement |
|---|---|---|
| Unit tests | Small calculations, conversions, parameter conventions, and input/error handling | Rust test modules beside the implementation |
| Public behavior and integration | Observable crate behavior and interactions between calculation layers | The owning crate's `tests/` directory |
| Independent reference comparisons | Matched formulations and inputs with qualified expected values | Reference data and generation/comparison tooling under `validation/`; small fixtures consumed by crate integration tests |
| Thermodynamic consistency | Identities, derivatives, limiting behavior, and invariance | Model and solver tests, with broader sweeps under `validation/` |
| Solver robustness | Correct phases, equilibrium, conservation, boundary handling, and explicit failure | Solver and property-package integration tests |
| Flowsheet and application integration | Unit-operation balances, coupled calculations, persistence, and stale-result handling | Flowsheet and application integration tests as those capabilities appear |
| Experimental validation | Deviations from measurements with source uncertainty and a defined domain | Separate qualified datasets and reports under `validation/` |
| Performance | Runtime, iteration counts, and memory for fixed workloads | Separate benchmarks introduced once real calculation paths exist |

Add tests with the capability they verify. Do not add numerical placeholders that pass without evaluating a real implementation. Demonstration inputs under `examples/` are not independent expected results unless separately qualified.

## 3. Reference cases and independence

Each reference case must record:

- A stable case identifier and purpose: implementation verification, model comparison, or experimental validation.
- Component identities, composition values and basis, temperature, absolute pressure, units, phase specification where applicable, and all other prescribed inputs.
- Exact model variant, alpha function, mixing rules, component data and interaction parameters, ideal contributions, and caloric reference states. Record any deliberately chosen zero-parameter approximation.
- Reference source, software version or commit where applicable, data version, generation method/command, usage terms, and known limitations. Experimental data also need reported uncertainty where available.
- Expected output quantities, expected phase behavior or failure, applicable domain, and per-quantity absolute and relative tolerances with a rationale.

Keep small qualified reference fixtures in Git so routine Rust tests can run without network access, commercial software, or external reference installations. Preserve sufficient output precision. Generation tools may use external software, including non-Rust tools; that does not add a runtime dependency to the Rust engine.

Expected values must come from an independent implementation, independently evaluated published formulation, or qualified measurements. A second wrapper around our own engine is not an independent reference. Check whether candidate references share code or parameter sources and document the remaining independence limits. Consistency and round-trip tests complement independent comparisons; they cannot replace them because shared errors may cancel.

Use NIST's [AGA8/GERG sample code and datasets](https://pages.nist.gov/AGA8/) as candidates for GERG implementation verification. Qualify the specific files and domain before adoption. For PR-BM and NRTL, select references that match the implemented variant and parameter conventions; ordinary PR is not a PR-BM verification reference.

Regenerate expected values only through an explicit reviewed change recording the reason, source/version changes, and numerical differences. Never update expected values automatically just to make a failing test pass. If no suitable independent reference is available, record that gap and limit the capability's verification claim.

## 4. Required numerical scenarios

### Types, data, and single-phase models — A1–A3

- Verify unit conversions, mass/mole composition conversion, component identification, and parameter ordering. Test the specified handling of invalid inputs, missing parameters, and nonfinite values.
- Compare PR and PR-BM separately against matched independent references for roots, density, compressibility, fugacity coefficients, and caloric properties in the supported domain.
- Check the BM alpha formulation and its derivatives, including behavior around any piecewise transition in the selected formulation.
- Check `h = u + p*v` with consistent molar or mass units and reference states. Check `cp = (∂h/∂T)p` in smooth single-phase regions.
- Compare analytical or automatic derivatives with independently evaluated numerical derivatives using multiple suitable step sizes; use higher precision for difficult diagnostic cases when available. Phase boundaries and singular regions need separate treatment.
- Verify that permuting components and all associated parameters leaves physical results unchanged; compare composition-dependent outputs after undoing the permutation.
- Verify pure-fluid limits and appropriate low-density limits within supported model domains. Test trace and zero components according to the calculation contract.
- For cubic models, verify that each pure component reproduces its own critical temperature and pressure through the vanishing first and second pressure derivatives at the critical point. Where the cubic has three real roots, verify that the selected root minimizes the Gibbs energy, and that phase-specified requests return the requested root or an explicit error.
- For the BM alpha function, verify continuity of the value and its first temperature derivative at the transition at `T = Tc`.

The following identities must hold for each implemented equation of state to within a tolerance justified by floating-point conditioning, not by reference quality. Evaluate them at the same grid of single-phase states used for reference comparisons.

| Identity | Checks |
|---|---|
| `Σ n_i ∂ln φ_i/∂n_j = 0` at constant T and p | Gibbs–Duhem consistency of composition derivatives |
| `∂ln φ_i/∂n_j = ∂ln φ_j/∂n_i` | Symmetry of the fugacity Jacobian used by flash solvers |
| `Σ x_i ln φ_i = ln φ_mix` | Partial molar quantities sum to the mixture value |
| `A_res(λn, λV, T) = λ A_res(n, V, T)` | Euler homogeneity of the residual Helmholtz energy |
| `ln φ_i → 0` and `Z → 1` as density → 0 | Ideal-gas limit |
| `cp − cv = −T (∂p/∂T)_V² / (∂p/∂V)_T` | Consistency of caloric derivatives |
| `(∂s/∂p)_T = −(∂v/∂T)_p` | Maxwell relation linking volumetric and entropic parts |
| Speed of sound from `cp/cv` and `(∂p/∂ρ)_T` | Cross-check of second derivatives where speed of sound is reported |

Use property-based tests (with the `proptest` crate, added as a development dependency when the first model lands) for component permutation, uniform scaling of amounts, the identities above, and composition normalization. Record the proptest seed for each failure and keep any failing case as a fixed regression fixture.

Derivative verification is particularly important for caloric properties and solver Jacobians; the [NIST teqp paper](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=933725) illustrates how independent derivative checks can expose precision errors.

### Stability, TP flash, and PH flash — A4–A5

- Test Rachford–Rice on its own with prescribed K-values: all `K_i > 1` and all `K_i < 1` (no root in `[0, 1]`, single-phase result), `K_i` close to 1, negative-flash ranges where the solver supports them, and monotonicity of the residual within the bracketing interval.
- Check pure-component saturation pressures against independent references through fugacity equality.
- Cover known single-phase and two-phase states with independently established phase behavior. Check normalized phase compositions, bounded phase fractions, component balances, and fugacity equality for active components.
- Test stability assessment separately from solving the phase split. A converged root or small residual alone does not establish a stable equilibrium; finite TPD trials do not prove global stability.
- Cover phase appearance/disappearance, near-critical states, nearly identical phases, trace components, and sensitive initial guesses. Classify cases as supported, diagnostic, or outside the declared domain before execution.
- Test TP→PH recovery only where the state is uniquely recoverable under the specified constraints, checking temperature, phases, compositions, and enthalpy as appropriate. Include single-phase and two-phase states.
- Test PH flash at bubble and dew points (vapor fraction exactly 0 and 1) and for a pure component and a narrow-boiling mixture inside the two-phase region, where temperature alone does not determine enthalpy.
- Exercise invalid inputs, unsupported domains, missing data, and forced iteration exhaustion. Return explicit diagnostics rather than successful values containing NaNs or stale results. Each iterative solver must have a documented iteration limit, and at least one test must reach it deliberately.

### Saturation points and phase envelope — A13

- Compare bubble- and dew-point pressure and temperature with independent references for the demonstration mixtures, including the hydrocarbon dew point of the natural-gas case.
- Verify that a saturation result satisfies fugacity equality with an incipient-phase fraction of zero.
- Verify that the phase-envelope trace reports its covered range and returns an explicit diagnostic near the critical point rather than unverified values.

### GERG, NRTL, and supported domains — A6–A7

- Qualify reference formulations, component inventories, ideal/residual contributions, and caloric reference states before comparison.
- For GERG, cover the declared component and state ranges with independent property and derivative comparisons. Do not infer flash validation from single-phase density agreement. The acceptance set is a fixed, versioned list of reference states; every state inside the declared domain must pass, and states outside it must return an explicit out-of-domain outcome rather than being skipped.
- For NRTL, test directed parameters, component permutations, temperature dependence, activity coefficients, vapor-model coupling, and caloric consistency.
- Maintain a supported-domain statement for each model/property/solver combination. Keep implementation error separate from differences against other models or experimental measurements.

### Connected flowsheet simulator and interface — A8–A9

- Mixer and separator cases must close component balances; check energy balances under the specified heat/work assumptions. A separator must use the qualified flash path.
- Heater/cooler duty must match the stream enthalpy change; an adiabatic valve with no shaft work and negligible kinetic/potential changes must preserve specific enthalpy while using PH flash for the outlet state.
- Scaling all prescribed feed flow rates and extensive specifications consistently must scale extensive outputs while preserving intensive states for applicable unit-operation specifications.
- Run the complete mix/cool/expand/separate demonstration without the GUI before accepting its graphical workflow.
- Save/load must preserve inputs, model/data selections, and results after recalculation within defined tolerances. Editing inputs must invalidate affected results; failures must remain visible and must not leave old values marked current.

## 5. Tolerances, reproducibility, and failure reports

For scalar comparisons, use `abs(actual - expected) <= atol + rtol * abs(expected)`, with documented quantity-specific `atol` and `rtol`. Absolute tolerances carry the quantity's units. Handle expected errors and nonfinite outputs explicitly rather than allowing them through numerical comparisons.

Retain the initial targets in the [development plan](development-plan.md): composition-sum error at most `1e-10`, normalized component-balance error at most `1e-8`, active-component log-fugacity residual at most `1e-7`, and scaled energy-balance error at most `1e-6` for regular, well-conditioned cases. Define balance normalization, scaling floors, and active-component rules with the relevant calculation contracts before testing. These residual targets are distinct from property accuracy against references or measurements.

Set property tolerances before implementation using reference quality and conditioning. Record them in the table below as part of A2; a model milestone cannot be accepted while its rows are still open. Tolerances may differ per model and reference; add rows rather than widening a shared one.

| Quantity | Model / reference | `atol` (unit) | `rtol` | State domain | Rationale |
|---|---|---|---|---|---|
| Mass density ρ | PR-exact and PR-exact-BM / teqp 0.23.2 | 1e-7 kg/m³ | 2e-9 | 30 homogeneous states/model, including 18 regular grid states | Double precision, matched constants and regular roots |
| Compressibility factor Z | same separate model fixtures | 2e-10 | 2e-9 | same, every admissible root | Independent NumPy roots plus teqp pressure refinement |
| ln φ_i | same | 2e-9 | 2e-9 | mechanically stable roots in same fixture | Independent residual Helmholtz derivatives |
| Residual and total enthalpy | same + Cantera 3.2 NASA7 for total | 2e-5 J/mol | 2e-9 | same | Matched R and ideal reference offsets |
| Residual and total entropy | same + Cantera 3.2 NASA7 for total | 2e-7 J/(mol K) | 2e-9 | same | Same-T,p ideal departures and mixing entropy |
| cp and cv | same + Cantera 3.2 NASA7 | 2e-5 J/(mol K) | 2e-8 | smooth stable roots; BM Tc excluded | Second-derivative and root conditioning budget |
| Speed of sound w | same | 2e-4 m/s | 2e-8 | same | Independent residual derivatives and cp/cv |
| Ideal gas h | NASA7 / Cantera 3.2.0 | 2e-6 J/mol | 1e-11 | 64 ideal states, 200..6000 K including fit join | Independent integration and reference offsets |
| Ideal gas s | same | 2e-8 J/(mol K) | 1e-11 | same | Independent cp/T integration |
| Ideal gas cp | same | 2e-8 J/(mol K) | 1e-11 | same | Independent retained polynomial evaluation |
| Pure saturation pressure | PR, PR-BM / teqp 0.23.2 | 1e-6 Pa | 2e-9 | 30 pure states, `flash-v1.json` | Independent superancillary + VLE refinement; regular coexistence |
| Saturation Z | same | 2e-10 | 2e-9 | same | Refined phase densities |
| TP vapor fraction | PR, PR-BM / teqp + SciPy | 2e-8 | 0 | 22 TP states, `flash-v1.json` | Independent Powell split and finite-start TPD |
| TP phase composition | same | 2e-8 | 1e-7 | same | Active-component equilibrium and trace budget |
| TP phase Z | same | 2e-9 | 2e-8 | same | Split/root error propagation |
| PH / process temperature | PR, PR-BM / teqp, Cantera + SciPy | 2e-6 K | 0 | 28 PH + 14 process states | Independent energy root/coexistence solve |
| PH / process h | same | 2e-5 J/mol | 2e-9 | same | Same ideal offsets and AD departures |
| PH / process fraction and composition | same | 2e-8 | 0 / 1e-7 | same | Fraction / composition budgets respectively |
| Bubble/dew temperature | PR, PR-BM / native Clapeyron 0.6.30 | 2e-6 K | 0 | 32 p/T cases + 30 envelope pressure samples, `a13-v1.json` | Bounded dry demo mixtures and propane; [A13 status](a13-status.md) |
| Bubble/dew pressure / sampled envelope pressure | same | 0.02 Pa | 2e-8 | same | Explicit branch bounds; verified prefix only, near-critical diagnostic stop |
| Saturation phase composition | same | 2e-8 | 1e-7 | same | Both incipient/material compositions; unchanged pre-comparison budget |
| Saturation phase molar volume | same | 1e-12 m^3/mol | 2e-8 | same | Both phase volumes; guards exclude coalescing states |
| Critical-point / retrograde envelope continuation | PR, PR-BM | deferred | deferred | A13 follow-up | Outside bounded qualification; no extrapolated values or accuracy claim |
| Homogeneous properties and derivatives | GERG-2008 / NIST AGA8 C++ 2.01 | 2e-8 in each reported SI unit | 2e-10 | 273 derivative states / 269 admissible property states, `gerg-v1.json` | Matched native ideal reference; fixed numerical comparison budget |
| Molar density | GERG / NIST + SciPy | 2e-7 mol/m³ | 2e-9 | 42 TP roots, `gerg-density-v1.json` | Independent bracket/refinement |
| TP pressure | GERG / NIST + SciPy | 1e-5 Pa | 2e-9 | 42 TP roots | Independent pressure refinement |
| TP Z | same | 2e-9 | 2e-9 | same | Root propagation |
| TP h | same | 1e-5 J/mol | 2e-9 | same | Matched native ideal reference |
| TP s | same | 2e-7 J/(mol K) | 2e-9 | same | Matched native ideal reference |
| TP cp, cv | same | 2e-6 J/(mol K) | 2e-8 | same | Independent derivative and root conditioning |
| TP sound | same | 2e-4 m/s | 2e-8 | same | Second derivative conditioning |
| TP dp/dρ, dp/dT | same | 2e-5 Pa m³/mol, Pa/K respectively | 2e-8 | same | Pressure derivatives at refined roots |
| ln φ | GERG / teqp AD | 2e-8 | 2e-9 | 45 states | Independent chemical-potential derivatives |
| Equilibrium phase composition | GERG / NIST + teqp + SciPy | 2e-8 | 1e-7 | 13 TP + 15 PH cases | Independent finite-start stability and split |
| Equilibrium vapor fraction | same | 2e-8 | 0 | same | Lever-rule and split conditioning |
| PH temperature | same | 2e-6 K | 0 | same | Independent energy solve |
| Equilibrium molar density | same | 2e-5 mol/m³ | 2e-8 | same | Phase root/split propagation |
| Equilibrium Z | same | 2e-9 | 2e-8 | same | Phase root/split propagation |
| Equilibrium h | same | 2e-5 J/mol | 2e-9 | same | Matched caloric reference |
| Equilibrium s | same | 2e-7 J/(mol K) | 2e-9 | same | Matched caloric reference |
| Equilibrium cp | same | 2e-5 J/(mol K) | 2e-8 | same | Second derivative conditioning |
| Molar volume | PR, PR-BM / Clapeyron 0.6.30 | 1e-12 m³/mol | 2e-9 | 60 homogeneous + phase states of 28 PH cases | Independent cubic roots, BigFloat at dilute bound |

The approved machine-readable [registry](../validation/tolerances.json) preserves
all exact cubic, TP, PH, process, Clapeyron, bounded A13 and GERG budgets from their generators
and fixtures, including GERG TP root, fugacity and equilibrium quantities.
`python validation/consolidate-tolerances.py` rejects registry/fixture drift.
The original proposals are retired to `tolerances.preparation-history.json`
for the historical preparation checker; they are not acceptance budgets.
Zero relative tolerances are intentional absolute-only budgets. Units follow
the fixture quantity name and SI output contract. Separate model/reference
rows are never combined by selecting the loosest tolerance.

Each row also states the minimum number of reference states and how they cover the domain, for example a temperature–pressure grid per mixture. Critical regions and trace components require separate justified criteria; do not silently relax tolerances after failures. Experimental comparisons must state measurement uncertainty and application accuracy requirements separately from numerical agreement limits.

The homogeneous fixture, exact state list, preselected budgets, reference
independence and computational/physical validity distinction are documented in
[A1–A3 status](a1-a3-status.md). These rows qualify implementation agreement for
the listed fixed fixtures only; process qualification is recorded in [A2 process qualification](a2-process-qualification.md).
The manual/monthly numerical-validation workflow runs the fixed independent
fixtures plus 2048 fixed-seed permutation/scaling/identity cases per platform in
release mode. It never regenerates expected values.

Generated cases must record a fixed random seed, generator version, domain, and sampling method. Preserve failing generated inputs as small regression fixtures. Reports identify the code commit, toolchain/platform, case, model/data versions, expected and actual results, residuals, tolerances, phase outcome, and solver status/iterations where applicable. Report failures and unsupported states explicitly; exclude neither silently from success rates.

A11's [fixed-mixture experiment](a11-status.md) in PR #85 reports interpolation
error against the already qualified direct Rust model, rather than an independent
engine reference. Its separate value/derivative targets, timing evidence, C^1
cell-edge discontinuities in second derivatives and finite-lattice qualification
remain scoped to exact PR, the fixed natural gas and 280-420 K / 1-80 bar.
A12's [DWSIM integration gates](a12-status.md) are cross-model transport/mapping
checks, not replacements for the tighter engine reference budgets above.

## 6. CI and development gates

The existing Windows/Linux workflow runs formatting, Clippy, workspace tests, and documentation builds; documentation warnings fail the build, and build artifacts are cached between runs. Numerical unit and integration tests will run through its existing `cargo test --workspace --locked` step as they are implemented. Larger numerical sweeps belong in the validation workflow and should run with `--release`.

Results may differ in the last bits between Linux and Windows, because of differing math-library implementations. Tests must therefore compare with tolerances, never with exact floating-point equality, and the same tolerance must pass on both platforms.

Once executable models exist, introduce a separate scheduled/manual validation workflow with reproducible state-space sweeps and reference reports. Keep a small representative subset in ordinary pull-request tests. Define the sweep workload and supported domains with the model implementation; do not schedule an empty validation job now.

Keep performance runs separate from correctness gates. Before A11 starts, set numeric targets for both the required speedup and the maximum interpolation error of values and derivatives over the approved domain. Record hardware/runner conditions, workloads, repetitions, toolchain, runtime distributions, and iteration counts; do not gate pull requests on noisy shared-runner timings. Numerical correctness remains a prerequisite for any speedup claim.

A capability is complete only when its intended behavior, domain, and failure outcomes are documented; meaningful unit/integration tests and applicable independent reference checks pass; evidence and limitations are recorded; and the relevant milestone's acceptance criteria are met. Each fixed bug needs a regression case when its observable behavior can be reproduced. Coverage percentages and test counts are supporting signals, not acceptance evidence by themselves.
