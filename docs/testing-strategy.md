# Testing Strategy

This is the testing guide for the thermodynamic engine and the connected flowsheet simulator. It defines ongoing requirements; the [calculation contract](calculation-contract.md) and [reference evidence](../validation/README.md) describe the current homogeneous implementation and executed checks. Phase equilibrium and flowsheet test requirements below remain future work.

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

Derivative verification is particularly important for caloric properties and solver Jacobians; the [NIST teqp paper](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=933725) illustrates how independent derivative checks can expose precision errors.

### Stability, TP flash, and PH flash — A4–A5

- Cover known single-phase and two-phase states with independently established phase behavior. Check normalized phase compositions, bounded phase fractions, component balances, and fugacity equality for active components.
- Test stability assessment separately from solving the phase split. A converged root or small residual alone does not establish a stable equilibrium; finite TPD trials do not prove global stability.
- Cover phase appearance/disappearance, near-critical states, nearly identical phases, trace components, and sensitive initial guesses. Classify cases as supported, diagnostic, or outside the declared domain before execution.
- Test TP→PH recovery only where the state is uniquely recoverable under the specified constraints, checking temperature, phases, compositions, and enthalpy as appropriate.
- Exercise invalid inputs, unsupported domains, missing data, and forced iteration exhaustion. Return explicit diagnostics rather than successful values containing NaNs or stale results.

### GERG, NRTL, and supported domains — A6–A7

- Qualify reference formulations, component inventories, ideal/residual contributions, and caloric reference states before comparison.
- For GERG, cover the declared component and state ranges with independent property and derivative comparisons. Do not infer flash validation from single-phase density agreement.
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

Set property tolerances before implementation using reference quality and conditioning. Critical regions and trace components require separate justified criteria; do not silently relax tolerances after failures. Experimental comparisons must state measurement uncertainty and application accuracy requirements separately from numerical agreement limits.

Generated cases must record a fixed random seed, generator version, domain, and sampling method. Preserve failing generated inputs as small regression fixtures. Reports identify the code commit, toolchain/platform, case, model/data versions, expected and actual results, residuals, tolerances, phase outcome, and solver status/iterations where applicable. Report failures and unsupported states explicitly; exclude neither silently from success rates.

## 6. CI and development gates

The existing Windows/Linux workflow runs formatting, Clippy, numerical unit/integration tests, and documentation builds. Its `cargo test --workspace --locked` step includes the small independent fixture comparisons as implementations are added. Passing it establishes only the evidence actually exercised by those tests, not general numerical or experimental accuracy.

A separate scheduled/manual workflow now runs the implemented deterministic 200-case homogeneous robustness grid and reports failures explicitly. The workload is recorded in its test source and the calculation contract; it does not constitute an independent accuracy sweep. Keep independent reference comparisons in ordinary pull-request tests and add broader reference sweeps as formulations and supported domains expand.

Keep performance runs separate from correctness gates. Record hardware/runner conditions, workloads, repetitions, toolchain, runtime distributions, and iteration counts; do not gate pull requests on noisy shared-runner timings. Numerical correctness remains a prerequisite for any speedup claim.

A capability is complete only when its intended behavior, domain, and failure outcomes are documented; meaningful unit/integration tests and applicable independent reference checks pass; evidence and limitations are recorded; and the relevant milestone's acceptance criteria are met. Each fixed bug needs a regression case when its observable behavior can be reproduced. Coverage percentages and test counts are supporting signals, not acceptance evidence by themselves.
