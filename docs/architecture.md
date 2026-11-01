# Architecture

The repository has six library crates and two binary crates. The first implementation connects `thermo-cli` → `property-packages` → `thermo-models` / `thermo-data` → `thermo-types`. Shared serde dependencies provide JSON input/output; reference-generation Python packages are optional development tools, not Rust runtime dependencies. `thermo-solvers`, `flowsheet-core`, and `thermo-desktop` remain scaffolds. No GUI framework is selected.

## Intended boundaries

- `thermo-types` defines shared calculation vocabulary without model or GUI dependencies.
- `thermo-data` owns parameter access and provenance; `thermo-models` implements thermodynamic formulations.
- `thermo-solvers` owns iteration and convergence diagnostics, keeping numerical algorithms separate from formulations.
- `property-packages` assembles supported combinations of data, models, and solvers.
- `flowsheet-core` consumes property capabilities and owns stream/unit-operation calculations.
- `thermo-cli` and `thermo-desktop` consume the calculation layers. Calculation libraries must never depend on either application.

The [calculation contract](calculation-contract.md) documents the implemented SI types, mixture validation, homogeneous PR/PR-BM evaluation, and JSON request/result boundary. Other responsibilities above remain development targets. The cubic-specific root bracketing is internal to `thermo-models`; future flash solvers can depend on models without introducing a model/flash-solver dependency cycle.

## Calculation contracts

The current contract makes SI units, absolute pressure, composition basis, model/data versions, homogeneous-root semantics, numerical guardrails, and distinct failure categories explicit. It exposes residual h/s relative to an ideal gas at the same T/P/composition; total caloric reference states remain to be established. Unsupported capabilities and nonconvergence must not be reported as successful values. Future stateful flowsheet/GUI code must invalidate derived results when inputs change.

## Data and verification

Keep parameter provenance in `data/`, reusable demonstration cases in `examples/`, and independent reference comparisons in `validation/`. Keep historical claims separate from newly executed checks. See the development plan for numerical acceptance targets and deferred capabilities.
