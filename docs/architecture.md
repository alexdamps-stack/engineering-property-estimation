# Architecture

The repository starts with six library crates and two binary crates. Library crates contain module documentation only. Binary crates report that implementation is pending. There are no public calculation functions, cross-crate dependencies, external dependencies, or selected GUI framework yet.

## Intended boundaries

- `thermo-types` defines shared calculation vocabulary without model or GUI dependencies.
- `thermo-data` owns parameter access and provenance; `thermo-models` implements thermodynamic formulations.
- `thermo-solvers` owns iteration and convergence diagnostics, keeping numerical algorithms separate from formulations.
- `property-packages` assembles supported combinations of data, models, and solvers.
- `flowsheet-core` consumes property capabilities and owns stream/unit-operation calculations.
- `thermo-cli` and `thermo-desktop` consume the calculation layers. Calculation libraries must never depend on either application.

These are intended responsibilities, not precommitted Rust trait signatures. Add dependencies only when real interfaces are introduced. The scope of this setup is infrastructure; the development plan explicitly schedules numerical API and GUI decisions later.

## Calculation contracts to establish in A1

SI units, absolute pressure, composition basis, reference states, model/data versions, validity domains, and distinct failure categories must be explicit. Unsupported capabilities and nonconvergence must not be reported as successful values. Input changes invalidate derived results.

## Data and verification

Keep parameter provenance in `data/`, reusable demonstration cases in `examples/`, and independent reference comparisons in `validation/`. Keep historical claims separate from newly executed checks. See the development plan for numerical acceptance targets and deferred capabilities.
