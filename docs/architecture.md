# Architecture

The repository has six library crates and two binary crates. `thermo-types`,
`thermo-data` and `thermo-models` provide validated inputs, embedded component
data and native homogeneous PR/PR-BM calculations; see [A1–A3 status](a1-a3-status.md).
`thermo-solvers` and `property-packages` remain scaffolds. `flowsheet-core` holds
a GUI-independent topology model: units, ports, streams, specifications, and
pre-calculation diagnostics. `thermo-desktop` is a GUI prototype built on
`egui`/`eframe` and depends on `flowsheet-core`. Using `egui` here is for the
prototype only and does not decide the final framework (see the
[flowsheet and interface proposal](flowsheet-and-ui.md), section 4.7).
Canvas layout, zoom, and selection stay in the desktop crate.

## Fundamental engine constraint

The engine performs its calculations in process in native Rust. It must not
delegate property or solver calculations to external APIs, services, CLIs,
subprocesses or foreign thermodynamic runtimes. Calculation requires no network
connection, Python installation, teqp, Cantera, CoolProp executable or other
external calculation program. Unsupported capabilities return explicit errors
rather than falling back to an external engine.

The repository's own `thermo-cli` is an application that calls the Rust library;
the engine never invokes or depends on it. Separately installed reference tools
under `validation/` may generate independent expected values for development
verification. They are outside the engine dependency/runtime path. Normal Rust
tests read committed fixtures offline and never launch those tools. Component
constants and polynomial coefficients are embedded data with retained provenance,
not calls to the software projects from which they were transcribed.

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
