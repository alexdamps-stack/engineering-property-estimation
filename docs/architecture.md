# Architecture

The workspace (`crates/*`) has eleven crates: nine library crates, the
`thermo-cli` binary, and `thermo-desktop`, which builds the editor executable
from a library target plus a thin `main.rs`.

- `thermo-types`, `thermo-data` and `thermo-models` provide validated inputs,
  embedded component data and native homogeneous calculations: PR/PR-BM
  ([A1–A3 status](a1-a3-status.md)), GERG-2008 ([A6](a6-status.md)) and the
  NRTL excess model with a gamma–phi package ([A7](a7-status.md)).
- `thermo-solvers` implements Rachford–Rice, TPD stability, two-phase TP/PH
  flash and pure saturation ([A4 status](a4-status.md)), native GERG density
  inversion and TP/PH flash, NRTL bubble/dew and flash, and bounded PR/PR-BM
  mixture bubble/dew points and envelope tracing ([A13](a13-status.md)), and
  multiphase (VLLE/LLE) stability and TP/PH flash for PR/PR-BM and NRTL-IG
  ([VLLE status](vlle-status.md)).
- `property-packages` supplies the immutable flowsheet packages behind the
  `FlowsheetPackage` interface: cubic (PR, PR-BM, SRK), NRTL-IG, GERG-2008 and
  EOS-CG-2021. Cubic and NRTL-IG also provide the three-phase `vlle_tp`/`vlle_ph`
  operations; the others return `UnsupportedCapability`. Each retains model, interaction, data and caloric-reference
  metadata. New models and packages follow the
  [property-model requirements](property-model-requirements.md).
- `flowsheet-core` holds a GUI-independent topology and sequential calculation
  engine: units, ports, material streams, specifications, balances and
  diagnostics ([A8 record](a8-status.md)), with tear-stream recycle solving
  ([A10](a10-status.md)).
- `flowsheet-project` holds the versioned project file, package selection,
  result snapshots and per-unit result invalidation for the editor
  ([A9](a9-status.md)). It is independent of any GUI framework.
- `property-tables` is the bounded A11 tabulation experiment for one fixed PR
  mixture ([A11](a11-status.md)). No other crate depends on it; it does not
  replace the direct model in the flowsheet or project format.
- `thermo-bridge` is the A12 external-integration boundary
  ([A12](a12-status.md)): a versioned JSON-over-C-ABI library (`cdylib` and
  `rlib`, library name `thermocore_bridge`) for the DWSIM client in
  `integrations/dwsim`. It only calls `CubicPackage` TP/PH flash and contains
  no model or solver code.
- `thermo-cli` runs reproducible JSON cases (`examples/`) against the
  libraries.
- `thermo-desktop` is the native flowsheet editor built on `egui`/`eframe`. It
  depends on `flowsheet-core`, `flowsheet-project`, `property-packages` and
  `thermo-types`. egui is the selected GUI framework (see the
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
- `flowsheet-project` owns persistence and result currency; it does not depend on a GUI framework.
- `thermo-bridge` translates between the foreign-host JSON contract and `property-packages`; it adds no calculations.
- `thermo-cli` and `thermo-desktop` consume the calculation layers. Calculation libraries must never depend on either application.

These are intended responsibilities, not precommitted Rust trait signatures. Add dependencies only when real interfaces are introduced. (Historical: this list was written for the initial infrastructure setup, before the numerical API and GUI decisions; the crate list above records the current state.)

## Calculation contracts (established in A1)

SI units, absolute pressure, composition basis, reference states, model/data versions, validity domains, and distinct failure categories must be explicit; see the [calculation contract](calculation-contract.md). Unsupported capabilities and nonconvergence must not be reported as successful values. Input changes invalidate derived results.

## Data and verification

Keep parameter provenance in `data/`, reusable demonstration cases in `examples/`, and independent reference comparisons in `validation/`. Keep historical claims separate from newly executed checks. See the development plan for numerical acceptance targets and deferred capabilities.
