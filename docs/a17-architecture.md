# A17 architecture: unit groups, ports, energy streams

Status: foundation step of work package A17 (#170), November 22, 2026.
Read with [a17-implementation-brief.md](a17-implementation-brief.md) (rules,
WSL build) and [unit-operation-gap-analysis.md](unit-operation-gap-analysis.md)
(the 48 items, #171–#218). Progress per group is recorded in
[a17-status.md](a17-status.md).

The foundation adds every missing palette kind (drawn only), a variable and
typed port model, energy streams (#171, implemented and tested), one module
per work-package group in each of the three crates, a grouped palette, and
project format version 4. Seven groups then work in parallel; each edits its
own files and touches the shared files only at the hook points below.

## 1. Groups and the files they own

| Group | Kinds (sub-issues) | `flowsheet-core/src/units/` | `flowsheet-project/src/records/` | `thermo-desktop/src/panels/spec/` |
|---|---|---|---|---|
| rotating | pump, compressor, expander (A14), multistage compressor #182; PS flash | `rotating.rs` | `rotating.rs` | `rotating.rs` |
| heat_transport | air cooler #179, fired heater #180, multi-stream exchanger #181, pipe segment #183, relief valve #184 | `heat_transport.rs` | `heat_transport.rs` | `heat_transport.rs` |
| simple_logical | component splitter #172, multiplier #173, selector #174, stream cutter #175, tank #178, adjust #212, set #213, spreadsheet #214, balance #215, recycle #216, control blocks #217, user unit #218 | `simple_logical.rs` | `simple_logical.rs` | `simple_logical.rs` |
| reactors | #191–#197, reaction sets, formation properties | `reactors.rs` | `reactors.rs` | `reactors.rs` |
| columns | #185, #186, #187, #190 | `columns.rs` | `columns.rs` | `columns.rs` |
| multiphase | VLLE flash, #176, #177, later #188, #189 | `multiphase.rs` | `multiphase.rs` | `multiphase.rs` |
| solids | solids substream model, #198–#211 | `solids.rs` | `solids.rs` | `solids.rs` |

A group may add **new** files freely: submodules of its own module (for
example `units/columns/stage_solver.rs` after turning `columns.rs` into
`columns/mod.rs`), its own integration tests
(`crates/flowsheet-core/tests/<group>*.rs`, `crates/flowsheet-project/tests/<group>*.rs`),
validation scripts under `validation/<group>/`, and docs
`docs/a17-<group>.md`. A new crate (for example a reaction-data or solids
property crate) is allowed when it has a single owner; list it in the report.

`UnitKind::group()` in `crates/flowsheet-core/src/unit.rs` is the
authoritative kind-to-group map.

## 2. What each group module contains

Every `crates/flowsheet-core/src/units/<group>.rs` follows the same pattern:

```rust
pub(crate) fn info(kind: UnitKind) -> KindInfo          // label, tag series, palette category,
                                                        // fixed material + energy ports, issue
pub(crate) fn status(kind: UnitKind) -> CalculationStatus // Calculated | DrawnOnly(reason)
pub enum <Group>Spec { <Kind>(<Kind>Spec), ... }        // one variant per kind
pub(crate) fn default_spec(kind: UnitKind) -> <Group>Spec
impl GroupSpec for <Group>Spec { kind, issues, variable_ports, consumes_energy }
pub(crate) fn calculate(kind, &<Group>Spec, &UnitInputs, &CalcContext)
    -> Result<UnitOutput, String>                       // stub: Err(not_implemented(kind))
```

`UnitSpec` has one variant per group wrapping the group's type
(`UnitSpec::Reactors(ReactorSpec)`, …). `UnitSpec` is no longer `Copy`, so a
spec may hold vectors (stoichiometry, per-component fractions, stage
profiles). Pump, compressor and expander keep the shared `UnitSpec::Pressurize`
/ `UnitSpec::Expand` variants; their metadata, status, calculation
(`rotating::calculate_core`) and form (`spec/rotating.rs::core_form`) are owned
by the rotating group, which may change those two variants and their records
in `records/mod.rs` (list it as a shared change).

## 3. Hook points in shared files

Only these edits to shared files are expected; keep them minimal and list
them in the report.

| File | Hook | Who |
|---|---|---|
| `flowsheet-core/src/units/<group>.rs` → `status()` | Return `CalculationStatus::Calculated` for a kind once it is implemented. This is the only switch; `UnitKind::calculation_planned` follows. | owner |
| `flowsheet-core/src/unit.rs` → `UnitKind::default_spec` | None needed (dispatches to `default_spec`). | — |
| `flowsheet-core/src/units/mod.rs` → `CalcContext` | Do **not** edit; add methods in an `impl CalcContext<'_>` block in your own file (section 5). Move a helper into `mod.rs` only through the integrator. | integrator |
| `flowsheet-core/src/engine.rs` | None for ordinary units. Engine-wide changes (new balance terms, solids in `MaterialStream`, nested solvers) are shared changes. | listed |
| `thermo-desktop/src/symbols.rs` → `TRACED` | Change your kinds' symbol or set `reviewed_ports: true` after reviewing nozzle positions (then the line-work test checks them). | owner |
| `thermo-desktop/src/view.rs` → `block_size`, `port_exit_in` | Add a match arm for your kind when the generic size or nozzle side is wrong (for example a column's distillate at the top). | owner, one arm per kind |
| `thermo-desktop/src/canvas.rs` → `name_anchor` | Optional arm for your kind. | owner |
| `thermo-desktop/src/panels/unit.rs` → `unit_results` | Optional rows for group results; prefer adding `UnitCalculation`/`UnitResult` fields through the integrator. | listed |
| `flowsheet-project/src/records/mod.rs` | None for group records. | — |
| `flowsheet-project/src/project.rs` → `FORMAT_VERSION`, `migrate` | **Integrator only.** Groups never bump the format version. | integrator |
| `flowsheet-project/src/lib.rs` → `ENGINE_REVISION` | Integrator only (stale-marks stored results). Mention in the report if your change alters existing numerics. | integrator |
| `docs/a17-status.md` | Your own section only. | owner |

`UnitKind`, `KindRecord` and the palette order already contain every A17 kind;
adding a further kind is an integrator change.

## 4. How to add a specification field, record, form and test

1. **Spec** (`units/<group>.rs`): add a field to `<Kind>Spec` (enums for
   mutually exclusive choices, `Option<f64>` for required values without a
   default, SI units, absolute pressure) and keep `Default` meaningful.
   Report missing/invalid values in `GroupSpec::issues` using
   `crate::check_optional` / `crate::check_value` with a `Sign`. A specified
   value that the kind cannot use gets `SpecProblem::NotAvailableForKind`;
   counts out of range get `SpecProblem::CountOutOfRange`.
2. **Record** (`records/<group>.rs`): add the field to the kind's struct
   variant with `#[serde(default)]` (or `default = "fn"`), map it in both
   `From<&Spec>` and `TryFrom<Record>`, and return every `f64` from
   `numbers()` (JSON cannot store NaN/inf; `Project::to_json` checks these).
   Records use `deny_unknown_fields`, so an older build refuses a newer file
   instead of dropping fields. **Do not change `FORMAT_VERSION`**; note
   "fields added behind format version 4" in your report.
3. **Form** (`panels/spec/<group>.rs`): edit the field with the helpers of
   `panels/spec/mod.rs` (`choice`, `opt_field`, `value_field`, `count_field`,
   `note_row`, `drawn_only_note`). Child modules can call these private
   helpers directly. Keep `drawn_only_note` while the kind is drawn only.
4. **Tests**: unit tests in your module (`#[cfg(test)]`), engine tests in
   `crates/flowsheet-core/tests/<group>.rs` (behaviour, failure paths,
   component and energy balance via `UnitCalculation::balance`, independent
   reference values from `validation/`), and a record round trip in
   `records/<group>.rs` tests. The shared contract test
   `units::testing::check_drawn_only_contract(UnitGroup::…)` stays in your
   module.

## 5. Calculation context API

`calculate(kind, spec, inputs, ctx)` receives:

- `inputs.material: &[MaterialStream]` — one per material inlet port in port
  order (an unconnected optional inlet is a zero-flow stream with no state);
- `inputs.energy: &[Option<f64>]` — one per **energy inlet** port: the duty in
  W delivered into the unit, or `None` if unconnected or if the stream's
  source is an energy source set to "calculated" (the unit then determines
  the value);
- `ctx: &CalcContext` with
  - `package()`, `model()`, `component_count()`, `component_ids()`, `feed()`;
  - `tp(z, T, p)`, `ph(z, p, h)` → `EquilibriumState` (two-phase VLE);
  - `ps(z, p, s)` → explicit error until the rotating group implements it
    (the stub lives in `units/rotating.rs`);
  - `homogeneous(x, T, p, RootSelection)` → `HomogeneousState` (h, s, cp,
    cv, ρ, Z, ln φ, speed of sound, …, no stability claim);
  - `ln_fugacity_coefficients(x, T, p, RootSelection)` → `FugacityState`;
  - `k_values(x_liquid, y_vapor, T, p)` → `K_i = φ_i^L / φ_i^V`.

All methods return `Result<_, String>` so calculations use `?` and fail
explicitly. Groups extend the context **additively** in their own file:

```rust
// crates/flowsheet-core/src/units/multiphase.rs
impl CalcContext<'_> {
    pub(crate) fn vlle(&self, z: &[f64], t: f64, p: f64) -> Result<ThreePhaseState, String> { … }
}
```

The calculation returns `UnitOutput { outlets, duty_w, transferred_duty_w,
material_exchange, warnings, energy }`: one `MaterialStream` per material outlet port in port
order; `duty_w` is the net external duty (positive into the unit) and enters
the energy balance that the engine checks (1e-8 relative component, 2e-9
scaled energy); `energy` gives one value per energy port (positive in the
port's direction) or stays empty for the default — a unit with exactly one
energy port puts `duty_w` through it (negated for an outlet).

The multiplier and selector supply `material_exchange` for their additional
flowsheet boundary. The engine subtracts its signed component and enthalpy
flows from the outlet-minus-inlet residuals. Positive flow enters the unit.
The result preserves the boundary contribution separately from heat duty.

## 6. Port model

- `PortSpec { name, number, direction, medium, optional }`;
  `medium` is `Material` or `Energy`; `number` marks a member of a numbered
  group (`feed 2`); `label()` renders it.
- `Unit::ports()` = `kind.ports()` (fixed material) ++ `kind.energy_ports()`
  (fixed, optional) ++ `spec.variable_ports()` (from `GroupSpec`). A
  `PortRef.port` indexes this list. Fixed indices therefore never change, and
  projects from format versions 1–3 keep their material port indices.
- Variable groups (column feeds and side draws, multi-stream exchanger
  passes, selector inlets) come from the specification.
  `Flowsheet::try_set_spec` keeps each connected stream on the port with the
  same identity (`PortSpec::same_port`) and re-indexes it; a change that
  would remove a connected port returns `SpecChangeError::ConnectedPortRemoved`
  and changes nothing (the desktop shows it in the status line).
- Optional ports produce no "not connected" diagnostic. Every energy port of
  an equipment unit is optional; the energy source/sink ports are required.
- `connect` rejects mixing media (`ConnectError::MediumMismatch`). Material
  streams are named `S<n>`, energy streams `Q<n>`.
- Desktop: material inlets on the left, outlets on the right, energy inlets
  on the top edge and energy outlets on the bottom edge, spread evenly
  (`view::unit_port_pos`). Energy streams are drawn dashed.

## 7. Energy streams (#171)

An energy stream carries a duty in W, **positive in the stream direction**:
into the destination unit, out of the source unit. Unit duties keep the
convention "positive into the unit".

| Kind | Energy port | Role |
|---|---|---|
| Energy source | `Q` out (required) | *Calculated* (default): reports the duty its destination determines. *Specified*: sets the duty. |
| Energy sink | `Q` in (required) | Records the duty it receives. |
| Heater | `Q` in (optional) | Determines (outlet T or duty), or consumes with the thermal choice *Energy stream*. |
| Cooler | `Q` out (optional) | Determines; the stream carries the removed heat (−duty). |
| Pump, compressor | `power` in | Rotating group (A14). |
| Expander | `power` out | Rotating group (A14). |
| A17 kinds | see `info()` | Owning group; default role *determines*. |

Each energy stream must be determined by exactly one end
(`UnitSpec::energy_role` → `EnergyRole::{Determines, Consumes}`):

- determines → consumes: the value follows the calculation order (heat
  integration cooler → heater, specified source → heater);
- consumes → determines: allowed only from an energy source (back-filled
  after the destination is calculated);
- both determine: `DiagnosticKind::EnergyOverspecified`;
- neither: `DiagnosticKind::EnergyUnderspecified`;
- a consuming unit without its stream: `DiagnosticKind::EnergyStreamMissing`.

The engine checks that a delivered duty equals the duty the unit used
(relative 1e-9), stores `Calculation::energy_streams`, and the project keeps
`results.energy_streams`. Energy streams are part of the calculation order
but are never torn: automatic tears choose material streams, and explicit
energy tears are rejected. A unit taking its duty from a stream within a
recycle loop is solved inside each pass. A source that consumes an energy
outlet (other than an energy source) is not supported in steady-state
sequential order (`EnergyAgainstFlow`).

Tests: `crates/flowsheet-core/tests/energy.rs`,
`crates/flowsheet-project/tests/project.rs`
(`energy_streams_and_every_unit_kind_round_trip_in_version_4`).

## 8. Project format

Format version **4** (`FORMAT_VERSION`): energy specifications, the heater
`energy_stream` choice, `results.energy_streams`, every A17 kind, and group
records nested as

```json
"spec": {"type": "columns", "unit": {"kind": "distillation_column", "feeds": 1, "side_draws": 0}}
```

Versions 1–3 load through `project::migrate` (one JSON step per version;
3 → 4 is the identity because version-3 content is a subset) and are written
back as version 4. **Only the integrator bumps `FORMAT_VERSION`** and adds a
`migrate` step; groups add record fields with serde defaults behind version 4.

## 9. Palette and symbols

The desktop add palette is grouped into collapsible categories (`PaletteCategory`):
Streams, Mixing/splitting, Separators, Heat transfer, Pressure change,
Columns, Reactors, Logical and Extensibility. The deferred Solids category
stays in the project registry for saved units. Kinds that are drawn only
are shown in grey with the reason and sub-issue in the hover text, in the
specification form, and as a warning message; a solve stops in preflight at
the first of them (`"<kind> is drawn only; … (#nnn)"`).

Symbols come from `docs/design-reference/iso-10628` (`symbols.rs::TRACED`,
attribution CC BY-SA 4.0 kept in each embedded SVG). Logical blocks (adjust,
set, spreadsheet, balance, recycle, control, user unit, multiplier) have no
ISO 10628 symbol and are labelled boxes (`fallback_code`); feed, product,
splitter and the energy boundaries keep project drawings. New symbols start
with `reviewed_ports: false`: their nozzles sit on the block edge until the
owning group reviews positions against the line work.

## 10. Checks

In WSL (see the brief): `cargo fmt --all`,
`cargo clippy --workspace --all-targets -- -D warnings`,
`cargo test --workspace` with your own `CARGO_TARGET_DIR`. UI snapshot
changes are reviewed and regenerated with
`UPDATE_SNAPSHOTS=1 cargo test -p thermo-desktop --test ui`. Model, caloric
or solver changes also run the Clapeyron check
(`julia --project=validation/julia validation/julia/check.jl`).
