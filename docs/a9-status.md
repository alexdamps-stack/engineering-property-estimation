# A9 native interface and flowsheet editor

Implemented November 8, 2026 on `alexdamps-stack/A9`, based on `main` commit
`3241cab` (A8 merged). The repository-pinned Rust is 1.99.1.

| Issue | Delivered evidence |
|---|---|
| #52: framework decision | Section 4: criteria assessed on the A9 editor, 100-unit frame-time measurement, recorded decision for **egui** in [`flowsheet-and-ui.md`](flowsheet-and-ui.md) §4.7 and §5; Slint not prototyped, reasons recorded |
| #53: create, calculate, modify | Section 2: the editor consumes the A8 engine through the GUI-independent `flowsheet-project` crate; background solve with progress, cancel and auto-solve; statuses, results, stream table, feed compositions, new-project dialog |
| #54: save/reload and invalidation | Section 3: versioned JSON project format, per-unit input fingerprints, failures and cancellations never leave old values current, save → reload → recalculate reproduces results exactly |
| #55: Windows packaging | Section 5: statically linked executable, Inno Setup per-user installer, portable ZIP and notices; CI builds, installs, self-tests the installed copy and uninstalls on `windows-latest` |
| #51: parent acceptance | Section 6: the acceptance script of `flowsheet-and-ui.md` §4.1 runs headless (`--self-test`, `run_acceptance`) and through the interface (`tests/ui.rs`) |

## 1. Architecture

```
thermo-desktop (egui)  ──►  flowsheet-project  ──►  flowsheet-core (A8 engine)
  canvas, panels, dialogs     package selection,        topology, specs,
  background SolveJob         project file, results,    sequential solve
                              fingerprints, acceptance  ──► property-packages (PR/PR-BM)
```

`flowsheet-project` has no GUI dependency. It holds everything that decides
calculation results, persistence and result currency. All of it is tested
headlessly. The desktop crate only converts layout and display settings and
renders status. `flowsheet-core` gained two non-numerical APIs:
`Flowsheet::restore` rebuilds stored units and streams through the ordinary
connection checks, and `Flowsheet::solve_with_progress` reports each unit before
it is calculated and can stop between units. `solve` is unchanged in behavior.
No model, solver, data, fixture or tolerance was changed.

## 2. Create, calculate, modify (#53)

- **Project creation.** *File ▸ New project* opens a dialog. It offers the
  embedded components in database order, PR or PR-BM, and the PH search
  interval (200..600 K by default). No interaction data is applied: the
  project records an explicit all-zero kij matrix with the source note
  "explicit synthetic zero-kij approximation selected in the desktop New
  project dialog; not measured mixture accuracy". The *Demonstration* button
  loads the qualified A8 process exactly: same compositions, 330/350 K at
  5 MPa, 100/20 mol/s, 240 K cooler, 1.5 MPa valve, PR-BM, the A2 zero-kij
  note and a 200..400 K PH interval. Its streams are named as in
  `validation/fixtures/process-v1.json`.
- **Inputs.** Unit specifications use the existing forms. Feeds get a mole
  fraction editor in package order, with its sum and a *Normalize* button. A
  missing composition or a sum outside `1e-10` makes the feed underspecified
  and appears in Messages. Changing the equation of state in the project
  panel invalidates every result.
- **Calculation.** *Solve* (F5) or auto-solve runs `solve_project` on a worker
  thread with a snapshot of the inputs and the cached immutable package. The
  toolbar shows `completed/total` units. *Cancel* (Esc) drops the job at once
  and the worker stops before its next unit. A property calculation already
  running is not interrupted, but its result is discarded. Auto-solve starts
  only when the inputs differ from the last attempt and no pointer drag is in
  progress.
- **Status (proposal §4.6).** Each unit has exactly one status, shown by
  colour *and* badge shape: underspecified `!`, ready (no badge), solved
  (green check), warning (amber triangle), failed (red cross) or stale (grey
  dots with a dashed outline). The properties panel shows duty, outlet T/p/flow
  per port, balance residuals and warnings. The stream panel shows flows, T, p,
  vapor fraction, h, s, a phase table (fraction, density, Z, compositions),
  provenance, iterations, residuals and solver diagnostics. The stream table
  has a status row, and *Copy* puts it on the clipboard as tab-separated text
  in display units. Stale values are grey italic and marked "stale". A failed
  unit shows the cause and no values. Nonfinite numbers are never displayed.

## 3. Persistence and invalidation (#54)

**Project format, version 1** (`*.flowsheet.json`, written by
`Project::to_json`) contains:

- `format_version: 1` and a `generator` string;
- `package`: variant (`PR`/`PR-BM`), ordered component identifiers, the full
  kij matrix, `interaction_source`, and the PH interval;
- `units`: id, kind, tag and a tagged specification record. Values are in SI
  with absolute pressure; `null` means unset;
- `streams`: id, name, and `from`/`to` unit and port;
- `feeds`: mole fractions per feed in package order;
- `layout` (world coordinates) and `display.unit_set`;
- `results`: the last successful snapshot with its provenance (engine
  revision, model identifier and version, data identifier and version), the
  completion time, the calculation order and per-unit and per-stream values,
  each with the input fingerprint it was computed for.

The version is read before the content. Any other or missing
`format_version` fails with a clear message, and there is no silent
migration. Unknown fields, mismatched specification variants, duplicate
identifiers, invalid connections, feed fractions for non-feeds or of the
wrong length, and layout entries for unknown units are all rejected. Opening
a bad file leaves the current project unchanged. Floats round-trip exactly
(`serde_json` `float_roundtrip`), so save → load → save is byte-identical.
New units after a reload never reuse an identifier or a generated tag.

**Invalidation.** A unit fingerprint (64-bit FNV-1a over a canonical
description) covers:

- the package selection together with the running program's provenance;
- the unit kind and specification;
- the feed composition, for feeds;
- for each inlet, the fingerprint of the upstream unit and the outlet port.

It is computed recursively, so an edit changes the fingerprints of exactly
the edited unit and everything downstream of it. A stream result belongs to
its source unit and port. Names and layout are not inputs. A stored result is
*current* only while its recorded fingerprint equals the current one. A
different engine revision, model or data version, package or input therefore
shows it as *stale*. `ENGINE_REVISION` must be changed whenever
`flowsheet-core` numerics change.

**Failures.** A failed or cancelled solve returns no results, and the
previous snapshot is kept unchanged. The failed unit is shown as failed while
its inputs are unchanged. Its own result and every downstream result are
stale by fingerprint. Upstream results whose inputs did not change remain
current: they are the values the current inputs produce. A cancelled solve
only reports the cancellation.

## 4. Framework evaluation and decision (#52)

The A9 editor was assessed against the §4.7 criteria. The decision and the
summary table are in [`flowsheet-and-ui.md`](flowsheet-and-ui.md) §4.7 and §5.

| Criterion | Evidence |
|---|---|
| Canvas performance | `canvas_frame_time_with_100_units` (ignored by default; run with `cargo test --release -p thermo-desktop --test ui -- --ignored --nocapture`): 150 units (50 feed → heater → product chains) and 100 streams, with pointer movement in every frame. Mean CPU time per frame was **0.78 ms**: full window layout, canvas hit-testing, routing, tag placement and tessellation, excluding GPU submission. Measured November 8, 2026 on an AMD Ryzen 7 PRO 4750U in WSL Ubuntu, release build. Not a gated metric |
| Forms and tables | §2 above; `egui::Grid` suffices for the demonstration scale |
| Keyboard and clipboard | Shortcuts listed in §2; TSV copy checked by `stream_table_copies_values_with_units_and_status` |
| High DPI | egui/winit scale-factor handling; the vector canvas scales with zoom. Not snapshot-tested at a scale factor other than 1 |
| Windows packaging | Section 5 |
| Accessibility | eframe's AccessKit integration is on. The UI tests query buttons, radio buttons and labels through the accessibility tree. Canvas items have no accessible nodes (limitation) |
| License | egui/eframe MIT or Apache-2.0. Slint is GPLv3, royalty-free with attribution, or commercial |

Slint was not prototyped; the reasons are recorded in `flowsheet-and-ui.md`
§4.7. **Decision: egui.**

## 5. Windows packaging (#55)

- `cargo build --profile dist --target x86_64-pc-windows-msvc` with
  `-C target-feature=+crt-static`. The `dist` profile is `release` plus thin
  LTO, one codegen unit and stripped symbols. The executable needs no
  Visual C++ redistributable; wgpu uses the system's Direct3D 12 or Vulkan
  drivers.
- `packaging/windows/build.ps1` stages the executable, `README.txt` and
  `notices/`. The notices are `NOTICE.txt` (ISO 10628 symbol attribution,
  CC BY-SA 4.0; embedded CoolProp/Cantera/thermo data; egui fonts), the three
  data license texts, and `THIRD-PARTY-CRATES.txt`, generated from the locked
  dependency graph with each crate's version, declared license and
  repository. The script then writes the portable ZIP and runs Inno Setup 6
  (`thermo-desktop.iss`).
- The installer installs per user by default, without administrator rights
  (all-users is selectable). It requires Windows 10 or later on x64 and
  creates Start-menu entries for the editor, the demonstration and the
  notices, with an optional desktop icon and an uninstaller.
- `packaging/windows/smoke-test.ps1` installs silently, checks the installed
  files and Start-menu shortcuts, runs the installed executable's `--version`
  and `--self-test`, then uninstalls and confirms removal. The CI job
  `windows-package` runs build, self-test of the built executable and the
  smoke test, and uploads the installer, ZIP and reports.

Sizes from the CI run on `windows-latest` (run 42511465567, commit `e410dfd`):

| Artifact | Size |
|---|---|
| `thermo-desktop.exe` (static CRT, stripped) | 17.0 MB |
| `PropertyEstimationFlowsheet-0.0.0-x64-setup.exe` | 7.06 MB |
| `PropertyEstimationFlowsheet-0.0.0-x64-portable.zip` | 7.20 MB |

On the Windows runner, both the built executable and the installed copy
passed `--self-test`. The installed copy reproduced the Linux values: separator
inlet 213.981 K and 219.063 K, and a reload difference of 0.

Release builds use the Windows GUI subsystem, so `--self-test` and
`--version` write to a file named on the command line.

## 6. Acceptance (#51)

The §4.1 acceptance script is automated at two levels:

1. `flowsheet_project::run_acceptance`, exposed as `thermo-desktop --self-test`:
   - create and solve the demonstration, then check all results current and
     the A8 balance limits;
   - change the valve to 20 bar(a): exactly the valve, separator, products
     and their streams turn stale;
   - recalculate;
   - save, reload (byte-identical), check the results current and
     recalculate: the largest relative difference is **0**;
   - a valve pressure above its inlet fails at the valve and leaves the old
     valve result stale;
   - cancellation returns nothing;
   - format version 2 is rejected.
2. `crates/thermo-desktop/tests/ui.rs` drives the same steps through the
   interface: *Solve* button, F5, the *Cancel* button, the New-project dialog
   with *Create* and *Enter composition*, and save and reopen with
   `FlowsheetApp::with_file`. Snapshots record the solved demonstration with
   the stream table and a failed-valve/stale-separator state.

The demonstration project reproduces the independent A2 process fixture
within its approved budgets for both PR and PR-BM
(`demonstration_project_reproduces_the_independent_process_fixture`). This
uses the project layer the desktop uses, not the A8 example code.

## 7. Verification

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | Passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | Passed |
| `cargo test --workspace --locked` (WSL Ubuntu, lavapipe) | 192 passed, 0 failed, 1 ignored (the frame-time measurement), including pixel comparison of all 6 UI snapshots |
| `cargo doc --workspace --no-deps --locked` with `RUSTDOCFLAGS=-D warnings` | Passed |
| Release `thermo-desktop --self-test` (Linux) | Passed, reload difference 0 |
| GitHub CI run 42511465567 | Workspace on Ubuntu and Windows, live Julia 1.11.7 Clapeyron check, and the Windows installer job (build, self-test, install, installed self-test, uninstall): all passed |

Snapshot baselines changed intentionally. All four existing images now show
the menu bar and readiness statuses, and two images were added. All were
regenerated in WSL with lavapipe and reviewed visually. Native Windows
compilation on the development machine is still blocked by Smart App Control
(OS error 4551, as in A8). Native Windows evidence comes from the CI runner.
Julia is not installed locally; the live Clapeyron comparison ran in CI.
No model, data, fixture or tolerance changed.

## 8. Limitations and follow-up

- No undo/redo, copy/paste of units, typed-unit input (`5 barg`), gauge
  pressure, CSV export, localization or locale-aware number formats
  (proposal §4.3–4.5).
- One flowsheet-wide cubic package; NRTL/GERG packages are not selectable in
  the editor. Interaction data cannot be entered in the dialog. The only
  offered choice is the explicit zero-kij approximation.
- Cancellation acts between units. A single long flash runs to completion in
  the background, and its result is discarded.
- The project format has no migration path yet; the next format version must
  add one explicitly.
- The installer is unsigned, so Windows SmartScreen and Smart App Control may
  warn about or block it until code signing is arranged. No file association
  is registered for `*.flowsheet.json`.
- `THIRD-PARTY-CRATES.txt` lists declared licenses. A full license-text
  bundle and a legal review are needed before public distribution.
- Recycles remain A10, saturation specifications A13.
