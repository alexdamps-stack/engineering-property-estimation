# Engineering Tools — Property Estimation

A Rust workspace for a planned standalone thermodynamic property engine and steady-state flowsheet simulator.

**Status: homogeneous PR/PR-BM implementation.** Validated SI types, eight pinned component records, NASA7 ideal contributions, cubic roots, fugacities and caloric properties are implemented with offline independent reference comparisons. The CLI evaluates reproducible JSON cases. A4 adds TPD stability assessment, two-phase TP flash and pure saturation pressure with independent reference comparisons ([A4 record](docs/a4-status.md)); A5 adds PH flash with a pure-fluid saturation lever rule, energy diagnostics and independent references ([A5 record](docs/a5-status.md)); A7 adds a directed, temperature-dependent NRTL liquid with an explicitly selected ideal-gas vapor, pinned methanol/water data, and gamma-phi bubble/dew points and TP/PH flash with independent Clapeyron references ([A7 record](docs/a7-status.md)). `flowsheet-core` supplies the A8 sequential engine with native unit operations, balances and the qualified headless demonstration ([A8 record](docs/a8-status.md)); A9 makes `thermo-desktop` a usable `egui` editor: it creates, solves in the background, modifies, saves and reloads projects through the GUI-independent `flowsheet-project` crate, marks exactly the affected results stale after edits, keeps failures visible, and ships as a Windows installer ([A9 record](docs/a9-status.md)). See the [A1–A3 verification record](docs/a1-a3-status.md) for formulation, tests, domains and remaining acceptance gates.

## Development plan

A10 adds deterministic tear-stream recycle solving, closure diagnostics,
bounded iteration and cancellation, plus desktop/project support and an
analytically checked recycle demonstration ([A10 record](docs/a10-status.md)).
Run `cargo run -p thermo-desktop --locked -- --recycle-demo`, open
`examples/recycle.flowsheet.json`, or run the headless
`cargo run -p flowsheet-project --example recycle --locked`.

A6 now supplies a [native 21-component GERG-2008 homogeneous path](docs/a6-status.md),
including native ideal/residual derivatives and caloric properties. Run
`cargo run -p thermo-cli --locked -- examples/gerg-natural-gas.json` for a
temperature/molar-density calculation, or use `examples/gerg-tp-natural-gas.json`
for pressure inversion with an explicit density bracket. Both return analytic
fugacity coefficients. Native GERG TP/PH flashes support the declared dry-gas
domain; examples are `gerg-tp-flash-binary.json` and
`gerg-ph-flash-natural-gas.json`. The A6 record gives bounds and verification.

Read the [English development plan](docs/development-plan.md) for priorities, dependencies, acceptance criteria, and indicative effort. The [original German plan](ENTWICKLUNGSPLAN.md) is preserved as the planning baseline.

Initial numerical priorities are PR/PR-BM, validated component data, stability assessment, and TP/PH flash. Native GERG-2008, NRTL, a flowsheet editor, and later acceleration are separate work packages. Excel integration is deferred.

## Workspace

| Crate | Intended role |
|---|---|
| `thermo-types` | Units, identities, composition bases, states, errors |
| `thermo-data` | Component data, interaction parameters, provenance |
| `thermo-models` | Equations of state, activity models, derivatives |
| `thermo-solvers` | Root finding, phase stability, flash calculations |
| `property-packages` | Compatible model/data combinations |
| `flowsheet-core` | Streams, unit operations, balances, calculation order |
| `flowsheet-project` | Versioned project files, package selection, result snapshots and invalidation |
| `thermo-cli` | Command-line entry point |
| `thermo-desktop` | Native flowsheet editor (`egui`) |

Additional directories: `packaging/` for the Windows installer, `data/` for pinned source transcriptions, `examples/` for reproducible cases, `validation/` for independent reference comparisons, and `docs/` for planning and architecture.

## Getting started

Install [Rust through rustup](https://www.rust-lang.org/tools/install). The repository pins its toolchain in `rust-toolchain.toml`. On Windows, install the MSVC C++ build tools and Windows SDK required by the Rust MSVC toolchain; on Ubuntu/WSL, install `build-essential` and `ca-certificates` first. GitHub CI checks Ubuntu and Windows.

From the repository root:

```sh
rustup show
cargo build --workspace --locked
cargo run -p thermo-cli --locked
cargo run -p thermo-cli --locked -- examples/homogeneous-binary.json
cargo run -p thermo-desktop --locked
```

The desktop command opens the flowsheet editor. Add `-- --demo` to start with the qualified demonstration project (mix, cool, expand, separate; press F5 to solve), or pass a `*.flowsheet.json` project file. `-- --self-test report.json` runs the headless acceptance workflow without a window. The Windows installer is built by `packaging/windows/build.ps1` (see the [A9 record](docs/a9-status.md)). The A8 headless engine runs the independently qualified process through native unit-operation APIs:

```sh
cargo run -p flowsheet-core --example demonstration --locked -- PR-BM
cargo run -p flowsheet-core --example demonstration --locked -- PR 100
```

The optional final argument scales all feed flows. JSON output includes streams,
phase properties, model/data/reference provenance, duties, balances and solver
diagnostics. See the [A8 implementation and verification record](docs/a8-status.md).

For future flowsheet/UI design, see the [process symbol reference](docs/design-reference/iso-10628/README.md):
296 editable SVG symbols, a searchable visual catalog, source references, and AI reuse guidance.

## Checks

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
cargo doc --workspace --no-deps --locked
julia --project=validation/julia validation/julia/check.jl
```

The Julia command is the standard independent PR/PR-BM comparison using Julia
1.11.7 and the locked Clapeyron.jl 0.6.30 environment. It runs on every push/PR
and in scheduled numerical validation. Ordinary Rust tests also compare the
committed Clapeyron reference values offline, without requiring Julia. See
[comparison scope and evidence](validation/clapeyron-checks.md).

### Desktop UI tests

`crates/thermo-desktop/tests/ui.rs` runs the whole editor headlessly with [`egui_kittest`](https://docs.rs/egui_kittest). The tests never open a window or move the real mouse:

- **Interaction tests** inject pointer and key events: placing units, connecting nozzles, moving and deleting.
- **Snapshot tests** render the UI off-screen and compare it with the PNG baselines in `tests/snapshots/`.

After an intended visual change, set `UPDATE_SNAPSHOTS=1` to re-record the baselines, then review the PNGs before committing. On CI the snapshots are rendered but not compared pixel by pixel.

On Windows with **Smart App Control** enabled, freshly built unsigned test executables can be blocked ("An application control policy has blocked this file"). Run the UI tests in WSL instead (needs Ubuntu with rustup and `mesa-vulkan-drivers`):

```powershell
wsl -d Ubuntu --cd "$PWD" -- bash scripts/wsl-ui-tests.sh
wsl -d Ubuntu --cd "$PWD" -- env UPDATE_SNAPSHOTS=1 bash scripts/wsl-ui-tests.sh
```

Numerical tests compare homogeneous PR and PR-BM separately against teqp and NASA7 against Cantera. They verify implementation agreement; A4 stability/flash tests compare against an independent SciPy/teqp equilibrium solution; global stability proof and experimental accuracy are not claimed. See the [testing strategy](docs/testing-strategy.md), [contributing](CONTRIBUTING.md), [architecture](docs/architecture.md), and [artifact inventory](docs/artifact-inventory.md).

## Access and licensing

This project is public and released under the [MIT License](LICENSE). Third-party data, symbols and fonts keep their own licenses, recorded next to them (for example `data/source-extracts/` and `docs/design-reference/iso-10628/`). Workspace packages remain `publish = false`. Access to linked historical ChatGPT conversations depends on the reader's account permissions.
