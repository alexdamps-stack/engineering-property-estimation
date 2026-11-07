# Engineering Tools — Property Estimation

A Rust workspace for a planned standalone thermodynamic property engine and steady-state flowsheet simulator.

**Status: initial scaffold.** Thermodynamic models, flash calculations, and component datasets are not implemented. `flowsheet-core` holds a GUI-independent flowsheet topology model, and `thermo-desktop` is an `egui` prototype of the flowsheet editor; it cannot solve or display results. The CLI prints a scaffold notice only. Passing CI verifies the workspace structure, not numerical accuracy.

## Development plan

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
| `thermo-cli` | Command-line entry point |
| `thermo-desktop` | Native flowsheet editor (`egui` prototype) |

Additional directories: `data/` for qualified data, `examples/` for reproducible cases, `validation/` for independent reference comparisons, and `docs/` for planning and architecture. Each currently contains guidance rather than numerical content.

## Getting started

Install [Rust through rustup](https://www.rust-lang.org/tools/install). The repository pins its toolchain in `rust-toolchain.toml`. On Windows, install the MSVC C++ build tools and Windows SDK required by the Rust MSVC toolchain; on Ubuntu/WSL, install `build-essential` and `ca-certificates` first. GitHub CI checks Ubuntu and Windows.

From the repository root:

```sh
rustup show
cargo build --workspace --locked
cargo run -p thermo-cli --locked
cargo run -p thermo-desktop --locked
```

The desktop command opens the flowsheet editor prototype. Add `-- --demo` to start with the demonstration flowsheet (mix, cool, expand, separate). You can place units, connect ports, and enter specifications; solving stays disabled until the flowsheet engine (A8) exists.

## Checks

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
cargo doc --workspace --no-deps --locked
```

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

No numerical tests or validation claims are present yet. Add meaningful reference-based tests as models are introduced, following the [testing strategy](docs/testing-strategy.md). See also [contributing](CONTRIBUTING.md), [architecture](docs/architecture.md), and the [artifact inventory](docs/artifact-inventory.md).

## Access and licensing

This project is intended for a private GitHub repository. No open-source license has been selected, and workspace packages are marked `publish = false`. Access to linked historical ChatGPT conversations depends on the reader's account permissions.
