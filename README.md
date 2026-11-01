# Engineering Tools — Property Estimation

A Rust workspace for a planned standalone thermodynamic property engine and steady-state flowsheet simulator.

**Status: first homogeneous engine slice.** PR and PR-BM now calculate explicitly selected roots, density, fugacity coefficients, and residual enthalpy/entropy through a JSON CLI. Eight component records are bundled with provenance. Sixteen independent teqp reference cases check this implementation. Phase stability, TP/PH flash, total caloric properties, GERG, NRTL, flowsheet solving, and the graphical interface remain planned. See the [calculation contract](docs/calculation-contract.md) for exact formulation and limits.

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
| `thermo-desktop` | Future native graphical application |

Additional directories: `data/` for versioned constants and provenance, `examples/` for runnable JSON requests, `validation/` for independent reference fixtures and their generator, and `docs/` for planning, calculation contracts, architecture, and testing.

## Getting started

Install [Rust through rustup](https://www.rust-lang.org/tools/install). The repository pins its toolchain in `rust-toolchain.toml`. On Windows, install the MSVC C++ build tools and Windows SDK required by the Rust MSVC toolchain; on Ubuntu/WSL, install `build-essential` and `ca-certificates` first. GitHub CI checks Ubuntu and Windows.

From the repository root:

```sh
rustup show
cargo build --workspace --locked
cargo run -p thermo-cli --locked -- components
cargo run -p thermo-cli --locked -- calculate --input examples/natural-gas.json
cargo run -p thermo-cli --locked -- calculate --input examples/binary-hydrocarbon.json
```

The examples use explicit illustrative interaction parameters. Results are homogeneous roots with `phase_stability: "not_assessed"`; they are not equilibrium flashes. Residual h/s are not total enthalpy/entropy. `thermo-desktop` is still a console placeholder and does not open a GUI.

## Checks

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
cargo doc --workspace --no-deps --locked
cargo test -p property-packages --test state_sweep --locked -- --ignored --nocapture
```

The ordinary suite covers 15 tests, including 16 independent reference cases, derivatives, identities, composition invariance, critical-root rejection, input failures, and CLI behavior. The separately invoked 200-case grid checks robustness, not independent accuracy; a scheduled/manual workflow runs it after that workflow reaches the default branch. These checks do not establish experimental accuracy or a fully validated operating domain. Follow the [testing strategy](docs/testing-strategy.md) when extending the implementation. See [reference evidence](validation/README.md), [contributing](CONTRIBUTING.md), [architecture](docs/architecture.md), and the [artifact inventory](docs/artifact-inventory.md).

## Access and licensing

This project is intended for a private GitHub repository. No open-source license has been selected, and workspace packages are marked `publish = false`. Access to linked historical ChatGPT conversations depends on the reader's account permissions.
