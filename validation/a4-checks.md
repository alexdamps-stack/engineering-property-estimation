# A4 checks — stability, TP flash, pure saturation

Date: 2026-11-07. Starting revision: `538d9f0`. The results apply to the local,
uncommitted A4 change. Toolchain: the repository's pinned Rust toolchain on
Windows and in WSL Ubuntu.

| Check | Outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| Windows `cargo clippy -p thermo-types -p thermo-data -p thermo-models -p thermo-solvers -p thermo-cli --all-targets --locked -- -D warnings` | passed |
| WSL `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| Windows `cargo test -p thermo-types -p thermo-data -p thermo-models -p thermo-solvers -p thermo-cli --locked` | passed, 40 tests |
| WSL `cargo test --workspace --locked` | passed, 89 tests including 13 headless UI tests |
| Windows `RUSTDOCFLAGS=-D warnings cargo doc` for the numerical crates | passed |
| WSL CLI `examples/tp-flash-binary.json` | two phases, β = 0.3942916347, ln f residual 6.1e-12 |
| WSL CLI `examples/tp-flash-natural-gas.json` | two phases, β = 0.9874922296, ln f residual 8.7e-11 |
| WSL CLI homogeneous examples | unchanged Z 0.9184554903 and 0.8809015715 |
| WSL CLI at 150 K | nonzero exit with `OutsideValidityDomain` |
| Reference generation (`generate-flash-references.py`) | 22 flash and 30 saturation states; every assertion passed; no state skipped |

Maximum deviations from the independent fixture were:

| Quantity | Deviation |
|---|---|
| β | 2.5e-11 |
| Phase mole fractions | 3.5e-11 |
| Phase Z | 3.3e-12 |
| Saturation pressure | 2.0e-14 relative |

All are well inside the pre-set budgets in the fixture. Phase counts agree for
all 22 states. The input data hash in the fixture
(`072cc688…b147`) matches `data/source-extracts/components.provisional.json`.

Windows application control blocked the build scripts of the desktop
dependencies and the standalone CLI executable (os error 4551). This matches the
A1–A3 record. Those checks ran in WSL instead, and no OS policy was changed.
Hosted CI has not been run, and no independent code review of A4 has been done
yet.
