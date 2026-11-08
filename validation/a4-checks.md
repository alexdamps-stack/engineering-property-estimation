# A4 checks — stability, TP flash, pure saturation

## Newton acceleration follow-up — 2026-11-08

Baseline: `68ef34a`. Rust 1.99.1 in WSL Ubuntu; Julia 1.11.7 with the locked
Clapeyron 0.6.30 environment. The [scoped review](../docs/a3-a4-acceptance.md)
records A3/A4 acceptance criteria and the eighteen-state iteration comparison.
The November 7 report below is retained as historical evidence.

| Check | Current outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| WSL `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| WSL `cargo test -p thermo-solvers --lib --test flash --locked` | passed: 9 unit tests and 10 flash tests |
| WSL `cargo test --workspace --locked` | passed, including cubic/PH, GERG, NRTL and headless UI tests |
| WSL `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| Julia 1.11.7 standard `validation/julia/check.jl` | passed: fresh 60 homogeneous / 28 PH PR/PR-BM states and all six NRTL comparison tests |
| WSL A3 homogeneous natural-gas and binary CLI examples | passed: Z = 0.9184554902895139 / 0.8809015715485475, with model/data/property-unit/diagnostic metadata |
| `python3 validation/consolidate-tolerances.py` | passed; registered budgets match fixtures |
| `git diff --check` | passed |

The 22 flash and 30 saturation comparisons pass with the original budgets.
All eighteen difficult binary states match the eigenvalue-only solution and use
6–20 split iterations instead of 63–599. Three focused unit tests cover coupled
row pivoting, backtracking and rejected acceleration. The committed reference
fixtures, generators, Julia manifest and tolerance registry are unchanged.

The extended 80/20 near-critical-region probe exposed drift toward the trivial
solution with residual-only Newton acceptance. Physical trial fractions and a
split-Gibbs guard resolved it; all twelve added 80/20 states are retained as
regressions, alongside the six 75/25 states. Gibbs-increasing Newton candidates
are also explicitly rejected by the unit checks.

Windows initially ran the flash suite successfully, then application control
blocked an updated test executable (OS error 4551). Windows workspace Clippy
also failed to load `bytemuck_derive` in desktop dependencies. WSL checks are
used for the complete local run; no Windows policy was changed. The first
live Julia attempt generated its candidate but could not launch Cargo because
the shell's PATH omitted it; the corrected launcher reruns the full standard
command. Neither failed attempt is counted as a passing check.

The final post-guard workspace and live Julia runs both pass. The live cubic
maximum deviations remain 4.66e-9 K in PH temperature, 3.06e-11 in phase
fraction, 5.71e-7 J/mol in PH enthalpy and 2.51e-12 J/(mol K) in homogeneous
cp. No fixture or tolerance was replaced to obtain agreement.

Hosted CI for baseline main passed on both platforms in
[run 42494429639](https://github.com/alexdamps-stack/engineering-property-estimation/actions/runs/42494429639).
Hosted results for this change and PR review remain recorded on its GitHub PR.

## Original A4 implementation checks

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
