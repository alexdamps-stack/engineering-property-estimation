# A1 and A2 formal acceptance review

Review date: 2026-11-08. Baseline: `b477020` (main, including PRs #72–#75).
This review addresses GitHub #8/#5 and #15/#16/#9. Historical worker reports
remain available; their pending foundation gates are superseded here.

## A1 decision

Accepted for the foundation scope in issue #5. A0 explicitly authorizes fresh
implementation, with no recovered source imported. The workspace, pinned Rust
toolchain, formatting, Clippy, documentation and cross-platform CI are present.
The calculation contract and implementation have been reviewed together:

- Finite SI quantities, absolute pressure/temperature, explicit composition and
  normalization policies, retained zeros/traces and identity/order checks.
- Hand-calculated mass/mole and flow conversions. Standard-volume flow requires
  named numeric conditions and an explicit ideal/Z assumption; no implicit Nm³.
- Stable input/calculation failure variants; nonfinite and range failures do
  not return successful NaNs or silently lose positive contributions.
- Untranslated EOS liquid density, separate model/data/caloric metadata and
  a shared exact SI gas constant. NASA7 source offsets are retained while ideal
  h/s are shifted to the explicit 298.15 K / 100000 Pa convention.
- Dedicated quantity types and f64 are justified for A1. Analytical derivatives
  and scalar root isolation are recorded with A3; no general matrix library
  is required by the types. Subsequent solver choices remain with their owners.

`crates/thermo-types/tests/public_contract.rs` covers public boundaries, unit
conversions, normalization, mass/mole conversion, ordering/permutation,
zero/trace behavior, standard-volume assumptions and arithmetic range errors.
These exercise observable behavior and hand-derived expectations. The A2 data
identities, bases and reference conventions are reconciled with this contract.

The historical Windows error 4551 is superseded by successful hosted
Linux/Windows checks, including [main CI run 42488612278](https://github.com/alexdamps-stack/engineering-property-estimation/actions/runs/42488612278)
on `b477020`. No numerical model accuracy follows from type/CI acceptance.

## A2 decision

Accepted for the data and fixed reference-case scope in issue #9. Eight
canonical components have source-pinned constants, retained NASA7 correlations
and fit ranges, notices and caloric-reference provenance. Exact PR and exact
PR-BM selections are kept separate from rounded-factor/default-data models.
Independent teqp 0.23.2 / Cantera 3.2.0 and locked Julia 1.11.7 / Clapeyron
0.6.30 references document formulation, ideal-reference and parameter matching.
Sources, software versions and numerical independence are recorded in the
A1–A3, A4, A5, A6 and Clapeyron verification records.

Issue #15 is satisfied by [the stream qualification](a2-process-qualification.md):
independent states for both feeds, mixer, cooler, valve, separator vapor and
separator liquid, for each cubic variant. These qualify the exact fixed inputs
and energy solutions, not a continuous phase diagram over the search envelope.
The new Rust regression also checks phase splits and the associated balances.

Issue #16 is satisfied by `validation/tolerances.json` and the updated testing
strategy: existing pre-comparison cubic, TP, PH, process, Clapeyron and GERG
budgets are consolidated without widening. Each entry identifies its fixed
reference set, domain and rationale. The registry/provenance checker runs in
CI. Preparation proposals are archived explicitly as historical material.
Bubble/dew and phase-envelope budgets remain deferred to A13 (#66), as required
by issue #16; they do not prevent this bounded A2 acceptance.

A2 acceptance is implementation verification with matched shared input data.
It supplies no experimental accuracy claim, global stability proof or solid/
multi-liquid/critical-region qualification. A8 must still demonstrate the
connected unit-operation process; A4's remaining acceleration gate is separate.

## Verification of this change

Executed locally in WSL Ubuntu on the pinned Rust toolchain, with the locked
Julia 1.11.7 environment and pinned Python reference dependencies:

| Check | Result |
|---|---|
| Independent process generator | Passed: 14 fixed stream states, no Rust calls |
| `cargo test --workspace --locked` | Passed, including the new process regression |
| `cargo test -p thermo-solvers --test process_reference --locked` | Passed after final test cleanup |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | Passed |
| `RUSTDOCFLAGS='-D warnings' cargo doc --workspace --no-deps --locked` | Passed |
| `cargo fmt --all -- --check` | Passed |
| `julia --project=validation/julia validation/julia/check.jl` | Passed: fresh PR/PR-BM and NRTL candidates, all comparisons |
| `python validation/consolidate-tolerances.py` | Passed: budgets match fixtures and process provenance hashes match |
| `powershell -NoProfile -File validation/check-preparation.ps1` | Passed: historical preparation artifacts retained |
| `git diff --check` | Passed |

No committed pre-existing numerical fixture or comparison budget was replaced
or widened. Clippy initially flagged a range-index loop in the new regression;
the final iterator-based balance check passes. The final formatting check was
run directly after a temporary shell script's CRLF argument error; Clippy and
documentation had already completed successfully in that script.
