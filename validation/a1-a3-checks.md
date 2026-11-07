# A1–A3 continuation checks and review

Date: 2026-11-07. Starting revision: `d6ab6d8`; results apply to the local
uncommitted continuation, not to that unchanged commit. Toolchain on both
platforms: rustc 1.99.1, cargo 1.99.1. Reference interpreter: CPython 3.12.12
in a separate WSL environment, dependencies pinned in
`reference-requirements.txt`. Application binaries have no Python dependency.

| Check | Outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| Windows `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| WSL scoped numerical Clippy with warnings denied | passed |
| Windows `cargo test -p thermo-types -p thermo-data -p thermo-models -p thermo-cli --locked` | passed, 23 tests |
| WSL same scoped numerical tests | passed, 23 tests |
| WSL `cargo test --workspace --locked`, software Vulkan | passed, 72 tests including 13 headless UI tests, plus doctests |
| WSL `PROPTEST_CASES=2048 cargo test -p thermo-models --release --locked` | passed; independent fixtures and fixed-seed sweep |
| Windows `cargo doc --workspace --no-deps --locked`, `RUSTDOCFLAGS=-D warnings` | passed |
| WSL CLI natural-gas and binary JSON examples | both passed; Z approximately 0.9184554903 and 0.8809015715 respectively |
| Explicit reference generation | passed, 60 cubic states and 64 NASA7 states; no failed states skipped |
| `powershell -NoProfile -File validation/check-preparation.ps1` | passed historical source/provisional artifact checks |
| `git diff --check` | passed |

Windows full workspace testing compiled successfully, then OS application
control blocked the headless UI test executable with error 4551 before it ran.
The Windows standalone CLI executable was likewise blocked. These are not
passing Windows execution results or assertion failures. The complete workspace
tests and both CLI cases were executed successfully in WSL instead. No OS policy
was changed. Hosted GitHub CI has not been launched by this local task.

## Independent code review

A separate reviewer inspected production code, reference generation, fixtures,
tests and the status document. It reported two P2 findings:

1. A global Cardano discriminant tolerance misclassified low-pressure vapor
   roots as critical. Replaced with monotone-interval root isolation, scaled
   stationary-root checks and bounded bisection. Regression references now
   include 400 K, 0.01 Pa for both PR and PR-BM; exact critical and forced
   iteration-exhaustion outcomes are tested.
2. Cantera's exact-1000-K branch selection made the adapter's offset correction
   discontinuous. The adapter now evaluates low and high intervals explicitly;
   all eight species have exact-join reference fixtures. The source coefficients
   and comparison budgets were not changed to hide a failure.

The new low-pressure reference also exposed an absolute SciPy density tolerance
too coarse for that density. It now scales with density and independently
asserts the pressure residual. Result provenance additionally includes the
ordered component identities and actual interaction matrix, addressing the
reviewer's optional suggestion.

The reviewer rechecked both fixes, independently probed Cantera at and around
the join, and executed the eight consistency/reference test functions. Final
review outcome: **no remaining material findings**. This is code review, not
experimental validation or full milestone approval.

See [implementation status](../docs/a1-a3-status.md) for the supported
homogeneous scope and the still-open demonstration-phase/PH and full A2/A3
acceptance gates. No commit, staging, push, or deployment was performed.

## Second independent code review

On the user's request, a fresh reviewer examined the current code, fixtures,
generator, dependencies and documentation. The existing nine model tests and
the 23 scoped numerical/CLI/type tests pass. Both earlier P2 fixes remain valid.
This review found **one new P2 issue**, resolved in the following section:

- `crates/thermo-models/src/cubic.rs:381–382`, `residual_helmholtz`: taking the
  logarithm of a ratio of nearly equal volumes loses the attractive term at
  extremely low density. Pure methane PR at 300 K, V=1e30 m³ and n=1e15 mol
  returns +0.06685068622043411 J; independent 80-digit Decimal evaluation
  gives −0.13535394993250723943 J. The public helper accepts these inputs,
  although equivalent pressure is below the separate TP evaluator's envelope.
  Recommended correction: evaluate the equivalent logarithm with `ln_1p`,
  or explicitly bound the supported density, and add a dilute-limit regression.
  The finding was left unfixed during that review-only pass; closure follows below.

The reviewer also explicitly audited the user's fundamental architecture rule:
the engine's production calculation path uses native Rust with embedded data;
there are no external calculation API, network, CLI, subprocess, Python, teqp
or Cantera runtime calls. The project's own CLI consumes the engine in process.
Independent reference tools are outside that path, and ordinary tests read
offline fixtures. The rule is now explicit in `docs/architecture.md`.

## Resolution of all review findings

At the user's request, the remaining P2 dilute-limit error is now fixed. A shared
`attraction_integral` evaluates the logarithm as
`ln_1p(2*sqrt(2)*b/(v+(1-sqrt(2))*b))/(2*sqrt(2)*b)` in Helmholtz, fugacity
and caloric paths. This preserves the attraction contribution when the original
ratio rounds to one. A finite V/n check prevents volume-division overflow from
becoming successful zero energy after this change.

`generate-helmholtz-reference.py` independently evaluates six methane PR/PR-BM
states at 80 decimal digits with Python's standard library only. It covers a
regular volume and both dilute review cases; expected values are committed in
`fixtures/helmholtz-dilute-v1.json`. The new regression test failed on the original
implementation before the fix and now passes with a relative tolerance of
5e-14 and no absolute floor. Existing fixtures and tolerances were not changed.
The generator is development-only; normal tests use the static fixture and the
engine gains no external runtime, API or CLI dependency.

Final checks after the fix:

- Windows scoped numerical/CLI/type tests: **24 passed**.
- WSL full workspace, including UI tests and doctests: **73 tests passed**.
- WSL model release tests with 2048 fixed-seed cases: passed.
- Workspace formatting, Clippy with warnings denied, warning-free documentation,
  and `git diff --check`: passed.

The independent reviewer rechecked the patch, regenerated the original probe
outside the repository and obtained −0.13535394993250718 J for the 300 K,
V=1e30 m³, n=1e15 mol PR case, agreeing with the 80-digit reference. All ten
model test functions passed. Review conclusion: **P2 resolved, no new material
findings**. All findings from both reviews are closed. The remaining A2/A3
phase/process acceptance gates are still open as documented; they are not
unfixed code-review findings.
