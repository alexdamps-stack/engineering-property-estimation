# A5 checks: PH flash and energy foundations

Date: 2026-11-07. Starting revision: `ccb3cdf`. These results apply to the local
A5 change on the pinned Rust toolchain in WSL Ubuntu. Windows application
control blocked the serde build script (OS error 4551); no OS policy changed.

| Check | Outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| `cargo test --workspace --locked` | 96 passed, including 13 headless UI tests |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| A5 solver tests | 6 passed, including all 28 independent PH reference states |
| CLI tests | 3 passed, including PH recovery and schema failures |
| `generate-ph-references.py` with pinned requirements | 28 reference states, all assertions passed |
| Fixture provenance | component source and all three adapter SHA-256 hashes match |
| Standalone PH CLI example | 225.416724931 K, two phases, scaled energy residual `7.08e-11` |
| `git diff --check` | passed |

The maximum deviations over the 28 independent states were:

| Quantity | Maximum absolute deviation |
|---|---|
| Temperature | `4.66e-9 K` |
| Total molar enthalpy | `5.71e-7 J/mol` |
| Vapor fraction | `3.06e-11` |

All states satisfy the default scaled energy tolerance `1e-10`. Tests also
compare individual phase enthalpies and compositions against the independent
fixture and check component balances. Expected values were generated using
teqp 0.23.2, Cantera 3.2.0, NumPy 2.3.4 and SciPy 1.16.3 from the pinned
requirements. No expected values were taken from the Rust engine.

Reference selection exposed a limitation of the existing independent A4
MINPACK adapter for a candidate 240 K/5 MPa valve inlet. The final independent
valve case uses the existing A4-qualified 250 K/3 MPa inlet, with an independent
outlet temperature solve at 1.5 MPa. No reference convergence checks were
relaxed, no existing fixture was replaced, and every selected case passes.
The additional Rust consistency test retains the 240 K/5 MPa expansion.

Hosted CI and independent code review were not run. These checks establish
local implementation agreement within the documented scope; they do not prove
global stability, experimental accuracy or complete flowsheet acceptance.
