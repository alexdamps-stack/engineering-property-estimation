# A7 checks: NRTL and liquid properties

Date: 2026-11-08. Starting revision: `cfdcb01`. These results apply to the local A7
change on the pinned Rust toolchain in WSL Ubuntu, as A5 did. Julia 1.11.7 with
the locked Clapeyron 0.6.30 environment and Python 3 (standard library only, for
source transcription) ran in the same WSL environment.

| Check | Outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| `cargo test --workspace --locked` | 126 passed: 12 NRTL model, 9 gamma–phi solver, 6 NRTL reference and 1 new CLI test |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| `julia --project=validation/julia validation/julia/check.jl` | passed in 206 s: PR/PR-BM (60 homogeneous, 28 PH, deviations unchanged from A5) and NRTL live candidates |
| `python3 -I validation/prepare-nrtl-source-data.py` | 2 pure-liquid and 1 NRTL transcription; thermo/ChemSep cross-check passed |
| `thermo-cli examples/nrtl-tp-flash-methanol-water.json` | two phases, β = 0.27685, ln-fugacity residual `3.2e-12` |

Maximum deviations against the independent Clapeyron references (committed
fixture and live candidate are identical):

| Quantity | Max abs. deviation | Largest budget fraction |
|---|---|---|
| ln γ (binary, ternary) | `2.2e-16` | `1.0e-4` |
| ∂ln γ/∂T | `1.2e-18 1/K` | `3.9e-6` |
| g^E, h^E | `3.4e-13 J/mol` | `6.2e-5` |
| c_p^E | `1.4e-15 J/(mol K)` | `8.2e-9` |
| p^sat | `6.1e-9 Pa` | `8.2e-3` |
| Pure-liquid h | `3.3e-11 J/mol` | `4.3e-7` |
| Saturated-liquid density | `7.7e-9 mol/m³` | `0.35` (1e-12 relative budget) |
| Bubble/dew T | `2.7e-10 K` | `2.7e-2` |
| Bubble/dew p | `7.0e-10 Pa` | `2.4e-5` |
| Flash compositions | `2.1e-12` | `3.9e-3` |
| Flash phase h | `1.4e-7 J/mol` | `2.7e-3` |
| PH temperature | `1.2e-9 K` | `0.12` |
| PH vapor fraction | `6.8e-11` | `6.8e-2` |

Generator tolerances were set before comparison. When the fixture was first
generated, a two-phase liquid molar volume differed by `1.6e-12` relative,
against an initial `1e-12` phase-volume budget. That budget was below the 1e-9
relative composition budget that the volume inherits through the solved
liquid composition. It was corrected to 1e-9 before the fixture was adopted.
The pure correlation density keeps its 1e-12 budget. No other budget changed,
and no reference value was taken from Rust.

The live check first ran the NRTL generator with `--compile=min`, like the
PR/PR-BM generator. Its 256-bit BigFloat solves had not finished after 15
minutes, which would exceed the CI timeout. `check.jl` now runs only the NRTL
generator compiled (about 1 minute). The algorithms are unchanged.

Hosted CI and independent code review were not run. These checks establish
local implementation agreement within the documented scope. They do not
establish experimental accuracy, global liquid stability, or flowsheet
integration.
