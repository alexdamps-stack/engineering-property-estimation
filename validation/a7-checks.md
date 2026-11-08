# A7 checks: NRTL and liquid properties

Date: 2026-11-08. Starting revision: `cfdcb01`. These results apply to the local A7
change on the pinned Rust toolchain in WSL Ubuntu, as A5 did. Julia 1.11.7 with
the locked Clapeyron 0.6.30 environment and Python 3 (standard library only, for
source transcription) ran in the same WSL environment.

| Check | Outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| `cargo test --workspace --locked` | 128 passed: 12 NRTL model, 11 gamma–phi solver, 6 NRTL reference and 1 new CLI test |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| `julia --project=validation/julia validation/julia/check.jl` | passed in 198 s with the restored 1e-12 volume budget and 1e-13 split tolerance: PR/PR-BM (60 homogeneous, 28 PH, deviations unchanged from A5) and NRTL live candidates |
| Hosted CI on PR #75 (pre-review commit) | Workspace ubuntu/windows and independent Clapeyron check passed |
| PH sweep through the release CLI (after review fixes) | 48 near-zero-enthalpy targets at 1 atm (8 feeds × 6 targets) and 168 TP→PH round trips over 100 Pa–500 kPa, 3 feeds, 274–445 K: no failures |
| `python3 -I validation/prepare-nrtl-source-data.py` | 2 pure-liquid and 1 NRTL transcription; thermo/ChemSep cross-check passed |
| `thermo-cli examples/nrtl-tp-flash-methanol-water.json` | two phases, β = 0.27685, ln-fugacity residual `8.7e-14` |
| Saturation stress test through the release CLI (tolerance 1e-13) | 540 bubble/dew pressure/temperature cases, 1e-6 to 1 − 1e-6 methanol, 700 Pa–500 kPa: none failed; max ln-fugacity residual `4.6e-14` |

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
| Bubble/dew T | `6.3e-13 K` | `6.3e-5` |
| Bubble/dew p | `7.0e-10 Pa` | `2.4e-5` |
| Flash compositions | `5.4e-14` | `9.9e-5` |
| Flash phase h | `3.5e-9 J/mol` | `7.0e-5` |
| Flash phase molar volume | `3.5e-18 m³/mol` | `4.1e-2` (1e-12 relative budget) |
| PH temperature | `5.1e-10 K` | `5.1e-2` |
| PH vapor fraction | `2.8e-11` | `2.8e-2` |

Generator tolerances were set before comparison and none is widened. No
reference value was taken from Rust.

**Phase molar volume (1e-12 relative).** The first comparison failed. The liquid
of the TP state at 345 K, 1 atm, z = (0.7, 0.3) differed by `1.6e-12` relative.
An intermediate commit widened that budget to 1e-9. The review flagged this
against AGENTS.md, and the 1e-9 budget was reverted. The cause was in the
engine, not the budget. Liquid volume is `v = Σ x_i/ρ_i(T)`, so a composition
error δx changes it by about `δx (v_1 − v_2)/v ≈ 0.9 δx`. The split stopped at
an ln-fugacity residual of 1e-11, which left about 2e-12 composition error
against the 256-bit reference. The default `ActivityOptions::tolerance` is now
1e-13. Composition deviations fell from `2.1e-12` to `5.4e-14`. The volume
passes the original budget at 4 % of it, and bubble/dew temperatures moved
from `2.7e-10 K` to `6.3e-13 K`. Regenerating the fixture changed only the
generator hash and that budget value; every reference value is bit-identical.
The tighter tolerance was stress-tested. 540 bubble/dew pressure and
temperature cases ran for feeds from 1e-6 to 1 − 1e-6 methanol and 700 Pa to
500 kPa; none failed, and the largest ln-fugacity residual was `4.6e-14`. The
PH sweep below was repeated.

The live check first ran the NRTL generator with `--compile=min`, like the
PR/PR-BM generator. Its 256-bit BigFloat solves had not finished after 15
minutes, which would exceed the CI timeout. `check.jl` now runs only the NRTL
generator compiled (about 1 minute). The algorithms are unchanged.

## Code review (PR #75)

A high-effort automated code review of `main...alexdamps-stack/a7-NRTL` produced
ten findings. It was not an independent human review. Resolution:

1. **Fixed: PH failed for two-phase targets near h = 0.** The `max(|h_t|, 1 J/mol)`
   scale gave a 1e-10 J/mol tolerance, below the split's attainable precision.
   The reviewer saw 14 of 48 failures at 1 atm. The two-phase mixture range
   now uses the two-phase enthalpy span as its scale floor (see
   `docs/a7-status.md`). The single-phase branches keep the A5 criterion.
   Covered by `ph_two_phase_targets_near_zero_enthalpy_converge`.
2. **Fixed: PH needed bubble and dew temperatures inside the bracket.** For
   example, z = (0.4, 0.6) at 1 kPa failed for a vapor target. Saturation
   temperatures outside the bracket now truncate or remove that branch.
   Covered by `ph_with_saturation_outside_the_bracket_solves_the_reachable_branches`.
3. **Fixed: the phase molar-volume budget had been widened (1e-12 → 1e-9).**
   This breaks the AGENTS.md rule. The budget is back at 1e-12, and the
   engine's split is converged further instead (described above).
4. **Fixed:** a liquid with no interacting pair reported `WithinDeclaredDomain`
   at any temperature. It now reports `Unassessed`.
5. **Fixed:** the CLI provenance now carries the component and liquid/NRTL data
   identifiers and versions. Bubble and dew outputs report the liquid's NRTL
   validity.
6. **Fixed:** the PH two-phase search runs split-only TP flashes, without TPD
   trials, using Illinois regula falsi instead of bisection. The accepted
   temperature is flashed again with the trials.
7. **Fixed:** pure-liquid states are computed once per temperature in the
   bubble, dew and TP solvers. In the follow-up after merge (see below),
   `ActivityModel::liquid` and the new `ActivityModel::pure_liquids` check the
   temperature domain once rather than once per component.
8. **Fixed:** the dew-temperature solve reuses the dew point's excess state and
   its last evaluation.
9. **Fixed:** the shared `stability::normalized` and
   `COMPOSITION_SUM_TOLERANCE` replace the duplicates.
10. **Fixed:** the unused `ActivityModel::ln_k_values` and a redundant TPD
    renormalization were removed.

The post-fix sweep found one more defect. At the exact dew temperature, a
split converged to β = 1 + 5e-15 and was rejected. A vapor fraction within
1e-10 outside (0, 1) is now reported as the saturated single phase
(`A7-SPLIT-AT-SATURATION-BOUNDARY`).

These checks establish
local implementation agreement within the documented scope. They do not
establish experimental accuracy, global liquid stability, or flowsheet
integration.

## Post-merge re-verification and closure

Date: 2026-11-08. Revision: `68ef34a` (main after A6 GERG-2008 #76 and A1/A2
acceptance #77 merged on top of A7 #75), plus the review finding 7 follow-up.
Same WSL Ubuntu environment, pinned Rust toolchain, Julia 1.11.7 and the locked
Clapeyron 0.6.30 environment.

| Check | Outcome |
|---|---|
| `cargo fmt --all -- --check` | passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| `cargo test --workspace --locked` | 147 passed (146 before the follow-up, plus `pure_liquids_match_per_component_states_and_reject_bad_index`) |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| `julia --project=validation/julia validation/julia/check.jl` | passed on the final tree (PR/PR-BM and NRTL live candidates); committed reference fixtures unchanged |

The NRTL deviations and budget fractions are unchanged from the table above
(flash compositions `5.4e-14`; the TP-flash vapor fraction, not tabulated
above, `9.2e-14`). No fixture was regenerated and no tolerance was changed.

Finding 7 follow-up: `ActivityModel::liquid` evaluated each component through
`pure_liquid`, so the temperature and correlation-range check ran once per
component after `check_domain` had already run it. A private evaluator now
assumes a checked temperature. `liquid` and the new
`ActivityModel::pure_liquids(T)`, which the solvers' per-temperature pure-state
cache uses, each validate the domain once. `pure_liquid` with an index outside
the component list now returns `InvalidInput(ComponentMismatch)` instead of
panicking. Numerical results are bit-identical; the offline reference
comparisons are unaffected.

The A7 work-package acceptance (#41) is met within the documented scope:
directed parameters, temperature dependence, the explicit vapor model and
caloric coupling pass formulation-matched Clapeyron references and the
consistency checks, and supported domains and missing-data failures are
explicit (`docs/a7-status.md`).
