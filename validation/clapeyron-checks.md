# Clapeyron standard numerical check

The independent Julia comparison is part of the standard PR/PR-BM correctness
checks. It supplements the existing teqp/Cantera references, thermodynamic
identities and solver robustness tests.

## Run and CI behavior

Install Julia **1.11.7** and the repository's Rust toolchain. From the repository
root, run:

```sh
julia --project=validation/julia validation/julia/check.jl
```

The script instantiates `validation/julia/Manifest.toml`, verifies the committed
fixture's input, generator and environment hashes, generates a temporary Julia
candidate and runs `cargo test -p thermo-solvers --test clapeyron --locked`
against that candidate. It inherits `CARGO_TARGET_DIR` and other environment
settings. Instantiation needs network access on the first run; the Rust engine
has no Julia dependency. Julia's first compilation can take several minutes.

The check runs in an Ubuntu job on every push/PR and in the monthly/manual
numerical workflow. Ordinary Windows/Linux workspace tests always compare
the committed `validation/fixtures/clapeyron-v1.json` offline. A live generation
error, changed provenance, nonfinite result or failed comparison fails the
job. No case is skipped, and no committed expected value is updated by CI.
Hosting these jobs does not by itself configure repository branch protection.

Reference regeneration is a separate, explicit action:

```sh
julia --project=validation/julia -e 'using Pkg; Pkg.instantiate()'
julia --project=validation/julia validation/julia/generate.jl --output candidate.json
```

Review numerical changes, provenance and supported scope before adopting a
candidate as the committed fixture. Never use Rust results as expected values
or widen tolerances simply to make a comparison pass.

## Matching formulations and references

The environment locks Clapeyron.jl **0.6.30**, source commit
`3440a15ba224664cac5b24441ac8adacbaf326f9`, and its transitive dependencies.
The generator uses native Clapeyron `PR`, `PRAlpha`/`BMAlpha`, `vdW1fRule`
and `NoTranslation`. It supplies all component constants and pair interactions
from the repository inputs, including explicit zero covolume interactions.
It checks exact PR critical factors and the gas constant before evaluation.

For the 200..600 K comparison range, Clapeyron's `ReidIdeal` polynomial exactly
represents the five low-range NASA7 heat-capacity coefficients multiplied by
R. Clapeyron integrates and differentiates this polynomial itself. Its
`:ideal_gas` reference sets each pure ideal enthalpy and entropy to zero at
298.15 K and 100000 Pa; assertions verify that anchor for every model.
The NASA formation/entropy integration constants cancel in this convention.
No fitted approximation to heat capacity is introduced.

Clapeyron's residual entropy at fixed T,V differs from the engine's departure
at fixed T,p. The adapter therefore subtracts Clapeyron ideal entropy at the
same T,p,x from its total entropy. This matches the engine's residual entropy
definition (equivalently, the T,V residual plus R ln Z). Enthalpy has no ideal
pressure dependence and needs no analogous correction.

The 28 PH targets come from the independent teqp/Cantera fixture. They provide
only pressure, feed composition, enthalpy target, variant and parameters to
Clapeyron's native `ph_flash`: no reference temperature or phase-composition
solution is supplied as an initial guess. The comparison includes vapor,
liquid, natural-gas splitting, pure propane/water coexistence endpoints and
interior fractions, a narrow-boiling binary and an isenthalpic valve outlet.
Phase ordering is normalized by molar volume before comparison.

At saturation endpoints Rust retains an incipient zero-fraction phase while
Clapeyron can return only the material phase. Material-phase counts and full
properties are compared after excluding fractions no larger than the declared
`2e-8` phase-fraction tolerance. The actual PH-returned volume selects the
property root; a new minimum-Gibbs request would be ambiguous at coexistence.

The 60 homogeneous inputs are the existing A3 grid and diagnostic cases. Both
mechanically stable outer roots are evaluated when present. The comparison
covers Z, molar volume, mass density, fugacity coefficients, total and residual
enthalpy/entropy, cp, cv and sound speed. It also checks phase states in the PH
results, together with temperature, phase count, fractions and compositions.

The two 0.01 Pa dilute-limit cases use Clapeyron's own generic arithmetic at
256-bit `BigFloat` precision, then round its outputs to Float64 for storage.
At that pressure its Float64 cubic volume calculation loses precision (the
initial PR comparison differed by `8.56e-7` in Z). All other cases use Float64.
The fixture records precision per homogeneous case. This changes reference
evaluation precision, not the EOS, input parameters, Rust engine or tolerances.

The Julia adapter never invokes Rust. Shared parameters and NASA coefficients
are deliberate; derivatives, thermodynamic evaluation and equilibrium solves
come from Clapeyron. Agreement establishes implementation verification for
these matched cases, not experimental accuracy or a proof of global stability.
Default Clapeyron database results and other model formulations are separate
model/data comparisons and are not covered by this check.

Upstream primary sources consulted:

- [PR constructor and parameters](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.30/src/models/cubic/PR/PR.jl)
- [Boston–Mathias alpha](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.30/src/models/cubic/alphas/BM.jl)
- [Reid polynomial ideal contribution](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.30/src/models/ideal/ReidIdeal.jl)
- [Reference-state implementation](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.30/src/database/params/ReferenceState.jl)
- [Native PH solver](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.30/src/methods/property_solvers/multicomponent/flash/PH.jl)

## Tolerances

Budgets were selected before inspecting Julia/Rust deviations, reusing the A3
and A5 regular-state budgets. Comparisons use `atol + rtol * abs(reference)`.

| Quantity | Absolute tolerance | Relative tolerance |
|---|---|---|
| PH temperature | 2e-6 K | 0 |
| Molar phase fraction | 2e-8 | 0 |
| Phase composition | 2e-8 | 1e-7 |
| Z | 2e-10 | 2e-9 |
| Molar volume | 1e-12 m³/mol | 2e-9 |
| Mass density | 1e-7 kg/m³ | 2e-9 |
| Total/residual enthalpy | 2e-5 J/mol | 2e-9 |
| Total/residual entropy | 2e-7 J/(mol K) | 2e-9 |
| cp, cv | 2e-5 J/(mol K) | 2e-8 |
| Sound speed | 2e-4 m/s | 2e-8 |
| ln fugacity coefficient | 2e-9 | 2e-9 |

The Rust PH energy residual must additionally meet `1e-10`, scaled by
`max(abs(h_target), 1 J/mol)`. Generation rejects Julia PH results whose scaled
energy residual exceeds `1e-8`; the actual comparisons still apply the tighter
quantity-specific budgets above. These are numerical implementation budgets,
not physical-model accuracy estimates.

## Local verification record, 2026-11-08

Environment: WSL Ubuntu, pinned Rust toolchain, Julia 1.11.7, Clapeyron.jl
0.6.30 and JSON3 1.14.4, with the committed Julia manifest. The Julia archive
was downloaded from the official Julia distribution and checked against the
SHA-256 published in its versions index.

All 28 PH and 60 homogeneous reference cases were generated without Rust
calls. The offline Rust comparisons pass with the original budgets. Maximum
absolute deviations are:

| Quantity | Maximum deviation |
|---|---|
| PH temperature | `4.66e-9 K` |
| PH phase fraction | `3.06e-11` |
| PH molar enthalpy | `5.71e-7 J/mol` |
| Homogeneous mass density | `4.78e-12 kg/m³` |
| Homogeneous molar enthalpy | `2.01e-10 J/mol` |
| Homogeneous cp | `2.51e-12 J/(mol K)` |

Individual phase compositions and the other listed caloric, fugacity and
volumetric quantities pass their budgets as well. The component source,
input fixtures, generator and Julia manifest hashes match the committed
reference provenance. A temporary candidate with a deliberate `+0.01 K`
error was rejected by the Rust PH comparison with a nonzero exit; the
committed fixture remained unchanged. Both workflow YAML files were parsed
and their push/PR/manual/scheduled triggers and live-check commands inspected.

The initial Float64 dilute-limit failure and the T,V versus T,p residual
entropy convention were resolved in reference evaluation as documented above.
No production Rust EOS/solver formulas, component data or comparison tolerances
were changed to obtain agreement.

| Final local check | Result |
|---|---|
| Standard `validation/julia/check.jl` command | passed: fresh Julia candidate, both Rust comparison tests, committed fixture unchanged |
| `cargo test --workspace --locked` | 98 tests passed, including 13 headless UI tests and both Clapeyron tests |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| `cargo fmt --all -- --check` | passed |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| Deliberate incorrect-reference guard | nonzero failure as expected for a +0.01 K deviation |
| Fixture input/generator/manifest SHA-256 checks | passed |
| Workflow YAML parsing and trigger/command inspection | passed |
| `git diff --check` | passed |

The live check reproduced the same maximum deviations as the committed
fixture. Hosted GitHub Actions and an independent code review have not been
executed in this local session.

## A7 NRTL and gamma–phi check

`check.jl` also verifies `fixtures/nrtl-v1.json` (input, generator and manifest
hashes), regenerates it live with `generate_nrtl.jl` and runs
`cargo test -p thermo-solvers --test nrtl_clapeyron` with
`THERMO_NRTL_REFERENCE`. Clapeyron's `NRTL` and `aspenNRTL` supply g^E and
activity coefficients. Its `PolExpSat`/`PolExpLiquid` evaluate the transcribed
CoolProp ancillaries, and `ReidIdeal` evaluates the NASA7 ideal gas.
Temperature derivatives use ForwardDiff through those functions.

Clapeyron 0.6.30's native gamma–phi VLE (`GammaPhi` with a `FluidCorrelation`
fluid) fails with a missing `a_res` method. The generator therefore solves the
modified-Raoult equations (`y_i p = x_i γ_i p_i^sat`, ideal-gas vapor, no
Poynting factor) in 256-bit BigFloat. Pure liquids use Clausius–Clapeyron,
matching the engine's declared formulation. Tolerances were fixed in the
generator before comparison. One budget was corrected before the fixture was
adopted: phase molar volume went from 1e-12 to 1e-9 relative, matching the
composition budget it depends on. The pure density correlation keeps 1e-12.
Details are in [A7 checks](a7-checks.md).
