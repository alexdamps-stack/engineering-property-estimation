# A6 executed checks

Date: November 8, 2026. Completion checks used the pinned Rust 1.99.1 toolchain
in WSL Ubuntu: all 115 workspace tests passed. Windows formatting and Clippy
also passed. Native Windows workspace execution passed 76 tests, including the
new GERG TP/PH CLI test and homogeneous/fugacity tests, then Windows Application
Control blocked the offline `clapeyron` test executable (OS error 4551). An
ordinary rerun encountered the same block; a targeted new density-test binary
was also blocked before execution. No OS policy changes were attempted. This
is an execution restriction, not an assertion failure or a complete Windows
workspace pass. Hosted CI and external technical review were not performed.

The prior homogeneous/density stage's record reported 109 tests passing natively
on Windows after a transient application-control block. That earlier result
does not supersede the current equilibrium-stage execution restriction. No
committed cubic/flash/PH/Clapeyron fixtures or their tolerances were replaced.

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | passed |
| `cargo test --workspace --locked` in WSL | 115 passed, including 13 headless UI tests |
| Native Windows workspace test execution | 76 passed; then OS error 4551 before the offline Clapeyron test executed |
| `RUSTDOCFLAGS=-Dwarnings cargo doc --workspace --no-deps --locked` | passed |
| GERG homogeneous/fugacity and bracketed-density integration tests | 9 passed |
| GERG equilibrium integration tests | 5 passed: 13 TP and 15 PH references, root policy, consistency and failures |
| CLI tests including GERG homogeneous, TP/PH and schema checks | 6 passed in WSL and Windows |
| `python3 validation/prepare-gerg.py` | exact parameter and 273-state fixture reproduction passed |
| `python validation/check-gerg-teqp.py` with pinned reference environment | 273 ideal/residual AD derivative states passed |
| `python validation/generate-gerg-density-references.py` | 42 NIST/teqp/SciPy density roots reproduced without writing |
| `python validation/generate-gerg-fugacity-references.py` | 45 teqp AD fugacity states reproduced without writing |
| Independent GERG equilibrium regeneration to a temporary candidate | exact reproduction of all 13 TP and 15 PH cases, with explicit trajectory outcomes |
| `cargo run -p thermo-cli --locked -- examples/gerg-natural-gas.json` | passed |
| `cargo run -p thermo-cli --locked -- examples/gerg-tp-natural-gas.json` | passed |
| Julia 1.11.7: `julia --project=validation/julia validation/julia/check.jl` | passed: 60 homogeneous and 28 PH cubic states |
| Julia 1.11.7: `julia --compile=min --project=validation/julia validation/julia/check-gerg.jl` | passed: native GERG EOS/AD checks on 13 TP and 15 PH states plus two fresh binary TP solves |
| Modified CI YAML files | parsed successfully with pinned ruamel.yaml |
| `git diff --check` | passed |

NIST reference execution used g++ 15.2.0 and Python 3.14. The teqp AD check used
Python 3.12 and teqp 0.23.2 with the existing pinned reference dependencies.
The generator also supports Python 3.12, used by scheduled/manual CI.
The C++ compiler zero-initializes automatic storage to avoid exporting
indeterminate unused cells from upstream `SetupGERG`. All 273 numerical states
were compared before/after this extraction fix and remained unchanged. Only the
new GERG outputs were written; existing cubic/flash/PH/Clapeyron fixtures and
their tolerances were preserved.

All 273 states pass native ideal/residual Helmholtz and reducing-quantity
comparisons against NIST. 269 pass the local property checks and compare
pressure, Z, both first pressure derivatives, u, h, s, cv, cp, sound, a and g.
The four rejected reference states are n-heptane, n-octane, n-nonane and
n-decane at 450 K and 1000 mol/m³; all have negative `dp/drho|T`. They remain
in the derivative fixture and explicitly test failure from the property API.

After the reviewed ideal-reference offset mapping, maximum absolute teqp/NIST
differences across all 273 states were:

| Dimensionless quantity | Maximum difference |
|---|---|
| residual Helmholtz energy | `2.44e-12` |
| any scaled residual first/second derivative | `5.82e-12` |
| ideal Helmholtz energy | `7.14e-10` |
| scaled ideal first temperature derivative | `3.17e-9` |
| scaled ideal second temperature derivative | `9.28e-9` |

teqp recalculates exact ideal-reference integration constants while NIST keeps
rounded values. The first raw cross-check exposed this mismatch at 200 K.
Source review confirmed the different reference conventions; only the validation
adapter applies an explicit `A+B/T` ideal-reference transformation. Neither
production parameters nor comparison budgets were adjusted. Characteristic
temperature rounding accounts for the remaining small ideal differences.

The CLI case at 300 K and 1000 mol/m³ produces pressure `2377948.03612 Pa`,
Z `0.953336959188`, mass density `17.963346 kg/m³`, cp `39.8944417473 J/(mol K)`
and sound speed `414.730924942 m/s`, with native model/data/reference metadata.

## Pressure inversion and fugacity continuation

The 42 independent density roots pass quantity-specific density, pressure,
caloric, acoustic and pressure-derivative budgets stored in
`fixtures/gerg-density-v1.json`. Maximum native/NIST density difference is
`5.67e-7 mol/m³`; the largest iteration count is 36 with the fixture's
asymmetric brackets. Every returned state satisfies the default scaled pressure
tolerance `1e-10`. Tests additionally cover wide intervals, exact endpoint
roots, permutation/retained zeros, distinct roots at identical T/P, invalid
options, unbracketed pressures, encountered negative stiffness and forced
iteration exhaustion. No phase identity or global stability is inferred.

Analytic composition derivatives and fugacity coefficients pass all 45 teqp AD
states in `fixtures/gerg-fugacity-v1.json`, including inactive and trace species.
Maximum absolute ln(phi) error is `1.20e-13`. The residual Gibbs/Euler relation,
pure limit, component ordering and finite differences of total residual
Helmholtz energy at fixed T,V with two step sizes also pass. Production
composition derivatives are analytic, including the reducing-function terms.

The TP CLI example at 300 K and 10 MPa returns density `4856.29336155 mol/m³`
and Z `0.825541967699` after 17 iterations from `[1000,7000] mol/m³`, with
scaled pressure residual `-6.59e-11`. Its output includes ln(phi), native caloric
reference, density brackets and unassessed phase validity. Schema tests reject
missing/mixed inputs, reversed brackets and unsupported pressure; qualified
GERG flash modes are handled separately without caller-supplied densities.

The dense propane reference at 250 K and 0.1 MPa exposed a limitation of the
NIST logarithmic-pressure density iteration: broad or negative-pressure initial
guesses failed. A close positive-pressure seed from an independent teqp root
search permits NIST comparison; a fresh SciPy bracket verifies the same root.
This seed is recorded explicitly. The Rust safeguarded bracket solver succeeds
from `[12000,13000] mol/m³`, including a negative-pressure lower endpoint.
No reference convergence tests, production parameters or tolerances were relaxed.
The added fixtures and their source/generator/input hashes reproduce read-only
with the pinned Python 3.12 environment (teqp 0.23.2, NumPy 2.3.4, SciPy 1.16.3).
Scheduled/manual numerical CI now runs both reproductions.

The live Julia 1.11.7 / Clapeyron 0.6.30 standard check was rerun after this
continuation and passed with temporary
independent values for all 60 homogeneous and 28 PH cubic states. Maximum
homogeneous enthalpy deviation was `2.01e-10 J/mol`; maximum PH temperature
deviation was `4.66e-9 K`. The committed Clapeyron fixture remained unchanged.
This checks existing cubic behavior, not GERG verification.

## A6 stability, TP and PH completion

The native equilibrium path passes all 13 independent TP cases and all 15 PH
targets in `fixtures/gerg-equilibrium-v1.json`. It covers single vapor/liquid,
regular binary/natural-gas splits, retained zeros, a narrow-boiling ethane/propane
mixture, pure-propane saturation at 280 and 300 K with phase fractions
0/0.25/0.75/1, and an isenthalpic natural-gas valve outlet. Every split meets the
default maximum absolute log-fugacity residual `1e-10`, normalized component
balance `1e-8`, Gibbs decrease and post-split finite TPD checks. Every PH result
meets scaled energy error `1e-10` with its explicit 1 J/mol floor. Temperature,
phase amounts, compositions, density and calorics also pass the fixture's
quantity-specific budgets.

The valve reference independently expands the dry 330 K / 5 MPa inlet to 3 MPa
at native molar enthalpy `424.567733332 J/mol`; the independently solved outlet
temperature is `322.303293256 K`. The executable PH example uses this target.

Additional regressions verify rejection of interior multiparameter roots,
permutation and retained-zero consistency, unsupported components/states,
invalid options, forced TP/PH/stability iteration failure, zero target enthalpy
and its absolute floor, pure TP coexistence without a phase-fraction specification,
and pure PH with retained inactive species. The low-temperature natural-gas
case exposed a pressure-convergence check being applied to an interior candidate
that was subsequently discarded. Convergence is now enforced on the selected
exterior roots; their budgets were not relaxed. Pure PH saturation requires
log-fugacity equality <=1e-12 and uses the native enthalpy lever rule.

The independently implemented reference minimizes modified TPD using SciPy
L-BFGS-B with 18 starts per search (Wilson plus 16 seeded Dirichlet starts), then
solves phase splitting with trust-region least squares and hybrid-Newton polishing.
It records 25 TPD searches / 450 attempted trajectories over the TP cases;
159 trajectories report optimizer or exterior-root failures, and every search
also has successful minimizations. Exceptions contribute no candidate. Valid
negative values from unconverged finite trajectories can show instability, but
failed stable termination is never a proof of stability. These records remain
in the fixture. Native Rust propagates failures from its own required TPD trials.
The full fixture reproduces exactly in a separate temporary candidate.

Fresh Clapeyron GERG EOS/AD calculations independently confirm the pressure,
calorics and active-component chemical-potential equality of all 28 equilibrium
states, plus balances and PH energy residuals. Two binary TP cases additionally
pass fresh native Clapeyron flash solves without supplied solutions. Matching
the ideal reference requires constant h/s shifts determined at 298.15 K and the
same native ideal molar density. Cp and all residual EOS derivatives remain
unchanged; no comparison budgets were widened.

Clapeyron 0.6.30's general GERG flash initialization returned nonfinite states
for some attempted cases (including its optional full-TPD initialization of the
first binary case). Those are backend qualification failures, not successful
independent flashes. The new GERG Julia check therefore verifies the physically
admitted independent states with native EOS/AD and retains only the two verified
fresh binary solves. Complete Rust TP/PH solver verification comes from the
independent NIST/teqp/SciPy solutions. The existing mandatory cubic Clapeyron
standard check continues to execute fresh homogeneous and PH calculations.

Push/PR and scheduled numerical workflows now run the separate live GERG
EOS/AD check. Scheduled/manual Python validation regenerates the equilibrium
reference only into a temporary candidate and compares it with the preserved
fixture. No release publication or hosted CI run was performed locally.

The A6 native model and its declared solver domains now have completed local
implementation and numerical-qualification evidence. Experimental accuracy,
global stability proof, near-critical/aqueous/heavier-fluid equilibrium,
expanded domains and phase envelopes remain outside this qualification. Hosted
CI and external technical review remain release/review gates. See
[A6 scope](../docs/a6-status.md).
