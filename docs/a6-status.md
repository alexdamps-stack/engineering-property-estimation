# A6 native GERG-2008 implementation and scope

A6 model implementation and local numerical qualification are complete for the
declared model/property/solver domains. The acceptance evidence and external
review status are recorded in [A6 checks](../validation/a6-checks.md).

The native Rust `thermo_models::GergModel` evaluates fixed-composition GERG-2008
Helmholtz energy at temperature and molar density. It includes all 21 GERG
components, pure-fluid residual terms, composition-dependent reducing functions,
all binary parameters and departure models, and the native ideal-gas functions.
The runtime has no C++, Python, Julia, teqp, REFPROP or FFI dependency.
Native pressure inversion, stability assessment, TP flash and PH flash now
connect this model to the solver and CLI layers. The supported component and
state domains differ by capability; they are explicit below.

The [pinned NIST source and terms](../data/source-extracts/gerg/README.md) are
reviewable alongside the derived tables. Native GERG uses `R = 8.314472 J/(mol K)`
and its own masses and ideal functions. It does not substitute the cubic
NASA7 contribution or the shared modern gas constant. Integration constants
retain NIST's rounded pure-ideal reference near h=s=0 at 298.15 K and 101325 Pa.
Its enthalpies must not be interchanged with the cubic package's reference.

## Properties and derivatives

`helmholtz` returns dimensionless ideal/residual Helmholtz quantities and
scaled derivatives `delta*ar_delta`, `delta²*ar_deltadelta`, `tau*ar_tau`,
`tau²*ar_tautau`, and `delta*tau*ar_deltatau`, plus ideal temperature derivatives.
Each residual term uses analytic logarithmic slopes and curvature. The ideal
path includes mixing entropy with the exact zero-fraction limit. Reducing
temperature/density and molar mass are available for inspection.

`evaluate` computes pressure, compressibility, mass density, internal energy,
enthalpy, entropy, Helmholtz and Gibbs energies, cv, cp, speed of sound,
`dp/dT|rho` and `dp/drho|T`. All are molar SI, with density in mol/m³.
States identify model/data versions, caloric reference and unassessed physical
validity. Fixed arrays and immutable per-model data replace upstream mutable
global calculation caches. Composition order is canonicalized by component ID.

`fugacity` adds analytic residual composition derivatives and returns ln(phi)
in the original input order, including retained zero components. Derivatives
include both the explicit pure/departure weights and the composition dependence
of the reducing temperature and volume. The chemical-potential transformation is
`mu_i^r/(RT) = ar + ar_d + ar_xi - sum(x_j*ar_xj)` and
`ln(phi_i) = mu_i^r/(RT) - ln(Z)`. Production composition derivatives use no
finite differences. The returned `GergFugacityState` retains the homogeneous
state's unassessed physical validity.

## Pressure inversion and explicit root selection

`thermo_solvers::gerg_density` accepts a GERG model, T, absolute pressure and
`GergDensityOptions` with explicit lower/upper molar-density bounds. There is no
default bracket or implicit phase choice. Bounds must be ordered and inside
1e-3..=30000 mol/m³; temperature and target pressure retain the model limits.
Endpoint residuals must bracket the target, allowing a converged endpoint.

The solver uses analytic pressure slopes in safeguarded Newton steps, with
bisection whenever the proposed step is outside the central 80% of the updated
bracket. Every visited density must have positive isothermal stiffness. Final
states additionally pass the model's pressure and cv checks. Trial pressures
may be negative or exceed 35 MPa while bracketing a supported target; this does
not expand the final-state pressure domain. A sampled unstable state fails
explicitly rather than being silently skipped.

Acceptance requires `abs(p_eos-p_target)/max(p_target,1 Pa)` to satisfy the
requested tolerance (default 1e-10, permitted positive values <=1e-7). A small
step or collapsed interval alone never establishes convergence. Default maximum
iterations are 100; invalid options, unbracketed targets, encountered unstable
states, interval collapse and exhausted iterations return explicit errors.
The result reports the signed pressure residual, iteration count, initial and
final density brackets, model/reference metadata and diagnostic.

A bracket selects a homogeneous root. Neither positive stiffness nor a pressure
root proves physical phase identity, root uniqueness, completeness of roots,
monotonicity over the entire interval or equilibrium stability. Separate narrow
brackets can select low- and high-density roots at identical T/P; metastable
roots are intentionally possible. Near-critical and general automatic root
enumeration remain unqualified.

## Stability and equilibrium integration

`gerg_outer_roots` samples a bounded isotherm on logarithmic and linear density
grids, refines pressure crossings, and keeps the exterior branches outside the
sampled negative-stiffness regions. Multiparameter EOS can contain additional
nonphysical loops inside their two-phase regions (also discussed in this
[NIST phase-equilibrium study](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=930138)); these interior candidates
are excluded before enforcing final exterior-root convergence. An exterior
root must pass the native pressure, mechanical and caloric checks. Minimum
Gibbs selection compares only the admitted exterior roots. Finite sampling
does not establish root completeness or critical-region validity.

`gerg_assess_stability` and `gerg_tp_flash` use the existing A4 TPD and phase-split
algorithms through a private EOS adapter. Public cubic APIs remain unchanged.
Wilson vapor/liquid and near-pure-component starts initialize the finite TPD
search; pinned A2 critical constants and acentric factors are used only for
these guesses and the single-phase label heuristic. Native GERG supplies every
pressure, chemical potential, fugacity and caloric value. Unresolved trials
propagate an error. A split must satisfy fugacity equality, component balances,
Gibbs decrease and post-split finite-trial stability. No PR calculation is used
as a fallback for GERG failure.

`gerg_ph_flash` uses native GERG enthalpy and a safeguarded secant/bisection
temperature solve around equilibrium TP evaluations. Its explicit temperature
bracket defaults to 230..350 K; inner TP failures propagate. Acceptance requires
scaled energy error `abs(h-h_target)/max(abs(h_target),1 J/mol)` <=1e-10 by
default, with permitted positive tolerances <=1e-6. Interval collapse alone
never establishes convergence. Pure coexistence below 0.95 Tc uses a bracketed
liquid/vapor log-fugacity equality solve and an enthalpy lever rule, retaining
incipient phases at beta=0/1. Pure TP coexistence without an energy specification
returns an explicit unsupported/underdetermined result rather than inventing a
phase fraction. Retained zero components are supported in this path.

## Independent verification

The offline [fixture](../validation/fixtures/gerg-v1.json) comes from unmodified
NIST C++ 2.01, executed independently of Rust. It covers 273 T-density states:
two per pure component, every one of the 210 binary pairs, 20 natural-gas
states and a full 21-component state. Every state checks ideal/residual first
and second derivatives and reducing variables. 269 states also pass local
admissibility and compare 12 properties; the four locally unstable states
require explicit failure from the property API.

Budgets fixed before Rust comparison are `2e-8 + 2e-10*abs(reference)` in the
reported SI units for each property, reducing quantity and dimensionless
derivative. This conservative double-precision reproduction budget is not an
experimental accuracy claim. Tests also cover finite-difference pressure and
energy derivatives, Maxwell's entropy-density relation, energy identities,
component-order invariance, retained-zero pure limits, dilute behavior,
unsupported components, mass-basis rejection and computational limits.

The additional [teqp AD check](../validation/check-gerg-teqp.py) independently
checks all 273 ideal/residual derivative states using teqp 0.23.2. Source review
found that teqp recalculates ideal integration constants for exact h0=s0=0,
whereas NIST retains rounded constants. The check maps only ideal offsets by
`A+B/T`, determined at 298.15 K from NIST constants, before comparison. This
changes neither cp/cv nor residual derivatives. Published rounded reduced
ideal characteristic temperatures also produce small differences. No Rust
parameters, committed cubic fixtures or tolerances were changed to match teqp.

See [executed checks](../validation/a6-checks.md) for commands and outcomes.

The separate [density fixture](../validation/fixtures/gerg-density-v1.json)
contains 42 independently solved TP roots: all 21 pure fluids, 16 natural-gas
states, low/high-density ethane and propane pairs, and a 21-component mixture.
The pinned NIST density routine supplies one reference root; teqp 0.23.2 AD
pressure and SciPy Brent solving independently verify the same root. Explicit
asymmetric brackets avoid initializing Rust at an already solved midpoint.
Quantity-specific budgets and generator/source hashes accompany the fixture.
The reference-only dense-propane solve needs a close positive-pressure initial
guess for NIST's logarithmic-pressure Newton routine; this is recorded in the
generator. Rust uses its independently implemented safeguarded bracket solve.

The [fugacity fixture](../validation/fixtures/gerg-fugacity-v1.json) contains
45 teqp AD states, including infinite dilution and trace components. Tests
also verify the residual Gibbs/Euler relation, pure-fluid limit, original
component ordering, retained zeros and finite differences of total residual
Helmholtz energy at fixed T,V with two step sizes. Both new generators default
to read-only reproduction; writing requires an explicit candidate output path.

The [equilibrium fixture](../validation/fixtures/gerg-equilibrium-v1.json) adds
13 TP cases and 15 PH targets: regular single and split binary/natural-gas
states, retained zero components, a narrow-boiling ethane/propane mixture, pure
propane coexistence at 280/300 K with beta 0/0.25/0.75/1, and a natural-gas
isenthalpic valve. NIST C++ supplies independent pressure and native calorics;
teqp supplies AD chemical potentials; SciPy supplies multi-start L-BFGS-B TPD
minimization, trust-region least squares and hybrid-Newton phase-split polishing.
Valve temperature is independently solved with Brent's method. Rust is never
called during generation. Pressure-root sampling is independently denser than
the native search. Quantity-specific budgets precede Rust comparisons.

The independent TPD search retains failed reference trajectories explicitly in
`reference_tpd_trajectories`. Exceptions produce no candidate. A finite negative
TPD evaluated on an unconverged trajectory can still demonstrate instability;
unconverged positive values do not establish stable termination. Each recorded
TPD search also has successful minimizations. Some sampled CO2-rich compositions
have no admissible exterior root inside the bounded search. This is finite
regular-branch evidence, not an asserted global minimum. Rust independently
fails if one of its required trials encounters an unsupported state.

`validation/julia/check-gerg.jl` additionally evaluates native Clapeyron 0.6.30
EOS derivatives on all 13 independent TP and 15 PH equilibrium states. It checks
pressure, calorics, active-component chemical-potential equality, balances and
PH energy residuals on the physically admitted exterior branches. Two binary
cases also pass fresh native Clapeyron TP solves with Wilson guesses and no
expected solution supplied. The package's general GERG flash initialization
produced nonfinite states for other attempted cases, so it is not represented
as an independently qualified general GERG flash reference. Complete native
Rust TP/PH solves are compared with the independently executed NIST/teqp/SciPy
solutions above.

The Julia validation adapter explicitly maps constant ideal h/s zeros at the
same native ideal reference density to NIST's rounded convention. It changes
no EOS, Cp, density or fugacity parameters and no comparison budgets. This
live GERG check supplements the mandatory existing PR/PR-BM Clapeyron check,
which still executes fresh native Clapeyron homogeneous and PH calculations.

## Supported domains and remaining gates

The computational envelope is deliberately bounded to 200..=450 K and
1e-3..=30000 mol/m³. Property evaluation additionally requires positive pressure
<=35 MPa, positive `dp/drho|T` and positive cv. These are arithmetic and local
admissibility checks, not a claim of physical validity for every composition
in this rectangle. Numerical evidence applies to the fixture states and tested
regular branches. The low-level derivative API intentionally permits unstable
homogeneous states for diagnostics.

All 21 component identities are listed by `GergModel::supported_components()`.
Methanol and unknown IDs return `MissingParameter`, even at zero fraction.
Mass input requires an explicit conversion; GERG molar masses can be obtained
from pure models. No implicit aliases or interaction defaults are provided.

The CLI supports `model: GERG-2008` with `calculation: homogeneous-trho`
([T-density case](../examples/gerg-natural-gas.json)) or `homogeneous-tp`
([pressure inversion case](../examples/gerg-tp-natural-gas.json)). The TP mode
requires `pressure_pa` and `density_bracket_mol_per_m3`; density is an output.
The T-density mode rejects pressure/bracket inputs. Both modes return ln(phi)
and native caloric-reference metadata; TP also reports convergence evidence.
The strict schema rejects mixed or ambiguous input specifications and cubic
interaction matrices. GERG `tp-flash` and `ph-flash` select exterior roots through
the native equilibrium adapter; they reject user density inputs. PH requires
`enthalpy_j_per_mol` on the native GERG reference and optionally accepts
`temperature_bracket_k`; temperature is an output. See the
[binary TP example](../examples/gerg-tp-flash-binary.json) and
[natural-gas PH valve example](../examples/gerg-ph-flash-natural-gas.json).

| Capability | Component domain | Computational bounds and evidence |
|---|---|---|
| Homogeneous T-density and analytic fugacity | All 21 native components | 200..450 K, 0.001..30000 mol/m³; positive p<=35 MPa and local checks for properties; 273 derivative / 269 property states, 45 fugacity states |
| Explicit-bracket pressure inversion | All 21 native components | Same final-state bounds; caller selects density interval; 42 independent density roots |
| Exterior roots, TPD, TP flash | Methane, ethane, propane, nitrogen, CO2 | 230..350 K, 0.1..10 MPa, same density limits; finite-search evidence at the 13 reference TP states and regular regression cases |
| PH flash | Same dry five-component fluid domain | Explicit 230..350 K bracket; native enthalpy; pure coexistence below 0.95 Tc; 15 independent PH targets |

The rectangular limits are computational envelopes, not a guarantee that every
state or every trial composition inside them is supported. Failed root searches,
stability trials, phase splits and energy solves return explicit errors. A
single-phase label uses reducing density/temperature and a pseudo-critical
pressure heuristic; split vapor/liquid labels use phase molar volume. These
labels do not establish experimental phase-boundary accuracy.

Water/aqueous hydrocarbon equilibrium, multiple-liquid/VLLE or solid phases,
heavier-component equilibrium, near-critical states, general automatic root
completeness, phase envelopes (A13), an expanded envelope and experimental
accuracy remain outside the qualified A6 solver scope. All 21 homogeneous
component models remain available. The native model and declared solver path
are implemented; local completion evidence is recorded in the check report.
Hosted Windows/Linux CI execution and external technical review are not claimed
by the local record and remain release/review gates.
