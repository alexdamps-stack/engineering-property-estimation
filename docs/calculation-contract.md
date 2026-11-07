# Calculation contract — provisional A1 foundation

Current implementation decisions and check evidence are recorded in
[A1–A3 status](a1-a3-status.md). In particular, A3 now selects exact critical
factors, analytical temperature derivatives, scalar root isolation without a
matrix dependency, and cp-integrated NASA7 reference offsets. The original
foundation and worker history below remain for audit.

This contract defines the dependency-free `thermo-types` API introduced after
scaffold baseline `bf86da9`. It implements validated inputs and conversions,
not thermodynamic models, flashes, property packages or numerical accuracy
claims. On November 7, 2026 the coordinator reported reviewing A0's
explicit-unavailability decision and authorized a fresh Rust foundation.
The candidate homogeneous Rust model code located during A0 is not copied here.
A1 acceptance still requires coordinator review of this contract, the types,
checks, and compatibility with A0 and A2; this document does not accept A1.

## Units and state inputs

All stored values are finite `f64` in SI. Private newtypes protect construction;
there are no unchecked arithmetic operator implementations or implicit defaults.

| Type / quantity | Unit | Accepted range |
|---|---|---|
| `Temperature` | K, absolute | strictly positive |
| `Pressure` | Pa, absolute | strictly positive |
| `MolarMass` | kg/mol | strictly positive |
| `MolarFlow` | mol/s | nonnegative |
| `MassFlow` | kg/s | nonnegative |
| `StandardVolumeFlow` | reference m³/s | nonnegative, with conditions and assumption |
| `PositiveDimensionless` | 1 | strictly positive |
| `FiniteSi` | explicitly documented by caller | signed, finite |

`FiniteSi` is a low-level finite scalar, not dimensional analysis: future public
property outputs must give quantities specific names and units. Molar energy
and entropy will use J/mol and J/(mol K); mass-specific counterparts J/kg and
J/(kg K); density kg/m³ and molar volume m³/mol. Do not mix these bases.
Dedicated SI newtypes are the A1 units decision; no units dependency is needed
for this scope. A later units-library adoption requires a reviewed migration.

Temperature conversion is explicit through `from_celsius`; pressure conversion
through `from_bar_absolute`. Gauge inputs must first provide their ambient
absolute reference upstream. Zero absolute pressure and temperature are rejected;
vacuum and zero-temperature limits are mathematical model limits, not valid
finite state requests. Negative zero is accepted for nonnegative flows/fractions
and rejected for positive quantities. Negative flow directions belong in future
network connectivity, rather than negative physical feed amounts.

NaNs, infinities, detected overflow and positive-to-zero underflow are errors.
IEEE subnormals are otherwise accepted. Extreme conversions can fail even when
an algebraically rearranged calculation might fit: this conservative API never
promises arbitrary dynamic range. Ordinary rounding of representable traces is
permitted; positive fractions must not silently disappear.

## Identities, ordering and composition

`ComponentId` is an opaque, stable, case-sensitive ASCII key, allowing letters,
digits, `_`, `-`, `.`, and `:`. Examples are `methane` and `cas:74-82-8`.
The data layer must choose one canonical key per physical species and qualify
aliases, isomers and pseudo-components. A display name or vector index is not
identity. The type does not resolve aliases or validate registry numbers.

`Composition` always carries `CompositionBasis::Mole` or `Mass`. Its nonempty
ordered `(ComponentId, fraction)` entries are immutable through the public API.
Duplicate identities, negatives, nonfinite values and all-zero mixtures fail.
Individual zero components remain in place; no sorting, pruning, trace cutoff
or implicit parameter reordering occurs. Every supplied molar mass, including
one for a zero entry, must match identity and order exactly. Length mismatch or
misordering returns `ComponentMismatch`. A joint permutation of fractions,
identities and parameters must leave physical results invariant.

Every constructor requires a `NormalizationPolicy`:

- `RequireUnitSum` accepts only an absolute sum error at most `1e-10`.
- `Normalize` explicitly treats entries as relative nonnegative amounts.

Both policies store fractions rescaled to unit sum, subject to floating-point
rounding. Requiring a unit sum therefore permits a documented small correction;
it never accepts arbitrary percentages. Normalization scales by the largest
entry first to avoid overflowing large finite relative amounts. If a positive
contribution becomes zero, construction fails with `NumericalRange`.

Conversion uses supplied positive molar masses `M_i` in kg/mol:

```
w_i = x_i M_i / sum_j(x_j M_j)
x_i = (w_i / M_i) / sum_j(w_j / M_j)
M_mix = sum_i(x_i M_i) = 1 / sum_i(w_i / M_i)
mass_flow = molar_flow M_mix
```

Internally basis conversion scales molar masses to reduce overflow. `molar_mass`
converts to mole basis and evaluates its arithmetic weighted sum. Passing a
mixture molar mass directly to a flow conversion is an explicit caller assertion
that it belongs to that composition. No dataset or molar mass is guessed.

All strictly positive components are active at the type boundary, including
traces. Zero components are retained for bookkeeping and excluded from future
logarithmic residuals; limiting fugacity coefficients for zeros require an
explicit model capability. A future solver's active-set threshold must be
declared and tested rather than silently changing the mixture.

## Named reference-volume flows

Neither Sm³ nor Nm³ has an implicit default. `ReferenceConditions` requires a
nonblank name, positive absolute temperature and positive absolute pressure.
`StandardVolumeFlow` stores these conditions and an explicit
`GasVolumeAssumption`: either `IdealGas` or `CompressibilityFactor(Z)` with a
supplied strictly positive finite Z at those reference conditions for that
mixture. The name labels the convention; the numeric conditions are authoritative.
The caller must qualify phase/applicability; no EOS evaluates or validates Z.

```
V_ref_dot = n_dot Z_ref R T_ref / p_ref
n_dot = V_ref_dot / (Z_ref R T_ref / p_ref)
R = 8.31446261815324 J/(mol K)
```

The constant is the rounded `f64` product of the exact SI defining molar and
Boltzmann constants. No external reference package is used. The explicit hourly
constructor converts m³/h to stored m³/s using 3600 s/h. Changing reference
conditions requires conversion through molar flow and a new explicit assumption.
This is gas reference-volume bookkeeping, not actual process volume, liquid
standard volume, or a phase-equilibrium calculation. For mass flow convert
through molar flow with the supplied mixture molar mass.

## Model, data, reference and failure vocabulary

`ModelMetadata` identifies implementation version and exact formulation choices;
PR and PR-BM must be separate selections, including alpha function, mixing rule
and any translation. `DataMetadata` identifies a dataset version and source;
`CaloricReference` identifies the enthalpy/entropy zero and ideal contributions.
Caloric reference states are separate from gas-volume reference conditions.
These are plain records, not source qualification or validated model selections;
future data/package boundaries must reject absent or inconsistent metadata.

`ValidityStatus` distinguishes within-declared-domain, extrapolated and
unassessed. Within-domain is a caller/model assessment against a separately
specified domain, never inferred from finite numbers. `Phase` labels individual
vapor, liquid or supercritical results without claiming stability. Diagnostics
carry a stable code, severity and human-readable explanation.

`InputError` distinguishes nonfinite, nonpositive, negative, malformed identity,
empty composition, zero total, duplicate component, non-normalized input,
component mismatch and numerical range. `CalculationError` distinguishes invalid
input, missing parameter, unsupported capability, outside validity domain,
nonconvergence and numerical failure. Errors implement `std::error::Error`;
diagnostic text is not a machine-readable substitute for their variants.
Missing interactions must not become zero by default. Unsupported states and
nonconverged iterations must not produce successful values with NaNs.

Future property results must include phase, units, model/data/reference identity,
validity and applicable convergence evidence. This package deliberately does not
introduce generic result containers or speculative solver traits. Iteration
limits, residual definitions, stable diagnostic codes and solver-specific
statistics must be defined alongside actual algorithms. Input edits must
invalidate derived results at the owning application/package layer; no cached
result or mutable input state is introduced here.

## Liquid density and caloric coupling

The initial cubic-model baseline is **untranslated PR/PR-BM density from its own
EOS volume**, with fugacity and caloric residuals derived from that same EOS.
No Peneloux translation or separate liquid-density correlation is silently
enabled. This choice establishes a consistent starting point, not a claim of
accurate liquid density.

A separately selected liquid-density correlation may later provide an explicitly
labelled engineering density only. It must not overwrite EOS volume in
fugacities, residual enthalpy/entropy or derivative identities. Any volume
translation intended to change thermodynamic volume must define its potential,
temperature/composition dependence, fugacity corrections, and caloric derivatives
consistently, with identity tests before acceptance. Selection of such a method
is deferred to A3 because no qualified correlation, parameters or formulation
has been accepted by A2. Report EOS and engineering density distinctly if both
are exposed.

## Numerical decisions and verification gates

Use `f64`, checked constructors and dedicated newtypes now. No linear-algebra
crate is selected in A1: these types require none. A3 must assess actual matrix
sizes, factorization needs and conditioning before selecting a crate; do not
hand-roll a general matrix library in this foundation.

Derivative implementation is deferred to A3 until the PR-BM alpha formulation
and caloric reference inputs are qualified. Choose analytical derivatives or
automatic differentiation there, recording the rationale and supported orders.
Finite differences serve independent derivative checks with multiple step sizes,
not an undocumented production default. Density, fugacity and caloric derivatives
must share the chosen thermodynamic potential. No derivative dependency is added.

Tests in this crate exercise known SI and gas-law values, hand-evaluated
mass/mole examples, explicit policies, malformed/nonfinite inputs, range failures,
zero/trace/pure limits and joint component permutation. Small scalar checks use
`abs(actual - expected) <= 1e-14 + 1e-12 abs(expected)` in the quantity's SI
unit; this covers floating-point arithmetic only, not model accuracy. Exact
comparisons are reserved for identity, enums, error classes and structure.
These checks use no external models and cannot qualify PR/PR-BM or dataset
accuracy. Composition unit-sum tests use tighter arithmetic checks than the
input acceptance tolerance.

The existing regular-case residual targets remain composition sum `1e-10`,
normalized component balance `1e-8`, active log-fugacity residual `1e-7`, scaled
energy balance `1e-6`. A4/A5 must define actual component/energy normalization
and scaling floors with their algorithms; A1 supplies no balance implementation
and makes no residual acceptance claim. A2 must qualify reference cases,
property-specific tolerances, domains and caloric inputs before model acceptance.
Model permutation/property sweeps and proptest remain with A3 as required by
the testing strategy; no dependency is added for deterministic A1 cases.

Acceptance requires review of A0 compatibility and A2 interfaces plus formatting,
package tests, Clippy and warning-free documentation on the pinned toolchain.
The worker's completion report records actual check outcomes and environment
limitations; the presence of this contract does not imply all milestone gates
have passed. Cross-platform workspace CI remains a coordinator/review gate.

## Reconciliation with concurrent A2 preparation

The eight provisional canonical IDs in
[reference qualification](reference-qualification.md) are `methane`, `ethane`,
`propane`, `nitrogen`, `carbon-dioxide`, `oxygen`, `water`, and `methanol`.
All satisfy `ComponentId`; their registry mappings remain data-layer records,
not a hardcoded list in this crate. The A2 examples' absolute Pa/K states,
mole fractions and mol/s flows agree with this contract. Every positive component
is active, zeros are retained, and the untranslated cubic baseline agrees.
Use the same `R = 8.31446261815324 J/(mol K)` for subsequent EOS and ideal
caloric terms; a reference with a different R needs an explicit matching override.

A1 agrees with A2's proposed caloric convention, subject to coordinator review:
pure-species ideal-gas molar h and s are zero at 298.15 K and 100000 Pa.
Preserve NASA7 source constants and explicitly transform their offsets to this
reference; do not discard the original formation/absolute-entropy provenance.
Mixture ideal entropy includes `-R sum(x_i ln(x_i))` and
`-R ln(p/100000 Pa)`, using the limiting value `0 ln(0) = 0` for zero entries.
Ideal h/s and consistent same-T,p,x EOS departures form the total properties;
this does not make total real-fluid h/s zero at the reference state.
No NASA7 evaluator, offset transformation or caloric calculation is implemented
here, and their independent verification remains an A2/A3 gate.

Exact critical-constraint cubic factors versus rounded historical PR factors
remain an explicit model/reference selection gate for A3, before implementation
or tight comparisons. Metadata must identify the selected variant; tolerances
must not hide differences. Independent expected outputs, demonstration phase
qualification, source fit/join handling, and final property-specific tolerances
remain open in A2. This contract does not accept the provisional demonstration
phases or source-extraction checks as numerical model validation.

## A1 worker check record — November 7, 2026

Environment: Windows PowerShell, `rustc 1.99.1 (58cf8f9c3 2026-10-13)` and
`cargo 1.99.1 (5f94df478 2026-08-27)`, matching the pinned toolchain.

| Scoped check | Final result |
|---|---|
| `cargo fmt -p thermo-types -- --check` | passed |
| `cargo test -p thermo-types --locked` | passed: 2 unit tests, 11 public-behavior tests, 0 doctests |
| `cargo clippy -p thermo-types --all-targets --locked -- -D warnings` | passed |
| `cargo doc -p thermo-types --no-deps --locked` with `RUSTDOCFLAGS=-D warnings` | passed |
| `git diff --check -- crates/thermo-types docs/calculation-contract.md` | passed for tracked changes; new files also reviewed |

The first Clippy run failed on `items_after_test_module`; moving the test module
after the implementation fixed it, and the complete final check set above then
passed. An attempted `rg` inspection failed because `rg` is unavailable on this
machine; `Get-Content` supplied the inspection instead. An initial orchestration
mail read used an unsupported `--from` flag and was retried successfully with
`--terminal`; neither operational issue blocks the package checks.

Changed files are this contract, `crates/thermo-types/src/lib.rs`, the new
`composition.rs`, `error.rs`, `flow.rs`, `metadata.rs`, `quantity.rs` modules,
and `crates/thermo-types/tests/public_contract.rs`. No dependency, manifest or
lockfile changes were required. Other workers' files were not edited, and no
staging, commit, reset or checkout was performed. Shared workspace checks and
Linux CI were not run by this worker. Coordinator contract/API review, A2
reconciliation and cross-platform acceptance remain open; A3 owns the explicitly
deferred matrix/derivative selection and any translated/correlated density path.
