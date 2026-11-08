# A1–A3 implementation and verification record

Date: 2026-11-07. Scope: homogeneous, fixed-composition PR and PR-BM.
Completed checks and the independent review are recorded in
[the verification report](../validation/a1-a3-checks.md).
This supersedes the scaffold-only descriptions for `thermo-data`,
`thermo-models`, and `thermo-cli`. The earlier preparation records are retained
as history. No recovered candidate source or numerical fixture was imported.

## Implemented contract

- A1: retain the existing validated SI types, explicit composition/flow bases,
  exact SI gas constant, error classes and untranslated EOS baseline.
- A2: load the eight pinned CoolProp 7.2 cubic/Cantera 3.2 NASA7 records through
  an immutable database; reject unknown identities. Preserve upstream constants,
  a6/a7, source hashes and notices in the original extract. No alias guessing.
- A3: separate `PR-exact-1976` and `PR-exact-BM-Clapeyron-0.6.30` selections;
  exact critical-constraint factors `OmegaA=0.4572355289213822`,
  `OmegaB=0.07779607390388846`, `R=8.31446261815324 J/(mol K)`.
  These are deliberately distinct from rounded-factor CoolProp PR.
- Linear covolume and quadratic attraction mixing with a complete, explicitly
  supplied symmetric, finite, temperature-independent interaction matrix,
  zero diagonal and `kij<1`. Missing interactions never become zero.
  Examples explicitly choose synthetic zeros; one reference case uses 0.035.
  Parameter order follows composition; metadata includes identities and matrix.
- All admissible real roots are isolated by partitioning the cubic at its
  derivative's real zeros and bisecting each monotone interval, up to 128
  iterations per root. Bracket width and scaled polynomial residual are checked.
  A forced-zero-iteration unit test verifies the exhaustion error. Coalescent
  admissible roots are diagnostic errors. There is no generic matrix solver.
- Root requests mean smallest volume, largest volume, or minimum Gibbs energy
  **among roots at the supplied composition**. They do not classify the phase
  or establish mixture stability. Mechanically unstable/spinodal and calorically
  unstable evaluations fail. Z roots may include the unstable middle root.
- Return EOS volume/density, all ln(phi), molar h/s/u, cp/cv, speed of sound,
  pressure derivatives, model/data/caloric identities and unassessed physical
  validity. CLI output names every property's SI unit and reports the limitation.

Production derivatives of alpha and attraction through second temperature order
are analytical. Finite differences at multiple steps are used only for tests.
No AD or linear algebra dependency is needed for the scalar cubic; matrix
selection for future flash Jacobians remains deferred to A4. Liquid density is
the untranslated EOS result, with no separate correlation or translation.

## Caloric convention and computational envelope

NASA7 cp/R uses the first five retained coefficients. Integrate cp and cp/T
from 298.15 K so pure ideal h=s=0 there at 100000 Pa. At 1000 K, the upper
interval owns cp; integrals cross the join continuously, removing the original
fit's a6/a7 offset difference rather than treating it as an energy jump. Those
original source offsets remain in the extract for audit. Ideal mixture entropy
includes `-R sum(x ln x) - R ln(p/100000)` with `0 ln 0=0`.

EOS departures use the same T,p,x ideal reference. If `q` denotes
`ln((v+(1+sqrt(2))*b)/(v+(1-sqrt(2))*b))/(2*sqrt(2)*b)`, then
`h_res=RT(Z-1)+(T*a_T-a)*q`, `s_res=R*ln(Z-B)+a_T*q`, and
`cv=cv0+T*a_TT*q`; `cp=cv-T*p_T^2/p_v` and
`w^2=-(cp/cv)*v^2*p_v/M`. All quantities are molar, with M in kg/mol.

The shared attraction integral uses the algebraically equivalent `ln_1p`
form to preserve the dilute limit. Six additional 80-digit Decimal references
cover regular and very dilute methane PR/PR-BM Helmholtz energies. The public
Helmholtz helper accepts finite V/n>b at 200–600 K independently of the TP
pressure bounds and reports a numerical error if V/n overflows.

NASA7 evaluation accepts 200–6000 K without extrapolation. Homogeneous total
property evaluation is bounded to 200–600 K and 0.01–100000000 Pa, with
diagnostic exclusions at singularities. Mechanical pressure derivatives alone
accept 100–1000 K to test all eight critical points, including methane and
nitrogen below the caloric range. These are computational bounds, not an
experimentally qualified physical validity domain. Water/methanol cubic values
are implementation checks, not an accepted polar-fluid or liquid model.

BM uses squared Soave alpha at Tr<=1 and
`exp(2*(1-1/d)*(1-Tr^d))`, `d=1+m/2`, above it. Alpha and its first derivative
are continuous at Tc; its second derivative generally is not. The low branch
owns the reported one-sided alpha second derivative at Tr=1. Total caloric
evaluation returns an explicit unsupported error exactly at any active
component's Tc rather than assigning a two-sided cp/cv there.

## Independent references and tolerances

`validation/generate-references.py` is a new adapter, using **teqp 0.23.2**
residual Helmholtz AD, **NumPy 2.3.4** polynomial roots and **SciPy 1.16.3**
pressure refinement. teqp Mathias–Copeman `[m,0,0]` implements the selected PR
alpha even for high acentric factors (avoiding PR78's different polynomial).
Twu `[2*(1-1/d),1,d]` algebraically implements BM above Tc. Each component's
branch is fixed at the requested T for local derivatives. This is not a Julia
execution of Clapeyron, and transition continuity is checked separately.
The pinned teqp header and Clapeyron BM source URLs and hashes are in the
fixture. teqp's runtime metadata and R are asserted during generation.

Cantera 3.2.0 evaluates the retained NASA7 coefficients with independently
computed reference offsets, including both fit intervals and the exact join.
Expected values never call Rust. Shared parameter sources remove parameter
differences intentionally; agreement verifies implementation, not measurements.
Software and numerical-library versions and input-file SHA-256 are retained.

teqp is reference tooling under its pinned [NIST permission and warranty
notice](https://github.com/usnistgov/teqp/blob/v0.23.2/LICENSE.md); it is not
redistributed as application code. Cantera/CoolProp data notices remain in
`data/source-extracts`. The Clapeyron source is inspected as formulation evidence
under its recorded MIT terms; no Julia implementation is copied.

The offline fixture contains **60 cubic states** (30 per model): the original
18-state regular grid per model, eight pure states, a three-root propane state,
a nonzero interaction case, a zero-component case, and the lower pressure bound.
The **64 ideal states** cover eight species at 200, 298.15, 600, 999.9, 1000,
1000.1, 3000 and 6000 K. No failed state is skipped. Mechanical middle roots
remain in the root comparison but do not receive stable caloric assertions.

Absolute-plus-relative budgets were fixed in the generator before Rust
comparisons. See its `TOLERANCES` and the fixture for exact machine-readable
values. Z: 2e-10/2e-9; mass density: 1e-7 kg/m³/2e-9; ln(phi): 2e-9/2e-9;
h: 2e-5 J/mol/2e-9; s: 2e-7 J/(mol K)/2e-9;
cp/cv: 2e-5 J/(mol K)/2e-8; sound: 2e-4 m/s/2e-8.
They cover double-precision evaluation and root refinement away from singular
states; they are not accuracy claims across the entire computational envelope.

Consistency checks cover h=u+pv, cp=d(h)/dT, Maxwell entropy/volume derivatives,
cp−cv, Gibbs–Duhem, fugacity-Jacobian symmetry, partial-molar Gibbs summation,
all eight pure critical constraints, BM transition and alpha derivatives,
pure/zero/trace limits, low density, Euler homogeneity and failures.
Property tests use a fixed seed `0xA1A2A3`, 64 cases by default, on the declared
smooth grid 420–550 K, 0.1–3 MPa and binary fractions 0.05–0.95; they check
permutation, amount scaling, normalization and Euler/energy identities.

## Reproduce and review

Reference dependencies are separate from the Rust application:

```sh
python3.12 -m venv /path/to/reference-venv
/path/to/reference-venv/bin/pip install -r validation/reference-requirements.txt
/path/to/reference-venv/bin/python validation/generate-references.py --output /tmp/cubic-candidate.json
```

Generation requires network access only to hash the two pinned source headers.
Compare the candidate with `validation/fixtures/cubic-v1.json`; replacing expected
values requires an explicit reviewed change. Do not regenerate on test failure.
This change added exact-1000-K fixtures after review found Cantera's join branch
needed explicit low/high evaluation in the adapter. Existing non-join values
were not relaxed and no tolerance was widened.

```sh
cargo run -p thermo-cli --locked -- examples/homogeneous-natural-gas.json
cargo run -p thermo-cli --locked -- examples/homogeneous-binary.json
cargo test -p thermo-types -p thermo-data -p thermo-models -p thermo-cli --locked
```

The CLI accepts one strict JSON case with molar fractions, absolute pressure,
model, root selection and complete interactions. It rejects unknown fields,
malformed quantities and missing data, returning a nonzero status on failure.

## Acceptance boundaries

A1 and A2 are accepted in the [formal review](a1-a2-acceptance.md). The
[demonstration qualification](a2-process-qualification.md) supplies independent
feed, mixer, cooler, PH valve and separated-stream evidence for both variants.
The [approved tolerance registry](../validation/tolerances.json) supersedes the
historical preparation proposals. A3 homogeneous implementation evidence and
its computational bounds above are retained; A4/A5 own stability and TP/PH
solver acceptance. A13 owns bubble/dew and envelope qualification. Acceptance
does not assert experimental accuracy, a continuous process phase map, global
mixture stability or critical-neighborhood accuracy. A8's connected flowsheet
execution remains a separate milestone.
