# A2 reference preparation and qualification

**Implementation update:** [A1–A3 status](a1-a3-status.md) records the implemented
homogeneous data/model path and executed independent teqp/Cantera references.
It resolves the earlier exact-factor, offset and homogeneous expected-value
gaps. [A2 process qualification](a2-process-qualification.md) now supplies the
demonstration phase/PH states, and [formal acceptance](a1-a2-acceptance.md)
records the A1/A2 decisions and approved tolerance registry.
The preparation-only statements below describe the earlier first-wave record.

Prepared 2026-11-07 against scaffold baseline `bf86da9`. **Preparation only:**
A1 integration, independent output generation, phase qualification and numerical
tolerance approval remain pending. No PR, PR-BM, flash or experimental acceptance
is claimed. Files here deliberately use provisional schemas and are not engine
input formats. Other workers own the calculation contract and recovery decision.

## Source register and permission evidence

The generator records exact source URLs and SHA-256 of UTF-8 response text in
`data/source-extracts/components.provisional.json`; these are text hashes, not
claims about compressed HTTP bytes. All numeric extraction uses invariant culture.

| Source | Pinned version / use | Terms / qualification |
|---|---|---|
| [CoolProp cubic data](https://github.com/CoolProp/CoolProp/blob/98b3523d5daa98454618d381d2ae53f7471d216b/dev/cubics/all_cubic_fluids.json) | v7.2.0, commit `98b3523d5daa98454618d381d2ae53f7471d216b`; eight CAS-matched constant records | [MIT](https://github.com/CoolProp/CoolProp/blob/98b3523d5daa98454618d381d2ae53f7471d216b/LICENSE), full notice retained locally; qualified as reproducible software parameter definitions, not a measurement dataset |
| [Cantera gas data](https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/data/nasa_gas.yaml) | v3.2.0, commit `4a8358eb80cfeb50474386b5f9ec0b3a83519889`; eight NASA7 gas records | [BSD terms](https://github.com/Cantera/cantera/blob/4a8358eb80cfeb50474386b5f9ec0b3a83519889/License.txt) retained; source notes preserved |
| [NASA-TM-4513](https://ntrs.nasa.gov/citations/19940013151) | McBride, Gordon, Reno (1993), document 19940013151; upstream publication identified by Cantera | NTRS explicitly marks US Government work, public use permitted; publication describes polynomial fits with formation-inclusive enthalpy constants; original per-species literature and uncertainty not independently re-audited |
| [Clapeyron BM source](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/3440a15ba224664cac5b24441ac8adacbaf326f9/src/models/cubic/alphas/BM.jl) | v0.6.30, commit `3440a15ba224664cac5b24441ac8adacbaf326f9`; proposed PR-BM reference | [MIT](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/3440a15ba224664cac5b24441ac8adacbaf326f9/LICENSE.md); inspected, no source copied and no Julia execution |
| [Peng and Robinson (1976)](https://doi.org/10.1021/i160057a011) | original PR formulation citation | No paywalled table or code copied; software source fixes the actual candidate convention |
| [NIST teqp cubic documentation](https://pages.nist.gov/teqp-docs/en/main/models/cubics.html) | displayed 0.22.0, candidate exact cubic constants and residual derivatives | Documentation is not a pin for recovered candidate 0.23.2; package source/version matching still needed |

Constants are `Tc` in K, `pc` in Pa absolute, acentric factor dimensionless,
molar mass in kg/mol. Identifiers are methane/74-82-8, ethane/74-84-0,
propane/74-98-6, nitrogen/7727-37-9, carbon-dioxide/124-38-9,
oxygen/7782-44-7, water/7732-18-5 and methanol/67-56-1. The source's
`N-PROPANE` is mapped by CAS, not name guessing. Parameters are computational
definitions without uncertainty estimates; matching them verifies implementation
agreement only. There are no fitted interaction tables. All example pairs use
an explicit synthetic, symmetric, temperature-independent `kij=0` approximation;
this is not a default for missing physical data. No volume translation is selected.

## Formulation matching: separate gates

**PR:** the pinned [CoolProp implementation](https://github.com/CoolProp/CoolProp/blob/98b3523d5daa98454618d381d2ae53f7471d216b/src/Backends/Cubics/GeneralizedCubic.cpp)
uses attraction factor 0.45724, covolume factor 0.07780, and
`m = 0.37464 + 1.54226*w - 0.26992*w^2` with the usual squared
`1+m*(1-sqrt(T/Tc))` alpha. Candidate comparison is quadratic attraction and
linear covolume mixing, shared constants, explicit interaction matrix, no
translation. Specify R explicitly during generation and in the engine.
CoolProp's PR backend alone is **not** a PR-BM verification reference.

**PR-BM:** proposed version is Clapeyron v0.6.30 `PR` + `BMAlpha` +
`vdW1fRule` + `NoTranslation`, with all critical/acentric/molar-mass parameters
overridden from the selected data and both interaction matrices explicitly zero.
In its PR branch BM uses the same m polynomial below; for `Tr<=1` it uses
the squared PR expression. For `Tr>1`, let `d=1+m/2`; the continuation is
`alpha=exp(2*(1-1/d)*(1-Tr^d))`. This identifies a specific candidate variant,
not every commercial method named PR-BM. Value and first derivative match at
`Tr=1`; second-derivative continuity is not assumed. Confirm derivative behavior
with independent execution before A3.

Important mismatch: CoolProp 7.2 uses rounded cubic factors whereas
[Clapeyron cubic equations](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/3440a15ba224664cac5b24441ac8adacbaf326f9/src/models/cubic/equations.jl)
derive factors from critical constraints. teqp's documentation reports
`OmegaA=0.4572355289213822`, `OmegaB=0.07779607390388846` and
`R=8.31446261815324 J/(mol K)`. Rounded and exact variants cannot share tight
implementation tolerances. A1/A3 must select/version this choice; then use the
matching reference or explicit documented parameter overrides, rather than
widening tolerances. Reference caloric defaults must also be overridden.

Clapeyron and CoolProp are independent projects from this scaffold; this is
architectural independence only, not proof of independent parameter research.
Using the same constants intentionally removes one source of discrepancy and
cannot validate their experimental accuracy. Audit any future source reuse for
shared algorithms; a copied reference evaluator would reduce independence.

## Ideal-gas contribution and caloric states

The transcribed NASA7 records each contain two seven-coefficient intervals:
200–1000 K and 1000–6000 K. They are ideal-gas functions, not liquid heat
capacities or a claim that stable gas exists at every such state. With T in K,
`cp0/R = a1+a2*T+a3*T^2+a4*T^3+a5*T^4`; a6 and a7 set enthalpy/entropy
integration constants. Coefficients retain the source ordering and species note.
See the [Cantera species documentation](https://www.cantera.org/3.2/reference/thermo/species-thermo.html)
for cp, h and s definitions. Gas constants and reference pressure must match
before comparing Cantera values with the cubic engine. No extrapolation below
200 K is authorized; source fit boundaries require separate join diagnostics.
The 200–400 K process search interval lies inside the fit intervals, but this
does not qualify cubic phase stability or PH behavior.

Propose engine ideal enthalpy and pure-species ideal entropy each zero at
298.15 K, 100000 Pa, via integration of cp0 and cp0/T from that temperature.
This is a process-simulation reference proposal, **not** the source's formation
enthalpy/absolute entropy convention and not zero total real-fluid h/s. Preserve
the original constants for audit, explicitly transform offsets when generating
references. Include ideal mixing entropy `-R*sum(x_i*ln(x_i))`, pressure term
`-R*ln(p/p_ref)`, and the zero-component limiting convention agreed with A1.
Total h/s combine these ideal contributions and consistent residual quantities
relative to the ideal mixture at the same T,p,x. Residual entropy at constant
T,density must not be confused with constant T,p departure. Do not apply arbitrary
real-fluid NBP/IIR offsets or silently combine CoolProp alpha0 with NASA7.

Status: exact gas-correlation transcription, reported fit range and redistribution
terms are established. Engine R, pressure/entropy convention, join ownership,
offset implementation, original fit uncertainty and independently evaluated
cp/h/s values remain open. No NASA7 evaluator or caloric engine was implemented.

## Inputs, phase claims, state counts and tolerances

`examples/cases.provisional.json` specifies a dry five-component natural gas
and a methane/propane binary, both on mole-fraction and mol/s bases with
absolute Pa. `examples/process.provisional.json` specifies two unequal feeds,
adiabatic mixing at 5 MPa, cooling to 240 K, an isenthalpic valve to 1.5 MPa,
and equilibrium separation. Pressure loss, heat/work assumptions, search
envelope, required stream properties, duties and diagnostics are explicit.
These are chosen engineering inputs, not literature measurements.

No exact feed or outlet phase has been established by an independent calculation,
so every expected phase is null, including the feeds. No experimental or matched
model result supports a two-phase separator assertion yet. Water and methanol
are excluded to keep the intended scope to hydrocarbon vapor/liquid behavior;
solid CO2, multiple liquids and phase stability still require domain screening.
If references predict a single-phase outlet, revise and review setpoints rather
than changing expected phases by assumption. PH must stop diagnostically outside
the proposed envelope. There are no standard-volume flows or hidden reference
conditions. Outputs identify model/data versions and basis.

`validation/tolerances.preparation-history.json` supplies unit-bearing atol/rtol proposals
for density, Z, each ln(phi), residual h/s, cp, sound speed and saturation,
with rationales and minimum counts. All are **pending**, not approved thresholds.
The regular grid has 3 temperatures x 3 pressures x 2 mixtures = 18 candidate
states per model, 36 across PR and PR-BM, with no silent skipping on failure.
There are 9 proposed pure-saturation states/model, 6 mixture bubble/dew
states/model (temperatures not selected), 24 BM transition diagnostics/variant,
40 NASA7 property states and 24 join diagnostics. Single-phase qualification
and conditioning screening precede assigning the regular property tolerance.
Near-critical, trace and phase-boundary tolerances need separate evidence.
GERG remains a separate future qualification using the
[NIST AGA8 candidate](https://pages.nist.gov/AGA8/); no GERG tolerances invented.

Retain project targets of composition sum 1e-10 (dimensionless), normalized
component balance 1e-8, active ln fugacity residual 1e-7 and scaled energy
balance 1e-6. Proposed normalization is max component residual in mol/s divided
by max(total inlet mol/s, 1 mol/s); energy residual W divided by
max(sum(abs(inlet enthalpy flow W))+abs(Q)+abs(Wshaft), 1 W).
Active-component cutoff and scaling agreement with A1 remain pending; these
residual targets are not property accuracy against measurements.

## Recovered local candidate: qualification outcome

On coordinator instruction, inspected read-only the original checkout at
`C:/Users/PCUser/Documents/ChatGPT/Engineering Tools - Property Estimation`,
commit `0985d1ac85cdbddc386467286044f20489964fc6`, including the complete
`validation/generate_cubic_references.py`, initial reference JSON, and notice.
The exact file hashes and origin are in [A0 recovery](artifact-recovery.md).
**No candidate data or code was imported.**

Its 16 cases use teqp 0.23.2 residual Helmholtz derivatives and independent
density bracketing, with Mathias–Copeman `[m,0,0]` for PR and Twu
`[2*(1-1/d),1,d]` above Tc for BM. Algebraically the latter gives the BM
continuation above; executing exact teqp source, derivative checks and both
branches is still necessary. Selecting a branch from fixed requested T may
be appropriate for local derivatives away from Tc, but does not establish
transition continuity. It supplies roots, fugacity, density and residual h/s;
it supplies no total h/s, cp, sound speed, stability or flash evidence.
A 12000-point density sign-change search can miss double roots or tangencies;
it cannot establish completeness near critical or spinodal conditions.

PyPI availability probes returned HTTP 200 for teqp 0.23.2 and CoolProp 8.0.0;
these verify published version existence, **not** installed versions, generator
execution or fixture integrity. The candidate uses CoolProp 8 fluid EOS reducing
states, whereas this preparation pins CoolProp 7.2 cubic data; do not equate
data versions or hashes. Its source note labels original local code unlicensed,
so code redistribution permission remains unresolved despite upstream CoolProp
MIT terms. Locally inspecting it is authorized; adopting or copying its generator
requires a separate provenance/permission decision. The candidate can reduce
reference design work, but **cannot currently close any expected-value gate**.

## Reproduction, checks and next actions

Run from repository root, using Windows PowerShell 5.1 or compatible PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File validation/prepare-source-data.ps1 -OutputDirectory data/source-extracts
powershell -NoProfile -ExecutionPolicy Bypass -File validation/check-preparation.ps1
git diff --check -- data validation examples docs/reference-qualification.md
```

The first command intentionally regenerates only source transcriptions and
upstream notices; it requires network access. Review diffs and source hashes.
It never generates thermodynamic expected values. The second is offline and
checks identities, source pins, units/shapes, finite values, ranges, normalized
dry compositions, state counts, pending phase/expected-value status and three
malformed compositions. These checks validate artifacts, not physical models.

Reference-availability probes `python --version` and the attempted Python import
probe both reported the Windows Store alias error (Python unavailable), and
`Get-Command julia,py,node` found no executables. No reference-generation command
was run or succeeded. Initial `rg` discovery failed with CommandNotFoundException;
bounded PowerShell file enumeration replaced it. No full-drive search, secret
access, model implementation, engine dependency change or root manifest edit
was made. The exact validation outcomes and changed-file list accompany the
worker report under `validation/preparation-checks.md`.

Late read-only comparison with the concurrent A1 draft found compatible canonical
IDs, absolute SI states, mol/s flows, untranslated cubic density and explicit
missing-parameter behavior. Its gas-volume R is 8.31446261815324 J/(mol K),
which is the proposed shared model R; engine/caloric adoption still needs review.
A1 retains zero components and treats every positive component as active; future
reference adapters must preserve order and that convention. This is an initial
compatibility observation, not completed integration or acceptance.

Before A2 acceptance: reconcile schemas and SI/bases with A1; decide exact versus
rounded cubic factors, R and caloric offsets; install an approved reference-only
runtime; pin and audit the selected reference source and dependency lock;
generate both variants separately with metadata and conditioning diagnostics;
qualify all demonstration phases and PH states, then review numerical tolerance
budgets from precision/conditioning evidence; select actual saturation states;
preserve small expected-value fixtures with reviewed regeneration commands;
resolve recovered-code permissions before import. No current expected-value or
phase fixture may be inferred from this preparation's passing integrity checks.
