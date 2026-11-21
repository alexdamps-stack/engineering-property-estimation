# Independent numerical validation

The [approved A2 tolerance registry](tolerances.json) consolidates the exact
cubic, flash, PH, process, Clapeyron, A13, NRTL, Helmholtz and GERG budgets.
Check registry/fixture agreement, registration of every committed fixture and
process input/adapter hashes with `python validation/consolidate-tolerances.py`.
The offline Rust Clapeyron, NRTL, A13 and A4 flash tests also recompute each
fixture's recorded `input_sha256` and fail when an input file has changed. The retired preparation proposals
in `tolerances.preparation-history.json` are historical only. Bounded PR/PR-BM
bubble/dew and envelope qualification is implemented by A13
([below](#a13-bubbledew-and-envelope-references)); critical-point and
retrograde envelope continuation remains deferred.

The [A2 process fixture](fixtures/process-v1.json) adds all seven demonstration
streams for each cubic variant, independently generated with the pinned Python
reference environment using `generate-process-references.py --output <candidate>`.
Review candidates before adopting them; the generator never calls Rust.
See [stream qualification](../docs/a2-process-qualification.md) for the fixed
state domain, phases, assumptions and separation from A8 flowsheet acceptance.

### Reference stability acceptance

The independent cubic and GERG Python TPD searches require a successful finite
minimization before using a nonnegative or near-zero minimum as stability
evidence. An unsuccessful trajectory may still witness instability when its
finite TPD is below the existing `stable_tm_floor` (-1e-8). If no successful
finite minimum or instability witness exists, generation fails as unresolved;
it does not emit a stable state. This rule also applies to phase stability after
a split and to the PH and process generators that reuse the cubic search.
Finite starts remain evidence, rather than proof of global stability.

Run acceptance regressions in the pinned Python reference environment with
`python -m unittest discover -s validation/tests -v`. They simulate unsuccessful
and successful optimizations, negative witnesses and nonfinite results for both
generators, independently of Rust and EOS evaluation.

The 21 November 2026 regeneration for #134 reviewed all 22 cubic flash cases,
30 pure saturation states, 28 PH cases, 14 process streams, and 13 GERG TP / 15
GERG PH cases against the preceding fixtures. Phase classifications, physical
states, compositions, caloric targets and tolerances were unchanged. Six cubic
and 17 GERG near-zero TPD diagnostic leaves changed, by at most 5.91e-14.
Provenance hashes were refreshed for the edited adapters, including the
Clapeyron fixture's PH input hash; the flash environment record now reflects
the pinned NumPy 2.3.4 and SciPy 1.16.3 actually used. The GERG Julia EOS/AD
check also passed all states and two fresh binary solves. No Rust values were
used to generate references, and no tolerance was widened.

For native GERG-2008, ordinary Rust tests compare 273 homogeneous derivative
states and 269 locally admissible property states against independently executed
NIST C++ 2.01. Reproduce without writing using `python3 validation/prepare-gerg.py`
(requires g++); `--write` explicitly replaces outputs and requires diff review.
The byte-exact reproduction assumes g++ targets the generic x86-64 baseline. On a
toolchain that defaults to AVX2/FMA (e.g. Ubuntu's amd64v3 variant, as on GitHub's
`ubuntu-latest`; check `g++ -Q --help=target | grep march`), put a `g++` link to
`validation/gxx-baseline.sh` first on `PATH`, as the numerical-validation workflow
does. The Python generators likewise assume the AVX2 kernel path; on AVX-512
machines set `OPENBLAS_CORETYPE=Haswell` and disable NumPy's AVX-512 dispatch with
`NPY_DISABLE_CPU_FEATURES` (values in `.github/workflows/numerical-validation.yml`).
`python validation/check-gerg-teqp.py` adds teqp 0.23.2 AD comparisons using the
existing pinned reference environment. See [A6 scope](../docs/a6-status.md) and
[executed checks](a6-checks.md), including ideal-reference compatibility.

The additional GERG density and fugacity generators default to read-only
reproduction of 42 TP roots and 45 chemical-potential states:

```sh
python validation/generate-gerg-density-references.py
python validation/generate-gerg-fugacity-references.py
```

The density check requires g++ and the pinned Python reference dependencies.
Both accept `--output <candidate.json>` for explicit candidate generation;
review any differences before replacing a fixture. They never call Rust.

GERG equilibrium references add 13 TP cases and 15 PH targets, including pure
coexistence, a narrow-boiling mixture and a valve. Generate an explicit temporary
candidate with the pinned Python requirements and g++:

```sh
python validation/generate-gerg-equilibrium-references.py --output /tmp/gerg-equilibrium-candidate.json
julia --compile=min --project=validation/julia validation/julia/check-gerg.jl
```

Generation uses NIST pressure/calorics, teqp AD and independent SciPy stability
and phase-split algorithms. The Julia command checks native Clapeyron EOS/AD
equilibrium residuals on those exterior-branch states and performs two fresh
binary TP solves. It does not claim a general Clapeyron GERG flash qualification;
nonfinite backend initialization results are excluded from evidence, not treated
as successful phases. Reference-state mappings and supported domains are in the
[A6 record](../docs/a6-status.md). Existing PR/PR-BM standard checks stay mandatory.

Clapeyron.jl is a mandatory independent standard check for supported PR/PR-BM
calculations. Run `julia --project=validation/julia validation/julia/check.jl`
with Julia 1.11.7. The locked Julia environment is validation-only. Its live
generator evaluates 28 PH and 60 homogeneous states using matched parameters,
and the Rust integration tests compare the current engine with those temporary
values. The same comparisons run offline against
[`fixtures/clapeyron-v1.json`](fixtures/clapeyron-v1.json) in ordinary workspace
tests. Push/PR and scheduled/manual numerical CI run the live check. See
[Clapeyron checks](clapeyron-checks.md) for provenance, tolerances and limits.

The A5 [PH fixture](fixtures/ph-v1.json) adds 28 independently generated
PR/PR-BM cases with caloric targets, phase fractions/compositions, pure
coexistence endpoints, a narrow-boiling mixture and isenthalpic valve outlets.
Generation reuses the reviewed A3 caloric and A4 equilibrium adapters, with
their pinned dependencies in `reference-requirements.txt`. Ordinary Rust tests
consume the fixture offline. See [A5 status](../docs/a5-status.md) and
[checks](a5-checks.md) for evidence and limits.

```sh
python validation/generate-ph-references.py --output /tmp/ph-candidate.json
```

Review candidate values and provenance before replacing the checked-in fixture.

This directory contains A2 preparation tooling, [archived preparation tolerance proposals](tolerances.preparation-history.json), and an independently generated [homogeneous PR/PR-BM fixture](fixtures/cubic-v1.json). The fixture contains 60 cubic and 64 NASA7 states generated by pinned teqp/NumPy/SciPy/Cantera. See [A1–A3 status](../docs/a1-a3-status.md) for exact formulation, accepted comparison budgets, reproduction and limits. The A4 [flash fixture](fixtures/flash-v1.json) from `generate-flash-references.py` adds 22 stability/TP-flash and 30 pure-saturation states; see [A4 status](../docs/a4-status.md) and [checks](a4-checks.md).

The separate [Helmholtz regression fixture](fixtures/helmholtz-dilute-v1.json)
contains six regular/dilute states evaluated independently at 80 decimal digits.
To explicitly generate a candidate using only Python's standard library:

```sh
python3 validation/generate-helmholtz-reference.py --output /tmp/helmholtz-candidate.json
```

Review differences before replacing expected values. Neither generator runs in
the engine or ordinary Rust tests. See [review closure](a1-a3-checks.md).

The [testing strategy](../docs/testing-strategy.md) defines the reference-case requirements, independence rules, numerical scenarios, tolerances, and reporting conventions. Keep small qualified fixtures in Git for offline crate integration tests. Reference-generation tools may use external software; ordinary Rust tests must not require it. Record and review every regeneration of expected values.

Record the reference implementation and version, model formulation, composition, parameters, units, reference states, domain, and per-property absolute/relative tolerances. Keep comparison code separate from the production implementation to avoid comparing a model with itself.

Report implementation verification, model-to-model comparisons, and experimental validation separately. Historical prototype test counts are not current repository results. Ordinary CI runs the offline numerical comparisons; the monthly/manual numerical-validation workflow adds a fixed-seed release sweep and re-runs the Python reference generators into temporary candidates, comparing each with its committed fixture via `compare-reference-candidate.py`. Neither replaces expected values.

Run `powershell -NoProfile -ExecutionPolicy Bypass -File validation/check-preparation.ps1` from the repository root for offline artifact-integrity checks. The explicitly invoked `prepare-source-data.ps1 -OutputDirectory <directory>` downloads pinned sources and transcribes selected records; it does not implement or evaluate any model. Review regeneration before adopting changes. Exact commands, hashes, checks and limitations are in [preparation checks](preparation-checks.md) and [reference qualification](../docs/reference-qualification.md). These files are provisional and do not finalize A1 integration or A2 acceptance.

## A7 NRTL / gamma–phi references

[`fixtures/nrtl-v1.json`](fixtures/nrtl-v1.json) comes from
[`julia/generate_nrtl.jl`](julia/generate_nrtl.jl) in the same locked Julia
environment. It holds 30 methanol/water and 9 directed-ternary activity states,
10 pure-liquid states, 36 bubble/dew points, 8 TP flashes and 16 PH targets.
`check.jl` verifies its provenance and compares a live candidate with the Rust
engine through `THERMO_NRTL_REFERENCE`. Ordinary tests compare the committed
fixture offline. Explicit regeneration:

```sh
julia --project=validation/julia validation/julia/generate_nrtl.jl --output candidate.json
```

The source data are transcribed only by `python3 -I validation/prepare-nrtl-source-data.py
--output-directory <dir>`. See [A7 status](../docs/a7-status.md) and [checks](a7-checks.md).

## A13 bubble/dew and envelope references

[`fixtures/a13-v1.json`](fixtures/a13-v1.json) comes from
[`julia/generate_a13.jl`](julia/generate_a13.jl) in the same locked Julia
environment, using native Clapeyron bubble/dew pressure and temperature solvers
on the inputs in [`a13-inputs.json`](a13-inputs.json). It holds 62 points: 32
PR/PR-BM bubble/dew pressure and temperature cases and 30 PR-BM envelope
samples. `check.jl` verifies its provenance, regenerates a live candidate,
checks that it reproduces the committed fixture and runs
`cargo test -p thermo-solvers --test a13_clapeyron` with `THERMO_A13_REFERENCE`.
Ordinary tests compare the committed fixture offline. See
[A13 status](../docs/a13-status.md) for the formulation, supported domain and
limits.
