# Parallel work plan

## A11-A13 reconciliation - 2026-11-09

The first-wave plans and scoped acceptance records below remain historical.
The current completion/merge gates are:

| Package | Delivered scope and record | Remaining gate / deferred extensions |
|---|---|---|
| A11, [PR #85](https://github.com/alexdamps-stack/engineering-property-estimation/pull/85) at `c54d3d2` | [Status and measurements](a11-status.md): one exact-PR fixed-natural-gas table, 280-420 K and 1-80 bar; separate value/derivative targets, speedup, explicit direct/error fallback | Merged in PR #85 after final CI/review; variable composition, PR-BM, other packages and flowsheet integration deferred; finite-lattice/TPD qualification and C^1 edge limits remain |
| A12, merged PR #84 | [Status](a12-status.md): versioned C ABI/.NET boundary and real DWSIM 9.0.4 Automation separator cross-check | DWSIM property-package/CAPE-OPEN/GUI hosts and Excel remain follow-ups; cross-model integration tolerances are not engine reference budgets |
| A13, PR #86 | [Status](a13-status.md): 62 independent Clapeyron points, dry PR/PR-BM bubble/dew p/T, natural-gas dew point and bounded sampled envelope with near-critical stop; pure-fraction correction independently cleared | Final-head CI/merge pending; no critical-point estimate, retrograde continuation, wet VLLE or global-stability proof |

The [testing strategy](testing-strategy.md) and [approved registry](../validation/tolerances.json)
now expose the unchanged A13 quantity budgets. A11 compares interpolation against
the qualified direct model; it does not generate an independent engine reference.
A13 does not enable A14-A16 unit specifications or project-format changes by itself.
A13 is rebased on merged A11/A12 main; both the `property-tables` and
`thermo-bridge` lock entries and the A13 package-test dependency are retained. No historical acceptance is broadened.

## A3/A4 scoped review — 2026-11-08

The [A3/A4 acceptance review](a3-a4-acceptance.md) records separate PR and
PR-BM homogeneous acceptance and reviews bounded two-phase stability/TP flash,
including safeguarded Newton acceleration for #30. It supersedes the pending
A4 coordinator-review note in the historical implementation report. Current
check outcomes are in [A4 checks](../validation/a4-checks.md); hosted checks
and PR review remain visible on GitHub. A13 phase boundaries/envelopes,
critical-point accuracy and global-stability proof are not included.

## Formal A1/A2 acceptance — 2026-11-08

A1 is accepted after review of the calculation contract, public type/conversion
and failure tests, A0's fresh-implementation decision and A2 compatibility.
SI, composition/flow bases, named standard-volume conditions, untranslated
liquid-density baseline, caloric convention and staged numerical decisions
are explicit. Current Linux/Windows CI has passed, superseding the historical
first-wave Windows execution restriction. A2 is accepted for the pinned data,
matched formulations/references, fixed demonstration states and approved
quantity budgets. The [formal acceptance record](a1-a2-acceptance.md) provides
evidence and scope; [process qualification](a2-process-qualification.md) records
each stream. A13 bubble/dew/envelope and A8 connected process execution remain
their own gates. The first-wave reviews below are historical.

Continuation on November 7: the [A1–A3 verification record](a1-a3-status.md)
documents the subsequent homogeneous implementation, independent references
and review. The first-wave history below is retained; no first-wave worker is
still assigned. Full milestone gates are not superseded by implementation tests.

Prepared November 7, 2026 against `bf86da9`, using the maintained
[development plan](development-plan.md), artifact inventory, architecture, and
testing strategy. All eight crates are still scaffolds. This dispatch starts the
foundation work; it does not authorize skipping numerical acceptance gates.

## First wave

Three Codex workers share the current checkout with exclusive file ownership.
They must not commit, reset, stage, or overwrite another worker's changes.

| Worker | Scope | Exclusive write ownership | Deliverable |
|---|---|---|---|
| Artifact inventory (A0) | Bounded search for earlier ThermoCore source and reports; assess reuse | `docs/artifact-inventory.md`, `docs/artifact-recovery.md` | Search evidence, recovered artifact provenance or explicit unavailability, and a port-versus-fresh implementation decision |
| Calculation foundations (A1) | Draft the calculation contract and implement shared validated types | `crates/thermo-types/**`, `docs/calculation-contract.md` | SI and basis conventions, explicit failures and metadata, meaningful boundary/conversion tests, and numerical-foundation decisions |
| Data and reference preparation (A2 preparation) | Qualify primary sources and define demonstration/reference inputs | `data/**`, `validation/**`, `examples/**`, `docs/reference-qualification.md` | Initial component identities, candidate data and reference provenance, PR versus PR-BM formulation matching, demonstration inputs, justified tolerance proposals, and explicit gaps |

The A1 worker may draft and test types while recovery runs, but A1 acceptance
requires reviewing A0's decision. The A2 worker prepares source material and
case specifications independently; schemas and datasets remain provisional
until the A1 contract is reviewed. Neither concurrency nor passing infrastructure
checks completes A0, A1, or A2 automatically. Unavailable reference evidence must
remain an open gap rather than invented expected values or accuracy claims.

Workers read the existing project guidance, stay within their ownership, and
report changed files, checks and outcomes, source evidence, and remaining gates
through Orca orchestration. Shared root manifests, lockfile, architecture, testing
strategy, and existing development plans are coordinator-owned. Dependency or
cross-boundary changes require a request to the coordinator.

## Review and subsequent waves

1. Review A0's decision, reconcile A1 types and contract, then reconcile A2
   schemas, demonstration domain, caloric reference states, independently
   generated fixtures, and quantity-specific tolerances. Apply agreed changes
   to shared documentation and run the workspace checks. Explicitly record
   remaining gaps; do not mark a milestone complete without its evidence.
2. After A2 acceptance, split A3 into formulation/caloric implementation and
   independent comparison/consistency checks with agreed APIs and distinct
   ownership. Accept PR and PR-BM separately.
3. After A3 acceptance, implement and verify A4 stability and TP flash. A GUI
   API prototype may then proceed separately; it must consume the engine API.
4. After A4, PH flash (A5) and saturation/envelope work (A13) can proceed as
   separate tracks. GERG (A6) can proceed after A2 with integration after A4.
5. After A5, begin the headless flowsheet (A8) and the qualified NRTL path (A7).
   Accept the graphical workflow (A9) only after the headless demonstration.
   Recycles, splines, and external integration remain deferred by their gates.

Only the first wave is launched by this dispatch. Later assignments depend on
reviewed results, rather than speculative implementation against moving APIs.

## Coordinator review: A0

On November 7, 2026, reviewed the completed artifact inventory and recovery
record. A0 is accepted on the explicit-unavailability branch: proceed with fresh
Rust work, retaining the documented missing historical artifacts as evidence
gaps. The separate local Rust candidate at commit
`0985d1ac85cdbddc386467286044f20489964fc6` may be inspected for reuse, but its
code, data, declared reference versions, and numerical claims are not accepted
by this review. A1 was informed that recovery no longer blocks its foundation
work; A2 was asked to assess the candidate's formulation, provenance, terms,
and reference independence before proposing any import.

## Coordinator review: A2 preparation

Reviewed `docs/reference-qualification.md` and
`validation/preparation-checks.md`. The assigned preparation deliverable is
complete: pinned eight-component constants and NASA7 transcriptions, retained
notices, reproducible source tooling, provisional inputs, and explicit pending
tolerance proposals. The worker reported passing source reproduction and offline
artifact checks, including malformed-input rejection. No independent numerical
expected values were generated; A2 itself remains unaccepted. A1 was sent the
contract reconciliation points. Cubic-factor selection, caloric offsets,
reference execution, process phase/domain qualification, and evidence-based
tolerances remain gates before model acceptance. The recovered local candidate
was inspected but not imported.

## Coordinator review: A1 and first-wave checks

All three first-wave workers reported successful completion of their assigned
deliverables. Reviewed the A1 contract, quantity/composition/flow APIs, metadata,
and public-behavior tests. The shared gas constant, canonical IDs, explicit
bases, retained zeros, untranslated EOS baseline, and proposed ideal caloric
reference are compatible with A2 preparation. This historical foundation review
is superseded by the formal A1 acceptance entry above.

The worker reported 2 unit and 11 public-behavior tests passing on the pinned
Windows toolchain. Coordinator combined-workspace formatting, Clippy with
warnings denied, documentation with warnings denied, and `git diff --check`
passed. The combined workspace test command compiled the crates but stopped
before executing the `property-packages` test binary: Windows application
control blocked it with OS error 4551. This is an unresolved environment
restriction, not a test assertion failure or a passing workspace test result.
Linux/cross-platform CI is still outstanding. No policy bypass was attempted.

Next work is A2 reconciliation and independently executed reference evidence,
including explicit cubic factors and caloric offsets, before A3 model work.
There are no workers still assigned to the first wave. No commit or import of
the recovered model candidate was made.
