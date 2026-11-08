# A3 and A4 acceptance review

Review date: 2026-11-08. Baseline: `68ef34a` (main, including the A7 work).
Scope: GitHub #18 (already closed), #26 and its remaining implementation task
#30. This is a local implementation and evidence review; hosted checks on the
change and PR review are recorded separately on GitHub.

## A3 decision

PR and PR-BM are accepted separately for homogeneous calculations in the
documented computational envelope. The [A1–A3 record](a1-a3-status.md),
[calculation contract](calculation-contract.md), offline teqp/Cantera
comparisons and [live Clapeyron check](../validation/clapeyron-checks.md)
identify the exact formulations, component sources and shared caloric
reference. Tests compare roots, fugacities and caloric properties and cover
BM temperature derivatives/transition behavior, thermodynamic identities,
component permutation, zero/trace species and pure-fluid limits. The CLI
examples exercise both reference mixtures with model/data/unit/diagnostic
metadata. No further A3 implementation gap was found; #18 and #19–#25 remain
closed.

## A4 decision

The implementation covers the bounded two-phase PR/PR-BM acceptance criteria
in #26. Separate prescribed-K Rachford–Rice tests, multi-start finite TPD
trials, Gibbs decrease, phase reassessment, component balances and fugacity
equality are documented in [A4 status](a4-status.md). The unchanged independent
fixture has 22 flash and 30 pure-saturation states. Boundary/trace/zero cases,
phase appearance/disappearance and forced iteration exhaustion have explicit
results or diagnostic errors.

Issue #30 adds Newton acceleration of the active ln K system after the first
three SS evaluations. The Jacobian uses central differences of the complete
residual, including the implicit RR phase fraction. This is a numerical
Jacobian, rather than a new analytical fugacity-derivative API. A dense solve
with partial pivoting, a maximum ln K step of one and a residual-decreasing
backtracking line search safeguard the update. A rejected probe, singular
system or exhausted line search uses the existing SS/eigenvalue fallback.
Newton starts only from a physical phase fraction, retains physical trial
fractions and rejects split Gibbs increases beyond 1e-14 in dimensionless
per-mole units. The Gibbs guard prevents residual descent toward K=1 from
stranding states for which the previous solver converged.
Only the fugacity residual terminates a split; small Newton steps alone do not.
Inactive components are excluded from the system and retained as zero in
reported phases. Existing physical-fraction, Gibbs and post-split checks
continue to apply after convergence.

The new regression compares default Newton with the previous eigenvalue-only
algorithm at 250 K, explicit zero kij, for both models:

| Model | Methane/propane | Pressure (MPa) | Newton iterations | Eigenvalue-only iterations |
|---|---|---|---|---|
| PR | 75/25 | 8.5 | 7 | 70 |
| PR | 75/25 | 8.9 | 7 | 104 |
| PR | 75/25 | 9.025 | 6 | 96 |
| PR-BM | 75/25 | 8.5 | 8 | 63 |
| PR-BM | 75/25 | 8.9 | 7 | 79 |
| PR-BM | 75/25 | 9.025 | 7 | 118 |
| PR | 80/20 | 9.0 | 16 | 110 |
| PR | 80/20 | 9.1 | 15 | 136 |
| PR | 80/20 | 9.15 | 16 | 151 |
| PR | 80/20 | 9.2 | 16 | 171 |
| PR | 80/20 | 9.25 | 14 | 405 |
| PR | 80/20 | 9.3 | 14 | 599 |
| PR-BM | 80/20 | 9.0 | 16 | 88 |
| PR-BM | 80/20 | 9.1 | 17 | 113 |
| PR-BM | 80/20 | 9.15 | 18 | 143 |
| PR-BM | 80/20 | 9.2 | 18 | 146 |
| PR-BM | 80/20 | 9.25 | 19 | 146 |
| PR-BM | 80/20 | 9.3 | 20 | 186 |

All eighteen require matching phase fractions/compositions, log-fugacity residual
at most 1e-10 and component-balance residual at most 1e-12. These are iteration
counts, excluding stability trials and the extra EOS calls for the Jacobian
and line search; they are not a wall-clock speed claim. Unit checks exercise
coupled row-pivoting solves, domain-limited backtracking, singular/nonfinite
Jacobians, failed probes and exhausted line searches. Existing fixtures and
tolerances are unchanged.

The local review accepts A4 for this documented scope with the passing checks
in the [verification report](../validation/a4-checks.md). Merging the reviewed
change can close #30 and #26. Historical pending A4 coordinator-review wording
is superseded by this scoped review, without claiming an independent review
of this change.

## Limits retained

Finite TPD trials are evidence, not proof of global stability. There is no
new claim of critical-point accuracy, more than two phases, experimental model
accuracy or fitted interactions. A13 owns bubble/dew/envelope qualification;
A8 owns connected flowsheet execution. The separately implemented A5 PH and
A7 NRTL paths retain their own acceptance evidence. A3/A4 acceptance does not
expand any model or solver validity domain.
