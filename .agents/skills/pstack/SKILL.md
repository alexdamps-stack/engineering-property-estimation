---
name: pstack
description: Investigate, sweep for issues, fix bugs, and implement or refactor features in the Engineering Tools Property Estimation repository using its adapted pstack workflow.
---

# Project pstack workflow

Use this repository's contracts and executable evidence to route work. Resolve
paths below from the repository root. Read `AGENTS.md` and `CONTRIBUTING.md`
before choosing checks. Keep small documentation or mechanical edits focused;
they need diff and link checks rather than a full development ceremony.

## Route the task

| Request | Workflow and expected result |
|---|---|
| Explain or investigate | Trace the relevant entry point, callers, data, and tests. Answer with source evidence and distinguish observed behavior from inference. |
| Sweep or identify issues | Map the requested scope, inspect high-risk behavior and boundaries, and run targeted reproductions. Report findings without changing product code unless fixes are requested. |
| Fix a bug | Reproduce, trace the first incorrect state, make the smallest justified fix, then verify the same reproduction and relevant regression checks. |
| Feature or refactor | Identify domain types and affected callers, define observable acceptance criteria, implement in units that can be checked, and verify the resulting behavior. |

For a sweep, prioritize model/data conventions, solver failure handling,
conservation, persistence, stale-result invalidation, and interface contracts
as relevant to the requested scope. Check supported-domain statements before
calling a numerical result wrong. A confirmed finding needs a concrete trigger,
expected and actual behavior, impact, and a precise file/line reference.
Deduplicate symptoms with the same cause. Separate unverified suspicions and
coverage gaps from confirmed issues; state what was inspected even if no issue
is found. Do not invent findings or treat a missing tool as a product defect.

## Implement and verify

Inspect the working tree before edits and preserve unrelated changes. Choose
data shapes and boundary contracts before adding logic; inspect callers when
changing a function or type. Trace failures to their cause rather than adding
speculative guards. Add behavior-based regression evidence for reproducible
bugs according to `docs/testing-strategy.md`.

Run targeted checks while iterating. For code changes, finish with the ordinary
workspace formatting, Clippy, test, and documentation commands in `README.md`.
Use the pinned Rust toolchain. Apply these additional checks by scope:

- **Numerical:** for PR/PR-BM, NRTL/gamma-phi, caloric-data, or equilibrium-solver
  changes, also run `julia --project=validation/julia validation/julia/check.jl`
  with Julia 1.11.7. Follow `validation/clapeyron-checks.md` and
  `docs/testing-strategy.md`. Keep independent reference generation separate
  from the Rust engine. Review fixture or tolerance changes explicitly; never
  replace expected results or widen tolerances merely to pass a failing check.
- **Flowsheet/UI:** read `docs/flowsheet-and-ui.md` and the ISO 10628 README and
  manifest specified in `AGENTS.md` before choosing graphics. Use the existing
  `egui_kittest` interaction/snapshot tests and relevant desktop self-test or
  demonstration from `README.md`. Review intentional snapshot changes visually.
  Headless tests establish their covered behavior; claim live interaction only
  when actually exercised.
- **Other surfaces:** use the owning integration or packaging checks when those
  artifacts change. Documentation-only work needs valid references and a clean
  diff, without rerunning numerical checks for unchanged calculations.

Report the outcome, evidence, and material limits. Record commands actually run
and their results; identify missing prerequisites or inconclusive checks. For
fixes, include before/after behavior. Review the final diff before declaring done.

## Session compatibility

This is a local adaptation, not the full Cursor plugin. Use tools actually
available in the session. Cursor's `Task`, `poteto-agent`, slash commands,
`cursor-team-kit`, model slugs, and model-configuration rules are not prerequisites.
Follow session rules for delegation and inherit the parent model unless an
explicitly authorized choice applies. This skill does not authorize external
messages, publishing, commits, PRs, or merges beyond the user's task.

## Provenance and updates

Adapted from [pstack](https://github.com/cursor/plugins/tree/792440e8a560a4a3f63748cdee70491ecd412982/pstack)
at revision `792440e8a560a4a3f63748cdee70491ecd412982`:
`skills/poteto-mode/SKILL.md`, its Investigation, Bug fix, and Feature playbooks,
and `skills/principle-prove-it-works/SKILL.md`. Upstream currently has no
`sweep` skill; the sweep route above is a project-local procedure.

The upstream MIT notice is preserved in [LICENSE](LICENSE). Review upstream
changes against repository contracts and session capabilities before updating
this adapter. Keep the reviewed revision and attribution current.
