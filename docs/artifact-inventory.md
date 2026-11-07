# Artifact inventory

Updated: November 7, 2026, after the bounded A0 local search.

| Artifact | Status | Next action |
|---|---|---|
| German development plan | Preserved at repository root | Retain as historical baseline |
| English development plan | Included in `docs/development-plan.md` | Maintain as current planning document |
| Rust workspace | Initial scaffold | Complete calculation contracts under A1 |
| ThermoCore v0.3.0 Python source | **Unavailable in the searched locations**; no archive or identifiable source tree recovered | Proceed with fresh Rust work; reopen recovery on receipt of a concrete artifact |
| Earlier formulations and parameter sources | **Unavailable as primary artifacts**; planning summaries alone remain | Independently qualify formulations, conventions, parameter origins, and usage terms under A2/A3 |
| Earlier test results | **Raw tests, inputs, outputs, and logs unavailable**; historical claim of 53 tests and 111 GERG comparisons over 37 gas states only | No present pass or accuracy claim; reproduce only with recovered executable evidence |
| Rust solver/research reports | **Full earlier reports unavailable**; local documents repeat historical references | No inferred solver completion; qualify algorithms and failure behavior under later gates |
| Separate local Rust implementation | Original project checkout at commit `0985d1ac85cdbddc386467286044f20489964fc6` contains homogeneous PR/PR-BM code and 16 declared teqp reference cases | Candidate for coordinator review only; not ThermoCore recovery, not copied, executed, or accepted here |
| Component datasets and numerical references | Not included | Assemble and qualify under A2 |

The historical `/workspace/scratch/2a408b01b5b5` location is not a path available in this local repository. This inventory does not assert that linked archives or conversation content remain downloadable.

## A0 decision and acceptance limits

Choose **fresh Rust implementation rather than a port of the unavailable Python
prototype**. The bounded search is closed for this dispatch; the missing source,
formulations, parameters, raw tests, and full solver reports remain explicit gaps.
The [recovery record](artifact-recovery.md) records locations, search bounds,
inspection evidence, the separate Rust candidate's provenance, and later-recovery
compatibility requirements. The candidate may reduce later effort after review;
it does not override this worktree's scaffold baseline or the A1/A2/A3 gates.

This supplies A0's explicit-unavailability and decision evidence for coordinator
review. It does not establish numerical validation or completion of subsequent
milestones. No third-party code or data was imported and no models were implemented.
