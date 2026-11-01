# Independent numerical validation

Reserved for independent reference inputs, expected outputs, provenance, and comparison tools. No numerical validation is implemented yet.

The [testing strategy](../docs/testing-strategy.md) defines the reference-case requirements, independence rules, numerical scenarios, tolerances, and reporting conventions. Keep small qualified fixtures in Git for offline crate integration tests. Reference-generation tools may use external software; ordinary Rust tests must not require it. Record and review every regeneration of expected values.

Record the reference implementation and version, model formulation, composition, parameters, units, reference states, domain, and per-property absolute/relative tolerances. Keep comparison code separate from the production implementation to avoid comparing a model with itself.

Report implementation verification, model-to-model comparisons, and experimental validation separately. Historical prototype test counts are not current repository results. The initial scaffold CI only checks buildability, formatting, linting, and documentation.
