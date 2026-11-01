# Independent numerical validation

`cubic-reference.json` contains 16 independent teqp 0.23.2 reference cases for homogeneous PR/PR-BM roots, density, fugacity coefficients, and residual enthalpy/entropy. They include subcritical/supercritical alpha branches, multiple roots, both root selections, and nonzero interaction parameters. This is implementation verification, not experimental validation or phase-stability evidence.

The [testing strategy](../docs/testing-strategy.md) defines the reference-case requirements, independence rules, numerical scenarios, tolerances, and reporting conventions. Keep small qualified fixtures in Git for offline crate integration tests. Reference-generation tools may use external software; ordinary Rust tests must not require it. Record and review every regeneration of expected values.

Record the reference implementation and version, model formulation, composition, parameters, units, reference states, domain, and per-property absolute/relative tolerances. Keep comparison code separate from the production implementation to avoid comparing a model with itself.

## Running comparisons

Ordinary `cargo test --workspace --locked` runs offline reference comparisons along with derivative, identity, input-validation, composition-invariance, and CLI tests. It does not require Python. Fixture metadata contains the exact reference model JSON, input values, versions, and property-specific tolerances. See the [calculation contract](../docs/calculation-contract.md) for the supported implementation slice and tolerance rationale.

The broader robustness grid is separate:

```sh
cargo test -p property-packages --test state_sweep --locked -- --ignored --nocapture
```

It covers a fixed Cartesian grid of 200 model/composition/T/P/root selections. Output lists every failure and input; success establishes pressure reconstruction and finite results on that grid, not independent property accuracy. The scheduled/manual validation workflow runs the same command and includes output in its job summary once the workflow is on the default branch.

## Reference regeneration

Use Python 3.14 with the pinned packages in `requirements-reference.txt` in an isolated environment:

```sh
python -m venv .venv-reference
# Activate .venv-reference using your platform's shell command, then:
python -m pip install -r validation/requirements-reference.txt
python validation/generate_cubic_references.py
```

Add `--refresh-data` only when intentionally regenerating the selected CoolProp constants too. Review changes to both the data and expected values; generation never calls the Rust engine. teqp obtains properties through Helmholtz derivatives and brackets density roots numerically. The source implementations are independent, but both use the deliberately matched formulation and input constants; shared formulation/input errors remain possible. BM uses an algebraically equivalent Twu alpha representation above Tc and a linear Mathias–Copeman representation below Tc. This is not a reference for every other simulator's PR-BM variant.

Initial checks on Ubuntu/WSL with Rust 1.99.1: 15 ordinary tests passed (including all 16 reference cases and explicit critical-root rejection); the separately invoked 200-case grid passed; Clippy and documentation builds passed. Cross-platform CI checks each proposed change. Report current runs separately from historical prototype counts; no experimental measurements are included.
