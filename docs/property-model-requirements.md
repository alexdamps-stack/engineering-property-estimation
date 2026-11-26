# Requirements for new property models and packages

Status: November 24, 2026. This document applies when anyone adds or extends:

- a thermodynamic model (equation of state, activity model, corresponding states);
- a property package that combines a model with data and solvers;
- a property method used beside the equilibrium model: transport properties,
  an engineering liquid density, or pseudo-component characterisation;
- a capability of an existing package, for example a PS flash for NRTL-IG.

It builds on the [calculation contract](calculation-contract.md), the
[testing strategy](testing-strategy.md), the
[reference qualification rules](reference-qualification.md) and the
[architecture](architecture.md). Where this document and those disagree, the
calculation contract governs units, identities, metadata and failures. The
planned capabilities are tracked in #230 (#224–#229).

## 1. Fixed constraints

1. **Native and in-process.** Every model, solver and correlation runs in
   Rust inside the engine. No external program, service, Python runtime or
   foreign thermodynamic library may be called at run time (see the
   architecture's *Fundamental engine constraint*). Reference tools stay
   under `validation/`.
2. **Layering.** Parameters and provenance belong in `thermo-data`
   (transcribed sources under `data/`). Formulations belong in
   `thermo-models` and iterations in `thermo-solvers`. The supported
   model/data/solver combination is assembled in `property-packages`, and the
   user's selection and its persistence live in `flowsheet-project`.
   `flowsheet-core` consumes properties only through
   `FlowsheetPackage` and `CalcContext`.
3. **No silent substitution.** A package never falls back to another model
   for an operation it lacks. That includes the cubic package and an ideal
   gas. It never fills missing interaction parameters with zeros, and it never
   extrapolates silently. Unsupported operations return
   `CalculationError::UnsupportedCapability`. States outside the domain return
   `OutsideValidityDomain`. Missing data returns `MissingParameter`.
4. **One package per flowsheet** until #175 (stream cutter) defines how two
   packages exchange streams and reference states.

## 2. The package interface

### 2.1 Required by every package

Every flowsheet package implements `property_packages::FlowsheetPackage`:

| Operation | Requirement |
|---|---|
| `identity` | A new identity per built package, so that results computed with another selection become stale. |
| `component_ids` | Database identifiers in package order. Aliases are never treated as different components. |
| `temperature_interval` | The search interval used by PH (and PS) solves. It lies inside the model's declared domain. |
| `validate_composition` | Rejects compositions that the package cannot evaluate, for example components outside a qualified set. |
| `tp` | Equilibrium at T, p and overall composition. |
| `ph` | Equilibrium at p and molar enthalpy on the package's caloric reference, including the pure-component and narrow-boiling saturation path (phase fraction at saturation temperature). |

Every unit in the core palette needs TP and PH: mixer, splitter, heater,
cooler, heat exchanger, valve, separator, tank, fired heater, air cooler,
multi-stream exchanger and component splitter.

### 2.2 Returned states

`EquilibriumState` and each `EquilibriumPhase` must carry:

- the overall and per-phase mole fractions in package order, the phase
  fractions, and a phase label (vapor, liquid, supercritical) that makes no
  claim about stability;
- for each phase, the required `PhaseProperties`: T, p, Z, molar volume, mass
  density, ln φ per component, residual enthalpy, molar enthalpy, molar
  entropy and cp;
- `ModelMetadata`, `DataMetadata`, `CaloricReference` and `ValidityStatus`, plus
  the flag for an inferred regression interval where it applies;
- convergence evidence (fugacity, component, energy or entropy residuals as
  applicable), iterations, and diagnostics with stable codes.

Optional phase properties (cv, speed of sound, viscosity, thermal
conductivity) are `None` when the model does not supply them. The interface
reports them as *Unavailable*, and projects store them as `null`. They are
never estimated from another model in order to fill the field.

### 2.3 Capabilities that unit operations need

A package declares which operations it supports, and returns
`UnsupportedCapability` for the others. A unit operation fails at the first
operation that the selected package does not support, with an error naming
the unit and the operation.

| Operation | Needed by | PR / PR-BM / SRK | NRTL-IG | GERG-2008 | EOS-CG-2021 | Issue |
|---|---|---|---|---|---|---|
| TP, PH flash | All units | ✔ | ✔ | ✔ | ✔ (homogeneous gas branch) | — |
| PS flash | Pump, compressor, expander, multistage compressor | ✔ | ✔ (no liquid pressure dependence) | ✔ | ✔ (gas branch) | #226 |
| Homogeneous phase at T, p, x with root selection | Rigorous column stage equations, reactor phase volumes | ✔ | ✔ (liquid / ideal gas) | ✔ (exterior roots) | ✔ (gas root only) | #226 |
| ln φ and K-values | K-values, equilibrium and Gibbs reactors (package basis) | ✔ | ✔ (γ p^sat / p) | ✔ | ✔ ln φ (gas); K-values — | #226 |
| Bubble and dew temperature and pressure | Shortcut column, condenser and reboiler states, vapor-fraction specifications | ✔ (A13 domain) | ✔ (modified Raoult) | ✔ (A13 algorithm, dry set, 0.1–10 MPa) | — | #226, #227 |
| Formation properties on the package reference | All reactors | ✔ | ✔ (same ideal-gas reference) | — (native reference, no offset) | — (embedded offsets) | #226 |
| Three-phase (VLLE/LLE) stability and flash | Three-phase separator, decanter, three-phase column, extractor | — | — | — | — | #224 |
| Viscosity, thermal conductivity, surface tension | Pipe segment, relief valve, exchanger rating, column hydraulics | ✔ (8 components, [domain](transport-properties.md)) | ✔ (same) | ✔ (same) | ✔ (CO2 only; H2 has no data) | #225 |
| Engineering liquid density | Liquid-phase kinetics, pump hydraulics, relief and pipe sizing | — | — | — | — | #228 |
| Pseudo-components from assay data | Petroleum fractionator | ✔ ([status](pseudo-components-status.md)) | — | — | — | #229 |

Since #226 these operations are `FlowsheetPackage` methods (`ps`,
`homogeneous`, `ln_fugacity_coefficients`, `k_values`,
`saturation_temperature_interval`, `saturation_pressure`,
`saturation_temperature`, `formation_components`). Their default
implementations return `UnsupportedCapability`, so a package supports only
what it overrides. No unit uses `CalcContext::cubic_package()` or the
`FlowsheetPackage::cubic()` escape hatch any more; both remain only for
compatibility. **New unit code must not add `cubic()` dependencies.** Add the
operation to `FlowsheetPackage` instead (#226), implement it for each package
that can support it, and return `UnsupportedCapability` from the others. When
a package gains an operation, update this table and the units' status
records in [a17-status.md](a17-status.md).

## 3. Thermodynamic consistency

1. **One caloric reference per package.** Enthalpy and entropy use one
   declared `CaloricReference` for all phases and operations. TP, PH and PS
   must be mutually consistent: PH at the enthalpy of a TP state, and PS at
   its entropy, reproduce that state.
2. **Formation basis.** Reactors add standard formation enthalpies to stream
   enthalpies, which requires that each pure ideal gas has h = s = 0 at
   298.15 K and 1 bar. A package with a different reference (for example
   GERG-2008's) must either supply the offset that converts it, with a test
   against `Component::standard_gibbs_energy`, or reject reactors explicitly.
3. **Derivative identities.** Fugacity coefficients, residual enthalpy and
   entropy, cp and cv, and speed of sound must come from the same potential,
   with tests against numerical derivatives (for example
   ln φᵢ from ∂(nA_res)/∂nᵢ, and h = g + Ts). Activity models also need a
   Gibbs–Duhem test and a check of the ∂(gᴱ/T)/∂(1/T) enthalpy identity.
4. **Volume and density.** EOS volume is the thermodynamic volume. A separate
   liquid-density correlation is an engineering property only. It is reported
   separately and never replaces EOS volume in fugacities or caloric
   properties. A volume translation must define its effect on fugacities and
   caloric properties consistently (see *Liquid density and caloric coupling*
   in the calculation contract and #228).
5. **Side properties.** Transport properties, surface tension and engineering
   densities never feed back into equilibrium or caloric results. They carry
   their own method identity, data source and validity domain on the phase
   result.

## 4. Data and provenance

- Components are identified by the embedded database ID. New components need
  CAS number, formula, molar mass, the constants the model uses, ideal-gas
  heat capacity, and formation properties if reactors may use them. Each
  value needs its source and licence notice under `data/source-extracts/` or
  `data/`.
- Binary parameters are an explicit, complete matrix with a named source, or
  an explicitly labelled approximation. Missing pairs are an error.
- Each package states its temperature and pressure domain, its component set
  and its qualified phase regions. `ValidityStatus` is assessed against that
  declared domain and is never inferred from finite results. A domain weaker
  than the source states is flagged, as for NRTL regression intervals.
- Model and dataset versions go into `ModelMetadata` and `DataMetadata`. A
  change to formulation, parameters or data bumps the version, so that
  `Provenance` and the result fingerprints make saved results stale.

## 5. Project and interface integration

- Add a `PackageChoice` variant with a stable serialized label in
  `flowsheet-project/src/package.rs`. Selection fields are stored with serde
  defaults behind the current format version, unknown fields are refused, and
  a format version is bumped only by the integrator (see
  [a17-architecture.md](a17-architecture.md), section 8).
- Building the package validates components, parameters and domain before
  any unit is calculated. A package that cannot be built is a project error,
  not a unit failure.
- The desktop package selector names the model, the component set and the
  domain. Stream results show *Unavailable* for optional properties the
  package does not supply.
- `thermo-cli` examples (`examples/*.json`) cover the new package's
  homogeneous and flash entry points.

## 6. Verification

The rules of the [testing strategy](testing-strategy.md) apply. In
particular:

1. **Independent references.** Expected values come from tools that share
   parameters but no code with the engine (teqp, Clapeyron.jl, Cantera,
   CoolProp/REFPROP where permitted, SciPy solves), generated by scripts
   under `validation/`. The generator records versions and input hashes, and
   the Rust tests read committed fixtures offline.
2. **Matched formulations.** Compare against a reference with the same
   formulation and the same reference state. Report implementation agreement,
   model-to-model comparison and experimental validation separately.
3. **Tolerances first.** Declare tolerances per quantity before comparing.
   Never widen a tolerance or regenerate a fixture to make a failure pass;
   review the formulation, data and reference state instead.
4. **Coverage per capability:**
   - homogeneous properties at gas, liquid and supercritical points;
   - TP flash, including stability and single- and two-phase states;
   - PH (and PS once it exists) round trips from TP states, plus pure and
     narrow-boiling saturation;
   - saturation points and envelopes inside the declared domain;
   - failure paths: missing parameters, domain violations, non-convergence,
     and unsupported operations.
5. **Flowsheet evidence.** At least one unit using the new package is solved
   through the engine. Component and energy balances must close (relative 1e-8
   for components, 2e-9 scaled energy), and the project must survive a
   save/reload round trip.
6. **Standard checks.** Changes to PR/PR-BM or NRTL models, caloric data or
   equilibrium solvers also run the Clapeyron standard check
   (`validation/clapeyron-checks.md`). Every change runs `cargo fmt`,
   `cargo clippy -D warnings` and `cargo test --workspace`, in WSL where
   Windows blocks local builds.

## 7. Definition of done

A new package or capability is done when:

- [ ] it implements `FlowsheetPackage`, or extends it additively, without new `cubic()` dependencies;
- [ ] every unsupported operation returns `UnsupportedCapability`, and the table in section 2.3 is updated;
- [ ] metadata, caloric reference, validity domain and provenance are complete and versioned;
- [ ] the consistency tests of section 3 pass;
- [ ] independent fixtures and offline tests cover the capabilities in section 6, including failure paths;
- [ ] a flowsheet case solves with closed balances and reloads from a saved project;
- [ ] the package selector, stream results and CLI examples are updated;
- [ ] a status record (`docs/<package>-status.md`, or the matching `aNN-status.md`) states the domain, references, tolerances and limits, and the [architecture](architecture.md) crate summary mentions the package;
- [ ] the installed desktop prototype is refreshed and verified (see `AGENTS.md`).

## 8. Notes for the planned capabilities

| Issue | Additional requirements |
|---|---|
| #224 VLLE/LLE | Multiphase stability with trial phases of aqueous and hydrocarbon type. Phase labels for two liquids (for example by density or by the key component). Treat plait points and vanishing phases explicitly. Unit outlets map phases by a documented rule. |
| #225 Transport | Report the method per phase (gas, liquid, dense). Mixture rules and the pure-component data source are recorded. Outside the domain, the property is *Unavailable* rather than extrapolated. Units that need a missing property fail explicitly. |
| #226 Operations for every package | PS uses the same structure as PH, including the saturation path. Fugacity and K-values come from the package's own phases (for NRTL-IG, γ and the vapor model). Saturation solvers are bounded by each package's domain. |
| #227 Bubble/dew domain | Each added component set and range needs new independent saturation references. Near-critical states fail explicitly. |
| #228 Liquid density | Choose between a volume translation and a separate correlation, following section 3.4. Report the EOS density and the engineering density separately. Reactor concentrations must say which density they use. |
| #229 Pseudo-components | Characterisation correlations and their inputs (Tb, SG, M) are stored with the project. Pseudo-component IDs are stable and cannot collide with database IDs. Interaction parameters remain explicit. |
