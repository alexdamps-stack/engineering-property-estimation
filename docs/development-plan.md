# Development Plan — Engineering Tools / Property Estimation

Date: November 1, 2026. Planning baseline derived from the accessible project conversations and the initial local inventory. This is the complete English translation of the German plan, with a separate repository setup note below.

## Repository setup note

The initial repository now supplies the Rust workspace and documentation structure. All eight crates are scaffolds; no thermodynamic models, solvers, GUI, datasets, or numerical validation results are included. The historical inventory below describes the state before this repository was initialized. Creating the scaffold does not complete the numerical contracts, data work, or validation required by A1 and subsequent packages.

## 1. Goal and guiding decisions

Build a standalone Rust application for physical properties, mixtures, and steady-state process simulation. Keep the thermodynamic engine independent of both the user interface and flowsheet. First establish reproducible models and equilibrium calculations; then build a small graphical process simulator on that foundation.

Requirements taken from the most recent user contributions:

- Implement the software in Rust.
- Defer Excel integration.
- Provide a flowsheet interface as well as a calculation engine.
- Include GERG or PR-BM for hydrocarbons.
- Follow the discussed DWSIM concept of separating property packages from the simulator.
- Investigate spline-based property evaluation as a later performance improvement.

Open integration question: The Rust discussion also mentioned using the software in DWSIM. Whether a concrete DWSIM integration is required in addition to the standalone interface remains unresolved. Define the internal Rust API first; design a foreign-language interface only after resolving that requirement.

## 2. Inventory and limits of the evidence

| Area | Evidence available | Implication |
|---|---|---|
| Local project directory | Empty when initially inspected; no Git repository | Local project infrastructure and version control had to be established |
| Attachments to the planning chat | None | No source archive directly available |
| Earlier prototype | Described as Python project ThermoCore v0.3.0 in a related project conversation | Treat as prior technical work; do not assume an existing Rust engine |
| Prototype models | Earlier reports describe SRK, PR, PR-BM, NRTL, liquid-volume correlations, and a GERG adapter through teqp | Review formulations, parameters, and interfaces before porting |
| Earlier verification | A report states 53 tests and 111 GERG comparisons over 37 gas states | Historical report, not rerun during the planning session |
| Phase equilibrium | Earlier reports identify missing complete flash calculations and global stability assurance | Central development work in the new plan |
| Rust solvers and research | User requests are readable; several assistant responses are available only as content references | Their completion status cannot be inferred |

Earlier source code was generated at `/workspace/scratch/2a408b01b5b5`; the conversation links `ThermoCore_v0_3_0_Source.zip`. That historical remote path is not a source tree available in this local workspace. The inventory covers accessible conversation text and an additional look at the related earlier execution chat. It is not a fresh source review or a complete inspection of all files held in the ChatGPT project.

## 3. Recommended initial product scope

### First technical version

- Component and mixture definitions with explicit composition basis, SI units, and absolute pressure.
- PR and PR-BM as separately selectable models.
- Single-phase density, compressibility factor, fugacity coefficients, and caloric properties with a documented reference state.
- TP flash for selected hydrocarbon mixtures, including stability assessment within the defined test domain.
- A command-line application with reproducible inputs and outputs.

### First usable version

- A native Rust interface with component selection, mixture definition, and property results.
- A flowsheet supporting feeds, mixers, heaters/coolers, valves, and two-phase separators.
- PH flash as a prerequisite for valves and energy-based calculations.
- Saving and loading a versioned project format.
- An end-to-end demonstration: mix a gas stream, cool it, reduce its pressure, and separate the phases.

GERG-2008 is an important separate development track. A first PR-BM demonstration may be delivered while GERG is being implemented. A release advertised as GERG-capable requires a verified native GERG calculation path.

Deferred: Excel, dynamics, reactions, electrolytes, distillation columns, full petroleum characterization, and general multiple-liquid-phase systems. Initially represent kerosene and diesel as explicit surrogate mixtures.

## 4. Proposed architecture

The original plan proposed the following Rust crates. The repository now contains minimal scaffolds with these names; the responsibilities below remain development targets.

| Module | Responsibility |
|---|---|
| `thermo-types` | Units, component identities, mixtures, states, composition bases, and errors |
| `thermo-data` | Versioned component and interaction parameters, sources, and validity domains |
| `thermo-models` | PR, PR-BM, later GERG and NRTL; properties and derivatives |
| `thermo-solvers` | Root finding, stability assessment, TP/PH/later PS flash, and diagnostics |
| `property-packages` | Compatible combinations of models, data, and calculation capabilities |
| `flowsheet-core` | Streams, unit operations, connections, balances, and calculation order |
| `thermo-cli` | Reproducible examples, batch calculations, and regression entry point |
| `thermo-desktop` | Native interface and flowsheet editor using the same calculation API |

A property package defines the equilibrium and caloric methods and, where appropriate, a separate liquid-density method. Supported combinations and reference states must be verified. A density correlation must not silently replace the derivatives or fugacities of another thermodynamic model.

Calculation results must identify model and data versions, phases, units, validity status, and solver diagnostics. Missing parameters, unsupported phase regions, and nonconvergence must be distinguishable outcomes. After an input changes, a previously successful state must not remain displayed as current.

## 5. Work packages and acceptance criteria

| ID | Priority | Work package | Dependency | Complete when… |
|---|---|---|---|---|
| A0 | P0 | Recover and inventory earlier artifacts | None | Source, formulations, parameter sources, and tests are preserved or explicitly marked unavailable |
| A1 | P0 | Rust infrastructure and calculation contract | A0 or documented fresh implementation | Workspace, formatting, linting, and CI exist; SI, composition-basis, and error conventions are defined |
| A2 | P0 | Component data and reference cases | A1 | Initial components are unambiguous; data provenance, reference states, and valid domains are recorded |
| A3 | P0 | PR and PR-BM | A2 | Roots, fugacities, and caloric properties pass independent comparisons; BM formulation and derivatives are documented |
| A4 | P0 | Stability assessment and TP flash | A3 | Single/two-phase cases, boundary cases, and nonconvergence are handled correctly; balance and equilibrium tolerances are met |
| A5 | P0 | PH flash and energy foundations | A4 | Enthalpy reference and ideal contributions are consistent; TP→PH round trips and valve cases pass |
| A6 | P1 | Native GERG-2008 path | A2; integration after A4 | Components, parameters, and derivatives are checked against an independent reference; supported states are explicit |
| A7 | P1 | NRTL and liquid properties | A2, A4, A5 | Directed parameters, temperature dependence, vapor model, and caloric coupling are verified; missing data are visibly rejected |
| A8 | P1 | Steady-state flowsheet engine | A5 | The demonstration process runs without a GUI and closes material and energy balances |
| A9 | P1 | Native interface and flowsheet editor | API prototype after A3; calculation acceptance after A8 | A user can create, calculate, modify, save, and reload the demonstration project |
| A10 | P2 | Recycles and larger networks | A8, A9 | A defined recycle case converges reproducibly; cancellation and failure paths are checked |
| A11 | P2 | Tables and splines | Stable direct models and profiling | Speedup and error are measured over an approved domain; outside it, evaluation uses the direct model or returns a clear error |
| A12 | P2 | External integration | Stable API and clarified target scope | A concrete DWSIM or later Excel use case is demonstrated through a suitable interface |

P0 is the critical path to a reliable technical version. P1 adds the usable application and intended fluid systems. P2 follows a functioning end-to-end process.

## 6. Solver strategy

1. Define scaled residuals, physical bounds, termination conditions, and reproducible initial guesses.
2. Provide robust scalar root finding for density and outer temperature problems, accounting for domains and phase transitions.
3. Test Rachford–Rice separately for prescribed K values. This alone is not a thermodynamic flash.
4. Begin cubic-model TP flash with fugacity equality and damped successive substitution.
5. Integrate TPD stability assessment with multiple trial compositions; reassess stability after phase splitting. Finite starting-point searches must not be described as a mathematical proof of global stability.
6. Add Newton methods with step limits or line search once derivatives and comparison cases are reliable.
7. Initially implement PH flash as a controlled outer temperature iteration around TP flash; test difficult phase transitions separately.
8. Explicitly handle near-critical states, vanishing components, and nearly identical phases. A small step alone does not establish convergence.

Compare automatic differentiation, trust-region methods, homotopy, and reduction approaches on concrete difficult cases. Expand a general solver framework only when justified by demonstrated needs.

The separation of stability analysis and phase splitting is also described in the [CoolProp mixture documentation](https://coolprop.org/fluid_properties/Mixtures).

## 7. Models and validation matrix

| Fluid system | Intended calculation path | Main verification |
|---|---|---|
| Natural gas | PR-BM for initial flash cases; GERG-2008 as a separate high-accuracy model path | Density and caloric properties with identical composition and reference states; qualify the phase domain separately |
| Hydrocarbon mixtures | PR-BM | Alpha function, interaction parameters, VLE, and derivatives |
| Kerosene/diesel surrogates | PR-BM; separate liquid density only with documented coupling | Distinguish model error from parameter and surrogate errors |
| Pure methanol | Verified pure-fluid path or independent reference | Liquid density, vapor pressure, and enthalpy |
| Methanol/water | NRTL with an explicitly selected vapor model | Parameter conventions, VLE, temperature derivatives, and caloric consistency |
| Wet flue gas | Ideal gas as a limited baseline; condensation modeled separately | Wet/dry composition basis, water separation, and balances |

Initial components: methane, ethane, propane, nitrogen, CO₂, oxygen, water, and methanol. Add heavier surrogate components only with verified data.

Store inputs, composition basis, model variant, parameters, reference software version, sources, and expected quantities for each reference case. Compare PR only with an appropriate PR reference; standard PR is not a direct reference for PR-BM. Report model comparisons, implementation verification, and experimental validation separately.

Proposed initial numerical acceptance limits for regular, well-conditioned cases:

- Composition sums: absolute deviation at most `1e-10`; normalized component balance: at most `1e-8`.
- Two-phase equilibrium for active components: maximum absolute `ln(f_i,L/f_i,V)` at most `1e-7`.
- Energy balance: scaled error at most `1e-6`, with a documented absolute floor for the scaling denominator.
- Property comparisons: choose absolute and relative tolerances per quantity, before implementation, based on reference quality.

These are development targets, not achieved results or universal accuracy guarantees. Trace components and critical states require separate criteria and diagnostic cases.

## 8. GERG, interface, and splines: decision gates

**GERG:** The earlier teqp adapter does not itself provide a fully native Rust implementation. teqp and, where available, other independently verified reference code can generate comparison values outside the distributed application. Before native implementation, record the formulation scope, ideal and residual contributions, parameter provenance, and usage terms. The official [teqp GERG documentation](https://pages.nist.gov/teqp-docs/en/stable/models/GERG.html) references GERG-2004 and GERG-2008.

**Interface:** A short Rust GUI prototype must assess the canvas, connections, property editor, background calculation, and Windows packaging. Select the framework based on that prototype. The calculation model must not depend on the GUI framework. Separation of simulation and property packages follows the documented [DWSIM simulation concept](https://dwsim.org/help/main/06-simulation.html).

**Splines:** Profile direct calculations first. Use a fixed mixture in a clearly bounded single-phase domain for the first tabulation experiment. Measure value and derivative errors separately; exclude phase boundaries and critical regions. Interpolated state properties must account for thermodynamic consistency. Variable composition, memory demand, and reliable fallback rules require separate work. The earlier question about Dyad does not establish an already selected technical approach.

## 9. Indicative schedule

Nonbinding effort estimate for an experienced full-time developer; revise after A0–A3 using actual results. Durations are development weeks from the real start, not promised calendar dates.

| Stage | Approximate effort | Outcome |
|---|---|---|
| A0–A2: inventory, infrastructure, data | 1–2 weeks | Reviewable baseline and reference cases |
| A3: PR/PR-BM | 2–4 weeks | Reproducible single-phase calculations |
| A4–A5: TP/PH flash | 4–8 weeks | Technical version within a defined fluid and state domain |
| A8–A9: flowsheet and interface | 4–7 weeks | Usable PR-BM demonstration |
| A6: native GERG | Initially an additional 4–10+ weeks | Separate estimate after data and formulation review |
| A7: NRTL and liquid packages | Initially an additional 3–6 weeks | Extension to defined polar systems |

The first PR-BM application therefore represents roughly 11–21 development weeks. GERG, NRTL, critical states, and integrations expand this scope. Reliable existing artifacts may reduce effort; missing parameters and additional phase regions may increase it.

## 10. Immediate next work block

1. Preserve the earlier ThermoCore version and available solver reports; mark missing artifacts in the inventory.
2. Specify one demonstration process with its components, state domain, and required output quantities.
3. Document the PR-BM formulation, component data, interaction parameters, and caloric reference state.
4. Establish the Rust workspace, data format, and CLI input format. The initial repository supplies the workspace scaffold only.
5. Specify and implement one natural-gas case and one binary hydrocarbon case with independent references.
6. Proceed to TP flash and stability analysis only after the single-phase comparisons pass.

First numerical acceptance outcome: a reproducible Rust CLI calculates both cases, identifying the model, data version, units, and diagnostics, and passes the documented reference comparisons.

## 11. Project sources

- [Thermodynamic Solver Research](https://chatgpt.com/c/6adc861e-4de0-83ed-b04b-688703c6190d): Rust requirement and deferral of Excel; some assistant responses were accessible only as content references.
- [Rust statt C# — project conversation](https://chatgpt.com/c/6ae3c138-1120-83ed-930c-4760a2fc3ccc): standalone application, property-package/flowsheet separation, desired interface, and spline question.
- [Excel Addin Entwicklung](https://chatgpt.com/c/6ad7db8d-3c14-83eb-9b63-e96821e6ebe2): older development discussion; visible content does not reliably reconstruct implementation details.
- [Excel Tool Structure Plan](https://chatgpt.com/c/6ad67a9d-dd54-83eb-980f-b2c48b98df6a): Python prototype, example mixtures, PR-BM/GERG requirements, and reported tests.
- Related earlier execution chat `01a155d0-d5a2-7574-9f7e-7a2de227c07e`: supplementary records of the earlier Python implementation.
- Local filesystem and chat attachment checks on November 1, 2026.

External documentation supports the technical orientation. Architecture, priorities, acceptance limits, and effort are proposals in this plan. The German original is preserved at [`../ENTWICKLUNGSPLAN.md`](../ENTWICKLUNGSPLAN.md).
