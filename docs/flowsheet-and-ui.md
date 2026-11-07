# Flowsheet and User Interface — Proposal

Date: November 7, 2026. Status: **proposal for review.** This document suggests how the flowsheet engine (A8) behaves and how the desktop interface (A9) looks and works. Decisions marked *proposed* need confirmation; the GUI framework is chosen only after the prototype described in the [development plan](development-plan.md), section 8.

## 1. Target users and design principles

The [process symbol reference](design-reference/iso-10628/README.md) supplies 296
editable SVG equipment, valve, and piping drawings with a searchable catalog and
AI reuse guidance. Use it when selecting canvas graphics; symbol availability
does not imply an implemented unit-operation model.

The users are **process engineers**. They know commercial and open simulators (Aspen HYSYS/Plus, DWSIM, UniSim) and expect a process flow diagram (PFD), forms for specifications, stream tables, and engineering units. Design principles that follow from this:

- **Familiar conventions over novelty.** Canvas with unit symbols and streams, a properties panel, and a stream table, as in established simulators.
- **Traceability.** Every result shows its property package, model and data version, and solver status. Values are never shown without a status.
- **No silent wrong numbers.** Stale, out-of-domain, and failed results are visibly different from valid results. `NaN` is never displayed.
- **Keyboard and clipboard friendly.** Engineers move numbers between the simulator, spreadsheets, and reports all the time.
- **Small, dependable scope first.** The demonstration process must work end to end before more unit operations are added.

## 2. Flowsheet engine contract (A8)

### 2.1 Solution approach — *proposed: sequential-modular*

Unit operations are solved one at a time in topological order of the stream connections. Each unit receives fully defined inlet streams and computes its outlets with the property package.

Reasons: easy to understand and diagnose (a failure points to one unit), matches what process engineers know, and is sufficient for the acyclic demonstration process. Recycles (A10) are added later with tear streams and direct substitution or Wegstein acceleration. An equation-oriented solver is out of scope until a concrete need appears.

### 2.2 Streams — *proposed*

A material stream holds:

| Group | Content |
|---|---|
| State | Temperature, absolute pressure, total molar flow, composition (mole fractions) |
| Phase result | Vapor fraction, phase compositions and flows, phase identification |
| Properties | Molar enthalpy and entropy, molar mass, mass flow, per-phase density, Z, and further properties the package supports |
| Provenance | Property package, model and data versions, calculation status, solver diagnostics |

Energy streams carry a duty in W with sign convention "positive into the unit". One property package applies to the whole flowsheet initially; mixing packages between streams is deferred, because it needs explicit conversion rules.

### 2.3 Unit operations and specifications — *proposed*

Each unit has a fixed number of required specifications. The defaults below cover the demonstration process.

| Unit | Ports | Required specification (choose one where listed) | Defaults | Calculation |
|---|---|---|---|---|
| Feed | 1 out | T and p; flow as molar, mass, or standard volume; composition as mole or mass fractions. After A13 also p with vapor fraction (bubble/dew) | — | TP flash |
| Mixer | n in, 1 out | Outlet pressure: lowest inlet pressure, or a value not above it | Lowest inlet pressure; adiabatic | Sum of flows and enthalpies, then PH flash |
| Heater / cooler | 1 in, 1 out, energy | One of: outlet T, duty, or (after A13) outlet vapor fraction; plus pressure drop | Pressure drop 0 | TP or PH flash; duty reported |
| Valve | 1 in, 1 out | One of: outlet pressure or pressure drop | — | Isenthalpic PH flash |
| Two-phase separator | 1 in, 2 out (vapor, liquid), optional energy | Adiabatic, or one of outlet T / duty; plus pressure drop | Adiabatic, pressure drop 0 | Flash at outlet conditions; if single phase, the empty outlet carries zero flow and a warning |
| Splitter (*suggested addition*) | 1 in, n out | n−1 split fractions | — | Same state in all outlets |

Pumps, compressors, and expanders need a PS flash and efficiency models; they are deferred.

### 2.4 Specification check

Before solving, every unit is checked for missing, duplicate, or contradictory specifications (for example both outlet temperature and duty on a heater, or a valve outlet pressure above the inlet pressure). The result is a per-unit status (section 4.6) and a message that names the unit and the field. The solver does not start on an underspecified flowsheet; it reports what is missing.

### 2.5 When calculation runs — *proposed*

- **Manual solve is the default** (toolbar button and `F5`). An **auto-solve** toggle recalculates after each confirmed edit.
- Calculation runs in a background thread. The interface stays responsive and shows a progress indicator with a **Cancel** button.
- Any edit marks the affected unit and all downstream streams as **stale** immediately (section 4.6), even before the next solve.

### 2.6 Project file — *proposed*

A text-based, versioned format (JSON or RON) with a top-level `format_version`. It stores components, property package selection, units, streams, specifications, canvas layout, and display settings. Results are stored together with the code, model, and data versions that produced them; on load, results from a different version are shown as stale until recalculated. Text format keeps project files reviewable in Git.

## 3. Main window

```
┌───────────────────────────────────────────────────────────────────────────────────┐
│ File  Edit  View  Flowsheet  Tools  Help                                          │
│ [New][Open][Save] │ [Undo][Redo] │ [▶ Solve F5] [■ Cancel] [Auto-solve ☐] │ Units: Process metric ▾ │
├──────────┬───────────────────────────────────────────────────────┬────────────────┤
│ PALETTE  │ CANVAS                                                │ PROPERTIES     │
│          │                                                       │ V-100 (Valve)  │
│ Streams  │  [Feed 1]─S1─┐                                        │ Status: Solved │
│  Material│              ├─[M-100]─S3─[E-100]─S4─[V-100]─S5─[D-100]─S6 (vapor)
│  Energy  │  [Feed 2]─S2─┘              │Q-1                     │   └─S7 (liquid)│
│ Units    │                                                       │ Spec:          │
│  Mixer   │                                                       │ ◉ Outlet p     │
│  Heater  │                                                       │   [ 20.0 ] bar(a)
│  Cooler  │                                                       │ ○ Pressure drop│
│  Valve   │                                                       │ Results:       │
│  Separat.│                                                       │ T out  -28.4 °C│
│  Splitter│                                                       │ Vap.fr. 0.912  │
├──────────┴───────────────────────────────────────────────────────┴────────────────┤
│ [Stream table] [Messages (2)] [Diagnostics]                                         │
│           S1        S2        S3        S4        S5        S6        S7           │
│ T  °C     25.0      40.0      31.2     -10.0     -28.4     -28.4     -28.4          │
│ p  bar(a) 50.0      50.0      50.0      49.5      20.0      20.0      20.0          │
│ ...                                                                               │
├───────────────────────────────────────────────────────────────────────────────────┤
│ Package: PR-BM (data v0.3) │ Solved 14:02:11, 6 units, 0 errors, 1 warning        │
└───────────────────────────────────────────────────────────────────────────────────┘
```

| Area | Content |
|---|---|
| Toolbar | File actions, undo/redo, solve/cancel, auto-solve toggle, active unit set |
| Palette | Material and energy streams, unit operations; drag onto the canvas or click to place |
| Canvas | PFD with standard-like unit symbols, auto-generated tags (`M-100`, `E-100`, `V-100`, `D-100`, `S1` …), zoom, pan, snap to grid |
| Properties panel | Specification form and key results for the selected unit or stream; radio buttons for mutually exclusive specifications |
| Bottom panel | Stream table, message list, solver diagnostics |
| Status bar | Property package and data version, last solve time, error and warning counts |

The numbers in the sketch are illustrative layout content, not calculated results.

## 4. Proposals for the open items

### 4.1 Demonstration workflow

1. **New project** → dialog: choose components and property package (PR-BM for the demonstration).
2. Place two **Feeds**, a **Mixer**, a **Cooler**, a **Valve**, and a **Separator**; connect ports by dragging from outlet to inlet.
3. Enter feed conditions and compositions in the properties panel; unspecified fields are highlighted.
4. Specify cooler outlet temperature and valve outlet pressure; the separator stays adiabatic.
5. **Solve** → units turn to "solved"; read the stream table; check messages.
6. Change the valve pressure → downstream streams turn stale → solve again.
7. **Save**, close, **reopen**: inputs and layout are restored; results are shown as valid only if versions match.

This workflow is also the A9 acceptance script.

### 4.2 Result presentation

- **Stream table** with streams as columns and properties as rows. Users choose and order rows; phase-specific rows (vapor/liquid) are grouped.
- **Unit summary** in the properties panel: duty, pressure drop, vapor fraction, outlet temperature.
- **Copy** of any table selection as tab-separated text, so it pastes directly into a spreadsheet. **CSV export** of the stream table.
- **Phase envelope plot** for a stream after A13, with the stream's operating point marked.
- A **printable report** (PDF or HTML) with inputs, results, versions, and messages is a later addition.

### 4.3 Units of measure

- All calculations run in SI internally; conversion happens only at input and display.
- Predefined unit sets: **SI** (K, Pa, mol/s, W) and **Process metric** (°C, bar, kg/h, kmol/h, Nm³/h, kW) as the default for process engineers. Users can override individual quantities.
- **Gauge versus absolute pressure is always explicit** (`bar(a)` / `bar(g)`); gauge values require a project-wide atmospheric pressure setting.
- **Standard volume units show their reference conditions** (for example `Nm³/h (0 °C, 1.01325 bar)`, `Sm³/h (15 °C, 1.01325 bar)`), following the A1 contract.
- Every input field shows its unit and accepts a typed unit (`5 barg`, `-10 C`).

### 4.4 Interaction basics

- **Undo/redo** for all model edits (specifications, placement, connections), implemented as commands on the flowsheet model, not in the GUI layer.
- Copy/paste and duplicate of units including specifications; delete with confirmation if streams would be disconnected.
- Keyboard shortcuts for common actions (`F5` solve, `Ctrl+Z/Y`, `Del`, `Ctrl+C/V`).
- Auto-generated, editable tags following common PFD conventions.

### 4.5 Language and number format — *proposed*

- Interface text in **English first**, with all strings externalized (for example with the `fluent` crate) so a **German** translation can follow without code changes.
- Number format follows the operating-system locale for display and input (decimal comma in German locales). Clipboard copy uses the same locale as the user's spreadsheet application, which normally matches the system locale; CSV export uses a selectable separator.

### 4.6 Status and error display

Each unit and stream has one status. Status is shown by **color and icon**, so it is readable without color perception.

| Status | Canvas | Values in tables |
|---|---|---|
| Underspecified | Grey outline, `!` badge | Missing inputs highlighted |
| Ready (not solved) | Normal outline | Empty |
| Solved | Normal outline, ✓ badge | Normal |
| Warning (for example out of declared domain, single phase in separator) | Amber outline, ⚠ badge | Value shown with ⚠; tooltip names the reason |
| Failed (non-convergence, missing parameters) | Red outline, ✕ badge | No value; message names the cause |
| Stale (input changed since last solve) | Dashed outline | Grey italic values marked "stale" |

The message list links each message to its unit or stream; clicking selects it on the canvas. Diagnostics show solver iterations, residuals, and the model/data versions for the selected item.

### 4.7 GUI framework candidates

The prototype in development plan section 8 compares at least two candidates on the same small task: canvas with five units and connections, a properties form, a stream table, a background solve with cancel, and an installer for Windows.

| Candidate | Strengths | Points to check |
|---|---|---|
| `egui` (with a node-graph crate such as `egui-snarl`) | Fast to build tool-style UIs, good canvas control, simple threading model | Immediate-mode redraw cost for large flowsheets; look and feel; accessibility |
| `iced` | Retained, Elm-style architecture; clean separation of state and view | Custom canvas effort; maturity of text input and tables |
| `Slint` | Declarative UI, native-looking widgets, good tooling | Licensing terms for the intended distribution; custom canvas effort |

Evaluation criteria: canvas performance with about 100 units, table and form quality, keyboard and clipboard support, high-DPI rendering, Windows packaging and installer size, accessibility, license, and project activity. Proposed starting point: prototype with `egui` first, compare with `Slint`.

## 5. Open decisions

| Decision | Proposal | Needed by |
|---|---|---|
| Solution approach | Sequential-modular; tear streams in A10 | A8 start |
| Unit specification sets | Table in section 2.3 | A8 start |
| Splitter in first version | Include | A8 start |
| Project file format | Versioned JSON or RON | A8 |
| Default unit set | Process metric | A9 prototype |
| Interface language | English first, German prepared | A9 prototype |
| GUI framework | After prototype comparison | Before A9 implementation |
