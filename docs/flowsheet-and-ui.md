# Flowsheet and User Interface — Proposal

Date: November 7, 2026; decisions updated November 8, 2026. Status: **A8 and A9 decisions recorded in section 5.** This document describes how the flowsheet engine (A8) behaves and how the desktop interface (A9) looks and works. Remaining items marked *proposed* (for example undo/redo, typed-unit input and localization) are not implemented yet. The GUI framework decision (section 4.7) is based on the prototype described in the [development plan](development-plan.md), section 8; the implementation record is [A9 status](a9-status.md).

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

Reasons: easy to understand and diagnose (a failure points to one unit), matches what process engineers know, and is sufficient for the acyclic demonstration process. [A10](a10-status.md) adds deterministic tear streams and bounded direct substitution, with closure checks, cancellation and project/desktop support. Wegstein acceleration and an equation-oriented solver remain outside the implemented scope.

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

Pumps, compressors, and expanders need a PS flash and efficiency models; they are scheduled as A14 in the [development plan](development-plan.md) (section 5.1).

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

**Prototype status (November 7, 2026).** A first `egui` prototype is in `crates/thermo-desktop`. It edits the GUI-independent topology model in `flowsheet-core`. It covers:

- The canvas: place, move, connect, pan, and zoom.
- A technical-drawing look for the canvas:
  - Orthogonal stream routing with closed, filled ISO 128 arrowheads, and stream numbers in diamonds.
  - Monochrome symbols redrawn from the ISO 10628-2:2012 symbol sheets (Wikimedia Commons, CC BY-SA 4.0; see the [process symbol reference](design-reference/iso-10628/README.md)): static mixer X2673, heat exchanger X8079 for heater and cooler, valve 2101, shell-and-tube exchanger 2511, pump 2301, compressor 2302, turbine 2571, and vessel with dished ends 2062. Proportions follow the reference. Additions: a medium arrow on heater and cooler, and project symbols for splitter and feed/product, which the sheets do not cover.
  - CAD lettering, a model-space grid with major lines, grips on the selection, and snap markers on free nozzles.
  - The dark scheme looks like CAD model space; the light scheme looks like a plotted sheet.
- The properties form, with the section 2.3 specification sets and radio buttons.
- The pre-solve statuses *Underspecified* and *Ready*.
- The message list, a stream table without values, and the status bar.
- The SI and Process metric unit sets.
- Additional units beyond the demonstration: splitter (two outlets, one split fraction), two-stream shell-and-tube heat exchanger, pump, compressor, and expander. Pumps, compressors, and expanders can be drawn and specified (pressure plus an editable 75 % default efficiency) but are marked as not yet calculable, because they need a PS flash (section 2.3).
- Stream numbers placed beside their lines, avoiding equipment, status badges, other numbers, and all stream lines.

**A9 status (November 8, 2026).** The editor now also covers background solving with progress and cancel, an auto-solve toggle, the section 4.6 statuses on canvas and in tables, results in the properties panel and stream table with tab-separated copy, feed compositions, a new-project dialog, versioned save/load, and a Windows installer. It does not cover energy streams, splitters with more than two outlets, undo/redo, typed-unit input, CSV export, or localization.

**Evaluation and decision (#52).** The criteria were assessed on the A9 editor, which is the prototype's full successor; details and measurements are in [A9 status](a9-status.md), section 4.

| Criterion | egui result |
|---|---|
| Canvas with about 100 units | 150 units and 100 streams: 0.78 ms CPU per frame (release, layout, hit-testing and tessellation, excluding GPU). Frames are drawn only on input or while a solve runs |
| Tables and forms | Grid-based specification forms, radio-button choices, composition editor and stream table are sufficient. `egui_extras` tables remain available for larger tables |
| Keyboard and clipboard | F5, Esc, Del, Ctrl+N/O/S/Shift+S; stream table copied as tab-separated text through the platform clipboard |
| High DPI | egui follows the operating-system scale factor; the canvas is vector line work that scales with zoom. Not snapshot-tested at a scale factor other than 1 |
| Windows packaging | One statically linked executable without runtime dependencies; 17 MB executable; 7.1 MB Inno Setup installer and 7.2 MB portable ZIP, built and smoke-tested on the Windows CI runner |
| Accessibility | AccessKit is enabled; the headless UI tests locate widgets through the accessibility tree. The custom-painted canvas exposes no per-unit accessible nodes; unit status is also available in the properties panel and message list |
| License and activity | MIT/Apache-2.0; actively maintained |

A Slint prototype was **not** built. The comparison was judged unnecessary for this gate because egui meets every criterion above. Slint offers no immediate-mode painter, so the PFD canvas (routing, tags, hit-testing and symbols) would have to be re-expressed as declarative path elements. Its royalty-free or GPLv3 licensing would also need a distribution decision. The geometry in `view.rs` and the whole calculation and project layer are framework-independent, so a later port stays bounded. **Decision: egui is the GUI framework.**

## 5. Decisions and remaining interface gate

A8 decisions recorded November 8, 2026 for #47. Sections 2.1–2.3 are
accepted for the headless engine with the explicit implementation scope below.
The desktop layout and framework evaluation remain A9 proposals.

| Decision | Selected behavior | Implementation / owner |
|---|---|---|
| Solution approach | Deterministic sequential-modular calculation; acyclic API retained, recycle-capable API uses tear streams and bounded direct substitution | A8 plus [A10](a10-status.md); desktop/project solves use automatic tears |
| Unit specification sets | Section 2.3 enums; TP feeds, adiabatic PH mixer, T/duty heater/cooler, PH valve, adiabatic/T/duty separator | A8; vapor-fraction specs remain A13 |
| Initial input basis | Explicit ordered mole fractions and mol/s; SI absolute pressure; one immutable cubic package for every stream | A8; existing thermo-types conversions may be applied explicitly by callers |
| Splitter in first version | Include one inlet, two outlets and one split fraction, including 0 and 1 | A8; arbitrary port counts are deferred |
| Additional prototype exchanger | Hot/cold outlet temperature or transferred duty; independent side pressure drops; no terminal crossover | A8; pumps/compressors/expanders require PS and remain unsupported |
| Separator phase endpoints | Preserve existing phases for unchanged adiabatic conditions; absent outlet has zero flow, no state, and a warning | A8; a single supercritical fluid uses the vapor port by routing convention and retains its actual label |
| Project file format | Versioned JSON with integer `format_version`, package/component/input/layout selections and result provenance; unknown versions fail; incompatible results load stale | Implemented in A9 (`flowsheet-project`, format version 1; version 2 since A10 adds recycle results, version 1 still read); A8 example JSON is a calculation report, not a project file |
| Result invalidation | Compare full topology/specification snapshot, feed compositions and immutable package identity; a failed solve returns no partially current result | A8 API; A9 adds per-unit input fingerprints (unit, feed composition, package/provenance and all upstream units) so only affected results turn stale, plus background solve and cancel |
| Default unit set | Process metric display, SI calculations | Existing prototype; A9 |
| Interface language | English first; German localization deferred | Existing prototype; A9 |
| GUI framework | **egui** (evaluation in section 4.7; Slint not prototyped, reasons recorded) | Decided in A9 (#52); calculation and project layers stay framework-independent |

The numerical API and supported domain are recorded in [A8 status](a8-status.md).
No new equipment graphics are introduced by the headless engine; the existing
ISO symbol assets, attribution, license and explicit port definitions are retained.
