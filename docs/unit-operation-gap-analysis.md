# Unit-operation palette: gap analysis against Aspen HYSYS and Aspen Plus

Status: November 22, 2026. Tracked as work package A17 ([#170](https://github.com/alexdamps-stack/engineering-property-estimation/issues/170); sub-issues #171–#218, one per row below).

This compares the editor palette (`UnitKind::ALL` in
`crates/flowsheet-core/src/unit.rs`) with the standard model palettes of
Aspen HYSYS (steady-state *Model Palette*) and Aspen Plus (*Model Library*).
Only the unit operations and palette objects are compared. Analysis tools
(case studies, sensitivity, optimization), property methods and utilities are
out of scope. The Aspen names are given so that users recognise the operation;
they do not imply compatibility with either product.

## Current palette (12 items)

| Palette item | Engine status | HYSYS equivalent | Aspen Plus equivalent |
|---|---|---|---|
| Feed | Calculated | Material stream | Material stream |
| Product | Calculated | Material stream | Material stream |
| Mixer (2 inlets) | Calculated | Mixer | Mixer |
| Splitter (2 outlets) | Calculated | Tee | FSplit |
| Heater | Calculated | Heater | Heater |
| Cooler | Calculated | Cooler | Heater |
| Heat exchanger | Calculated | Heat Exchanger | HeatX |
| Valve | Calculated | Valve | Valve |
| Separator (V/L) | Calculated | Separator | Flash2 |
| Pump | Drawn only (A14) | Pump | Pump |
| Compressor | Drawn only (A14) | Compressor | Compr |
| Expander | Drawn only (A14) | Expander | Compr (turbine mode) |

Extensions of these items (n-port mixer/splitter, vapor-fraction
specifications, UA/approach exchanger, compressor curves) are covered by A14
and A15 and are not listed again below.

## Missing palette items

Each row is a sub-issue of A17. "Prerequisite" names the capability the
project does not yet have; section 3 of the [development plan](development-plan.md)
defers reactions, distillation columns, multiple liquid phases, electrolytes
and dynamics, so most items below are blocked by one of those decisions.

### Streams and connectivity

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Energy stream | Energy stream | Heat / work stream | Energy-port model (flowsheet-and-ui.md §2.2) |

### Mixing, splitting and manipulators

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Component splitter | Component Splitter | Sep, Sep2 | — |
| Stream multiplier / duplicator | — | Mult, Dupl | — |
| Stream selector | Selector (logical) | Selector | — |
| Property-package / stream cutter | Stream Cutter | ClChng | Multiple property packages per flowsheet |

### Separators and vessels

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Three-phase separator | 3-Phase Separator | Flash3 | VLLE flash |
| Decanter (liquid–liquid) | (3-Phase Separator) | Decanter | LLE flash |
| Tank / storage vessel | Tank | (Flash2) | — |

### Heat transfer

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Air cooler | Air Cooler | (HeatX / EDR air-cooled) | Air-side model, fan power |
| Fired heater / furnace | Fired Heater | (Heater + RGibbs) | Combustion or duty-only model |
| Multi-stream (LNG/plate-fin) exchanger | LNG Exchanger | MHeatX | Zone/pinch analysis (A15) |

### Pressure change and transport

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Multistage compressor with intercooling | (Compressor + Cooler) | MCompr | A14 |
| Pipe segment / pipeline | Pipe Segment, Compressible Gas Pipe | Pipe, Pipeline | Pressure-drop correlations, A14 for elevation/work terms |
| Relief valve | Relief Valve | (Pressure relief analysis) | Relief sizing basis |

### Columns

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Shortcut distillation | Short Cut Column | DSTWU, Distl | Distillation scope decision |
| Rigorous distillation column | Distillation Column | RadFrac | Equilibrium-stage solver |
| Absorber / stripper | Absorber, Reboiled/Refluxed Absorber | RadFrac (no condenser/reboiler) | Equilibrium-stage solver |
| Three-phase distillation | Three Phase Distillation | RadFrac (3-phase) | Stage solver + VLLE |
| Liquid–liquid extractor | Liquid-Liquid Extractor | Extract | Stage solver + LLE |
| Petroleum / multi-column fractionator | Column sub-flowsheets | PetroFrac, MultiFrac, SCFrac | Stage solver + petroleum characterisation |

### Reactors

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Conversion (stoichiometric) reactor | Conversion Reactor | RStoic | Reaction sets, heats of formation |
| Yield reactor | Yield Shift Reactor | RYield | Reaction-free mass/element balance |
| Equilibrium reactor | Equilibrium Reactor | REquil | Reaction sets, equilibrium constants |
| Gibbs reactor | Gibbs Reactor | RGibbs | Gibbs minimisation, formation properties |
| CSTR | CSTR | RCSTR | Kinetics |
| Plug-flow reactor | PFR | RPlug | Kinetics, ODE integration |
| Batch reactor | (dynamics) | RBatch | Kinetics, time integration |

### Solids handling

All require a solids substream / particle size distribution model.

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Cyclone | Cyclone | Cyclone | Solids model |
| Hydrocyclone | Hydrocyclone | HyCyc | Solids model |
| Baghouse / fabric filter | Baghouse Filter | FabFl | Solids model |
| Rotary vacuum / cake filter | Rotary Vacuum Filter | Filter | Solids model |
| Simple solids separator | Simple Solid Separator | SSplit, SSplit-style Sep | Solids model |
| Crystallizer | Crystallizer | Crystallizer | Solids model, solubility |
| Centrifuge | — | CFuge | Solids model |
| Electrostatic precipitator | — | ESP | Solids model |
| Venturi scrubber | — | VScrub | Solids model |
| Crusher | — | Crusher | Solids model, PSD |
| Screen / classifier | — | Screen, Classifier | Solids model, PSD |
| Dryer | — | Dryer | Solids model, moisture basis |
| Granulator / fluidised bed | — | Granulator, Fluidbed | Solids model, PSD |
| Solids washer / CCD | — | SWash, CCD | Solids model |

### Logical operations and flowsheet control

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| Adjust / design specification | Adjust | Design Spec | Outer solver loop |
| Set (linked variable) | Set | Transfer / Calculator | — |
| Spreadsheet / calculator | Spreadsheet | Calculator | Expression evaluation |
| Balance | Balance | (Mult / Calculator) | — |
| Explicit recycle block | Recycle | (Tear streams) | A10 tear streams exist; palette object only |
| Control blocks (PID, ratio, surge, transfer function) | PID / Ratio / Surge Controller, Transfer Function | (Aspen Plus Dynamics) | Dynamic simulation |

### Extensibility

| Missing item | HYSYS | Aspen Plus | Prerequisite |
|---|---|---|---|
| User / CAPE-OPEN unit operation | User Unit Op, CAPE-OPEN | User, User2, User3, CAPE-OPEN | Stable unit API (A12 bridge) |

Symbol choices for any new palette item follow
[`design-reference/iso-10628`](design-reference/iso-10628/README.md); a
symbol in that library does not imply the unit is calculated.
