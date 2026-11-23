# A17 status

Work package A17 (#170). One section per group; each group edits only its
own section. Conventions: [a17-architecture.md](a17-architecture.md).
Status values: *deferred*, *drawn only*, *calculated*, *calculated (partial)*, *not
feasible in steady state* (with the reason).

## Foundation (integrator)

| Item | Issue | Status | Notes |
|---|---|---|---|
| Energy streams, energy source/sink, heater/cooler energy ports | #171 | calculated | Tests in `crates/flowsheet-core/tests/energy.rs`; project format 4. |
| All A17 kinds in `UnitKind`, grouped palette, ISO 10628 symbols | #170 | done | Individual groups enable calculation as implementations become available. |
| Variable/numbered ports, group modules, format version 4 | #170 | done | See the architecture document. |

## Rotating

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Cubic pressure–entropy flash | A14 | calculated | Overall mole fractions, absolute pressure, molar entropy on the cubic caloric reference | 28 matched PR/PR-BM states in `validation/rotating/ps-v1.json`, consumed by `thermo-solvers/tests/ps.rs` | Available through `CubicPackage::ps`. Pump, compressor, expander and multistage compressor remain drawable. The fixture also holds outlet and power references for subsequent rotating-unit work. |

## Heat transport

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Air cooler | #179 | calculated (partial) | Specify process outlet temperature or nonpositive process duty, pressure drop, ambient air temperature, dry-air mass flow, and optional nonnegative fan power. | Analytic sensible-air and fan energy balance, TP/PH agreement and failure cases. | One material inlet and one outlet. ISO asset `s2-x2505-heat-exchanger-with-finned-tube-and-fan`. |
| Fired heater | #180 | calculated | Specify process outlet temperature, nonnegative process duty, or duty from an energy stream, with an assigned pressure drop. | Reverse of the independent process cooler reference for PR and PR-BM, plus sign, connection and zero-flow checks. | Duty-only model with one material inlet and one outlet. ISO asset `s2-2533-furnace`. |

Both units carry the process molar flow and overall composition through one
equilibrium outlet. Assigned outlet temperature uses TP and derives process duty
from the enthalpy change. Assigned duty uses PH at outlet pressure. Pressure
loss must leave positive absolute outlet pressure. The engine closes component
and process enthalpy balances using the signed process duty.

The air cooler models a dry, near-atmospheric air stream with constant
\(c_p = 1005\) J/(kg K). The air-side flow is imposed. Positive heat rejection
requires the process inlet to be warmer than ambient. The model accepts zero
terminal approach at either end. The process outlet must be at least as warm
as the air inlet, and the air leaving the exchanger must be at most as warm
as the process inlet. Sizing or predicting operating duty
from fan speed would require exchanger conductance, geometry and fan curves.
Humidity and evaporative cooling are outside the model.
[WMO dry-air heat capacity, summarized in Table 1](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=929976).

The process duty \(Q\) is negative for heat rejection. Fan power \(W_f\) is
assigned independently and dissipates into the air downstream of the
exchanger. The discharge air gains \(-Q + W_f\). Fan power stays outside
the process enthalpy balance, and a calculated energy source connected to
the existing fan power port reports only \(W_f\). A stopped process has an
empty outlet and zero process duty, with independent fan operation.

The fired heater's optional duty port reports process duty with a calculated
energy source in temperature or assigned-duty mode. In energy-stream mode
the connected stream supplies that duty. The model reports heat absorbed by
the process fluid and rejects heat removal. Applying a nonzero assigned
heat input to an empty process stream fails. Fuel use, combustion efficiency,
flue gas, emissions and furnace design require additional models.

Existing version 4 records retain their kinds and port identities. Old
records load with the outlet temperature unset and zero pressure drop.
The air cooler also requires assigned ambient temperature and airflow.

## Simple and logical

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Component splitter | #172 | calculated | Overhead recovery from 0 to 1 for each package component. Bottoms receives the remainder. | Analytic component split, empty outlet, signed duty and failure cases through `flowsheet-project/tests/component_splitter.rs`. Desktop form exercised headlessly. | The isothermal, zero-pressure-drop mode reuses the cubic TP flash at each flowing outlet. It calculates the net duty and can report it on its optional energy inlet. Recoveries are stored by component ID in project format 4. |
| Stream multiplier | #173 | calculated | Finite flow factor ≥ 0. Factor 2 doubles the stream into one outlet. | Analytic flow and enthalpy accounting in `flowsheet-project/tests/logical_streams.rs`. | Logical block with one inlet and one outlet. |
| Stream selector | #174 | calculated | Choose one of 2–12 connected inlets by its one-based port number. | Analytic forwarding and boundary accounting in `flowsheet-project/tests/logical_streams.rs`. | Retains ISO asset `s6-2103-valve-three-way-type-general`. |
| Tank | #178 | calculated | Adiabatic, specified outlet temperature, or signed duty, with a specified pressure drop. | Matched isenthalpic outlet reference, closure, failure and v4 reload checks in `flowsheet-project/tests/tank.rs`. | One inlet and one bulk outlet. Retains ISO asset `s1-x2063-tank-with-conical-roof-and-flat-bottom`. |
| Balance | #215 | drawn only | General envelope around selected existing streams. | Requires the prerequisites below. | The engine's unit balance checks remain available on all calculable units. |

The splitter prescribes component recoveries. Both outlets can contain vapor
and liquid depending on their composition at the inlet temperature and pressure.
An empty outlet has no thermodynamic state. Adiabatic and prescribed-duty
splitting would require a separate outlet-temperature solve. The existing
palette drawing retains ISO reference asset `s3-x2616-impact-separator`.

The multiplier and selector can exchange material with the modeled boundary.
The multiplier's external material is the inlet scaled by `factor − 1` at its existing state.
The selector forwards its chosen inlet as-is, including zero flow, and records
the other inlets as external outflow. Each unselected inlet carries its own
enthalpy across the boundary, so no mixing or flash is needed. These units
have zero heat duty. Results record the signed total and component molar
flows and their carried enthalpy, with positive flow into the flowsheet.
The same component and energy closure checks include these terms in both
acyclic solves and final recycle closure.

Both specifications and the optional result field are stored in project format
4. Older records preserve their ports and load with the factor or selection
unset. Specify the missing value before solving. The engine revision marks
results from older balance rules stale.

The steady-state tank carries the inlet flow and overall composition to a
single equilibrium outlet, including both phases when the fluid splits.
It uses TP for an assigned outlet temperature and PH for an assigned duty
or adiabatic pressure loss. Adiabatic service at unchanged pressure preserves
the inlet phase state. The optional energy inlet reports the calculated or
specified duty to a calculated energy source. Empty inlets give empty
outlets and a zero-duty warning. Inventory, level and time-dependent
operations are outside the steady-state model. Existing v4 tank records
default to adiabatic operation and zero pressure drop.

A Balance block can form an envelope around streams already connected
elsewhere and solve for an unknown stream variable. An audit-only first
step needs signed, non-consuming material and energy stream references
and result fingerprints for those references. It would report residuals
after recycle closure and account for heat duty and logical material exchange
exactly once. Solving unknowns additionally needs degree-of-freedom and rank
checks, explicit ownership of each calculated variable, and a reflash and
atomic update path that works within recycle sequencing.
[Aspen Plus User Guide, chapter 25](https://web.ist.utl.pt/ist11038/acad/Aspen/AspUserGuide10.pdf).

## Reactors

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|

## Columns

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|

## Multiphase

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|

## Solids

Solids handling (#198–#211) is deferred. The desktop add palette omits the
Solids category and its 14 units. Projects containing those units still
load, display and save with their existing kinds, ports and specifications.
