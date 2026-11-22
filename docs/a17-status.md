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

## Simple and logical

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Component splitter | #172 | calculated | Overhead recovery from 0 to 1 for each package component. Bottoms receives the remainder. | Analytic component split, empty outlet, signed duty and failure cases through `flowsheet-project/tests/component_splitter.rs`. Desktop form exercised headlessly. | The isothermal, zero-pressure-drop mode reuses the cubic TP flash at each flowing outlet. It calculates the net duty and can report it on its optional energy inlet. Recoveries are stored by component ID in project format 4. |
| Stream multiplier | #173 | calculated | Finite flow factor ≥ 0. Factor 2 doubles the stream into one outlet. | Analytic flow and enthalpy accounting in `flowsheet-project/tests/logical_streams.rs`. | Logical block with one inlet and one outlet. |
| Stream selector | #174 | calculated | Choose one of 2–12 connected inlets by its one-based port number. | Analytic forwarding and boundary accounting in `flowsheet-project/tests/logical_streams.rs`. | Retains ISO asset `s6-2103-valve-three-way-type-general`. |

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
