# A17 status

Work package A17 (#170). One section per group; each group edits only its
own section. Conventions: [a17-architecture.md](a17-architecture.md).
Status values: *drawn only*, *calculated*, *calculated (partial)*, *not
feasible in steady state* (with the reason).

## Foundation (integrator)

| Item | Issue | Status | Notes |
|---|---|---|---|
| Energy streams, energy source/sink, heater/cooler energy ports | #171 | calculated | Tests in `crates/flowsheet-core/tests/energy.rs`; project format 4. |
| All A17 kinds in `UnitKind`, grouped palette, ISO 10628 symbols | #170 | drawn only | Every group kind fails explicitly in preflight. |
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

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
