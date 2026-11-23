# A17 reactors

Status: implemented:
- conversion reactor (#191)
- yield reactor (#192)
- equilibrium reactor (#193)
- Gibbs reactor (#194)
- CSTR (#195)
- plug-flow reactor (#196)
- batch reactor (#197)

All seven are part of the A17 palette extension (#170).

Every reactor has one material inlet and one material outlet. Energy streams
are not modelled yet, so the duty is reported on the unit result. The outlet
is one equilibrium stream at the outlet conditions; route it to a separator
for a vapor/liquid split. New reactors start without an adiabatic default,
because an adiabatic exothermic reaction readily leaves the package
temperature domain.

| Unit | Aspen Plus / HYSYS / DWSIM | Chemistry | Thermal | Other inputs |
|---|---|---|---|---|
| Conversion reactor | RStoic / Conversion | Reactions, each with **one of**: fractional conversion of a key reactant, or molar extent (mol/s); applied in order | Adiabatic, outlet T, or duty | Pressure drop |
| Yield reactor | RYield / Yield Shift | Relative mass yields; inert components | as above | Pressure drop |
| Equilibrium reactor | REquil / Equilibrium | Independent reactions, each with K **from Gibbs energies** or a **correlation** `ln K = a + b/T + c ln T + d T` | as above | Fugacity basis, temperature approach, pressure drop |
| Gibbs reactor | RGibbs / Gibbs | All non-inert components; element conservation | as above | Fugacity basis, temperature approach, inerts, pressure drop |
| CSTR | RCSTR / CSTR | Kinetic reactions | as above | Volume (m³), pressure drop |
| Plug-flow reactor | RPlug / PFR | Kinetic reactions | **One of**: isothermal T, adiabatic, coolant (T_c and U) | Length and diameter (m), linear pressure drop |
| Batch reactor | RBatch | Kinetic reactions | **One of**: isothermal T, adiabatic | Reaction time (s); constant pressure |

## Shared foundations

### Data

Element counts come from the source formulas. Standard formation enthalpy and
absolute entropy come from the retained NASA7 constants a6 and a7 of the
lower interval:
- `Component::standard_formation_enthalpy`;
- `Component::standard_entropy`;
- `Component::standard_gibbs_energy(T)`, which is `dHf + h(T) − T (S0 + s(T))`
  on the cp integration the package already uses.

H₂ and CO were added for synthesis-gas chemistry in a separate extract,
`data/source-extracts/components-reacting.provisional.json`. It is produced
by `validation/prepare-source-data.ps1 -Set reacting` from the same pinned
CoolProp and Cantera sources and hashes. The script still reproduces the base
extract byte for byte (apart from line endings).

The base extract is unchanged, so every fixture that hashes it stays valid. A
package that uses an extension component reports the data version
`2026-11-07.1+reacting-2026-11-23.1`. All other packages keep their
provenance.

### Energy balance

Stream enthalpies are zero for each pure ideal gas at 298.15 K and 1 bar, so
reactors add the standard formation enthalpy of their generation:

    n_out h_out − n_in h_in + Σ_i g_i ΔHf°_i = Q

Adiabatic and duty specifications of the equilibrium, Gibbs and CSTR reactors
find the outlet temperature with a bracket that starts at the inlet
temperature and expands in doubling steps (10, 20, 40 K …), followed by Brent's
method. The energy balance is then closed exactly with a PH flash at the
converged composition. If the solution lies outside the package interval, the
message says whether the outlet would be too hot or too cold. The plug-flow
and batch reactors carry the formation-basis enthalpy along the integration
instead.

### Balances

The unit balance subtracts the generation from the component residuals and
adds `Σ g ΔHf°` to the energy residual. Recycle closure re-evaluates the same
terms.

Every reactor except the yield reactor enforces element conservation at
1e-8 relative per element. The mass residual is reported, because the pinned
CoolProp molar masses are not atom-consistent beyond about 1e-6 relative. The
yield reactor enforces mass and warns when elements do not balance.

Equilibrium reactors start from a trace floor of 1e-14 of the feed on every
species, which keeps the minimization interior. Species that can neither form
nor react leave at that trace level, which is far inside the 1e-8 element
budget.

### Specification checks

Before any property calculation, the engine checks:
- missing values, duplicate components, zero coefficients, reactions without a
  reactant or product, and a conversion key that is not a reactant;
- order components and the catalyst density whenever a rate is per catalyst
  mass;
- that every named component exists in the package;
- that every reaction conserves every element;
- that equilibrium reactions are linearly independent.

## Conversion and yield reactors

**Conversion reactor.**
- Reactions are applied in order. A conversion refers to the amount of the key
  reactant that enters that step (the HYSYS rank convention).
- A reaction that would consume more of a component than is present fails with
  the reaction number and the component.

**Yield reactor.**
- The non-inert feed mass is distributed over the yields, and inerts pass
  through.
- Mass closes by construction. When elements differ by more than 1e-6
  relative, the result carries a per-element warning, as in RYield.

## Equilibrium and Gibbs reactors

Both reactors minimize the Gibbs energy of one fluid phase,

    G/RT = Σ n_i (μ°_i/RT + ln(p/p°) + ln x_i + ln φ_i),

subject to linear invariants `B n = B n0`:

| Reactor | Invariants `B` | Standard potentials `μ°` |
|---|---|---|
| Gibbs | One row per element over the non-inert components, plus one unit row per inert | `g°_i(T)` from the data |
| Equilibrium | The null space of the stoichiometry of the listed reactions | The least-norm vector that reproduces every `ln K_j`: `μ°/RT = −ν (νᵀν)⁻¹ ln K` |

`ln K` comes from Gibbs energies (`−Σ ν g°/RT`) or from the user correlation.
Conditions are on the fugacity basis, with the ideal gas at 1 bar as the
standard state.

The minimizer works on the element potentials and ln N. It uses Newton's
method with row equilibration and backtracking and converges the scaled
constraint residuals to 1e-13. The outer loop updates the fugacity
coefficients by successive substitution until they change by less than 1e-11.

Options:
- **Fugacity basis:** *Package* uses the fugacity coefficients of the selected
  PR/PR-BM model for the reacting fluid, from a single homogeneous phase on the
  minimum-Gibbs root. *Ideal gas* uses partial pressures.
- **Temperature approach:** equilibrium is evaluated at the outlet
  temperature plus the approach, as in Aspen. The outlet state stays at the
  outlet temperature.
- **Inerts (Gibbs):** excluded components pass through unchanged.

The equilibrium reactor reports one extent per reaction. The Gibbs reactor
reports none.

Limits:
- Only one fluid phase reacts. When the outlet flash splits into vapor and
  liquid, the result carries a warning that equilibrium was computed for a
  single phase. Multiphase chemical equilibrium is not included.
- Equilibrium reactions must be independent; redundant ones are rejected with a
  message.

## Kinetic reactors

### Rate expressions

Each kinetic reaction has:
- stoichiometric coefficients;
- a reacting phase (vapor or liquid; a single supercritical phase counts as
  either);
- a concentration basis;
- a rate basis;
- a rate law.

The options follow Aspen Plus and DWSIM.

| Option | Choices |
|---|---|
| Concentration basis | Molar concentration (mol/m³), mass concentration (kg/m³), mole fraction, mass fraction, partial pressure `x p` (Pa), fugacity `x φ p` (Pa), each in the reacting phase |
| Rate basis | Per m³ of reacting phase (mol/(m³ s)), or per kg of catalyst (mol/(kg s)) with a catalyst bulk density (kg catalyst per m³ reactor) |
| Rate law | **Power law** `r = k_f Π q^a − k_r Π q^b` (reverse optional), or **LHHW** `r = k (K₁ Π q^a − K₂ Π q^b) / (Σ_i K_i Π q^c)^m` |
| Rate constants | Arrhenius `k = k0 (T/T_ref)^n exp(−E/R (1/T − 1/T_ref))`, or `k0 T^n exp(−E/RT)` without a reference temperature; LHHW constants `exp(a + b/T + c ln T + d T)` |

Rates are evaluated on a TP equilibrium state (isothermal) or a PH state
(adiabatic or coolant) of the local contents.

Phase volume uses no slip: a volume-based rate is multiplied by the reacting
phase's share of the contents' volume, computed from the phase molar volumes.
A reaction whose phase is absent has zero rate. A zero concentration with a
positive order gives zero rate; with a negative order it is an error.

Units of k0 follow from the basis, the orders and the rate basis; they are
entered in SI and are the user's responsibility.

### CSTR

The CSTR solves `ξ_j = V R_j(T, p_out, x_out)` for the reaction extents by
damped Newton iteration with a finite-difference Jacobian, keeping all flows
nonnegative. It starts from no reaction. If the direct solve fails, it
continues in volume from zero over 20 steps, which follows the low-conversion
branch deterministically.

Adiabatic and duty specifications use the temperature bracket described
above. If several steady states exist (adiabatic exothermic CSTR), the one
nearest the inlet temperature is found. Other branches are not searched.

### Plug-flow reactor

The plug-flow reactor integrates component flows, extents and the
formation-basis enthalpy flow along the length:
- `dn/dz = A Σ ν R`;
- `dH/dz = U π D (T_c − T)` for a coolant, 0 when adiabatic;
- the temperature is fixed when isothermal.

The pressure falls linearly with the specified drop. Integration uses adaptive
Dormand–Prince 5(4) at relative tolerance 1e-9 with at most 20 000 accepted
steps.

The explicit integrator does not handle stiff problems: very fast kinetics, or
a very large wall coefficient relative to the heat-capacity flow. Such cases
end with a "more than 20000 steps" or "step underflow" message rather than a
result. An implicit integrator is a follow-up.

### Batch reactor

The batch reactor runs at constant pressure (the inlet pressure) and integrates
the contents of one mole charged over the reaction time. The reacting volume is
`N v(T, p, x)` from the package.

In the steady-state flowsheet, the outlet is the time-averaged product. The
outlet molar flow is the inlet flow times the final amount per mole charged.
Cycle and down times do not change the averaged composition and are not
inputs.

The isothermal duty includes bringing the charge to temperature. Constant
volume operation would need a TV flash and is not included.

## Supported domain and limits

- **Temperature.** The PR/PR-BM package covers 200–600 K (0.01 Pa to 100 MPa),
  and every reactor outlet, every temperature at which an equilibrium is
  evaluated (including the approach) and every integration state must lie in
  it. Many adiabatic exothermic reactions leave it, and so does steam reforming
  at its usual temperatures. These fail visibly.
- **Components.** The embedded set is methane, ethane, propane, nitrogen,
  carbon dioxide, oxygen, water, methanol, hydrogen and carbon monoxide.
- **Hydrogen.** PR with zero kij is a known weak representation of hydrogen's
  volumetric behavior. There is no component-specific accuracy fixture for H₂
  or CO; they use the same verified cubic implementation as the other
  components.
- **Formation data.** Ideal-gas NASA7 formation properties agree with tabulated
  CODATA/JANAF/ATcT values within 0.6 kJ/mol. This is source agreement, not
  experimental qualification.
- **Not included.**
  - Multiphase chemical equilibrium.
  - Liquid activity bases (γ x) for kinetics.
  - Constant-volume batch.
  - Implicit (stiff) integration.
  - Ergun pressure drop.
  - Profiles in the result.
  - Energy ports.
  - Multiple steady-state search.
  - Selectivity specifications.

## Verification

| Fixture | Generator | Covers | Budgets |
|---|---|---|---|
| `validation/fixtures/reactions-v1.json` | `generate-reaction-references.py` | Conversion and yield cases; 8 formation enthalpies | Flows 1e-12 relative; duties 1e-9 relative |
| `validation/fixtures/reactors-v1.json` | `generate-reactor-references.py` | Equilibrium, Gibbs and kinetic cases; standard Gibbs energies of 7 components | Flows 1e-8; duties 1e-7; Gibbs energies 1e-6 J/mol |

Neither generator calls Rust, and CI regenerates both fixtures and compares
them with the committed copies.

`reactors-v1.json` has three equilibrium cases, all on the ideal-gas fugacity
basis and checked against Cantera `equilibrate("TP")` at rtol 1e-14 over the
same species:
- a Gibbs reactor with shift and methanation at 560 K;
- the water-gas shift with a 20 K approach;
- methanol synthesis with a correlation K, solved by SciPy `brentq`.

It also has three kinetic cases:
- a CSTR with a reversible power law on partial pressures (`brentq`);
- a catalytic plug-flow reactor with a power-law shift on mole fractions and
  LHHW methanation on partial pressures, with pressure drop (`solve_ivp`
  DOP853 at rtol 1e-13);
- a constant-pressure batch reactor on molar concentration, with the reacting
  volume from a separate Python PR implementation whose exact critical
  constants are derived independently.

The Rust tests remove the PR residual enthalpy of the solved inlet and outlet
and compare the remaining ideal-gas duty.

Further tests in `crates/flowsheet-core/tests/reactor_models.rs`:
- the equilibrium condition `Σν ln f = ln K` with package fugacities, and a
  measurable nonideal shift at 19 bar;
- duty round trips and adiabatic closure for the equilibrium, Gibbs and CSTR
  reactors;
- plug-flow thermal modes: U = 0 equals adiabatic, and a hot coolant heats the
  gas and increases conversion;
- dependent reactions, a missing catalyst density, unknown order components,
  and the required inputs of each default.

`crates/flowsheet-project/tests/reactors.rs` round-trips every option of every
reactor through format 4 and stores the extents of a solved CSTR.

## Project format

Format version 4 adds:
- the seven reactor kinds and their specification records:
  - reaction `coefficients`;
  - conversion or molar extents;
  - mass yields and inerts;
  - equilibrium `constant` (`gibbs_energy` or a `correlation`), `fugacity` and
    `temperature_approach_k`;
  - kinetics `reactions`, with `phase`, `basis`, `rate_basis` and `rate` (a
    power law with Arrhenius terms and orders, or LHHW), plus
    `catalyst_bulk_density_kg_m3`;
  - plug-flow and batch thermal choices;
- the optional `results.units.<id>.reaction` record.

Versions 1–3 are read unchanged. A reactor in a file that declares an older
version is rejected.

## Interface

The palette lists all seven reactors with traced ISO 10628 symbols.
ISO 10628-2 has no general reactor symbol, so vessels stand in, with
attribution retained:

| Reactor | Symbol |
|---|---|
| Conversion | Fixed-bed column `s1-x8015` |
| Yield | Fluidized-bed column `s1-x8014` |
| Equilibrium | Jacketed vessel `s1-x2069` |
| Gibbs | Spherical vessel `s1-2063` |
| CSTR | Jacketed vessel with agitator `s1-x8006` |
| Plug flow | Jacketed pipeline `s6-x409` |
| Batch | Vessel with heating/cooling coil `s1-x8000` |

The CSTR and batch nozzles sit on the vessel wall inside the jacket or coil.

The properties panel edits:
- each reaction as one coefficient per package component, with an equation
  preview;
- the conversion or extent;
- the K source and its correlation;
- inert components;
- for kinetics: the phase, basis, rate basis and rate law; Arrhenius
  parameters, with constants in scientific notation; orders; reverse and
  adsorption terms; and the catalyst density.

Results show the duty, the standard heat of reaction and the extent of each
reaction.
