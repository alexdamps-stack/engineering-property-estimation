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
| Pressure–entropy flash | A14, #226 | calculated | Overall mole fractions, absolute pressure, molar entropy on the package caloric reference | Cubic: 28 matched PR/PR-BM states in `validation/rotating/ps-v1.json` (`thermo-solvers/tests/ps.rs`). NRTL-IG and GERG-2008: TP→PS round trips and pure entropy lever rule against `nrtl-v1.json` and `gerg-equilibrium-v1.json` (`property-packages/tests/operations.rs`). EOS-CG-2021: gas-branch TP→PS round trips. | Available through `FlowsheetPackage::ps` for every package. |
| Pump, compressor | A14 | calculated | Outlet pressure or pressure increase; isentropic efficiency (0, 1] | `ng-compressor` and `water-pump` references in `ps-v1.json` (PR and PR-BM), consumed by `flowsheet-core/tests/rotating.rs` | Shaft power on the optional `power` inlet. |
| Expander | A14 | calculated | Outlet pressure or pressure drop; isentropic efficiency (0, 1] | `ng-expander` and `propane-two-phase-expander` references (lever-rule outlet) | Recovered power on the optional `power` outlet. Two-phase outlets are allowed and warned. |
| Multistage compressor | #182 | calculated | Stages 1–8, discharge pressure, stage isentropic efficiency, intercooler outlet temperature (two or more stages) | `ng-two-stage-intercooled` reference: stage outlets, stage powers, intercooler duty | Total power on `power`, heat removed on `intercooler duty`. |

All four machines use the isentropic-efficiency model. A PS flash at the
outlet pressure and the inlet molar entropy gives the ideal outlet and
\(\Delta h_s\). Pumps and compressors take \(\Delta h = \Delta h_s/\eta\),
expanders \(\Delta h = \eta\,\Delta h_s\), and a PH flash at the outlet
pressure gives the actual outlet. At \(\eta = 1\) the isentropic state is the
outlet. The shaft power \(W = \dot n\,\Delta h\) is the unit's external duty,
positive into the fluid, and closes the engine's component and energy
balances. Results report, per stage, the inlet and outlet pressures, the
isentropic and actual outlet temperatures, and the isentropic and actual
power (`UnitResult::rotating`, added behind format version 4).

Rules and failures:

- A pump inlet with any vapor fails, as does a compressor inlet with any
  liquid. A single supercritical phase is accepted by both. A pump outlet with
  vapor or a compressor discharge with liquid is a warning.
- Pumps and compressors reject an outlet pressure below the inlet pressure;
  expanders reject one above it. Equal pressures pass the inlet state through
  with zero power.
- Zero inlet flow gives an empty outlet, zero power and a warning.
- The multistage compressor uses an equal pressure ratio
  \((p_{out}/p_{in})^{1/N}\) and cools to the same temperature at constant
  pressure after every stage except the last (no aftercooler). An intercooler
  that would heat the gas or condense liquid fails; add a separator for
  condensate knockout.

Not modelled: polytropic efficiency, performance curves, speed, mechanical
and driver efficiency, pressure ratio as a specification, intercooler
pressure drop, power taken from an energy stream (for example an expander
driving a compressor). The machines call the package PS flash (`FlowsheetPackage::ps`, #226), available for PR/PR-BM/SRK, NRTL-IG, GERG-2008 and the EOS-CG-2021 gas branch. Under NRTL-IG the liquid enthalpy and entropy do not depend on pressure (no Poynting term), so a pump does no work in that model; the unit reports this as a warning. An engineering liquid density is #228.

## Heat transport

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Air cooler | #179 | calculated (partial) | Specify process outlet temperature, nonpositive process duty, or a fixed terminal approach. Assign pressure drop, ambient air temperature, dry-air mass flow, and optional nonnegative fan power. | Analytic sensible-air and fan energy balance, TP/PH agreement, approach targets and failure cases. | One material inlet and one outlet. ISO asset `s2-x2505-heat-exchanger-with-finned-tube-and-fan`. |
| Fired heater | #180 | calculated | Specify process outlet temperature, nonnegative process duty, or duty from an energy stream, with an assigned pressure drop. | Reverse of the independent process cooler reference for PR and PR-BM, plus sign, connection and zero-flow checks. | Duty-only model with one material inlet and one outlet. ISO asset `s2-2533-furnace`. |
| Multi-stream exchanger | #181 | calculated (partial) | For 2–12 material passes, assign each direction and pressure drop. Specify outlet temperature, signed duty or a fixed terminal approach on all but one pass. The last degree of freedom closes the heat balance. | Independent dilute-gas heat-capacity integral, endpoint feasibility and closure tests, and fixed-approach solve and v4 reload checks. | Full-span countercurrent core with two checked terminals. Retains ISO asset `s2-2516-heat-exchanger-of-plate-type`. |

The air cooler and fired heater carry the process molar flow and overall composition through one
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

### Fixed approach in the two-stream exchanger and air cooler

Both units offer a choice of hot-end or cold-end approach. The specified
temperature difference is a nonnegative equality at the chosen terminal.
Each setting replaces the outlet-temperature or duty specification.

| Unit | Fixed hot end | Fixed cold end |
|---|---|---|
| Two-stream exchanger | Hot inlet temperature minus cold outlet temperature. The cold outlet setpoint and pressure give duty by TP. PH at the hot outlet pressure closes the hot stream. | Hot outlet temperature minus cold inlet temperature. The hot outlet setpoint and pressure give duty by TP. PH at the cold outlet pressure closes the cold stream. |
| Air cooler | Process inlet temperature minus air outlet temperature before fan heat. The required rise from ambient determines process duty using the assigned dry-air heat-capacity rate. PH determines the process outlet. | Process outlet temperature minus ambient air inlet temperature. TP determines the process outlet and duty. The dry-air balance gives the air outlet before fan heat. |

The two-stream exchanger still takes both pressure drops as inputs. The air
cooler still takes its process pressure drop, ambient air temperature, dry-air
mass flow and optional fan power. Fan power raises the discharge temperature
after the hot-end reference point and reports on the existing fan port.

Fixed approach needs flowing process streams on both sides of the two-stream
exchanger and flowing process fluid through the air cooler. The solver rejects
an approach that reverses heat transfer, crosses the other terminal, or cannot
produce a valid outlet state. A zero approach is the limiting terminal case.
These checks cover the two terminal planes. Internal pinch, area and
predictive performance require further equipment models.

The results show the achieved approach and solved duty. Version 4 records
store the terminal choice and difference in kelvins. Earlier records keep
their thermal choice and connection identities.

### Multi-stream exchanger

Each `in k` connects through its own pass to `out k` with the same molar
flow and composition. The hot streams share one flow direction and the cold
streams run countercurrent across the same two terminal planes. Both directions
are required. Each pass takes a nonnegative pressure drop that must leave
positive absolute outlet pressure when flowing. Pressure drop is an input.

With $N$ passes, there are $N-1$ independent thermal specifications. In
outlet-temperature mode, TP determines enthalpy and duty from the setpoint.
In fixed-duty mode, PH determines the outlet. Positive signed duty heats the
pass. The one **Heat balance** pass uses PH with duty equal to minus the sum
of the other duties. The sum of duties must close to zero and the external
duty is zero. The results report the total internal transfer, the signed
duty into each pass, and the attained difference for fixed-approach passes.

**Fixed approach** replaces an outlet-temperature or duty input on the
selected pass. Choose the flowing inlet of an opposite-direction pass at the
same terminal. For a hot outlet the equality is $T_{h,out}-T_{c,in}=\Delta T$.
For a cold outlet it is $T_{h,in}-T_{c,out}=\Delta T$. The chosen gap fixes
the selected outlet temperature, TP gives its duty, and the heat balance
pass absorbs the remaining duty. Several independent fixed approaches can
reference known inlets. They cannot also prescribe duty or outlet temperature
on those passes.

The separate **Minimum terminal approach** is a feasibility bound. For
every flowing hot/cold pair the engine checks both $T_{h,out}-T_{c,in}$
and $T_{h,in}-T_{c,out}$ against it within $10^{-9}$ K. A zero-flow pass
requires explicitly zero duty and produces an empty outlet. It does not take
part in terminal comparisons or provide a reference for a fixed approach.
A zero-flow balancing pass cannot determine an outlet state.

Passing the terminal checks establishes endpoint feasibility for the
full-span core. Interior pinch, intermediate entries and exits, detailed
pressure profiles, UA, area, and plate-fin rating require a zone or
composite-curve model with heat-transfer and hydraulic data.
[Simulation and Design Methods for Multiphase Multistream Heat Exchangers](https://skoge.folk.ntnu.no/prost/proceedings/dycops-cab-2016/proceedings/media/papers/0172.pdf).

Legacy v4 exchanger records keep their count, numbering and connections.
They load with zero drops and minimum approach, one hot pass first, cold
passes after it, and the last pass selected for heat balance. Specify the
other thermal settings before solving.

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

Conversion, yield, equilibrium, Gibbs, stirred-tank, plug-flow and batch reactors
(#191–#197) run through the flowsheet engine and desktop. See the
[reactor model and validation record](a17-reactors.md) for specifications,
formation-basis balances, kinetic rate bases and supported limits.

## Columns

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Shortcut column | #185 | calculated | Light and heavy key, light-key recovery to distillate and heavy-key recovery to bottoms (each strictly 0–1, sum above 1), column pressure, total or partial condenser, and one of reflux ratio, multiple of minimum reflux or theoretical stages. | Seven PR/PR-BM cases from the independent teqp/Cantera/SciPy generator `validation/columns/` (binary, ternary with distributing methane, high-pressure let-down feed, partial condenser, specified stage count, both together), checked within the fixture's declared tolerances by `flowsheet-core/tests/shortcut_column.rs`. Project round trip, legacy reload and a headless desktop test. | Fenske–Underwood–Gilliland (Molokanov fit) with Kirkbride feed location. One feed, total or partial condenser, partial reboiler. Reports $N_{min}$, $R_{min}$, $R$, $N$, feed stage and the distillate recovery of every component. Retains ISO asset `s1-x8100-column-general`. |

The [rigorous distillation column](a17-columns-stages.md) (#186) solves
an equilibrium-stage cascade with a total or partial condenser and a partial
reboiler, a condenser and stage pressure-drop profile, and two Aspen-style
column specifications. Choose ideal trays or one
[ideal packed bed sized by HETP](a17-columns-packing.md), with numbered
feeds and liquid and vapor side draws. The [absorber or stripper](a17-columns-stages.md#absorber-and-stripper)
(#187) solves countercurrent adiabatic contact with ideal trays or one HETP
bed, yielding the outlet phase flows and the stage temperature profile.

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Rigorous distillation column | #186 | calculated (partial) | Condenser pressure plus condenser drop and a per-stage or whole-column drop; total, partial-vapor or partial vapor-and-liquid condenser; top specification (reflux ratio, reflux rate, distillate rate, condenser duty, distillate mole fraction or recovery) and bottom specification (distillate fraction, bottoms rate, reboiler duty, boilup ratio, bottoms mole fraction or recovery); 1–24 ideal trays or one HETP bed; numbered feeds, liquid and vapor side draws. | `distillation-v1.json` (PR/PR-BM ethane/propane), `extended-v1.json` (PR methanol/water at 1 atm, PR propane/methanol to 430 K) and `specs-v1.json` (partial condensers, vapor draw, pressure drop, reboiler duty, purity/boilup, condenser duty/recovery) from independent teqp/Cantera/SciPy generators; `flowsheet-project/tests/{distillation_column,columns_domain,column_specifications}.rs`. NRTL-IG methanol/water and GERG-2008 behaviour tests; liquid–liquid split, unsupported-package and infeasible-specification failures; headless desktop form test. | Any package component; package bubble/dew domain. Stage liquid–liquid check with the package VLLE flash. No tray efficiency, side heaters, subcooling or stage-temperature specifications. Specification fields added behind format version 4. |
| Absorber / stripper | #187 | calculated (partial) | Column pressure; 1–24 ideal trays or one HETP bed. | `absorber-v1.json` (PR/PR-BM ethane/propane) and the methane/CO2/methanol-solvent case of `extended-v1.json`; `flowsheet-project/tests/{absorber,columns_domain}.rs`, NRTL-IG behaviour test. | Any package component; stage temperatures in the package bubble/dew interval, pressure bounded by the package's TP/PH operations. Same liquid–liquid check. |
| Petroleum fractionator | #190 | calculated (partial) | One main column composed on the rigorous MESH solve: numbered feeds carrying characterised pseudo-components, total / partial-vapor / partial vapor-and-liquid condenser, the distillation top specification, a reboiler (with a distillation bottom specification) or stripping steam to an adiabatic bottom stage, liquid side products by molar flow or standard liquid volume flow at 60 °F (M/(SG ρw)), pressure profile, 1–24 trays or one HETP bed. | `petroleum-v1.json` from the independent teqp/NumPy/SciPy generator (8-cut crude-like TBP curve; two molar draws with reboiler; two volume draws plus a molar draw with partial condenser and pressure drop; steam-stripped with water kij = 0.5 approximation); `flowsheet-project/tests/petroleum_fractionator.rs` (reference agreement, balance closure, infeasible draws, missing characterisation, unsupported packages, round trip, drawn-only legacy record). | [a17-columns-petroleum.md](a17-columns-petroleum.md). **Not yet modelled:** side strippers and pump-arounds (need coupled cascades / stage energy ports in the MESH solve), ASTM D86/TBP cut-point specifications (no product distillation-curve calculation exists), three-phase condenser with free-water decant (realistic steam stripping with a cold drum fails as a condenser liquid split), MultiFrac/SCFrac multi-column coupling, furnace/overflash specifications and refinery correlations beyond the pseudo-component characterisation. Pseudo-components exist in the cubic packages only (250–600 K). Fields added behind format version 4. |

### Shortcut column

The column takes one feed and produces a distillate and a bottoms product.
It has one uniform pressure with no pressure drop. The condenser and
partial reboiler both operate at that pressure. The bottoms leave as
saturated liquid at their bubble point. A total condenser (the default)
returns the distillate as saturated liquid at its bubble point. A partial
condenser returns it as saturated vapour at its dew point. The column
pressure may not exceed the feed pressure: the feed is let down
isenthalpically and keeps its own enthalpy. The inputs are independent. The
operating specification is one of reflux ratio, multiple of minimum reflux
or number of theoretical stages. The keys must be two different
package components, and both must be present in the feed. The light key must
be more volatile than the heavy key at column conditions.

- **Feed condition.** The method takes
  $q=(h_{dew}-h_F)/(h_{dew}-h_{bub})$, using the saturated enthalpies of the
  feed composition at column pressure. Subcooled feeds give $q>1$ and
  superheated feeds give $q<0$.
- **Volatility and split.** Relative volatilities are package K-value ratios
  to the heavy key. Each is the geometric mean of its value at the overhead
  dew point and at the bottoms bubble point. These are iterated together with
  the Fenske total-reflux split to a fixed point, with $\max|\Delta\ln\alpha|<10^{-10}$
  in at most 100 iterations. Fenske (1932) gives $N_{min}$. The
  Hengstebeck–Geddes form gives the non-key distribution, and the key
  recoveries are met exactly.
- **Minimum reflux.** Underwood (1948) gives the minimum reflux from the root
  between the key volatilities. The keys must be adjacent: the model rejects
  a feed component whose volatility lies between them.
- **Stages and feed location.** Gilliland (1940) gives the stage count in the
  Molokanov et al. (1972) form
  $Y=1-\exp[(1+54.4X)/(11+117.2X)\cdot(X-1)/\sqrt X]$, with
  $N=(N_{min}+Y)/(1-Y)$, where $X=(R-R_{min})/(R+1)$. When the stage count
  is specified instead, $Y=(N-N_{min})/(N+1)$, the fit is inverted for $X$
  by bisection and $R=(X+R_{min})/(1-X)$. Stage counts include the reboiler
  and a partial condenser, and exclude a total condenser. Kirkbride (1944)
  gives the number of rectifying stages. The feed stage, counted from the
  top, is that number plus one; a partial condenser is stage 1.
- **Duties.** With a total condenser, $Q_C=(R+1)D(h_{V,top}-h_D)$, where the
  overhead vapour is at the dew point of the distillate. With a partial
  condenser, the condenser equilibrium is the dew point of the distillate
  composition $x_D$. The reflux is the incipient liquid $x_R$ of that dew
  point, so it sits at its own bubble point in equilibrium with the vapour
  distillate. The vapour entering the condenser has the composition
  $y_1=(x_D+Rx_R)/(R+1)$ from the constant-molar-overflow top balance, at
  its dew point. Then $Q_C=D[(R+1)h_V(y_1)-h_D-Rh_R]$. The relative
  volatilities, split and $R_{min}$ are the same for both condensers. The
  reboiler duty $Q_R$ closes the overall enthalpy balance. Both duties must
  be positive. The condenser-duty port reports $Q_C$ and the reboiler-duty
  port reports $Q_R$ to a calculated source or a sink. The net external duty
  is $Q_R-Q_C$.
- **Reported estimates.** $N_{min}$, the Underwood root, $R_{min}$, $R$,
  Gilliland $X$, $N$, the rectifying stages and feed stage, $q$, the key
  relative volatility, both duties, the temperature of the vapour entering
  the condenser, and the fraction of each feed component recovered in the
  distillate (the key and non-key distribution).

The solver fails explicitly when:

- the feed has zero flow or no state;
- a key is unknown, absent from the feed or the same as the other key;
- the column pressure exceeds the feed pressure;
- a bubble or dew point lies outside the package's bubble/dew domain, is near critical, or is a three-phase bubble point (the feed liquid splits);
- the keys are inverted or not adjacent;
- the volatilities do not converge;
- $R_{min}\le 0$, $R\le R_{min}$ or the multiple of minimum reflux is at most 1;
- a specified stage count is missing, not positive, or at most $N_{min}$;
- either duty is not positive or the bottoms are not hotter than the overhead;
- any result is not finite.

It warns when the Gilliland abscissa lies outside 0.01–0.9 or when a non-key
component distributes between the products.

**Domain and limits.** Bubble and dew points use the package's bounded
saturation solver. Under the cubic packages this covers every package component
(including water, methanol and pseudo-components) over the package temperature
interval and 1 kPa–10 MPa ([A13 status, #227](a13-status.md#227-extended-domain)).
A feed whose liquid splits into two liquids is rejected, but the stage
model does not test internal stages for a second liquid, so water–hydrocarbon
systems remain outside the validated scope. The results are design estimates.
Several things are not modelled:

- tray hydraulics, column diameter and efficiencies;
- side draws and multiple feeds;
- a mixed (vapour and liquid) distillate, or a subcooled reflux;
- non-adjacent keys, or more than one Underwood root;
- liquid–liquid phases;
- temperature approaches at the condenser and reboiler utilities.

Constant relative volatility and constant molar overflow are inherent
assumptions of the method.

**Persistence.** Version 4 records store the keys, recoveries, pressure,
operating choice (`ratio`, `minimum_multiple` or `theoretical_stages`) and
condenser. These fields were added behind format version 4 without a
version bump. A record without `condenser` loads as a total condenser, and a
total condenser is not written. An earlier empty `shortcut_column` record
loads with every input unset. Unit results carry the optional
`shortcut_column` estimate. Its `distillate_recoveries` list is empty for
results saved before it was added. Unknown fields and choices are refused.

**Rigorous columns.** The shortcut estimates for $N$, feed stage, $R$ and
$D$ can guide specifications for the [equilibrium-stage distillation
column](a17-columns-stages.md). That column and the absorber both support
ideal trays and HETP beds. Rate-based packing is a later batch; the
[petroleum fractionator](a17-columns-petroleum.md) (#190) reuses the
rigorous column with pseudo-component feeds, side products and stripping steam. Both rigorous
columns take their domain from the package: any package component and the
package's bubble/dew temperature interval, with stage properties from the
package's homogeneous-phase operation (PR/PR-BM/SRK, NRTL-IG and GERG-2008;
EOS-CG-2021 fails explicitly). After convergence every stage is tested for a
second liquid with the package's VLLE flash; a split fails the column, and a
package without VLLE (GERG-2008) publishes a warning that the check was
skipped. A detailed condenser or reboiler heat-transfer model should expose a fixed
temperature approach to a defined utility reference.

## Multiphase

| Kind | Issue | Status | Specification | Validation / reference | Notes |
|---|---|---|---|---|---|
| Liquid-liquid extractor | #189 | Calculated (PR, PR-BM, NRTL-IG; GERG-2008 and EOS-CG-2021 fail with `UnsupportedCapability`) | 1–50 ideal stages or a packed section (HETP, height or stage basis); column pressure; isothermal (temperature) or adiabatic; feed phase detect/heavy/light | `validation/multiphase/extractor-v1.json` (Python standard library, NRTL water/methanol/1-butanol, 1, 3 and 8 stages); single stage = VLLE TP flash; adiabatic single stage = VLLE PH flash; balances; project round trip | [a17-multiphase-extractor.md](a17-multiphase-extractor.md). Records added behind format version 4. |
| Three-phase distillation | #188 | Calculated (PR, PR-BM, NRTL-IG; GERG-2008 and EOS-CG-2021 fail explicitly) | 1–24 ideal trays or one packed bed (HETP); feed stage 2..N; column pressure; decanter temperature; reflux liquid light / heavy (whole liquid or reflux ratio) or fractions of both; feed-to-bottoms fraction | `validation/multiphase/three-phase-column-v1.json` (Python standard library, NRTL-IG water/1-butanol, four decanter-only cases covering all reflux modes) and `three-phase-column-vlle-v1.json` (Python standard library, two NRTL-IG water/1-butanol and two PR water/propane cases with two liquids on trays and reboiler); decanter = VLLE TP flash of the condensate; three-phase stage = VLLE PH flash of its contents; balances; energy ports; failure paths; project round trip and drawn-only record | [a17-multiphase-three-phase-column.md](a17-multiphase-three-phase-column.md). Full VLLE stage model: any tray and the reboiler carry one or two liquids, found automatically (phase-count oscillation, nonconvergence and unsupported packages fail explicitly). Per-stage phase counts and both liquids in the results and the desktop. PR water–hydrocarbon kij = 0.5 is the documented VLLE approximation. Records and results added behind format version 4. |

## Solids

Solids handling (#198–#211) is deferred. The desktop add palette omits the
Solids category and its 14 units. Projects containing those units still
load, display and save with their existing kinds, ports and specifications.
