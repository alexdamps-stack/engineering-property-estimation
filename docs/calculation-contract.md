# Initial Calculation Contract

Version 1 covers homogeneous PR76 and PR76 with Boston–Mathias supercritical continuation. It is the first implementation slice of A1–A3. This contract deliberately does not define flash, total caloric properties, or flowsheet APIs yet.

## Inputs and validation

`Temperature`, `Pressure`, and `MolarMass` have validated constructors; their stored values are finite and strictly positive. Units are K, Pa absolute, and kg/mol respectively. Component identity is an exact catalog ID; aliases and duplicate IDs are rejected.

A mixture contains 1–64 components, a declared mole or mass basis, nonnegative finite fractions summing to one within `1e-10`, and a complete, symmetric, temperature-independent binary interaction matrix. Only roundoff is normalized; arbitrary unnormalized amounts are rejected. Convert mass fractions using `x_i = (w_i/M_i)/sum(w_j/M_j)`. Zero fractions are retained, including infinite-dilution fugacity coefficients for absent components. Negative, nonfinite, missing, asymmetric, nonzero-diagonal, or `kij >= 1` interaction entries are rejected as applicable. There are no default interaction parameters.

The current evaluator imposes numerical guardrails of 80–1000 K, 1e3–1e8 Pa absolute, and acentric factors from -0.5 to 1.0. These guardrails are not a validated physical operating envelope. Near-multiple roots and mechanically unstable selected roots are explicitly unsupported. No full-domain convergence guarantee is made.

## Formulation

The equation is `p = RT/(v-b) - a(T)/(v²+2bv-b²)`. This implementation uses `R = 8.31446261815324 J/(mol K)`, `OmegaA = 0.4572355289213822`, and `OmegaB = 0.07779607390388846`, matching the exact critical constants used by teqp. They differ slightly from rounded textbook coefficients `0.45724` and `0.07780`.

Both variants use the original PR76 quadratic `kappa = 0.37464 + 1.54226*omega - 0.26992*omega²`; the PR78 high-acentric-factor extension is not implemented. Therefore the name PR-BM here does not assert exact compatibility with another simulator's complete property package or database.

For PR, `alpha = [1+kappa*(1-sqrt(T/Tc))]²` on its positive square-root branch. PR-BM uses the same expression for `T <= Tc`; above Tc it uses `alpha = exp[2*c*(1-(T/Tc)^d)]`, where `d = 1+kappa/2` and `c = 1-1/d`. The first temperature derivative is continuous at Tc; the second generally is not. Second derivatives and heat capacities are not exposed yet.

Use classical mixing: `a = sum_i sum_j x_i*x_j*sqrt(a_i*a_j)*(1-kij)` and `b = sum_i x_i*b_i`. The analytical temperature derivative includes both pure-component contributions; kij is constant.

Roots satisfy `Z³ + (B-1)Z² + (A-3B²-2B)Z - AB+B²+B³ = 0`, where `A=aP/(RT)²` and `B=bP/(RT)`. Bracket roots between cubic turning points and bisect. Only real roots with `Z > B` are returned. The caller selects `smallest` or `largest`; these names do not identify equilibrium liquid/vapor phases. A selected root must have negative `(∂p/∂v)T` and reproduce pressure with relative residual at most `1e-8`.

For `L = ln[(Z+(1+sqrt(2))B)/(Z+(1-sqrt(2))B)]`, compute standard homogeneous PR log fugacity coefficients, and departures:

- `hR = RT*(Z-1) + (T*da/dT-a)*L/(2*sqrt(2)*b)` in J/mol.
- `sR = R*ln(Z-B) + (da/dT)*L/(2*sqrt(2)*b)` in J/(mol K).

Departures are relative to an ideal gas at the same T, P, and composition, so they do not require an absolute caloric reference choice. An ideal caloric model and a common absolute enthalpy/entropy reference must be defined before total caloric properties, PH flash, or energy balances are implemented. Residual values alone cannot support those calculations.

Formulation sources: [teqp generalized cubics](https://pages.nist.gov/teqp-docs/en/latest/models/cubics.html); Peng and Robinson (1976), [A New Two-Constant Equation of State](https://doi.org/10.1021/i160057a011); Boston and Mathias (1980), *Phase Equilibria in a Third-Generation Process Simulator*, proceedings pp. 823–849. The BM continuation is also described in the Aspen-authored *Physical Property Methods and Models*, alpha-function chapter, [university-hosted manual](https://sites.chemengr.ucsb.edu/~ceweb/courses/che184b/aspenplus/PhysPropMeth%26Models.pdf).

## CLI and output

```sh
cargo run -p thermo-cli -- components
cargo run -p thermo-cli -- calculate --input examples/natural-gas.json
cargo run -p thermo-cli -- calculate --input examples/binary-hydrocarbon.json
```

Requests require all fields shown in the examples: `model` (`PR` or `PR-BM`), component IDs, fractions, `basis` (`mole` or `mass`), `temperature_k`, `pressure_pa`, `binary_interactions`, and `root` (`smallest` or `largest`). Unknown fields are rejected to catch unit/name mistakes. The two example matrices are illustrative explicit assumptions, not fitted physical data.

Successful output is JSON on stdout with schema/model/data versions, inputs converted to mole fractions, matrix, units in quantity names, all physical roots, selected properties, pressure residual, and limitations. `phase_stability` is always `not_assessed`; no equilibrium phase is reported. `validity` limits implementation evidence to reference cases, not arbitrary requested states. CLI failures produce a JSON error on stderr, no success output, and a nonzero exit code. The Rust API distinguishes invalid input, missing data, unsupported calculations, and numerical failures.

The library boundary is `property_packages::calculate(HomogeneousRequest)`; the CLI calls that API. Future flowsheet and desktop code can consume calculation capabilities without depending on the CLI. The initial API is pre-release and will be extended for A4–A5 rather than misrepresenting homogeneous calculations as a flash.

## Evidence and remaining work

The offline [reference fixture](../validation/cubic-reference.json) contains 16 cases generated by teqp, including both models, single and multiple roots, liquid-like and gas-like root selection, supercritical/subcritical alpha branches, and nonzero kij. It checks density, Z roots, log fugacity coefficients, and residual h/s. The independent generator uses Helmholtz-energy automatic derivatives and numerical density bracketing rather than Rust's root/departure code. BM is represented by teqp's Twu function with equivalent coefficients above Tc and Mathias–Copeman with only its linear coefficient below Tc. This verifies the chosen formulation, not a complete commercial PR-BM package.

Tolerances were fixed with fixture generation before comparison: relative `2e-8`; absolute Z `1e-10`, density `1e-7` in the respective output units, log phi `1e-9`, hR `1e-5 J/mol`, sR `1e-7 J/(mol K)`. These target double-precision implementation agreement for the selected well-conditioned cases and are not physical accuracy limits. Identity, derivative, composition conversion/permutation, pure-limit, failure-path, and CLI tests supplement these comparisons. The manually/scheduled 200-case state grid checks robustness and pressure reconstruction, not independent accuracy.

No experimental measurements, GERG, NRTL, stability assessment, full TP/PH flash, total h/s, cp, transport properties, flowsheet solving, or GUI are implemented in this slice. A1 is substantially established for homogeneous inputs; A2–A3 remain partial until broader data qualification, total caloric conventions, and additional coverage are complete. A4 remains gated on that work.
