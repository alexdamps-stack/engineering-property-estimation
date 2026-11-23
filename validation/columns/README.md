# A17 shortcut column references (#185)

`shortcut-v1.json` holds independent reference values for the
Fenske–Underwood–Gilliland(Molokanov)–Kirkbride shortcut distillation column
(A17 #185 method contract: alpha by successive substitution on the geometric
mean of the top dew-point and bottom bubble-point relative volatilities,
Hengstebeck–Geddes Fenske distribution, Underwood root between the key
volatilities, Molokanov fit of Gilliland, Kirkbride feed location, total
condenser, partial reboiler counted in N, reboiler duty by overall enthalpy
balance). The generator, `generate-shortcut-references.py`, cites the five
correlation sources in its docstring.

Cases (100 mol/s feed, explicit zero kij, components in the A13 set, all
temperatures 200–350 K):

| id | variant | feed (T, p) | column p | keys / recoveries | reflux |
|---|---|---|---|---|---|
| `c2c3-binary-PR-BM` | PR-BM | 295 K, 2.0 MPa (subcooled, q≈1.03) | 2.0 MPa | C2 0.98 / C3 0.97 | 1.3 R_min |
| `c2c3-binary-PR` | PR | same | 2.0 MPa | same | 1.3 R_min |
| `c1c2c3-ternary-PR-BM` | PR-BM | 260 K, 2.5 MPa (two-phase) | 2.5 MPa | C2 0.95 / C3 0.95 | R = 1.1 (≈1.54 R_min) |
| `c1c2c3-letdown-PR-BM` | PR-BM | 270 K, 4.0 MPa (liquid) | 2.5 MPa | C2 0.97 / C3 0.96 | 1.4 R_min |

Each case records its inputs, the feed's own equilibrium state, q, initial and
converged alpha (and the top/bottom alphas), iteration count, N_min, theta,
R_min, R, X, Y, N, Kirkbride ratio, N_R/N_S and feed stage, D, B, component
product flows, x_D, x_B, T_top (dew of x_D), T_D (bubble of x_D), T_bot
(bubble of x_B), h_F, h_D, h_B, h_V,top, Q_C, Q_R and the net duty. Budgets
(`tolerances`, `[absolute, relative]`) were fixed with their propagation
argument in the generator before any Rust comparison; do not widen them to
resolve a failure.

## Regenerate

In the pinned reference environment (`validation/reference-requirements.txt`;
on Windows run it in WSL2 Ubuntu):

```sh
python validation/columns/generate-shortcut-references.py --output /tmp/shortcut-candidate.json
```

Review the difference before replacing `shortcut-v1.json`. The run takes
about 15 s and is deterministic.

## Independence

Independent of the Rust engine: no Rust code is called. Fugacity coefficients
and residual properties come from teqp 0.23.2 automatic differentiation;
ideal-gas enthalpies from Cantera NASA7 via the reviewed A3 adapter
(`validation/generate-references.py`, same absolute caloric basis as the other
fixtures); the feed state from the A4 multi-start TPD flash adapter
(`generate-flash-references.py`) through the A14 phase-weighted TP state
(`rotating/generate-ps-references.py`). Bubble and dew temperatures are solved
here by inner ln K successive substitution inside an outer SciPy brentq on
temperature (a 200–350 K scan brackets the root); liquid and vapor use the
smallest and largest mechanically stable roots, as the contract specifies.
As a spot check this solver reproduced the Clapeyron-generated A13 bubble/dew
temperatures (`validation/fixtures/a13-v1.json`) to within 1e-8 K for all
bracketable mixture cases.

Not independent: the component constants and NASA7 coefficients
(`data/source-extracts/components.provisional.json`), the PR/PR-BM
formulation, the zero interaction matrix and the shortcut method contract
itself are shared with the engine. The fixture therefore checks implementation
agreement, not the accuracy of shortcut methods against rigorous columns or
experiment.
