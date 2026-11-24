# EOS-CG-2021 coefficients

The binary reducing and departure coefficients come from Neumann, Herrig,
Bell, Beckmüller, Lemmon, Thol and Span, *International Journal of
Thermophysics* 44, 178 (2023), and the referenced Beckmüller et al. hydrogen
mixture fit. Neumann et al. and the corresponding supplement are distributed
under CC BY 4.0. [Article and supplement](https://doi.org/10.1007/s10765-023-03263-6)
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

The pure-fluid references are listed with the coefficients. The original
CO₂ and normal hydrogen reference pair uses Span and Wagner (1996) and
Leachman et al. (2009).
[CO₂](https://doi.org/10.1063/1.555991)
[Normal hydrogen](https://doi.org/10.1063/1.3160306)
[Hydrogen mixtures](https://doi.org/10.1063/5.0040533)

The embedded JSON extracts all 16 pure-fluid coefficient sets and 120 binary
interactions from the cited torch-flash inventory at revision
`a65b099b891df2c8f315cfbfc57ba820d6d4ca0f`. Component names in pair records
were replaced by indices and departure blocks were stored as arrays. The
evaluator is implemented in Rust.
[Pinned inventory](https://github.com/ThermoPhase-FCSRG/torch-flash/blob/a65b099b891df2c8f315cfbfc57ba820d6d4ca0f/src/torch_flash/data/models/multiparameter/eos-cg-2021.yaml)
[Upstream data provenance and notices](https://github.com/ThermoPhase-FCSRG/torch-flash/blob/a65b099b891df2c8f315cfbfc57ba820d6d4ca0f/THIRD_PARTY_NOTICES.md)

The offline reference CSV holds independent teqp 0.23.2 outputs from the same
pinned repository. It supplies pressure and residual Helmholtz energy for 36
fixed-density points.
[Reference dataset](https://github.com/ThermoPhase-FCSRG/torch-flash/blob/a65b099b891df2c8f315cfbfc57ba820d6d4ca0f/tests/data/eoscg2021_co2_h2_teqp_reference.csv)

The flowsheet package evaluates mechanically admissible gas roots at 330–350 K
and 0.1–2 MPa. Some mixtures have no gas root in that range. Enthalpy and
entropy retain the embedded pure-fluid offsets. Condensable mixtures can
return metastable states. Phase equilibrium and transport are unavailable.
