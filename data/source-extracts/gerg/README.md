# GERG source and parameter provenance

Unmodified NIST AGA8 C++ GERG-2008 implementation, version 2.01 (April 2017),
pinned to [revision 3bdb9ab8ff317c618b0b59d1b704c2c86ddc5fce](https://github.com/usnistgov/AGA8/tree/3bdb9ab8ff317c618b0b59d1b704c2c86ddc5fce),
downloaded November 8, 2026. Files are `AGA8CODE/C/GERG2008.cpp`,
`AGA8CODE/C/GERG2008.h`, and the root `LICENSE` at that revision.
The retained [notice and disclaimer](LICENSE) permit use, modification and
distribution with the notice retained. These terms cover the retained upstream
material and derived parameter tables; they do not license original project code.

Kunz and Wagner, *J. Chem. Eng. Data* 57 (2012), 3032–3091,
[doi:10.1021/je300655b](https://doi.org/10.1021/je300655b), defines GERG-2008.
The source comments identify authors and contributing work. No REFPROP data or
binary library is used or redistributed.

`validation/prepare-gerg.py` compiles this source independently with g++ and
exports initialized tables after `SetupGERG` applies the reducing-parameter,
departure-function and ideal-reference transformations. `data/gerg-2008.json`
retains upstream array names and one-based indexing, all 21 molar masses, pure
terms, all pair parameters, departure models, and ideal terms. The source and
generator SHA-256 hashes are recorded in both the dataset and reference fixture.
The reference compiler uses `-ftrivial-auto-var-init=zero`: upstream setup
transforms unused local `bijk` cells as well as assigned cells. Zero-initializing
automatic storage makes unused exported cells deterministic. Review confirmed
that this did not change any of the 273 property/derivative reference states.

The Rust runtime embeds JSON and evaluates analytic expressions directly. It
does not compile, load, invoke or link this C++ source. Upstream mutable caches,
density root heuristics and its `1e-15` component cutoff are not ported. Positive
Rust mole fractions contribute without a cutoff; inactive components are skipped
exactly. Source availability and numerical agreement do not establish experimental
accuracy or global phase stability.
