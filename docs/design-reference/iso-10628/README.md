# Process symbol reference for AI design

This is the project's reusable visual reference for process equipment, piping,
valves, and related flowsheet symbols. It contains **all 296 symbols** from the
[Wikimedia Commons reference PDF](https://commons.wikimedia.org/wiki/File:ISO_10628-2_2012_Symbols.pdf).

Open [the searchable catalog](catalog.html) in a browser, or inspect the
[seven PNG contact sheets](previews/). The catalog works directly from disk and
does not need a server or network connection.

## Files

- `symbols/`: 296 standalone, transparent, editable SVG drawings. Open these in
  Inkscape or another vector editor; paths, curves, fills, and markers remain vectors.
- `manifest.json`: stable asset IDs, registration numbers, descriptions, sheet
  numbers, source drawing/label IDs, padded view boxes, and relative file paths.
- `catalog.html`: visual search by registration number or description and filtering
  by sheet. Click a drawing to open its SVG.
- `previews/sheet-1.png` through `sheet-7.png`: rendered catalog contact sheets for
  image-capable AI tools and human review.
- `sources/original.pdf` and `sources/sheet-*.svg`: unmodified upstream references.
- `sources/provenance.json`: download URLs and SHA-256 checksums.
- `sources/font/`: upstream osifont, its notices and license, and provenance.
- `verification.json`: render checks for every symbol and catalog checks.

| Sheet | Subject families | Symbols |
|---|---|---:|
| 1 | Vessels, tanks, columns with internals | 37 |
| 2 | Heat transfer, steam generators, furnaces, cooling towers | 33 |
| 3 | Filters, separators, centrifuges | 43 |
| 4 | Driers, size reduction, mixers, shaping machines, pumps, compressors | 50 |
| 5 | Conveyors, transport, feeders, motors, drives | 35 |
| 6 | Valves, safety fittings, piping | 60 |
| 7 | Equipment components, agitators, miscellaneous | 38 |
| **Total** | | **296** |

## Instructions for future AI design work

1. Read `manifest.json` and select a symbol by its description and registration.
   Inspect the SVG or contact sheet before using it. For example, search `pump`,
   `heat exchanger`, `column`, `valve`, or `mixer` in the catalog.
2. Reuse the selected SVG geometry. Scale uniformly and preserve distinguishing
   internals, nozzle stubs, line styles, and orientation. Do not replace equipment
   with generic pictograms or invent new internals for visual decoration.
3. Keep unit tags, operating values, selection highlights, and solver status
   badges outside the equipment drawing. Apply the project's
   [flowsheet/UI proposal](../../flowsheet-and-ui.md) when composing an interface.
4. Reference the asset ID in design notes or UI data so the choice is traceable.
   IDs include the sheet, upstream registration, and a description slug.
5. Preserve the registration strings supplied by the reference, including `C`,
   `X`, and `XC` prefixes. These are source registrations, not unit tags. The source
   includes apparent inconsistencies (such as `8063`); do not silently renumber them.
6. Connection ports are **not** defined by this library. Inspect nozzle endpoints
   and define explicit, reviewed port coordinates when implementing the editor.
   Canvas spacing, padding, and a symbol's displayed 240 px dimensions are not
   physical equipment dimensions.
7. This is a drawing reference, not evidence that a corresponding unit-operation
   model exists. Follow the implementation scope in the development plan.
8. Keep attribution and the asset license when reusing or adapting these drawings.

## Method and reproducibility

The original author already supplied editable SVGs made in Inkscape. A Node CLI
script using headless Chrome through Playwright re-authors each grouped drawing
as a separate SVG, preserving its source geometry, transforms, patterns, and
markers. It adds a padded view box, accessible title, provenance, and metadata;
sheet borders and adjacent descriptions are excluded. This is a vector extraction
and normalization, not a new freehand tracing or Excalidraw conversion.

Two source groups are out of sequence: centrifugal compressor `X8179` maps to
`g35204` on sheet 4, and support leg `C2005` to `g4862` on sheet 7. The duplicate
suction-filter label and repeated sheet-legend example are excluded. Fonts used
inside ten symbols are embedded as osifont so they render offline without relying
on an installed font. SVG-aware editors can still edit the internal text.

From the repository root, with Node.js and Google Chrome installed:

```powershell
npm ci --prefix tools/symbol-reference
node tools/symbol-reference/build.mjs
node tools/symbol-reference/verify.mjs
```

The build uses the saved sources and does not need a download. To fetch missing
upstream sources, run `node tools/symbol-reference/download.mjs` (requires curl;
uses `curl.exe` on Windows). Existing source files are retained. Check saved hashes
before accepting any intentional upstream update.

Verification checks all 296 SVGs for XML parsing, unresolved local references,
nonblank rasterization, ink at the canvas edges, source hashes, catalog search,
and sheet-filter counts. It generates the seven previews. These checks establish
rendering and inventory coverage, not independent ISO certification. The project
has not independently audited the Commons drawings against the published standard.

## Attribution and license

Original symbol sheets: **Steven Baltakatei Sandoval (Baltakatei)**,
2020, updated 2023. Source:
[ISO 10628-2 2012 Symbols on Wikimedia Commons](https://commons.wikimedia.org/wiki/File:ISO_10628-2_2012_Symbols.pdf).
Licensed **[Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/)**.

Changes in this project: individual symbol extraction, standalone SVG wrapping,
view-box padding, font embedding, searchable catalog, metadata, and contact-sheet
rendering. The adapted drawings and contact sheets are distributed under the same
CC BY-SA 4.0 license. Credit the original author, link the source and license,
indicate adaptations, and retain share-alike terms for adaptations of the assets.
This notice applies to these assets; the repository's licensing remains otherwise
as documented at its root.

osifont is by the upstream osifont contributors, maintained by hikikomori82. Its
unmodified font is embedded where needed under GPL version 3 with the upstream
font exception; see `sources/font/MANUAL`, `GPL-3.0.txt`, and its
[upstream project](https://github.com/hikikomori82/osifont). Font sources are
available upstream. The font exception permits independent documents under their
own terms.
