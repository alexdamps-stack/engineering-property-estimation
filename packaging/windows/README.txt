Property Estimation Flowsheet
=============================

A steady-state flowsheet editor for process engineers: draw a process,
specify feeds and units, calculate it with the Peng-Robinson (PR) or
Peng-Robinson/Boston-Mathias (PR-BM) property package, and save the project.

Start
-----
- Start menu: "Property Estimation Flowsheet", or "... (demonstration)" to
  open the mix / cool / expand / separate example.
- File > New project: choose components and PR or PR-BM.
- F5 solves in the background; Cancel (Esc) stops a running calculation.
- Ctrl+S saves a versioned *.flowsheet.json project file; Ctrl+O opens one.

Result status
-------------
Every unit and stream shows one status: underspecified (!), ready, solved
(green check), warning (amber triangle), failed (red cross), or stale
(grey dots, dashed outline: inputs changed since the shown values were
calculated). Stale values are grey and italic and are never reported as
current. Results stored in a project file are shown as current only if the
inputs, package, engine, model and data versions all match.

Scope and limits
----------------
- Acyclic flowsheets only; recycles are not supported yet.
- Pumps, compressors and expanders can be drawn but not calculated.
- The binary interaction parameters of new projects are an explicit
  zero-kij approximation; no measured mixture accuracy is claimed.

Command line
------------
  thermo-desktop.exe [project.flowsheet.json]
  thermo-desktop.exe --demo
  thermo-desktop.exe --self-test report.json   (headless check; exit code 0 = passed)
  thermo-desktop.exe --version version.txt

Third-party notices are in the "notices" folder.
