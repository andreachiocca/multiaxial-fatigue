# Handoff: automatic Abaqus and Ansys importers for multifat (new session prompt)

Paste everything below the line into a new session opened in
`C:\Users\tgros\PycharmProjects\multiaxial-fatigue`.

---

You are continuing work on `multifat` (multiaxial fatigue criteria), branch
`feat/fe-and-time-history-inputs` (local, not pushed; do not push, do not open a PR, do not add
co-author or AI attribution lines to commits). Read `AGENTS.md`, `docs/audit-checkpoint.md` (top
section "Milestone, FE CSV exchange format"), `docs/fe-csv-format.md`, `docs/fe-input-format.md` and
decisions `d086`-`d092` in `docs/decision-log.md` first. Language: English, direct, no emoji.

## Task

Add automatic importers that read solver output directly and write the FE CSV exchange format
(`multifat-fe-csv 1.0`, written with `multifat.fe.csvwriter.FeCsvWriter`, then importable with
`multifat.fe.csvio.import_fe_csv` and assessable with `multifat.history.nodes.assess_nodes`):

1. **Abaqus text exports**: Field Output Reports (`.rpt`, written by Abaqus/CAE "Report > Field
   Output") and, if worthwhile, `.dat` printouts, with node coordinates from the `.inp` (`*NODE`) or
   from `COORD` columns in the report.
2. **Ansys text exports**: `PRNSOL`/`NLIST` listings (`.out`, `.lis`, `.txt`) and Workbench
   "Export" text tables.
3. **Native Ansys `.rst`** through an OPTIONAL extra (for example `ansys` in `pyproject.toml`) using
   the MIT-licensed `ansys-mapdl-reader` (pyansys). Check first that it installs and imports on this
   machine (Python 3.12, numpy 2.x); if it does not, try an isolated environment, report what fails,
   and stop that sub-task rather than pinning old numpy for the whole project. Verify, with a
   public `.rst`, in which coordinate system `nodal_stress` is returned (global Cartesian is
   required; rotate or refuse otherwise) and the order of its components.

Native Abaqus `.odb` needs Abaqus and is out of scope; say so in the docs.

## Hard requirements

* Nothing is guessed: units (`stress_unit`, `length_unit`), `strain_convention`/`strain_unit`,
  `analysis_type`, `step_role`, `cycle` are arguments of the converter when the file does not state
  them (Abaqus reports carry no units). What the file does state (step, frame, increment, step time,
  "Unique Nodal" vs "Element Nodal", averaging criteria, "Not averaged across region boundaries",
  the source ODB name) is parsed, recorded in the CSV header (`notes`, `averaging`, `source`) and
  checked against the arguments.
* The CSV key is `(instance, node_id, material, step)`. An **Element Nodal** report has one row per
  element-node: collapse to one row per node only if the values of a node agree across elements
  within a tight tolerance (that is what the averaged "Unique Nodal" output looks like); otherwise
  refuse with advice to report Unique Nodal values or to give the element label (do not average
  silently). A report that mixes materials at an interface node must not be averaged globally
  (`averaging` values exist for that).
* Components: Abaqus `S.S11 S.S22 S.S33 S.S12 S.S13 S.S23` (note the report order 12, 13, 23);
  strain `E.E11...`, plastic `PE.`, elastic `EE.`; Ansys `SX SY SZ SXY SYZ SXZ`. The alias table
  and the sentinel tests in `tests/test_fecsv_header.py` and `tests/test_fecsv_import.py` already
  fix the component mapping; reuse `COLUMN_ALIASES`, do not duplicate it. Distinct values per
  component in every test (a swapped `S13`/`S23` must fail).
* Multi-frame and multi-step: one report may hold several frames, or there may be one report per
  frame; accept both; steps ascend; `time` from "Step Time" or "Total Time" as the file says and
  declared in `time_unit`.
* Third-party data is untrusted. Fetch each file into its own new empty directory (scratchpad), never
  run anything from it, read it as text, run Python that reads it with `-I`, check the repository
  LICENSE before using a file. Commit a data file as a fixture only if its licence allows it and it is
  small (a few kB), with the licence text and source URL beside it; otherwise make the test fetch-free
  by committing a tiny hand-written file in the same format.

## Public data already located (verify licence and content yourself)

* `renato-yuzup/axis-fem` (MIT): `Model Input Files/hexaedro_validacao/problema_1/abaqus/` has
  `problema_1_strain_stress.rpt` (4.6 kB, one frame, 8 element-nodes of one C3D8, columns
  `Element Label, Node Label, E.E11..E.E23, S.S11..S.S23`, "Averaged at nodes"; the expected physics
  is uniform axial stress `S33 = 50e3`, all other components about zero, which is a check), plus
  `problema_1.inp` (node coordinates), `.dat`, `.fil`, `.odb`. Other cases are in sibling folders.
* `icirauqui/cppfea` `data/c3d6_2/abaqus/abaqus.rpt` and `c3d6_3` (licence not checked).
* `JuliaFEM/Materials.jl` `examples/one_elem_disp_chaboche/unitelement_results.rpt` (52 kB, elastic-
  plastic unit element, licence not checked): useful for plastic strain `PE` and several frames.
* `lucas-bublitz/PHILLIPO.jl` `verification/*/*.rpt`; `Junyi88/SubroutineVault` `abaqus.rpt`.
* Ansys APDL text: `luchete80/MITC_python` `results/APDL/*/mitc.out` (27 kB, shell element, licence
  not checked). Look for a public `.rst` and for `PRNSOL` listings; if none is found, say so and build
  the parser from the documented listing format with a hand-written sample, and mark it untested on
  real output.
* `LallyLabTCD/localBasisAbaqus` (MIT) `.dat` files contain element tables for hyperelastic examples
  but no stress report; use only to test that `.dat` tables are recognised.

Use GitHub search (`gh` or the GitHub tools) for more: e.g. `"THE FOLLOWING TABLE IS PRINTED"`,
`"Field Output Report"`, `PRNSOL`, `extension:rst`.

## Design and structure

New module `src/multifat/fe/solvers/` (or `abaqus.py`, `ansys.py` in `fe/`), importing only
`multifat.fe`; no new import of `models`, `eval`, `imp`, `history`. Public functions such as
`abaqus_report_to_csv(report_paths, out_csv, *, stress_unit, length_unit, inp=None, ...)` and
`ansys_listing_to_csv(...)`, `ansys_rst_to_csv(...)`, plus a thin `multifat-fe-convert` console script
if cheap. Each returns the header written. The converters never write HDF5 directly: output is the CSV,
so that a user can inspect it and the one importer stays the single gate. Streaming where the file
allows (an `.rpt` is read frame by frame).

## Method

Use `superpowers:brainstorming` briefly (the scope above is already approved, so only the open design
questions), then `superpowers:writing-plans`, then execute natively test-first (failing test, minimal
code, pass), one commit per task, Conventional Commits, no attribution lines. Run the notation guard
(`tests/test_notation.py`) after every task: retired synonyms such as `mat`, `spec`, `dataset`,
`ratio`, `tolerance`, `young` fail it (use `material`, `tensor_spec`, `field`, `rank_ratio`, `tol`,
`E`). black line length 100, ruff correctness only, module docstring on line 1. Full `pytest -q` takes
16-19 minutes on this laptop: run targeted tests per task and the full suite once at the end (in the
background; do not edit files while it runs). Heavy runs belong on the workstations (see the global
CLAUDE.md); the regression snapshot (`tools/regression_snapshot.py run --skip-calibrate`, UniTN) is
needed only if an existing source file changes, which it should not.

Tests must include: the real `.rpt` files converted and re-imported (values against the file's own
numbers and the physical expectation), unit and convention arguments honoured, every refusal above
(non-unique element-nodal data, a missing unit, a frame that does not ascend, a report with no
stress, a changed column set between frames), multi-frame composition (build it by concatenating real
frames), and a mutation check (swap `S.S13`/`S.S23`; the test must fail). Deliver the documentation in
`docs/fe-csv-format.md` (section "Solver importers": what is read, what is refused, how the data was
tested, which files), decision-log entries `d093` onward, `AGENTS.md`, `README.md` and the checkpoint.

State plainly in the docs and the final report which parts were tested on real solver output and
which only on hand-written samples. Nothing here validates a fatigue criterion.

## Finish

Final gate: full pytest, `black --check .`, `ruff check .`, `tools/check_docs.py`,
`tools/check_scientific_reviews.py`, `tools/generate_method_docs.py --check`,
`tools/generate_workbook_docs.py --check`; `git diff origin/main --stat -- src/multifat ':!src/multifat/fe'
':!src/multifat/history'` must show nothing. Then a fresh whole-branch review by a separate reviewer
subagent (read the diff, not your summary), fix pass with failing tests first, and an update of
`docs/audit-checkpoint.md`. Commit locally and report; do not push.
