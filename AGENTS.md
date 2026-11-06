# Working on this repository

`multifat`: multiaxial high-cycle fatigue criteria (critical-plane, invariant, integral and
energy-based) evaluated on experimental datasets. Python 3.10+, `src/` layout, `uv` for the
environment. MIT licence (see `LICENSE`).

## Working protocol

- Open each session by reading `docs/audit-checkpoint.md`, checking Git status, recent commits and
  any open PR, and confirming that the checkpoint matches the source. Continue from its first
  unfinished task. Do not repeat completed audits or keep a second progress file.
- At each milestone, run the focused checks, update the checkpoint with the result, the modified
  files, the remaining work, blockers and the exact next action, then commit verified progress
  (push when authorised). Record untested changes explicitly before an interruption.
- Keep software regression and scientific validation separate. A passing regression snapshot shows
  that numbers did not change; it does not validate a method against experiments.
- Preserve numerics: any change that could move a number needs the regression snapshot before the
  commit. Preserve uncommitted work, use the existing branch and PR where appropriate, and do not
  force-push or merge without the user's instruction.

## Path map (from the upstream `fatigue/` layout)

`fatigue/` is now `src/multifat/`: `fatigue/models` -> `src/multifat/models`, `fatigue/eval` ->
`src/multifat/eval`, `fatigue/imp` -> `src/multifat/imp`, `fatigue/report` -> `src/multifat/report`,
`fatigue/planes` -> `src/multifat/planes`. The three numbered scripts (`0_calibrate...`,
`1_run_material.py`, `2_compare_results.py`) are console scripts (`multifat-calibrate`, `-run`,
`-compare`, `-pipeline`); `scripts/` keeps thin shims.

## Layout

- `src/multifat/models/{cp_methods,invariant_methods,integral_methods,energy_based_methods}`:
  one module per published criterion; `registry.py` wires the case-level models, `references.py`
  holds the citations, `base.py` the interfaces.
- `src/multifat/fe`: local FE histories (`model.py`, `tensors.py`, `store.py` HDF5 file version 1.1,
  `csvio.py` CSV header, columns and `import_fe_csv`, `csvwriter.py` stdlib-only writer, `vtk.py`
  importer, `ranges.py`, `bridge.py`, `synthetic.py` tube oracle and CSV generator); imported by
  nothing else; specs in `docs/fe-input-format.md` and `docs/fe-csv-format.md`; `h5py` is the optional
  extra `fe`.
- `src/multifat/fe/solvers`: strict Abaqus RPT and Ansys PRNSOL/NLIST/complete Workbench table
  converters; optional native RST uses `ansys-mapdl-reader` (extra `ansys`, inspected 0.56 API).
  Output is atomic FE CSV, never direct HDF5; aliases come from `fe.csvio`. Preserve explicit units,
  global frame, averaging and state metadata. Native RST is stress-only, default solids, one material,
  real static/transient. No ODB reader. MIT real Abaqus fixtures: `tests/data/solvers/axis-fem`.
  Ansys text samples are hand-written. Real RST test is opt-in via `MULTIFAT_TEST_RST`; ordinary
  tests never download or execute third-party content.
- `src/multifat/history`: general load histories (`signal.py` contract, `channels.py`, `proportionality.py`,
  `rainflow.py`, `events.py`, `accumulate.py` Basquin and Miner, `eligibility.py`, `assess.py`,
  `nodes.py` node histories, unit cases and `assess_nodes`); imported by nothing else; spec in
  `docs/load-histories.md`.
- `src/multifat/eval`: plane preparation and search (`evaluator.py`, `plane_search.py`), life and
  fatigue-index fitting (`calibration.py`, `fitting.py`, `material_fits.py`, `endurance.py`),
  Findley amplitude options (`findley_calibration.py`).
- `src/multifat/imp`: workbook importer (schema, tensors, harmonics) and the declarative overlays
  (`overlays.py`). `planes/`: plane frames. `report/`: CSV and figure output. `cli/`: the four
  console entry points.
- `tests/`: `methods/` reproduce published tables from `tests/data/papers/`; `test_notation.py`
  guards the vocabulary. `tools/` holds the doc generators and `regression_snapshot.py`.
- `docs/scientific-reviews/`: scientific dossiers (LaTeX, PDF, BibTeX, benchmarks), checked by
  `tools/check_scientific_reviews.py`.

## Conventions

- **Notation is data.** `docs/notation.md` holds the single symbol table (paper symbol, Python name,
  workbook key, unit). A function implementing a published equation uses that equation's symbols;
  `tests/test_notation.py` enforces the retired-synonym table and requires every public parameter
  name to appear in the table, or in `tests/data/notation_allow.txt` with its reason. Add the row
  before introducing a name.
- **Docstrings**: one-sentence summary; module docstring on the first line of the file (before
  `from __future__`, no leading newline) with the model, citation, equations in a literal block and
  a `Repository choices:` paragraph; numpydoc section names (`Parameters`, `Returns`, `Notes`,
  `References`, never `Remarks`). Details in `docs/adding-methods.md`.
- **Spelling**: `Boehme` in Python source, its docstrings and the tables generated from it;
  hand-written documentation keeps the author's `Böhme` (method catalogue and what is generated
  from it, `docs/bibliography.md`, the README). `tests/test_notation.py` enforces both directions
  and rejects any departures heading other than `Repository choices:`.
- **Adding a method**: `docs/adding-methods.md` (declared parameters, paper-reproduction test,
  catalogue update, `docs/method-audit.md` for departures from the source).
- No format churn: black owns layout (line length 100); ruff checks correctness only.

## Commands

```bash
uv run pytest -q                                   # full suite (about 4 minutes)
uv run pytest tests/methods -q                     # paper reproductions only
uv run black --check . && uv run ruff check .
uv run python tools/generate_workbook_docs.py      # after optional_params changes
uv run python tools/generate_method_docs.py        # after catalogue changes
uv run python tools/check_docs.py                  # local documentation links
uv run python tools/check_scientific_reviews.py    # scientific dossiers
uv run python tools/regression_snapshot.py         # numerical gate before a commit
```

`regression_snapshot.py` runs every model on every workbook plus the calibration and the
comparison, and hashes the outputs. Run it on the affected models when a change could touch
numbers, and in full before a release. Never commit a numerical change without it.

## Outputs and data

Workbooks live in `Experimental_Data/`; `multifat-run` writes to `Results/current/`, and the
pre-audit CSVs are archived in `Results/archive/pre-audit/`. `docs/workbook-format.md` is generated,
never edited by hand. Finding IDs (F, U, P, N) in comments refer to the audit ledgers; confirmed
paper misprints stay as strict `xfail`s with their reason.

Supplied workbooks are never edited in place. Their corrections are declared in
`src/multifat/imp/overlay_files/`, matched to a workbook by content hash and applied in memory when
the file is read (`--no-overlay` returns the supplied data). Add or change an overlay operation
instead of editing the workbook.
