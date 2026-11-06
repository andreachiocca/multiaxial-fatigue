# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

`multifat` is an **open-source platform for multiaxial fatigue assessment**. It began as the
comparison tool behind the FABER WG4.3 review article and stays public afterwards, for anyone who
needs to assess components under multiaxial loading. Its users make engineering decisions with
its numbers, so three principles take precedence over speed or convenience.

### 1. Scientific correctness comes first

Every method must compute what its original publication defines, or say exactly where it
departs.

- **Work from the source.** Implement a criterion only from its original paper, never from
  secondary summaries. If the paper is not available, ask the maintainer for it.
- **Check the whole method, not just the equation.** A method is more than its damage
  parameter. Verify each of these against the source:
  - the critical-plane rule and tie handling;
  - the shear amplitude for non-proportional paths;
  - mean-stress treatment;
  - how its constants are calibrated, and on which tests: the design set;
  - its life law;
  - its loading scope.
- **Prove it with numbers from the paper.** Reproduce its tables or figures in
  `tests/methods/`. A mismatch is investigated until it is explained. A misprint in the paper,
  or a table that contradicts its own equation, becomes a strict `xfail` with the reason. Never
  loosen a tolerance to hide a gap.
- **Make every adaptation visible.** Each repository choice goes in the module docstring
  (`Repository choices:`), the method guide and `docs/method-audit.md`.
- **Never extend a method silently beyond its scope.** Out-of-scope loading raises
  `LoadingScopeError`: the model is reported as skipped, not run on an invented formula.
- **Check physical consistency, not only the code.** Units, dimensions, the sign of mean
  stresses, the behaviour at the R = -1 tension and torsion limits, and closed-form special
  cases.
- **Keep the two kinds of evidence apart.** A regression snapshot shows that numbers did not
  change. A paper reproduction shows the method is implemented as published. An error
  statistic on calibration tests measures the fit, not predictive accuracy.

### 2. Easy to extend

New methods, life laws, data sources and comparison features will keep arriving. Extend through
the existing extension points instead of special-casing:

- the plane-level and case-level model interfaces, the registry, and declared parameters
  (`required_params`, `optional_params`);
- capability attributes such as `np_shear_rule`, `design_loading`, `native_life` and
  `comparison_capabilities`;
- the typed skip errors;
- the overlays, for data corrections;
- the generated catalogue and docs.

When a new need does not fit, add a small, documented extension point that the other models
inherit by default, rather than a branch on a model name. Keep outputs backward compatible:
append CSV columns, and prove existing numbers unchanged with the regression snapshot.

The project is expected to merge into FatPy, the FABER (COST CA23109) Python package. Move
towards it without restructuring:

- Keep the existing names. The *FatPy nomenclature* table of `docs/notation.md` records the
  FatPy key of each multifat name.
- Build new names from the FatPy abbreviation list (`fat`, `lim`, `ten`, `tor`, `ben`, `coef`,
  `exp`, `amp`, `eq`, …) when no paper symbol applies.
- Keep MPa internally; FatPy v01 gives E in GPa.
- Open questions about the FatPy draft are in `docs/research/fatpy-nomenclature-v01-feedback.md`.

### 3. Document so a user can trust and reproduce the result

Each method needs:

- a guide with the original equations in LaTeX, the repository choices, the calibration, and a
  worked example that the tests execute;
- a catalogue record and an audit row;
- a CHANGELOG entry.

A result is only useful if its provenance is clear: the source, the constants and the tests
behind them.

## Ask before assuming

If a source, an equation, a convention or the intended behaviour is unclear, **stop and ask the
maintainer**. Interview them with focused questions rather than guessing. This covers, for
example:

- which formulation, version or plane rule a paper intends;
- whether a departure from the source is acceptable;
- how a new feature should behave;
- which papers or data are needed.

Briefly state what is known and what is ambiguous, propose options with their scientific
consequences, and wait for the decision. A wrong assumption in a fatigue criterion produces
plausible-looking but unsafe numbers.

@AGENTS.md

The working protocol, conventions, command list and path map are in `AGENTS.md` (imported above).
Contribution rules and the "fully implemented" checklist are in `CONTRIBUTING.md`, and the
method-adding steps in `docs/adding-methods.md`. This file adds the architecture that spans
several modules.

## Running tests

```bash
uv run pytest tests/methods/test_ninic.py -q                        # one paper-reproduction module
uv run pytest tests/methods/test_ninic.py::test_table1_eq4_k -q     # one test
uv run pytest tests/test_pipeline.py -k scope -q                    # by keyword
```

CI (`.github/workflows`) runs, in its lint job:

- `ruff check`;
- `black --check`;
- `tools/generate_method_docs.py --check` and `tools/generate_workbook_docs.py --check`;
- `tools/check_docs.py`;
- `tools/check_scientific_reviews.py`.

It runs `pytest` on Python 3.10–3.13 and on Windows. A stale generated file fails lint, so rerun
the generators after touching the catalogue or a model's declared parameters.

## Data flow (multifat-pipeline)

1. **Import** (`imp/`). A workbook in `Experimental_Data/` is read with
   `experimental_data.load_material_workbook`.
   - Overlays (`imp/overlay_files/`, matched by content hash) correct the data in memory.
   - Each case carries two load-step tensors plus the harmonic tensors `S_mean`, `S_sin`,
     `S_cos` (stresses) and `E_*` (elastic strains from `E` and `nu`).
   - Slots: radial is `sigma_xx`, axial is `sigma_yy`, torsion is `tau_xy`.
2. **Classify** (`imp/classification.py`). The dataset is one of `fatigue_life`,
   `fatigue_limit`, `ambiguous` or `out_of_scope`, and the result is written to
   `<MAT>_classification.json`.
   - Life datasets go through the master-curve path.
   - Limit datasets go through `eval/endurance.py`, which reports the fatigue indices `FI_axial`,
     `FI_torsion` and `FI_load`.
3. **Calibrate** (`cli/calibrate.py`, called by `cli/pipeline.py`). This fills missing plane
   constants (`k_FS`, `k_FI`, `L_LI`, `k_MGSE`, `a_MSWT`) by minimising the life scatter of the
   R = -1 tension and torsion failures.
   - The record is written to `<MAT>_calibration.json` and archived under `calibrations/`
     with an id.
   - Constants come from, in order: the workbook, then `--constants`, then this calibration.
4. **Evaluate** (`cli/run.py::evaluate_model`).
   - The design set is the R = -1 uniaxial tension failures (`extract_uniaxial_tension`).
     A model that declares `design_loading = "torsion"` uses `extract_pure_torsion` instead
     (MWHS).
   - The master curve `DP = A N^b` is fitted on that set (`eval/fitting.py`).
   - Every other case is then evaluated.
5. **Write** (`run.write_variant`). One CSV per model and variant.
   - `Flag`: 0 = design, 1 = assessed, 2 = run-out.
   - Error columns: `Error_ln`, `Error_log10_dp` and `Error_load`, the load-equivalent error
     from `eval/load_error.py`.
   - Life predictions are suppressed when the curve is too flat (k > 15); only DP errors remain.
6. **Compare** (`cli/compare.py`). Writes `Comparison/`: statistics, plots, LaTeX tables and the
   eligibility matrix.

## Model interfaces (`models/`)

### Plane-level critical-plane methods (`cp_methods/`)

- They implement `evaluate_on_planes` (batched) and `evaluate_on_plane`, and return damage plus
  a selection metric on each candidate plane. The local `z` axis is the plane normal; `xz` and
  `yz` are the shear components.
- `CPModelAdapter` wraps them into case-level models. The plane search is `eval/evaluator.py`
  plus `eval/plane_search.py`: a 522-plane grid, seed planes, refinement, and ties resolved
  towards the larger damage.
- The `""` variant is the damage on the max-metric plane; `_ext` is the maximum damage.
- With harmonics, `prepare_planes` rewrites the two load steps on each plane:
  - shear: the endpoints of the shear range chosen by the method's `np_shear_rule` (`"rss"` by
    default, or `"longest_chord"`), placed around the mean;
  - normal: mean ± harmonic amplitude.

  So `(S0r + S1r)/2` is the plane mean, and half their difference is the amplitude.
- The adapter forwards only these attributes from the plane method: `native_life`,
  `life_fit`, `expected_dp_at_life` and `design_loading`. Add new ones to its list in
  `cp_adapter.py`.

### Direct models (`invariant_methods/`, `integral_methods/`, `energy_based_methods/`)

They implement `evaluate_case` and return a `CaseResult`.

### Comparison capabilities

Native-life models (MOREL, MDC, MMCCM, TVEIT2025, LIUMAHADEVAN) declare `comparison_capabilities`
(`models/comparison.py`). All other models use the shared master curve.

### Errors that mean "skip this model", not "fail"

| Error | Raised when | Reported as |
| --- | --- | --- |
| `MissingParamsError` | constants are missing | skipped |
| `LoadingScopeError` | the load is outside the criterion's declared scope | skipped, eligibility `out_of_scope` |
| `DesignSetError` | the material has no design tests for the model | skipped |

In load-scaling searches (`endurance.scaled_damage`), a `LoadingScopeError` at a trial load
counts as +inf damage, so the search backtracks.

### Stress-only models

`registry.stress_only_models()` finds, by probing, which models run without E and nu.

## When adding or changing a registered method

Besides the steps in `docs/adding-methods.md`, these tests hold hand-maintained lists:

- the registered-model count in `tests/test_tveit_reporting.py`;
- the analytical plane-value table in `tests/test_method_audit.py` (critical-plane methods);
- the longest-chord set `CHORD` in `tests/test_np_shear_rule.py`.

Also:

- **Paper reproduction.** Implement only from the original paper. Paper-reproduction tables
  live in `tests/data/papers/*.csv`, with DOI, table and page in `#` header lines. A confirmed
  misprint, or a table that contradicts its own equation, becomes a strict `xfail` with the
  reason.
- **Notation.** A public function parameter name must be in `docs/notation.md` or
  `tests/data/notation_allow.txt`. `kappa` is reserved for f_-1/t_-1, and `ratio` is a retired
  synonym.
- **Catalogue.** For a new catalogue record, also add its row to
  `docs/research/Multiaxial_Fatigue_Master_Catalogue.csv` and its query to
  `docs/research/external-code-search.json`.
- **Generated files.** `METHODS.md`, `docs/methods/*`, `docs/catalogue/methods.csv` and
  `docs/workbook-format.md` are generated: regenerate them, never edit them by hand.
- **Importer schema.** Do not add an optional model constant to `imp/schema.py` just to read it
  from the Summary. A Summary column named like a declared `optional_params` key is already
  passed to the model. Changing the schema alters every dataset fingerprint
  (`tests/test_dataset_overlays.py`).
- **Homepage figures.** Refresh the two README comparison figures, as described at the end of
  the method steps in `docs/adding-methods.md`. Redraw only those two images:
  `generate_gallery.py` re-renders all six with different fonts.
- **Regression snapshot.** `tools/regression_snapshot.py compare` matches CSV columns by
  position, so append new output columns at the end. To show that existing outputs are
  unchanged, run `run <dir>` on `main` and on the branch, then `compare`.

## Collaboration constraints

### FE input and history adapters

`src/multifat/fe/solvers` converts Abaqus nodal RPT and Ansys PRNSOL/NLIST or complete Workbench
tensor tables to atomic FE CSV; native solid RST is optional (`ansys` extra, inspected reader 0.56).
Use the shared CSV aliases and require explicit units, global frame, averaging and state metadata.
Native RST is stress-only, one material, default solids, real static/transient. ODB and DAT are
unsupported. Real MIT Abaqus fixtures are in `tests/data/solvers/axis-fem`; Ansys text samples are
hand-written; the real RST test uses `MULTIFAT_TEST_RST` and never downloads during ordinary tests.
See `docs/fe-csv-format.md` and `docs/load-histories.md`. CSV/PVD conversion must preserve source
files and existing destinations on refusal. Strain-reading criteria require a shared scalar
stress/strain path for exact counting. Continue from `docs/audit-checkpoint.md`.

- Do not modify the Tveit model (`TVEIT2025`, `invariant_methods/tveit2025.py`): another
  contributor owns it.
- Do not force-push; pushed branches are updated with merges, not rebases.
- Use one topic branch and one PR per method or topic.
- Stacked PRs are merged in order, retargeting each to `main`.
