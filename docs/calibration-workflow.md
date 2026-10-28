# Calibration history and experimental run-outs

> **Scientific review proposal:** The planned source-prescribed defaults, optional
> all-eligible/subset calibration, and independent-validation safeguards are defined in
> [the Issue #14 protocol](calibration-validation-protocol.md). The procedures proposed
> there are not yet implemented in production.

> **In this repository.** This document tracks upstream PRs #5-#6 (commits dc4fe6d-f0df0a7, and the
> Findley amplitude section of PR #9, 36ff0a1, 0105b69) and
> describes the Script 0 / Script 1 workflow that produced the records in `Calibrations/`. That versioned
> store is not used here: `multifat-pipeline` calibrates missing constants in one run and writes them, with
> their provenance, to `<MATERIAL>_calibration.json`; the records are read as inputs with
> `multifat-pipeline --constants` (see [`Calibrations/README.md`](../Calibrations/README.md)). Upstream file
> paths (`fatigue/...`) and script options refer to upstream commit 71eaac3. The run-out conventions
> (failure-only fits, censored rows in plots, no errors for run-outs) hold here too.

Excel supplies the material and experimental inputs. Script 0 determines method
coefficients and master curves, stores them in JSON and preserves earlier records.
Equivalent reruns reuse the existing record instead of creating a duplicate.
Script 1 loads a compatible record, evaluates cases and exports results. Script 2
compares results and produces plots. No coefficient transfer through Excel is required.
This workflow retains the present FABER WG4.3 scope of constant-amplitude loading
and smooth/plain geometries. It does not introduce a new fatigue formulation or
censored-data regression.

## Calibration stages

The previous workflow had two distinct fitting stages. Script 0 optimised five
method coefficients on uniaxial tension and torsion failures, with the existing
load-ratio fallback rules. Script 1 separately fitted the selected damage–life
curve using reversed axial failures. Both stages are now initiated by Script 0.

Symbols (`DP`, `N`, `A`, `b`, the plane coefficients) are the ones of
[notation.md](notation.md); the table below names the stored quantities, not new ones.

| Quantity | Meaning and provenance | Storage |
|---|---|---|
| Elastic constants, yield/ultimate strengths, supplied fatigue strengths and strain-life properties | Material inputs, including existing optional-property estimates where explicitly reported | Excel |
| `K` | Existing ratio of fatigue strengths; not an optimisation output | Excel |
| `k_FS`, `k_FI`, `L_LI`, `k_MGSE`, `a_MSWT` | Dimensionless coefficients of the existing FS, Findley, Li, MGSE_ZHU and scalar MSWT implementations | JSON |
| `A`, `b`, `A_surv` | Master curve `DP = A*N**b` and its existing damage-space shift | JSON, separately for standard and `_ext` variants |
| `STD`, `STD_ln_DP`, `R2_log`, `n` | Fit diagnostics and number of usable failure observations | JSON |
| Plane grid, search bounds, refinement counts, selection tolerances and scatter settings | Numerical calibration settings, distinct from material properties | JSON metadata and versioned procedure |
| `alpha_w`, `xi_ZHU` | Retired-model constants transferred from the former Summary; original provenance/units unverified | JSON legacy metadata; unused by active methods |

Coefficient-optimisation diagnostic curves and the downstream master curves are
different objects. They retain their own specimen populations and scatter settings.
The optimisation default is `stdnum=2`, while the analysis master-curve default is
`stdnum=0`, as before. The legacy optimisation report calls Findley `FI`; the registry
and master curves use `FIN`. `_ext` is a separate stored curve, not a new material.

```bash
# Fit the requested coefficient and both FS curve variants.
python scripts/0_calibrate_material_params.py --material 42CrMo4_QT --models FS
python scripts/1_run_material.py --material 42CrMo4_QT --model FS_ext --results-dir Results/current --plot-dir Results/current/plots

# BP/CAIM require no fitted method coefficients, but still need master curves.
python scripts/0_calibrate_material_params.py --material 42CrMo4_QT --models BP CAIM
python scripts/1_run_material.py --material 42CrMo4_QT --model BP --results-dir Results/current --no-plots
python scripts/2_compare_results.py --results_dir Results/current --out_dir Results/Comparison
```

Omitting `--models` calibrates all active methods; omitting `--material` selects
all workbooks. Searches can be expensive. `--n1`, `--n2`, `--plane-stride`,
`--n-grid`, `--n-refine` and the two scatter settings are explicit CLI options.
Changing them is a numerical choice that must be reported. An optimisation without a finite objective, or an invalid required curve,
prevents a successful record from being written. No default coefficient is silently substituted.

> **In this repository.** There is no store record to withhold: `multifat-pipeline` calibrates in
> one run and writes `<MATERIAL>_calibration.json` whatever the outcome, so the sentence above
> becomes "no default coefficient is substituted silently". A block with no finite objective keeps
> running on the constant's default value and the run report says so; the section below defines the
> vocabulary that records it.

## Calibration status

> **In this repository.** `multifat-calibrate` prints, and `<MATERIAL>_calibration.json` stores as
> `status` (with a one-line `reason`), exactly one outcome per material. The same value is copied to
> the `calibration_status` column of the eligibility matrix and to
> `Comparison/calibration_status.csv`. Upstream's `status` key means `"success"` on a stored record;
> the two vocabularies are unrelated.

| Status | Meaning | `multifat-calibrate` |
|---|---|---|
| `calibrated` | Every requested block was fitted | exit code 0 |
| `not_applicable` | No plain failure rows to calibrate on: a fatigue-limit dataset takes its constants from the R = -1 limits, a notched-only dataset has none to calibrate | exit code 0 |
| `infeasible` | A block has no finite objective over its search range, or the workbook has no usable uniaxial subset to fit on. `reason` keeps the failure message; `infeasible` records the block, the parameter and a one-line diagnostic (for example `DP increases with life for all L_LI in [0, 5] (slope +0.09 to +0.41)`), computed only when the fit fails | exit code 1 |
| `error` | An unexpected exception, kept with its traceback in the log | exit code 1 |

Calibration stays all-or-nothing per material: an `infeasible` material yields exactly the constants
it yielded before the status existed, so `1_Ti6Al4V` keeps its Basquin and strain-life estimates and
loses the plane constants of every block, not only the one that failed. Models that read a lost
constant with a fallback still run, and the run report marks them `ok (default constant)` with a note
naming the block and the default value; models that require it stay `skipped` with the same pointer
appended to their reason. Block isolation (`partial`) is reserved for a later change and is not
produced.

Each material is independent: one `infeasible` or `error` material does not stop the batch, the
materials after it are calibrated as usual, and the process exits 1 so a script can detect the
failure. Per material, `multifat-pipeline --material NAME` still completes and records the outcome.
A workbook the pipeline rejects before calibration (`10_En3B`, `11_GJS600` here: no usable test rows)
has no sidecar and does not reach one of these four statuses; the eligibility matrix reports it as
`out_of_scope`.

## Record structure and units

Schema version 1 is checked by `validate_record` in
[`calibration_store.py`](https://github.com/andreachiocca/multiaxial-fatigue/blob/71eaac3/fatigue/eval/calibration_store.py). Files are written
through a temporary file and published atomically without replacing an existing ID.
Floating-point values retain Python/JSON round-trip precision. Required coefficients
and curve values must be finite. Undefined optional scatter or correlation diagnostics
are stored as `null`, with that convention recorded explicitly.

| Field | Content |
|---|---|
| `calibration_id`, `created_at`, `status`, `schema_version` | Unique ID, UTC timestamp, successful completion and schema identity |
| `material_id`, `material_name`, `material_condition` | Workbook-stem identity; unknown condition remains `null` |
| `source` | Workbook, Summary/Dataset sheets, workbook hash and original-archive link where applicable |
| `material_inputs`, `material_input_units` | Material inputs used by the model, separated from fitted coefficients |
| `parameters`, `parameter_units`, `parameter_origin` | Coefficients, units and whether they were optimised, explicitly reused or transferred |
| `curves` | Per-model standard/ext series, fits and references to selected observations |
| `design_observations` | Common specimen/source-row metadata stored once across methods and variants |
| `curve_units`, `damage_parameter_units` | Cycles, criterion-specific DP units and coefficient/exponent conventions |
| `optimisation` | Coefficient-search results, diagnostic fits, selected cases, fallback information and settings |
| `counts` | Available failures/run-outs/unknowns, invalid exclusions and actual observations used for fitting |
| `input_fingerprint`, `formulation_fingerprints` | Compatibility identities for material/failure inputs and numerical source |
| `amplitude_settings` | Optional method-scoped Findley algorithm, definition version, evaluation path and implementation fingerprint |
| `procedure_version`, `implementation_revision`, `references`, `software` | Procedure and scientific/software traceability, including Git state |

The following excerpt comes from the initial AISI316L migration record. It is a
schema illustration, not a complete loadable record or a new optimisation result.
The [complete record](../Calibrations/AISI316L/20261016T112312876379Z_cee3997362e2.json)
includes all required fields, all 15 active methods and independent provenance.

```json
{
  "schema_version": 1,
  "status": "success",
  "material_id": "AISI316L",
  "parameter_origin": "Transferred legacy workbook values; original fitting settings unverified",
  "parameters": {"k_FS": 0.62, "k_FI": 0.28, "L_LI": 2.13, "k_MGSE": 0.5, "a_MSWT": 0.36, "alpha_w": 0.0, "xi_ZHU": 0.0},
  "counts": {"failure": 31, "runout": 13, "unknown": 0, "invalid_excluded": 0}
}
```

Stress inputs use MPa; strains are dimensionless with the existing tensorial versus
engineering-shear conventions. FS, MKBM, LIU2021, GSA, BP and CAIM return dimensionless
repository parameters. The other active methods return MPa-valued parameters under
these conventions. `b` is dimensionless, `N` is in cycles and `A`/`A_surv` have units
`DP * cycle**(-b)`. These are implementation units, not evidence that a repository
adaptation reproduces every component of the original publication.

## Selection and compatibility

Script 1 selects the most recent successful compatible record by declared UTC time,
then calibration ID. It does not use modification times or run Script 0 automatically.
`--calibrations-dir` selects another archive; `--calibration-id ID` selects a specific
historical record and still enforces compatibility. Script 0 can explicitly reuse
coefficients with `--reuse-calibration ID` and records the source ID. Its requested
model and settings must be present in that record.

Compatibility requires the material identity, material properties, failed-specimen
histories, model, plane grid, curve scatter setting, procedure version, implementation
revision and numerical source fingerprints to agree. Source hashes normalise line
endings and path separators. They include the selected criterion and shared numerical
helpers. Maintainers must increment the procedure version if fitting-selection rules
change. A corrupt record produces an error instead of silently selecting another file.

The input fingerprint excludes run-outs, unknown-status observations, specimen/row
IDs, workbook formatting and file timestamps. It conservatively includes all failed
observations because coefficient selection can broaden its load-ratio tolerance.
Consequently, changing a non-fitting failure may also require recalibration. Reordering
rows or adding/removing run-outs alone does not invalidate an existing record. The full
workbook hash is provenance only and is not the compatibility key.

The initial four migration records preserve the previous coefficients without claiming
their original fitting settings are known. Fresh optimisation can produce different
coefficients; that is a new calibration and must remain a separate historical record.
The transferred `alpha_w` and `xi_ZHU` values do not restore their retired models.

## Equivalent reruns

`calibration_signature()` computes a deterministic identity from the material and
failure-input fingerprint, numerical options, formulation fingerprints, fitted
values, fitting provenance and NumPy version. Timestamps, Git working state,
workbook formatting, source paths, row identifiers and censored-observation counts
do not distinguish equivalent results. Sequential reruns with equivalent content
return the existing file and its original ID. Changed numerical results or settings
create a new historical record. Existing IDs remain protected against overwriting.
No extra database or index file is used. Concurrent independent writers can still
produce equivalent records, so calibration exports should be run sequentially.

Reusing a record preserves its original source metadata and counts. Script 0
separately reports the current observation counts and those in the stored record,
so adding run-outs does not change the historical provenance or hide their current
exclusion. Compare the human-readable JSON metadata directly when reviewing two
historical calibrations. The signature is computed on demand, allowing existing
schema-version-1 files to remain readable without rewriting them.

## Findley amplitude selection

> **In this repository.** The option is `multifat-calibrate --shear-amplitude-algorithm`; it sets the
> Findley shear amplitude (MCE or MCC) in the k_FI calibration through the plane shear rule of the
> established search and returns `amplitude_settings` with the result. No record is written, so the
> store-specific text below (Script 0/1, JSON signatures, compatibility table, fingerprints) describes
> upstream only; `multifat-pipeline --constants` rejects a record whose Findley amplitude is not MCE
> and does not expose the option. See [harmonic shear amplitudes](shear-amplitude.md#findley-calibration-selection).

Upstream's Findley calibration accepts the following Script 0 option. The mathematical
definitions and amplitude/range conventions are given in
[harmonic shear amplitudes](shear-amplitude.md). Selection affects both optimisation
of `k_FI` on the established axial/torsional failure subset and master-curve fitting
on the established reversed axial failure subset. The complete Findley expression
is maximised on every candidate plane with the same selected amplitude.

```bash
python scripts/0_calibrate_material_params.py --material Al7075_T6 --models FIN \
  --shear-amplitude-algorithm mcc_harmonic
python scripts/0_calibrate_material_params.py --material Al7075_T6 --models FIN \
  --shear-amplitude-algorithm mce_harmonic
```

`method_default` is the default and resolves to MCE, retaining the historical
computational path exactly. Explicit `mce_harmonic` and `rss` use the shared
harmonic implementation; `rss` is stored as `mce_harmonic`. Explicit overrides
for other methods are unsupported. In a mixed batch, those pairs are reported
and skipped while eligible Findley calibrations proceed. Default batches retain
ordinary calibration for all methods. Run-outs remain excluded from both stages.

Explicit MCC/MCE requests require original stored `S_sin` and `S_cos` tensors for
all selected failures. Missing or invalid histories identify the material, FIN,
case and coefficient. The extractor's synthetic in-phase endpoint fallback is
not accepted for these requests. Default evaluation retains that fallback. Missing
material properties and run-out-only datasets retain existing eligibility rules.

The following excerpt is the exact amplitude metadata from the coarse MCC example
below. This is an optional extension of schema 1; historical files are not rewritten.
The definition version comes from the numerical library, not a second constant.

```json
{
  "amplitude_settings": {
    "FIN": {
      "algorithm": "mcc_harmonic",
      "definition_version": "harmonic-shear-v1",
      "evaluation_path": "harmonic",
      "implementation_fingerprint": "61d9e0b2d7a9d39939ae76f96fda437ad3c3219816ee725410e5ae4badc73dff"
    }
  }
}
```

`evaluation_path` distinguishes the historical adapter (`legacy`, necessarily MCE)
from explicit Phase 2 evaluation (`harmonic`). Its implementation fingerprint is
also stored under `formulation_fingerprints.FIN`. The legacy path retains the
existing production source identity. Harmonic evaluation additionally hashes
the Findley calibration adapter, Phase 2 evaluator and shared amplitude library,
along with the canonical algorithm and definition version. Source line endings
are normalised as in existing formulation hashes. These added dependencies affect
Findley harmonic records only. Other methods retain their existing hashes.

| Requested Findley algorithm | Record | Behaviour, subject to all other checks |
|---|---|---|
| Default, MCE or RSS | Verified historical MCE without amplitude metadata | Accept |
| MCC | Historical MCE | Reject |
| MCC | Compatible new MCC | Accept |
| Default, MCE or RSS | Compatible new MCE | Accept |
| MCE | MCC, even if newer or explicitly requested by ID | Reject |
| MCC or MCE | Different definition or relevant implementation | Reject |

Missing metadata alone does not establish MCE. The store recognises the verified
historical FIN fingerprint `f65f1330ca96ba1fca84513f61ddaf59c12b7b1af8e57c802de746e5953acb7b`
and also requires it to match the current legacy source, with unchanged material,
inputs, procedure, implementation revision and curve settings. Unknown provenance
or a changed formulation is not accepted. Alternative records have distinct source
identities, so removing their metadata cannot make them verified legacy MCE.

The signature includes method-scoped amplitude metadata. Equivalent sequential
reruns reuse their original ID; MCC/MCE identities remain distinct even if their
numerical fits coincide. Explicit coefficient reuse through `--reuse-calibration`
also requires a matching amplitude definition. Optimisation metadata and downstream
curve metadata must agree. Historical files and IDs remain immutable.

Script 1, eligibility and batch/review assessments call the loader without an
override and therefore require MCE for FIN. Algorithm compatibility is checked
before creation time or ID. A newer MCC record cannot displace MCE; an explicit
MCC ID or an MCC-only archive gives a controlled compatibility error. Script 1
performs no automatic recalibration and offers no MCC prediction option. New MCC
damage-life curves require a separately approved Phase 3B assessment integration.

### Complete numerical example

On the current `Al7075_T6` workbook, the commands above were run with the following
additional settings and a temporary `--calibrations-dir`:

```text
--n1 3 --n2 6 --n-grid 5 --n-refine 1 --plane-stride 1
--parameter-stdnum 2 --stdnum 0
```

The coarse grid is a reproducibility check. It retains the existing `[0, 3]`
coefficient bounds and selects 23 coefficient failures (11 axial and 12 torsional),
11 axial master-curve failures and no run-outs from the 14 available run-outs.

| Configuration | Optimised k_FI | Master A = A_surv | Master b | Scatter STD |
|---|---:|---:|---:|---:|
| Default MCE | 0.75 | 724.2113871168674 | -0.10692522368903569 | 0.5061625852437134 |
| Explicit MCE | 0.75 | 724.2113871168674 | -0.10692522368903569 | 0.5061625852437134 |
| Explicit MCC | 0.75 | 724.2113871168674 | -0.10692522368903569 | 0.5061625852437134 |

The coefficient-stage fit uses its own 23-point population and gives `A=634.6482916208016`,
`b=-0.10977780510589384`, `A_surv=445.01621466900093` and `STD=1.6167214068576232`
for the default. It is distinct from the final axial master curve. The original
main implementation was executed on these same inputs/settings; its coefficient
result, selected observations, master curve and damage series match the new default
exactly. MCC/MCE damage series also coincide here because the calibration paths
are proportional. Their record identities remain distinct. These new optimisation
results do not replace the historically transferred workbook coefficient.

Non-proportional diagnostic cases in the automated tests compare the actual
coefficient-objective callable against independent SVD calculations on every plane.
The existing circular and ranking-reversal examples distinguish MCC from MCE.
These checks verify numerical consistency and calibration reproducibility, not
independent predictive accuracy. Explicit MCE can differ from historical endpoint
encoding at float64 rounding or very small shear scales, as documented in the
[numerical limitations](shear-amplitude.md#controlled-findley-evaluation).

Upstream's `tests/test_findley_calibration.py` covers the calibration/JSON/production-safety
integration checks; here `uv run pytest tests/test_findley_calibration.py` checks the equivalent
parts that exist (selection, rules against the harmonic evaluator, calibration, pipeline guard).

## Run-out interpretation and plots

Only failures contribute to coefficient searches, regressions and master curves.
Central selectors and a defensive check in the generic fitting collector enforce this
policy. Run-outs retain the original loading tensors, specimen identity, source row
and recorded cycle count. Unknown status is not inferred from a large cycle count.
See [experimental inputs](../Experimental_Data/README.md) for the workbook convention.

Result CSVs add `Test_status`, `Specimen_id`, `Calibration_id` and `Calibration_file`.
The existing `Nf_exp` column contains recorded cycles; its run-out values are lower
bounds. Their ordinary life/damage errors and `DP_fit` are undefined. Script 2 also
masks these fields after reconstructing errors and excludes censored rows from
statistics, boxplots, distributions and damage-fit comparisons. CSVs without status
retain their historical failure-only interpretation; use explicit status for new data.

Run-outs appear as open circles with right-pointing arrows in damage–life and life
parity plots. The arrow indicates that the unknown failure life exceeds the recorded
abscissa. PNG and TikZ comparison exports retain this convention. These are qualitative
observations, not ordinary prediction-error measurements. `--hide-runouts` changes only
visualisation. Existing suppression of life inversion for flat, positive or excessively
shallow curves remains unchanged. Such cases retain the damage–life plot; saved parity
plots are omitted. No censored-data regression or survival-analysis model is added.

## Model output and fatigue-life comparison

The computational family does not determine the output quantity. In particular,
the current non-critical-plane BP and CAIM implementations return dimensionless
damage parameters instead of fatigue lives. All 15 currently enabled base
methods retain their existing damage-parameter workflow.

| Model output | Pipeline behaviour | Comparison |
|---|---|---|
| `damage_parameter` (existing default) | Load the compatible JSON curve and retain the existing slope safeguard | Life parity when inversion is permitted; damage-space comparison otherwise |
| `life_cycles` | Use the predicted cycles directly, without fitting a master curve or applying its inverse | Life parity and failure-only life errors |
| `index` | Reject generic curve fitting and automatic life conversion | A separately documented relation is required before life comparison |

A cycle-valued case model declares `output_quantity = "life_cycles"` and returns
`CaseResult(values={"": N}, output_quantity="life_cycles")`, where `N` is finite,
positive and expressed in cycles, not reversals. Model and result declarations
must agree. Script 1 uses measured material inputs and does not require a damage
master-curve archive on this path. A future method requiring fitted coefficients
must first provide its own explicit, traceable calibration procedure. Those
coefficients must not be estimated by fitting predicted lives to observed lives.

CSV rows use `Metric_mode=Nf_direct`. `Nf_expected` contains the native prediction,
while `CP_value`, `DP_fit` and damage errors remain undefined. No specimens are
labelled as master-curve design points. Run-outs remain lower-bound observations
in the parity plot and contribute no ordinary error. Script 2 accepts these CSVs
alongside the existing results. This interface is tested with a synthetic model;
it does not introduce or scientifically validate a new fatigue criterion.

## Verification boundary

The test suite covers JSON precision/history, required fields, incompatible inputs,
failure-only optimisation, generic curve fitting, explicit status import, error masking,
plots and Script 0 → JSON → Script 1 integration. The pre-change optimisation fixture
uses three specified planes with separated grid minima for a bounded software regression, not a recommended scientific
plane search. The 64-run comparison reproduces all 3,936 original failed-specimen rows
with zero observed numerical difference and adds 672 run-out rows. Calibration and
software regression remain distinct from independent scientific validation.

## Revised collection

Script 0 now checks and processes each selected method independently, reporting
skipped and failed combinations and optionally writing `--status-file` JSON.
Script 3 exports method-material eligibility and can execute compatible pairs.
The [dataset workflow](dataset-workflow.md) describes missing-property checks,
source-identity restrictions and common failure-observation comparisons.
