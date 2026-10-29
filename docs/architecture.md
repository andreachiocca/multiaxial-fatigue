# Data flow and architecture

The production entry points are `multifat-calibrate`, `multifat-run`, `multifat-compare` and `multifat-pipeline` in `src/multifat/cli/`. The numbered files in `scripts/` are compatibility shims. The pipeline validates and imports one or more workbooks, supplies or calibrates eligible constants, runs registered models and writes a per-model report before invoking comparison (see [workbook format](workbook-format.md)). `multifat-run` can execute one selected model. A fitted damage-parameter (DP) curve is used only where the method's comparison route calls for one; a native-life method can predict cycles without that curve. Neither entry point edits the source workbooks.

| Path | Role |
| --- | --- |
| `src/multifat/imp/experimental_data.py` | Workbook loading, material properties, stress/strain tensors and harmonics |
| `src/multifat/imp/schema.py` | Workbook input contract: column rules, value parsing, `WorkbookError`, validation report |
| `src/multifat/cli/pipeline.py` | Pipeline run by a single command across every registered model, with a report for each model run |
| `src/multifat/cli/run.py` | Case evaluation, model-specific life handling and per-method CSV export |
| `src/multifat/cli/compare.py` | Capability-filtered life and DP statistics, figures and native-life counts |
| `src/multifat/planes/sampling.py` | Candidate rotations of planes; the third column of each rotation matrix gives that plane's normal |
| `src/multifat/models/base.py` | Interfaces for models and result types |
| `src/multifat/models/cp_adapter.py` | Wrapper at case level for critical-plane methods |
| `src/multifat/models/registry.py` | Names made visible to the runner |
| `src/multifat/models/comparison.py` | Per-method declarations for native or fitted life and the availability of a DP reference |
| `src/multifat/models/utils/tensors.py` | Shared conversion of 3×3 inputs and shape checks |
| `src/multifat/models/integral_methods/` | Integral BP methods and hybrid critical-amplitude/integral-mean CAIM criteria |
| `src/multifat/eval/` | Evaluation, calibration and fitting on critical planes |
| `src/multifat/eval/calibration.py` | Reversed-axial design series and a DP fit only where the declared comparison route uses one |
| `src/multifat/eval/endurance.py` | Fatigue-index error on fatigue-limit datasets |
| `src/multifat/report/export_csv.py` | Per-method output in CSV format |
| `src/multifat/report/subsets.py` | Observation identities, dataset fingerprint and common subsets (from upstream 904c619) |
| `src/multifat/models/utils/shear_amplitude.py` | Harmonic MCC/MCE shear amplitudes that the criteria do not use (from upstream 88229e9) |

Material properties are obtained by the importer from the `Summary` sheet, whereas cases come from the `Dataset` sheet; there is one workbook per material ID, the material ID being the file name with its ordering number removed. Plain failures and run-outs are kept separate. Notched rows are discarded, and in the revised format rows with unknown loading are excluded (see [dataset-workflow.md](dataset-workflow.md)). Strains remain undefined when a workbook lacks E or nu, in which case only stress-only models are run. If a workbook holds fatigue-limit states only, it is routed to the fatigue index error (`eval/endurance.py`) instead of a life curve. Optional radial and torsional phase columns generate stress and strain harmonics (`mean`, `sin`, `cos`) for single-frequency sinusoidal loading. The critical-plane evaluator relies on them for non-proportional loading while retaining the legacy two-loadstep interface (`S0`/`E0` and `S1`/`E1`). Direct methods are given the harmonics, but they are required to use them explicitly. The symbols are listed in [notation.md](notation.md).

## Method and comparison flow

The registry currently exposes 17 critical-plane methods through `CPModelAdapter` and 18 direct case-level models. Both enter the runner through `get_models()` and return `CaseResult`. Conventional methods use a fitted DP-versus-life curve when eligible. Four methods have a declared native-life route: `TVEIT2025`, `MOREL`, `MDC` and `MMCCM`. Tveit's virgin endurance index is diagnostic and has no valid DP comparison. Morel has a separate empirical DP fit; MDC and MMCCM have model-law DP references. The [comparison capability register](comparison-capabilities.md) states the exact interpretation and limitations of each route.

The importer supplies physical harmonic tensors (`S_mean`, `S_sin`, `S_cos`) to methods that need phase information. The evaluator obtains reversed-axial design observations where available. Native life is computed by a method's own law; Morel additionally fits its life coefficients on suitable design cases, while Tveit uses supplied coefficients and can predict without axial design cases. `multifat-compare --plot-kind` selects life, damage or both output families using declared method capabilities and the result mode. A method's registration is not evidence of independent experimental validation.

Main CSV fields include `Material`, `Method`, `Method_family`, `Loading_type`, `Flag`, `Metric_mode`, `Nf_exp`, `Nf_expected`, `CP_value`, `DP_fit`, `Error_ln` and `Error_ln_dp`. `Flag=0` identifies design observations, `Flag=1` other failed specimens and `Flag=2` censored run-outs. Conventional `Metric_mode=Nf` inverts an eligible fitted curve; `DP_DIFF` suppresses that inversion when the slope safeguard applies. For `Metric_mode=Nf_direct`, `Nf_expected` is the native prediction also recorded in `Nf_native`, not a re-fitted Basquin estimate. `DP_fit` is defined only when a separate DP comparison is supported. Run-outs contribute no ordinary uncensored error. These field names remain stable for existing readers.

Exports also contain `Error_log10_dp`, `Implementation_revision` and `Method_reference`, and end with the observation identity columns `Test_status`, `Specimen_id`, `Case_id`, `Observation_id` and `Dataset_fingerprint`. Where a DP reference exists, the error is `log10(DP_fit/CP_value)`. The native-life summary separates design observations, assessment failures and run-outs, including predictions of infinite life. Neither these counts nor fitted design residuals are independent validation statistics. The harmonic conventions, corrected methods and excluded historical results are covered in [method-audit.md](method-audit.md).

## Development conventions and possible interoperability with FatPy

The [FatPy contributing](https://faberorg.github.io/FatPy/development/contributing/), [testing](https://faberorg.github.io/FatPy/development/testing/) and [code-style](https://faberorg.github.io/FatPy/development/code_style/) guides were reviewed by upstream commit 0105b69 (PR #9); the new calibration uses typed interfaces and one function per step. Numerical theory and usage remain in their technical documents. Here, the pytest suite, the `multifat` src layout and the black/ruff checks are retained; FatPy imports, packaging changes and repository-wide typing migrations are outside this change.

Plane-based methods and structural-mechanics utilities are documented as separate elements in the current [FatPy repository](https://github.com/faberorg/FatPy). The independent harmonic amplitude function (`multifat.models.utils.shear_amplitude`), tensor transformations and the plane/method evaluation interfaces are potential components for a future WG6 exchange. The similarity is architectural; API compatibility has not been confirmed. An exchange would require agreement on tensor/component order, tangent bases, harmonic phase conventions, units, amplitude versus range, output quantities and calibration provenance. No FatPy source or dependency is introduced.
