# User Guide

This walkthrough goes from a fresh checkout to a comparison of every criterion on one material.

---

## 1. Core Concepts

Three broad computational families of multiaxial high-cycle fatigue criteria are handled by `multifat`:

- **Critical-plane criteria** (e.g. Findley, Fatemi–Socie, Matake, McDiarmid, MWCM, MMCCM): the search targets the physical plane on which the combined damage from shear and normal stress/strain reaches its maximum.
- **Stress-invariant criteria** (e.g. Crossland, Sines, Dang Van, MMMH, Vu–Halm–Nadot): failure limits are expressed through the hydrostatic stress $\sigma_h$ together with deviatoric stress invariants ($J_2$, or the octahedral shear stress amplitude $\tau_{\mathrm{oct},a}$).
- **Integral and hybrid criteria** (e.g. Böhme–Papuga, CAIM, Morel, Robert–Fogue–Bahuaud, Papadopoulos): shear and normal stress amplitudes are integrated over a spatial reference sphere, or critical-plane amplitudes are combined with integral mean-stress terms.

Every symbol used below keeps one name across code and documentation: [notation.md](notation.md) gives the paper symbol, the Python name and the workbook key of each quantity.

### Damage Parameter ($DP$) vs. Fatigue Life ($N_f$)
The [`FatigueModel`](../src/multifat/models/base.py) interface is implemented by every criterion. A scalar **Damage Parameter ($DP$)** is what an evaluation of a load case returns:
- **Stress-based criteria**: an equivalent stress amplitude in $\mathrm{MPa}$ (e.g. Findley, Crossland).
- **Strain- and energy-based criteria**: an equivalent strain amplitude, or a strain energy density, in $\mathrm{MJ/m^3}$ (e.g. SWT, Ellyin–Kujawski).
- **Endurance criteria**: a dimensionless fatigue utilization index (e.g. Böhme–Papuga, CAIM).

Converting $DP$ into a predicted fatigue life $N_f$ relies on the calibrated uniaxial Basquin relation of the material ($DP = A \cdot N_f^b$). Where a method carries its own explicit life formulation (MDC, MMCCM, MOREL), the native predicted cycles ($N_{f,\mathrm{native}}$) are additionally obtained by direct computation.

---

## 2. Installation and Verification

### Recommended: Astral `uv`
`multifat` uses [`uv`](https://docs.astral.sh/uv/), which locks the dependencies and runs every command:

```bash
# Clone the repository
git clone https://github.com/andreachiocca/multiaxial-fatigue.git
cd multiaxial-fatigue

# Install dependencies and sync virtual environment
uv sync

# Run the test suite (about 1,700 tests, including paper-reproduction checks)
uv run pytest
```

### Alternative: Standard `pip`
`multifat` also installs into an existing Python 3.10+ environment:

```bash
python -m pip install -e .
pytest
```

---

## 3. Primary Workflow: The One-Command Pipeline

`multifat-pipeline` is the recommended route: it runs the whole analysis lifecycle in one deterministic pass.

### Running a Single Material
Every registered model on one material workbook:

```bash
uv run multifat-pipeline --material 42CrMo4_QT --jobs 4
```

### Running All Materials in Batch
To process the complete experimental library:

```bash
uv run multifat-pipeline --all --jobs 4
```

### What the Pipeline Does Automatically

```text
  Experimental Workbook (.xlsx)
                │
                ▼
  ┌───────────────────────────┐     Outputs: <MATERIAL>_validation.txt
  │   1. Schema Validation    │ ──► (schema integrity & stress-ratio audit)
  └─────────────┬─────────────┘
                ▼
  ┌───────────────────────────┐     Outputs: <MATERIAL>_calibration.json
  │ 2. Parameter Calibration  │ ──► (fitted Basquin curves & sidecar provenance)
  └─────────────┬─────────────┘
                ▼
  ┌───────────────────────────┐     Outputs: <METHOD>_<MATERIAL>.csv
  │  3. Multi-Model Execution │ ──► (DP & Nf across all 34 registered criteria)
  └─────────────┬─────────────┘
                ▼
  ┌───────────────────────────┐     Outputs: Comparison/ & run_report.tsv
  │ 4. Comparative Assessment │ ──► (scatter plots, LaTeX tables, eligibility)
  └───────────────────────────┘
```

1. **Validation**: the workbook is checked against the schema contract ([`docs/workbook-format.md`](workbook-format.md)). Execution stops at once when strict errors are found; warnings, by contrast, end up in `<MATERIAL>_validation.txt`.
2. **Calibration**: which constants the models require but the workbook does not provide is identified first. Uniaxial Basquin curves are then fitted on $R=-1$ tests, the interaction coefficients that are needed get calculated, and all of this is stored in `<MATERIAL>_calibration.json`, together with a cryptographic fingerprint of the dataset. **No modification is ever applied to the original Excel file.**
3. **Execution**: every one of the 34 registered criteria (along with their `_ext` variants) is evaluated over all failure and run-out observations. Should a model be missing required parameters, it is marked as skipped in `<MATERIAL>_run_report.tsv`, and the remaining models keep running.
4. **Comparison**: experimental lives and predicted lives are set against each other, log-error statistics ($\mu_{\log_{10}}, \sigma_{\log_{10}}$) are computed, and comparison plots plus summary tables are written out. On `--all` runs, the global eligibility matrix `Comparison/eligibility_matrix.csv` is produced as well.

Every output goes to `Results/current/`.

---

## 4. Output Artifacts

The pipeline writes the files below into `Results/current/` for every material it runs:

| Artifact | File | Description |
| :--- | :--- | :--- |
| **Validation Report** | `<MATERIAL>_validation.txt` | Detailed account of the column checks, of stress-ratio consistency, and of parsing notes. |
| **Calibration Sidecar** | `<MATERIAL>_calibration.json` | Calibrated constants, Basquin slopes, the calibration `status` and its reason, provenance metadata, and the dataset hash. |
| **Execution Summary** | `<MATERIAL>_run_report.tsv` | Tab-delimited status per model: execution status, reasons for skipping, row counts. |
| **Model CSVs** | `<METHOD>_<MATERIAL>.csv` | Complete row-by-row predictions: $DP$, predicted life $N_f$, log error, critical-plane angles. |
| **Comparison Plots** | `Comparison/<MATERIAL>/<MATERIAL>_expected_vs_experimental.png` | Scatter plots carrying Factor-of-2 ($F2$) and Factor-of-3 ($F3$) bands, next to the box plots, the error distributions and the damage-parameter parities of the same material. |
| **LaTeX Table** | `Comparison/summary_error_stats_table.tex` | Error statistics of every method and material in the run, with the style map, ready for direct inclusion in a publication; the same table is next to it as `summary_error_stats.csv`. |

---

## 5. Focused Single-Step Commands

Each pipeline stage is also a console command, for targeted experiments and development work:

### 1. Model Runner (`multifat-run`)
To run one chosen criterion on a single material, without invoking the full pipeline:

```bash
# Run the Böhme-Papuga criterion on 42CrMo4_QT
uv run multifat-run --model BP --material 42CrMo4_QT

# Run a stress-invariant criterion without generating individual plots
uv run multifat-run --model SINES --material 42CrMo4_QT --no-plots

# Run Findley with explicit shear amplitude algorithm
uv run multifat-run --model FIN --material 42CrMo4_QT --shear-amplitude-algorithm mce_harmonic

# Master curve with a fatigue limit (Kohout & Věchet 2001) instead of Basquin (also on multifat-pipeline)
uv run multifat-run --model SINES --material 42CrMo4_QT --sn-model kohout_vechet

# Design curve with a 90 % survival probability instead of the median (also on multifat-pipeline)
uv run multifat-run --model SINES --material 42CrMo4_QT --survival-probability 0.9
```

`--survival-probability` shifts the fitted master curve by $\Phi^{-1}(P_s)$ standard deviations of
its residuals, assuming normal scatter; it is a survival probability, not a confidence level, and it
does not apply to native life laws. See the [method audit](method-audit.md#damage-space-error-for-every-slope).

### 2. Parameter Calibrator (`multifat-calibrate`)
To compute the material constants that the models need, without evaluating them:

```bash
# Calibrate every workbook in Experimental_Data and print the constants
uv run multifat-calibrate --data-folder Experimental_Data

# Choose the shear-amplitude algorithm of the Findley fit
uv run multifat-calibrate --data-folder Experimental_Data --shear-amplitude-algorithm mce_harmonic
```

The batch calibrates every material independently: a material whose plane-constant calibration has no
finite objective (`1_Ti6Al4V` in this collection: its torsion/tension mix orders the Li damage
parameter against life) is reported as `infeasible` with the block and a measured diagnostic, the
materials after it are still calibrated, and the command exits with status 1. A summary table at the
end lists the status of every material.

### 3. Comparison Tool (`multifat-compare`)
To compare model outputs computed earlier, or to restrict a comparison to strictly matched failure subsets:

```bash
# Generate global comparative scatter plots and error tables
uv run multifat-compare --results_dir Results/current --out_dir Results/Comparison

# Compare models only on the exact subset of observations that both assess
uv run multifat-compare --results_dir Results/current --out_dir Results/Comparison/Common \
    --common-methods BP SWT --comparison-metric life --flag_filter other
```

### Direct Script Invocation
Standalone execution shims also exist in `scripts/` matching upstream naming (`0_calibrate_material_params.py`, `1_run_material.py`, `2_compare_results.py`, and `3_pipeline.py`). See [`scripts/README.md`](../scripts/README.md).

---

## 6. Verification & Quality Gates

Before contributing code, running sweeps or committing changes, verify the implementation with the test suites, style checks, documentation validators, and numerical regression gates:

```bash
# Run paper-reproduction unit tests only
uv run pytest tests/methods -q

# Run full test suite (about 1,700 tests)
uv run pytest -q

# Code formatting & lint checks (Black owns layout, Ruff checks correctness)
uv run black --check . && uv run ruff check .

# Check local documentation links & scientific reviews (LaTeX/dossiers)
uv run python tools/check_docs.py
uv run python tools/check_scientific_reviews.py

# Numerical regression snapshot (gate before any commit that could touch numbers)
uv run python tools/regression_snapshot.py
```

---

## 7. Adding Your Own Experimental Dataset

Testing your own fatigue data with `multifat` requires the following steps:

1. **Place the workbook** in `Experimental_Data/` as `<MATERIAL_NAME>.xlsx` (e.g. `Experimental_Data/Inconel718.xlsx`); the material name is the file stem without its ordering prefix (`6_42CrMo4_QT.xlsx` is `42CrMo4_QT`).
2. **Fill the two sheets the importer reads**:
   - `Summary`: one header row; the material values are read from the row labelled `Average`. Young's modulus and Poisson's ratio are required; the fatigue limits, the two $R = 0$ strengths, the yield and ultimate strengths and the method constants (`k FS`, `k FI`, `k MGSE`, `a MSWT`, `L LI`, ...) are optional, and a model whose constant is missing is skipped.
   - `Dataset`: one row per test, with columns found by header tokens: the radial, axial and shear mean and amplitude stresses, the `Phase AXIAL load` column of the revised format, the cycles to failure, `Notched` and `Run out (Yes/No)`.
   A life calibration needs at least two uniaxial $R = -1$ failures; a workbook of run-outs only is assessed through its fatigue limits instead.
3. **Run the pipeline**:
   ```bash
   uv run multifat-pipeline --material Inconel718
   ```

[`docs/workbook-format.md`](workbook-format.md) documents the full set of workbook column names, the header aliases, and the tolerance rules.

---

## 8. Python API Recipes

`multifat` imports from a script or a notebook:

### Recipe 1: Evaluate a Specific Criterion on Stress Tensors

```python
import numpy as np
from multifat.models.registry import get_models

# Retrieve the Böhme-Papuga model
model = get_models(["BP"])[0]

# Define 3x3 harmonic stress tensors (in MPa)
zero = np.zeros((3, 3))
axial_amplitude = np.diag([250.0, 0.0, 0.0])

# Evaluate case
result = model.evaluate_case(
    S0=axial_amplitude,
    S1=-axial_amplitude,
    E0=zero,
    E1=zero,
    params={"Sigm1": 240.0, "Taum1": 145.0, "Sig0": 400.0, "Tau0": 265.0},
    harmonics={"S_mean": zero, "S_sin": axial_amplitude, "S_cos": zero},
)

print(f"Calculated damage parameter (DP): {result.values['']:.4f}")
```

### Recipe 2: Build the Candidate Plane Set

A critical-plane criterion (e.g. Findley, Fatemi–Socie) is evaluated on a finite set of candidate planes. `build_R_list` returns the production set, as rotation matrices whose third column is the plane normal:

```python
from multifat.planes.sampling import build_R_list

# Candidate planes of the production grid
planes = build_R_list(n1=20, n2=40)
print(f"Sampled {planes.shape[0]} candidate critical planes.")
```

---

## 9. Frequently Asked Questions & Troubleshooting

### Why is a model reported as "skipped" in the run report?
A model is skipped when:
1. A parameter it needs is missing from the workbook and cannot be calibrated (a fatigue limit, a method coefficient, or the ultimate strength the estimate is built from).
2. Fewer than 5 uniaxial failure points are available, or they span less than one decade of life, so the Basquin fit is not identifiable.
3. The loading regime falls outside the mathematical scope of the method (e.g. non-zero mean stress applied to a criterion strictly formulated for zero mean).

### What do the calibration statuses mean?
`<MATERIAL>_calibration.json` carries one `status` per material, with a one-line `reason`; the batch calibrator prints the same table and the eligibility matrix repeats the value in its `calibration_status` column.

| Status | Meaning | `multifat-calibrate` |
| :--- | :--- | :--- |
| `calibrated` | Every requested constant was fitted | exit code 0 |
| `not_applicable` | No plain failures to calibrate on: a fatigue-limit dataset takes its constants from the R = -1 limits, a notched-only dataset has none to calibrate | exit code 0 |
| `infeasible` | A search found no finite objective in its bounds (or the workbook has no usable uniaxial subset); the block, the parameter and a one-line diagnostic are recorded | exit code 1 |
| `error` | An unexpected exception, kept with its traceback in the log | exit code 1 |

An `infeasible` material keeps the constants of every other step (for example its Basquin and strain-life estimates) and writes the same numbers as a run that ignores the failure; only the models that needed the lost constant change, and they say so: `ok (default constant)` in the run report when the model runs on the constant's default, a `skipped` reason pointing at the same block when it cannot run without it. Block isolation between fitted constants is not implemented, so `partial` is never written.

### How are run-out tests handled?
Specimens that reached run-out (`Run out = Yes` in the `Dataset` sheet) are **censored observations**:
- Every criterion evaluates them, so their damage parameter is computed.
- They are shown on the scatter plots with distinctive arrow markers, which indicate unfailed lower bounds.
- **Every parameter calibration, Basquin fit and error statistic excludes them**, so predictions are not skewed.

### Where can I find mathematical details for a specific criterion?
Start from the [method guides](method-guides/README.md): one page per criterion with the published equations, the repository implementation and its departures from the source, the calibration and a worked example. For every implemented criterion, the primary DOI and the implementation evidence are also cross-referenced in [`METHODS.md`](../METHODS.md) and [`docs/method-audit.md`](method-audit.md). Reproducible paper-validation test cases are located in `tests/methods/`.