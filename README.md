# multifat: a traceable library of multiaxial fatigue criteria

[![Quality checks](https://github.com/andreachiocca/multiaxial-fatigue/actions/workflows/quality.yml/badge.svg)](https://github.com/andreachiocca/multiaxial-fatigue/actions/workflows/quality.yml)
![Python 3.10-3.13](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

`multifat` implements published multiaxial high-cycle fatigue criteria (critical-plane, invariant,
integral and energy-based), checks each against its original paper, and runs all of them on the same
experimental datasets. One command on a workbook of fatigue tests returns every model's damage
parameter, predicted life and error statistics.

The project supports **FABER Working Group 4.3**
([WG4.3 activities](https://faber-cost.eu/wg4-fatigue-estimation-categories/wg4-3-multiaxial-high-cycle-fatigue/))
on its review article in multiaxial fatigue assessment, still in preparation. Scope:
constant-amplitude, mainly high-cycle fatigue of plain (smooth) specimens under proportional and
non-proportional loading, without notches, cracks or significant stress gradients.

## Quick start

> [!TIP]
> For an end-to-end walkthrough covering input data preparation, command-line workflows, and Python API usage, see the **[User Guide](docs/user-guide.md)**.

Python 3.10+ and [uv](https://docs.astral.sh/uv/), from the repository root:

```bash
uv sync
uv run multifat-pipeline --material AISI316L --jobs 4
```

`multifat-pipeline` takes a workbook path, `--material NAME` (the file name without its ordering
number, e.g. `AISI316L` for `8_AISI316L.xlsx`) or `--all`. For each material it

1. validates the workbook against the [workbook format](docs/workbook-format.md) and writes
   `<MATERIAL>_validation.txt`;
2. calibrates only the constants that the models need and the workbook lacks, writing them to
   `<MATERIAL>_calibration.json` with the run's provenance (workbook hash, failure-only
   `failure_fingerprint`, code revision, settings). The workbook is never modified; `--constants
   FILE` reads them from a `Calibrations/` record or an earlier sidecar instead;
3. runs every registered model and its `_ext` variant, one CSV per model. Run-outs stay censored
   observations (`Flag = 2`), outside every fit and statistic, and every row keeps its specimen
   and observation identity;
4. compares the results (statistics, plots and LaTeX tables) and, across workbooks, writes the
   method-by-material eligibility matrix.

A workbook of fatigue-limit states only (every plain row a run-out at the endurance basis, both
R = -1 limits in the Summary) is assessed by each model's fatigue index error instead of a life
curve (`FatigueLimit/`). `multifat-compare --common-methods A B` compares methods on the failures
they all assess. See [dataset workflow](docs/dataset-workflow.md).

Outputs go to `Results/current/`. A model that cannot run on a material (missing constants,
unidentifiable fit) is listed as skipped in `<MATERIAL>_run_report.tsv` with its reason. The others
still run.

Run single steps with `multifat-run --model BP --material AISI316L` (workbook constants only),
`multifat-calibrate` or `multifat-compare`. Without uv, `python -m pip install -e .` installs the
same console scripts. The numbered scripts in `scripts/` (from
`scripts/0_calibrate_material_params.py` to `scripts/3_pipeline.py`) are thin shims for them.

### Python API

```python
import numpy as np
from multifat.models.registry import get_models

model = get_models(["BP"])[0]
zero = np.zeros((3, 3))
axial_amplitude = np.diag([240.0, 0.0, 0.0])
result = model.evaluate_case(
    S0=axial_amplitude, S1=-axial_amplitude, E0=zero, E1=zero,
    params={"Sigm1": 240.0, "Taum1": 145.0, "Sig0": 400.0, "Tau0": 265.0},
    harmonics={"S_mean": zero, "S_sin": axial_amplitude, "S_cos": zero},
)
print(f"BP utilization: {result.values['']:.6f}")  # 1.000000
```

Stresses are in MPa. Reversed strengths are amplitudes and R = 0 strengths are maximum stresses.
`model.evaluate_case(...)` returns a `CaseResult` holding the scalar **Damage Parameter ($DP$)**:
an equivalent stress in MPa for stress models such as Findley and Crossland, an energy density in
MJ/m³ for SWT, a dimensionless fatigue index for endurance criteria such as BP and CAIM. Fatigue
life $N_f$ comes from mapping $DP$ through the material's calibrated uniaxial Basquin curve
($DP = A \cdot N_f^b$), or from the method's native life law where the source gives one (MDC,
MMCCM, MOREL). See [BP/CAIM usage](docs/bp-caim-usage.md) and
[architecture](docs/architecture.md).

## Visual results gallery

`examples/` holds a self-contained visualization suite. Every figure uses publication formatting
(Times New Roman, Computer Modern math, factor-of-2 scatter bounds) and is reproduced with one
command:

```bash
uv run python examples/generate_gallery.py
```

### 1. Critical-plane damage landscapes
Orientation-dependent damage distributions $\mathrm{DP}(\theta)$ under in-phase ($\delta = 0^\circ$)
versus non-proportional $90^\circ$ out-of-phase biaxial loading ($\sigma_a = 240\,\mathrm{MPa},
\tau_a = 150\,\mathrm{MPa}$; illustrative loading, not a material benchmark; $k_{FI} = 0.28$, $k_{FS} = 0.50$,
$S_y = 750\,\mathrm{MPa}$, $E = 206\,\mathrm{GPa}$, $\nu = 0.30$), showing how the phase shift moves
the critical plane and changes each criterion's damage peaks:

![Critical Plane Damage Landscapes](docs/images/critical_plane_polar.png)

### 2. Multiaxial fatigue life prediction correlation
Experimental versus predicted fatigue life ($N_{f,\mathrm{exp}}$ vs $N_{f,\mathrm{pred}}$) for the
plain axial and torsion tests of 42CrMo4-QT (`Experimental_Data/6_42CrMo4_QT.xlsx`; run-outs
excluded; no combined-loading tests in the workbook), against the Factor-of-2 ($F2$) and Factor-of-3
($F3$) scatter bands. Predictions come from the current pipeline
(`examples/build_predictions_csv.py`):

![Fatigue Life Prediction Correlation](docs/images/life_prediction_correlation.png)

### 3. Model benchmark & error dispersion
Cross-criterion comparison ranking models by log-error standard deviation
$s_{\ln} = \mathrm{std}(\ln(N_{f,\mathrm{exp}} / N_{f,\mathrm{pred}}))$ alongside systematic bias
$\mu_{\ln}$ (positive = conservative). Every test here is axial or torsional and all criteria are
calibrated on the same limits, so the ranking says nothing about multiaxial or non-proportional
performance:

![Comparative Benchmark of Multiaxial Fatigue Criteria](docs/images/model_benchmark_error.png)

### 4. Biaxial fatigue limit failure envelopes
Normalized failure boundaries ($\sigma_a / \sigma_{-1}$ vs $\tau_a / \tau_{-1}$) for fully reversed
proportional loading, calibrated on the 42CrMo4-QT fatigue limits ($\sigma_{-1} = 433$ MPa,
$\tau_{-1} = 344$ MPa): Gough-Pollard ellipse, maximum shear stress (Tresca, calibrated on
$\sigma_{-1}$), Crossland and Findley. Only the two calibration points are shown: the dataset has no
combined-loading fatigue-limit tests.

![Multiaxial Fatigue Limit Envelopes](docs/images/fatigue_limit_envelopes.png)

## Method coverage

<!-- METHOD_NAVIGATION_START -->

**81 missing method candidates** are available to explore: 69 have preliminary A/B feasibility (a new criterion or manageable extensions). Pick by scientific interest and available references/data; an off-roadmap proposal is welcome.

Browse contribution opportunities by family: [Critical plane](docs/METHOD_INDEX.md#family-critical-plane) · [Energy based](docs/METHOD_INDEX.md#family-energy-based) · [Integral and hybrid](docs/METHOD_INDEX.md#family-integral-and-hybrid) · [Invariant based](docs/METHOD_INDEX.md#family-invariant-based) · [Mesoscopic and damage](docs/METHOD_INDEX.md#family-mesoscopic-and-damage) · [Nonlocal and gradient](docs/METHOD_INDEX.md#family-nonlocal-and-gradient) · [Probabilistic](docs/METHOD_INDEX.md#family-probabilistic) · [Specialized](docs/METHOD_INDEX.md#family-specialized) · [Strain based](docs/METHOD_INDEX.md#family-strain-based) · [Stress based](docs/METHOD_INDEX.md#family-stress-based) · [Tensor and path based](docs/METHOD_INDEX.md#family-tensor-and-path-based).

**Available in the current API:** [FS](METHODS.md#i_fs), [FIN](METHODS.md#i_fin), [SWT](METHODS.md#i_swt), [MATAKE](METHODS.md#i_matake), [MCD](METHODS.md#i_mcd), [MGSE_YU](METHODS.md#i_mgse_yu), [MGSE_ZHU](METHODS.md#i_mgse_zhu), [MSWT](METHODS.md#i_mswt), [MKBM](METHODS.md#i_mkbm), [LIU2021](METHODS.md#i_liu2021), [LI](METHODS.md#i_li), [GSE](METHODS.md#i_gse), [GSA](METHODS.md#i_gsa), [BP](METHODS.md#i_bp), [CAIM](METHODS.md#i_caim), [DANGVAN](METHODS.md#i_dangvan), [OTT](METHODS.md#i_ott), [SWTD](METHODS.md#i_swtd), [CROSSLAND](METHODS.md#n01), [SINES](METHODS.md#n02), [PAPADOPOULOS](METHODS.md#n03), [MWCM](METHODS.md#n04), [VHN](METHODS.md#n20), [PBP](METHODS.md#n37), [TVEIT2025](METHODS.md#c02), [PCR](METHODS.md#g01), [PI](METHODS.md#g02), [MDC](METHODS.md#g06), [MMMH](METHODS.md#g07), [COIN](METHODS.md#g08), [VSE](METHODS.md#g15), [EK](METHODS.md#g16), [MMCCM](METHODS.md#g22), [MOREL](METHODS.md#g23), [RFB](METHODS.md#g41). Each link states the implemented scope and evidence; active does not automatically mean source completeness is established.

<!-- METHOD_NAVIGATION_END -->

Each registered method has a paper-reproduction test in `tests/methods/` and is checked against its
original publication (transcribed tables in `tests/data/papers/`, each citing DOI, table and page).
Departures from the source (shear-amplitude rule, elastic strains, repository adaptations) are
listed per method in the [method audit](docs/method-audit.md) and the
[complete catalogue](METHODS.md).

<!-- METHOD_PROGRESS_START -->

| Primary family | Fully implemented¹ | Partial | Missing | Requires review | Total |
| --- | --- | --- | --- | --- | --- |
| Critical plane | 6 | 1 | 11 | 0 | 18 |
| Energy based | 6 | 1 | 13 | 3 | 23 |
| Integral and hybrid | 5 | 0 | 4 | 0 | 9 |
| Invariant based | 5 | 0 | 8 | 0 | 13 |
| Mesoscopic and damage | 2 | 0 | 6 | 0 | 8 |
| Nonlocal and gradient | 0 | 0 | 7 | 0 | 7 |
| Probabilistic | 0 | 0 | 5 | 0 | 5 |
| Specialized | 0 | 0 | 3 | 0 | 3 |
| Strain based | 3 | 0 | 10 | 2 | 15 |
| Stress based | 1 | 0 | 12 | 0 | 13 |
| Tensor and path based | 1 | 0 | 2 | 0 | 3 |
| **TOTAL** | 29 | 2 | 81 | 5 | 117 |

<!-- METHOD_PROGRESS_END -->

**Complete** means the essential components are implemented and tested within the stated repository
scope, which may include documented adaptations. **Partial** means known essential gaps, and those
methods are disabled. **Requires review** means active code whose source completeness is not
established. A passing test is numerical verification, not independent validation against
experimental data.

<!-- METHOD_RECONCILIATION_START -->

The catalogue contains **81 missing method/scoped candidates**; the roadmap selects **6** for an initial, dependency-aware programme. The other **75 remain available for proposals**. The complete index retains **161 research records**, including **117 identifiable method/scoped records** and **44 bundles, related sources, unresolved overlaps/leads and wider workstreams** outside that progress denominator. These are not 161 independent fatigue theories.

<!-- METHOD_RECONCILIATION_END -->

The [count reconciliation](docs/COUNT_RECONCILIATION.md) accounts for every catalogue ID. The
[method index](docs/METHOD_INDEX.md) lists missing methods by family, with references and loading
scope. The [roadmap](ROADMAP.md) gives the suggested order. `docs/catalogue/methods.json` is the
source of truth.

## Add a method

1. Pick an ID in the [catalogue](METHODS.md) or the [roadmap](ROADMAP.md) (off-roadmap proposals are
   welcome) and read its original paper. Check the
   [external implementations](docs/EXTERNAL_IMPLEMENTATIONS.md) first, keeping their licences in
   mind.
2. Put a critical-plane criterion in `src/multifat/models/cp_methods/`, a direct one in
   `invariant_methods/`, `energy_based_methods/` or `integral_methods/`, and register it in
   `src/multifat/models/registry.py`. The interfaces are in
   [adding methods](docs/adding-methods.md).
3. Add a paper-reproduction test in `tests/methods/` and the citation in
   `src/multifat/models/references.py`, then update the catalogue and regenerate the docs.
4. Follow [CONTRIBUTING](CONTRIBUTING.md) for the checks and the evidence a PR needs.

## Add a workbook

1. Copy an existing workbook in `Experimental_Data/` and name it after the material (an ordering
   number such as `16_` may precede the name).
2. Fill the `Summary` sheet and one `Dataset` row per test. Only Young's modulus and Poisson's
   ratio are required. Without them only stress-only models run, and each other constant is needed
   only by the models that use it. The [workbook format](docs/workbook-format.md) lists the headers,
   units and required constants. The design set needs at least two uniaxial R = -1 tests with
   distinct lives.
3. Run `uv run multifat-pipeline --material <MATERIAL>` and read `<MATERIAL>_validation.txt` and
   `<MATERIAL>_run_report.tsv`.
4. Record where every value comes from (paper, table, page) in
   [Experimental_Data/PROVENANCE.md](Experimental_Data/PROVENANCE.md).

## Software map

| Component | Responsibility |
|---|---|
| `src/multifat/cli/` | `multifat-pipeline`, `-calibrate`, `-run`, `-compare` (shims `scripts/0_`-`scripts/3_*.py`) |
| `src/multifat/imp/`, `src/multifat/planes/` | Workbook contract and import (`imp/schema.py`); candidate planes |
| `src/multifat/models/` | Criteria, plane adapter, registry and source metadata |
| `src/multifat/eval/` | Plane search, calibration, material fits, life curves, fatigue index error |
| `src/multifat/report/` | CSV export, observation identities and common subsets |
| `tests/` | Paper reproductions, analytical and regression tests, optional material sweep |
| `tools/` | Doc generators, documentation check, regression snapshot, dataset install and audit |

A harmonic load cycle is represented by mean, sine and cosine 3×3 stress tensors. Elastic strains
follow from the workbook E and ν. Critical-plane methods search the plane with a grid plus
analytical seeds and local refinement. See [architecture](docs/architecture.md) for the structure
and [notation](docs/notation.md) for the symbols.

## Scope and maturity

Single-frequency harmonic cycles, elastic strains, damage-parameter assessment, a fitted life curve
per method and native life laws where the source gives one (MDC, MMCCM, MOREL). General periodic
paths, variable amplitude, plasticity, notches and gradients, and probabilistic methods are future
work, and the corresponding catalogue entries stay listed. Out-of-phase shear ranges follow each
source where it prescribes one (longest chord for FS, MKBM, LI, LIU2021). The other critical-plane
methods keep the RSS convention.

The CSVs in `Results/archive/pre-audit/` predate the numerical audit and are kept as history only
(see [Results/README.md](Results/README.md)). The experimental collection is the revised
October 2026 set (upstream 71eaac3, 15 workbooks, of which 7 fatigue-limit sets and 2 notched-only).
The provenance of each workbook, with its open questions, is in
[Experimental_Data/PROVENANCE.md](Experimental_Data/PROVENANCE.md).

## Development

```bash
uv run pytest
uv run ruff check .
uv run black --check .
uv run python tools/generate_method_docs.py --check
uv run python tools/generate_workbook_docs.py --check
uv run python tools/check_docs.py
```

CI runs these on Python 3.10-3.13 (plus Windows on 3.12). A change that can alter numerical outputs
also needs `uv run python tools/regression_snapshot.py run/compare` (see
[CONTRIBUTING](CONTRIBUTING.md)). The research provenance behind the catalogue (51-publication
review, assessment, saved bibliographic records) is frozen in
[docs/research](docs/research/README.md). Changes are listed in the [CHANGELOG](CHANGELOG.md), and project decisions in the [decision log](docs/decision-log.md).

## Licence and citation

Copyright (c) 2026 Andrea Chiocca, Michele Sgamma and Tommaso Grossi. `multifat` is released under
the [MIT License](LICENSE). No copyrighted article PDFs are distributed here. When referencing
calculations, cite the software as described in [CITATION.cff](CITATION.cff) with the exact commit,
together with each criterion's original publication from the catalogue.
