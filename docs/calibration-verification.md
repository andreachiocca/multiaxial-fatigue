# Calibration workflow verification

> **In this repository.** This document tracks upstream PRs #5-#6 (commits dc4fe6d-f0df0a7) and
> describes the Script 0 / Script 1 workflow that produced the records in `Calibrations/`. That versioned
> store is not used here: `multifat-pipeline` calibrates missing constants in one run and writes them, with
> their provenance, to `<MATERIAL>_calibration.json`; the records are read as inputs with
> `multifat-pipeline --constants` (see [`Calibrations/README.md`](../Calibrations/README.md)). Upstream file
> paths (`fatigue/...`) and script options refer to upstream commit 71eaac3. The run-out conventions
> (failure-only fits, censored rows in plots, no errors for run-outs) hold here too.

The comparison uses Git baseline `e70a8af3a7a06cfb49626bff16c04605c884830b`.
Its source and original workbooks were extracted into a separate temporary
directory and executed in the same environment as the revised implementation.
The initial unit baseline contained 20 passing tests. A bounded five-coefficient
optimisation baseline was also captured before changing the calibration code.

## Observed results

| Check | Result |
|---|---|
| Existing and new unit/integration tests | 39 passed |
| GitHub Actions | Push and PR checks passed on `f0df0a7`; the earlier tied-minimum fixture failure was diagnosed and corrected without a production-code change |
| Pre-change material/model sweep | 64 runs, 3,936 failure rows |
| JSON-based material/model sweep | 64 runs, 4,608 rows, including 672 run-out rows |
| Failure-row comparison | All original columns reproduced, maximum observed numerical difference 0 |
| Error-summary comparison | All summary statistics and sample counts exactly equal |
| Coefficient optimisation regression | All five coefficients and fit diagnostics reproduce the pre-change bounded baseline |
| Failed-only versus failed-plus-run-out data | Optimised coefficients and master curves agree exactly; fingerprints agree |
| Workbook preservation | All four original files archived byte-for-byte; all Dataset and retained Summary cells unchanged |
| JSON history | Equivalent exports reuse the existing file and ID; changed settings or results create distinct records; explicit IDs are retrievable; existing IDs cannot be overwritten |
| Compatibility | Changed failure data, properties, settings or numerical source are rejected; row reordering and run-out addition remain compatible |
| Plotting | Open run-out markers, lower-bound arrows and legends checked in PNG and TikZ pathways |
| End-to-end execution | Script 0 → JSON → Script 1 → CSV/plot completed on AISI316L; Script 2 comparison completed for the 64-run outputs |
| Native life outputs | Synthetic cycle-valued model bypasses master-curve fitting; errors, run-out masking, parity export and rejection of mismatched units or invalid lives checked |

The run-out count depends on the selected plain-specimen scope, not on cycle
thresholds. The workbooks themselves retain notched and other original observations.
The [workbook evidence](workbook-migration-evidence.json) records original public
blob IDs and the two changed ZIP members in each migrated file. No Dataset entry,
worksheet or other archive member was added or changed.

| Material | Imported failures | Imported run-outs |
|---|---:|---:|
| AISI316L | 31 | 13 |
| 42CrMo4_Santus | 70 | 15 |
| Al7075_T6 | 87 | 14 |
| Ti6Al4V_andrea_(AXIAL) | 58 | 0 |
| Total | 246 | 42 |

The sweep repeats those observations for 16 method/variant runs. Run-outs have
finite evaluated damage parameters but no ordinary life/damage errors. Their
predicted lives remain unavailable when the existing slope safeguard suppresses
life inversion. Adding them does not improve or worsen an error statistic.

## Reproduction

```bash
MPLBACKEND=Agg uv run pytest -q
MPLBACKEND=Agg python tests/run_material_validation.py /tmp/wg43-validation
python -m ruff check .
python -m pip check
python tools/generate_method_docs.py --check
python tools/check_docs.py
```

`tools/migrate_calibrations.py` and `tests/compare_material_regression.py` were one-off scripts of
that verification and are not kept in the repository; the pytest suite covers the same checks.

For an exact comparison against the previous pipeline, extract baseline commit
`e70a8af3a7a06cfb49626bff16c04605c884830b` into a separate checkout. Run its
`tests/run_material_validation.py /tmp/wg43-baseline-results` from that root,
then run the revised sweep from the current root and compare all baseline columns.

```bash
MPLBACKEND=Agg python tests/run_material_validation.py /tmp/wg43-revised-results
```

The comparator requires exact agreement, including undefined values, rather than
a numerical tolerance. The fresh rerun on 16 October 2026 uses the recovered five
local commits plus the native-life routing extension. Original `Results/*.csv`
files remain unchanged. Existing empirical material-property estimation warnings
occur in both baseline and revised sweeps and retain their previous interpretation.

The last command reports that no legacy columns remain on the migrated repository.
For a new calibration and saved plot, use the commands in the
[workflow](calibration-workflow.md), selecting a separate calibration directory
when an isolated example is preferred.

The active coefficient regression fixture uses a deliberately small three-plane
set, including a 45-degree plane that retains all 11 axial and torsional failure
observations. It was generated by the unchanged baseline commit in a separate
checkout. Its finite-grid objective minima are separated by at least 0.2697,
and the tests recompute those objectives rather than asserting only stored values.
This is not a converged scientific calibration. The earlier two-plane fixture is
retained for historical provenance but is not used to assert a unique MSWT
coefficient. It kept only five torsion points, making the MSWT objective invariant
to a scale factor within floating-point rounding. Different CI runners selected
different tied coefficients. The test fixture was corrected without changing
the optimiser, material data or fatigue formulations. The
four shipped migration records retain the former coefficients and identify their
original optimisation provenance as unverified. Their metadata was completed
before publication without changing the transferred numerical values or curves.

The numerical differences introduced by this task are confined to preserving
censored observations and leaving their error fields undefined. No fatigue
criterion, stress/strain calculation, critical-plane search or life-suppression
rule was changed. These results establish software regression consistency under
the investigated inputs, not independent predictive validation.

For a reproducible plotting example, use Script 0 and Script 1 with AISI316L and
BP as described in the workflow. Six axial failure observations define the curve;
the remaining failures and 13 run-outs provide qualitative context. Open circles
and arrows denote lower-bound lives. The shallow fitted slope suppresses life
inversion, so no life-parity plot is saved. This example illustrates data handling
rather than independent validation.
