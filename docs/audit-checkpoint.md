# Continuation checkpoint — 7 October 2026

## Completed

- Resumed the existing working tree and preserved the previously published Böhme–Papuga and CAIM implementations and their tests. BP and CAIM are distinct formulations and remain separately selectable.
- Reviewed all 20 registered base methods. Corrected engineering shear factors, Findley plane selection, Li mean shear, Liu's removable zero-shear singularity, harmonic strain decomposition, zero-shear cancellation, and critical-plane ties. See [method-audit.md](method-audit.md) for equation-by-equation scope and limitations.
- Disabled CARSPA, DANGVAN, OTT, SWTD and ZHU_EDP, including direct evaluation guards, with reasons and restoration requirements. Removed obsolete ZHU_EDP calibration.
- Added source metadata and a [bibliography](bibliography.md) for every method. CSVs identify the method reference and implementation revision; comparisons exclude disabled and pre-audit results by default.
- Added `Error_log10_dp = log10(DP_e/DP)` to export, statistics and plots for all slopes, including exactly flat curves. Fitting no longer inverts nearly flat curves. Invalid fits/calibration objectives fail explicitly.
- Preserved experimental workbooks and historical result CSVs. Corrected results must be regenerated; fitted method constants need recalibration before quantitative model rankings.

## Verification completed

- `MPLBACKEND=Agg python -m unittest discover -s tests -q`: **20 tests passed**, including the existing BP/CAIM reference cases and new numerical/reporting regressions.
- `MPLBACKEND=Agg python tests/run_material_validation.py /tmp/fatigue-audit`: **64 runs / 3,936 rows passed**, covering 15 retained base models plus FS_ext on all four supplied materials. Checks include finite damage/log ratios, source metadata, ratio recomputation and suppressed life predictions in DP_DIFF mode.
- `MPLBACKEND=Agg python 2_compare_results.py --results_dir /tmp/fatigue-audit --out_dir /tmp/fatigue-audit-plots`: passed, **21 PNG plots**, plus TikZ and statistics. All four materials have log10 damage statistics for 16 variants; the three with usable negative slopes also have life statistics. AISI316L has 496 damage-space rows and no life predictions. Inspected its damage boxplot visually.
- An exactly-flat synthetic curve passed CSV export, plotting and statistics; legacy CSVs missing DP_fit passed the explicit opt-in path without crashing.
- Five retained parameter calibrations passed a coarse-grid smoke check on 42CrMo4_QT: FS, FIN, LI, MGSE_ZHU, MSWT. This is an execution check, not final recalibration of workbook constants.
- Python parsing/compilation and `git diff --check` passed. Final diff reviewed; no raw data or archived Results changes.
- Non-failing warnings: empirical steel-property fallbacks where material coefficients are missing; existing Matplotlib deprecation/marker warnings.

## Unresolved evidence and exact next action

Full original-equation/plane-rule verification is still incomplete for **LIU2021, MGSE_ZHU, MKBM and LI**. MSWT is associated with the reported scalar modified-SWT form and is not claimed to implement Jiang's complete incremental model. These limitations are explicit in the audit; passing implementation tests is not independent scientific validation.

**Next action:** obtain an accessible full original Liu et al. article (DOI `10.1016/j.ijpvp.2021.104532`), inspect Eq. (27) and the definitions of its normal strain and shear ranges, and compare them with `fatigue/models/cp_methods/liu2021.py`. Add an independently calculated published load case, then correct or disable the model if necessary. Continue in order with Zhu 2018 (MGSE_ZHU/MSWT), Li 2011 (MKBM), and Li 2021 (LI), using the DOI mapping already saved in the bibliography. Do not repeat completed BP/CAIM integration or the general audit. Then recalibrate affected parameters and regenerate comparisons before interpreting rankings.

No pending test failure is known. The unresolved item is source completeness, not a failed test. The shared elastic, harmonic and finite-plane-grid limitations are documented in the audit.

## Files changed in this audit

- `0_calibrate_material_params.py`
- `1_run_material.py`
- `2_compare_results.py`
- `README.md`
- `docs/adding-methods.md`
- `docs/architecture.md`
- `docs/audit-checkpoint.md`
- `docs/bibliography.md`
- `docs/method-audit.md`
- `fatigue/eval/evaluator.py`
- `fatigue/eval/fitting.py`
- `fatigue/eval/metrics.py`
- `fatigue/models/cp_methods/car_spa.py`
- `fatigue/models/cp_methods/dangvan.py`
- `fatigue/models/cp_methods/fin.py`
- `fatigue/models/cp_methods/fs.py`
- `fatigue/models/cp_methods/gsa.py`
- `fatigue/models/cp_methods/gse.py`
- `fatigue/models/cp_methods/li.py`
- `fatigue/models/cp_methods/liu2021.py`
- `fatigue/models/cp_methods/mgse_yu.py`
- `fatigue/models/cp_methods/mgse_zhu.py`
- `fatigue/models/cp_methods/mkbm.py`
- `fatigue/models/cp_methods/mswt.py`
- `fatigue/models/energy_based_methods/swtd.py`
- `fatigue/models/energy_based_methods/zhu_edp.py`
- `fatigue/models/invariant_methods/ottosen.py`
- `fatigue/models/references.py`
- `fatigue/models/registry.py`
- `fatigue/models/utils/plane_history.py`
- `fatigue/report/export_csv.py`
- `tests/run_material_validation.py`
- `tests/test_method_audit.py`
