# Calibration records

Each material directory contains complete JSON records. Script 0 creates a new
record for each successful execution. Records are selected by their declared
UTC creation time and ID, with compatibility checks, instead of file timestamps.

The four initial records transfer the former workbook coefficients without
re-optimising them and reproduce all 15 active methods' master curves. Their
original coefficient-fitting settings are unverified. These records preserve
the previous analysis results and do not establish independent validation.

See the [calibration workflow](../docs/calibration-workflow.md) for the schema,
units, selection rules, provenance and commands. Keep useful new records under
version control, and inspect their provenance before committing them. Routine
plots and analysis exports belong in ignored result directories.

## Use with multifat-pipeline

`multifat-pipeline <workbook> --constants Calibrations/<material>/<id>.json` takes the plane
constants (k_FS, k_FI, L_LI, k_MGSE, a_MSWT) from a record instead of calibrating them; the
record must belong to the workbook's material ID; a record made on a different number of failures is
used but flagged in the sidecar (`constants_source.warning`). Its a_MSWT is in the weighting of the former
workbooks and is converted to the printed one, a' = a/(2 - a). The calibration sidecar records
the file, its SHA-256 and the calibration ID. The pipeline does not write new records here; its
own results go to `<MATERIAL>_calibration.json` in the results folder, which `--constants` also
accepts. The four records belong to the previous collection (material IDs `42CrMo4_QT`,
`AISI316L`, `Al7075_T6`, `Ti6Al4V_andrea_(AXIAL)`); the revised workbooks change the material
inputs, so a record reproduces former constants, not a calibration of the new data.

## Hand-entered TVEIT2025 constants

`<material>/TVEIT2025_stage1.json` holds the constants of the Tveit et al. (2025) model, entered by
hand (stage 1: the model runs with supplied constants; their calibration is not implemented yet).
These files use the sidecar layout (`params` and `provenance.workbook`), not the record schema.
`provenance` states how the constants were obtained. Each file sets:

- `A_TVEIT`, `C_TVEIT`, `K_TVEIT` and `L_TVEIT`;
- `Sigm1_TVEIT` and `Taum1_TVEIT`, the fatigue limits these constants were calibrated with. They
  replace the workbook Sigm1 and Taum1 for TVEIT2025 only, because the constants hold only
  together with them;
- the numerical settings of the time integration (`npts_TVEIT`, `min_blocks_TVEIT`,
  `max_blocks_TVEIT`, `tol_alpha_TVEIT`, `tol_damage_TVEIT`), at their default values.

`notes` explains each entry. Edit a copy to change a run:

    multifat-pipeline Experimental_Data/6_42CrMo4_QT.xlsx --models TVEIT2025 --no-calibrate --constants Calibrations/42CrMo4_QT/TVEIT2025_stage1.json

## Provisional TVEIT2025 records

The stage-1 TVEIT2025 constants are source-specific, manually supplied coefficients; they are neither independently validated nor transferable between experimental campaigns. The AISI316L record proposed in the contributor branch is deliberately withheld from this integration: its PSI load-cell-corrected axial stress scale (x1.1558) and reference fatigue limit (271.6 MPa) do not match the current Janssens workbook/reference (278.5 MPa). A source-compatible workbook and reviewed calibration are needed before enabling this record. Do not change the reference limit without re-estimating A, C, K, L.

For TVEIT2025, the model declares `ComparisonCapabilities(life_prediction="native", dp_comparison=False)`. The pipeline exports `Metric_mode=Nf_direct`, with `Nf_expected=Nf_native` for compatibility with fatigue-life plots; the diagnostic fatigue index remains in `CP_value` but `DP_fit` and its error are undefined. Other methods retain the normal calibrated master-curve comparison. Select `multifat-compare --results_dir Results/current --out_dir Results/current/Comparison --plot-kind life` for fatigue-life plots, `--plot-kind damage` for DP-only plots, or `--plot-kind both` (default). Tveit is included only in life plots. A `native_life_prediction_summary.csv` counts infinite-life predictions for failed specimens and run-outs, which are not converted into finite logarithmic residuals. The method is still at stage 1 and its constants are not independently validated.