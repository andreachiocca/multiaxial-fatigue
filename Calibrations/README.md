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
accepts. The four records belong to the previous collection (material IDs `42CrMo4_Santus`,
`AISI316L`, `Al7075_T6`, `Ti6Al4V_andrea_(AXIAL)`); the revised workbooks change the material
inputs, so a record reproduces former constants, not a calibration of the new data.

## Hand-entered TVEIT2025 constants

`<material>/TVEIT2025_stage1.json` holds the constants of the Tveit et al. (2025) model, entered by
hand (stage 1). They take precedence over the calibration of the pipeline (next section).
These files use the `hand_entered_constants` kind (#60): the constants under
`hand_entered_constants`, the batch under `material_id`, and in `provenance` how the constants were
obtained. They are not pipeline calibrations, so they carry no implementation revision or dataset
fingerprint and skip that check; the run sidecar labels their source "hand-entered constants" and
copies the provenance. Each file sets:

- `A_TVEIT`, `C_TVEIT`, `K_TVEIT` and `L_TVEIT`;
- `Sigm1_TVEIT` and `Taum1_TVEIT`, the fatigue limits these constants were calibrated with. They
  replace the workbook Sigm1 and Taum1 for TVEIT2025 only, because the constants hold only
  together with them;
- the numerical settings of the time integration (`npts_TVEIT`, `min_blocks_TVEIT`,
  `max_blocks_TVEIT`, `tol_alpha_TVEIT`, `tol_damage_TVEIT`), at their default values.

`notes` explains each entry. Edit a copy to change a run:

    multifat-pipeline Experimental_Data/6_42CrMo4_Santus.xlsx --models TVEIT2025 --no-calibrate --constants Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json

## TVEIT2025 calibration inputs and records

Without supplied constants, `multifat-pipeline` calibrates TVEIT2025 by the authors' procedure on
the workbook's own data (see the [method guide](../docs/method-guides/TVEIT2025.md)).
`<material>/TVEIT2025_inputs.json` holds, by hand, only what a workbook cannot: a fatigue limit at
a tensile mean stress (`mean_stress_fatigue_limit`, Al7075_T6), the material type for the FKM
estimate of A (`material_type`: steel for 42CrMo4_Santus and 42CrMo4_Reis; titanium alloy, which
FKM does not cover, for Ti6Al4V_Meneghetti and Ti6Al4V_Berto), a hand-entered A, the source of the
calibration tests and the time-integration settings. The calibration only reads it. A material
without the file is calibrated with the Goodman estimate of A.

    multifat-pipeline --material Al7075_T6 --models TVEIT2025 --calibration-archive Calibrations

calibrates and archives the record here as `<material>/<calibration_id>.json`. A later run with
the same workbook, overlay, inputs file, procedure version and implementation revision reuses that
record instead of calibrating again (a calibration takes up to several minutes). Editing the
inputs or the data gives a new calibration and a new record; old records stay.

The records committed for Al7075_T6, 42CrMo4_Santus, 42CrMo4_Reis, Ti6Al4V_Meneghetti and
Ti6Al4V_Berto were made this way from a clean checkout; each holds the constants, the three
calibration tests, the C range of the equivalent sets and the calibration confidence with its
findings. AISI316L has none until its fatigue limit is settled (#102).

## Provisional TVEIT2025 records

The stage-1 TVEIT2025 constants are source-specific, manually supplied coefficients; they are neither independently validated nor transferable between experimental campaigns. The AISI316L record proposed in the contributor branch is deliberately withheld from this integration: its PSI load-cell-corrected axial stress scale (x1.1558) and reference fatigue limit (271.6 MPa) do not match the current Janssens workbook/reference (278.5 MPa). A source-compatible workbook and reviewed calibration are needed before enabling this record. Do not change the reference limit without re-estimating A, C, K, L.

For TVEIT2025, the model declares `ComparisonCapabilities(life_prediction="native", dp_comparison=False)`. The pipeline exports `Metric_mode=Nf_direct`, with `Nf_expected=Nf_native` for compatibility with fatigue-life plots; the diagnostic fatigue index remains in `CP_value` but `DP_fit` and its error are undefined. Other methods retain the normal calibrated master-curve comparison. Select `multifat-compare --results_dir Results/current --out_dir Results/current/Comparison --plot-kind life` for fatigue-life plots, `--plot-kind damage` for DP-only plots, or `--plot-kind both` (default). Tveit is included only in life plots. A `native_life_prediction_summary.csv` counts infinite-life predictions for failed specimens and run-outs, which are not converted into finite logarithmic residuals. The stage-1 constants are not independently validated; the calibrated records above replace them only when no constants are supplied.
## Published constants

`<material>/<MODEL>_published.json` holds constants printed in a paper, entered by hand. The
layout is `published_constants`, `material_id` (checked against the workbook) and `reference`;
`notes` states scope and caveats. These are not calibrations of the workbook data, so they carry
no failure fingerprint. `42CrMo4_Reis/SSF_published.json` is the stress scale factor surface of
Anes et al. (2014), Eq. (31), which was fitted on the tests of that workbook (in-sample).

    multifat-pipeline --material 42CrMo4_Reis --models SSF --constants Calibrations/42CrMo4_Reis/SSF_published.json

A cross-batch transfer with 42CrMo4_Reis as the design batch accepts the same record (decision
d103, issue #71); it is the only supplied record a transfer accepts:

    multifat-pipeline --material 42CrMo4_FABEST --design-material 42CrMo4_Reis --models SSF --constants Calibrations/42CrMo4_Reis/SSF_published.json
