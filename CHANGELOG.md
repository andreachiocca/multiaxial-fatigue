# Changelog

All notable changes to this project. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project has no tagged release yet.
Finding IDs (F, U, P, N) refer to the audit ledgers kept with the carry-forward review.

## [Unreleased] (carry-forward of October 2026)

### TVEIT2025 calibration by the authors' procedure (stage 2)

- Without supplied constants, `multifat-pipeline` calibrates `A_TVEIT`, `C_TVEIT`, `K_TVEIT` and
  `L_TVEIT` on any fatigue-life workbook (`eval/tveit2025_calibration.py`; Tveit et al. 2025,
  Sect. 4; Ottosen et al. 2008, Sect. 5): the workbook's own `Sigm1`, `Taum1`; `A` in the authors'
  order (fatigue limit at a tensile mean stress, a hand-entered value, FKM for a declared material
  type, Goodman); `C`, `K`, `L` by Eq. (25) on three axial R = -1 failures of the source of the
  Summary limits, with Eq. (25) profiled over `C` because uniaxial data do not identify it.
  `--constants` takes precedence; `--no-calibrate` turns it off.
- Values entered by hand (`--constants`, or a Summary column named like a declared parameter)
  are not used by a model's own calibration, so it runs only when none is entered. With only some
  entered (e.g. `Sigm1_TVEIT`, `npts_TVEIT` or `A_TVEIT` alone), the model is not calibrated and
  is skipped with the reason (status `partly_supplied`), instead of being fitted with other
  inputs than those of its evaluation (review of #114).
- A repeated cross-batch transfer reuses its calibration: the record archived under
  `<target>_from_<design>` is searched too (review of #114).
- New extension point for models that calibrate their constants by their own procedure
  (`calibrated_params`, `calibrate_constants`): the result is kept under
  `model_calibrations.<MODEL>` in the sidecar and the archived record, and its tests become
  design rows (Flag 0). No model name in `cli/`.
- Calibration confidence (`eval/calibration_confidence.py`): named findings, major or minor, each
  with what to try (estimated A, failures below `Sigm1`, limits and tests of different sources,
  points close together, a constant on its search bound, the unidentified `C`, a calibration
  point whose integration reached `max_blocks_TVEIT`, ...); level high, medium or low. Printed,
  recorded and noted on the run-report row. A failure below `Sigm1` is major only when it failed
  before the reference life of `Sigm1`; after it, it agrees with a limit defined at that life
  (minor).
- The record names the inputs file relative to the repository, not by a local path.
- An identical calibration (workbook bytes, overlay, inputs file, procedure version, implementation
  revision) is reused from `Calibrations/<material>/` or the run's archive.
- Optional inputs `Calibrations/<material>/TVEIT2025_inputs.json` (schema 2) for Al7075_T6,
  42CrMo4_Santus, 42CrMo4_Reis, Ti6Al4V_Meneghetti and Ti6Al4V_Berto, and their first records;
  a record for AISI316L (Goodman A, confidence medium) once its fatigue limit was settled (#120).
  The records are regenerated for implementation revision 2026-11-18, with identical constants.
- The TVEIT2025 method sheet (`docs/method-sheets/sheets/TVEIT2025.tex`) describes this
  calibration.
- Calibration guide `docs/calibration-guides/TVEIT2025/` (LaTeX and PDF, with its figure
  script): the six steps, the confidence findings, and worked examples from the committed
  records (Al7075_T6, both Ti6Al4V and both 42CrMo4 batches) and from Al7075_T651_Wang.
- Changed numerically: TVEIT2025 now runs, with calibrated constants, in pipeline runs without
  `--constants` (it was skipped for missing constants); no other model changes.

### AISI316L tensile fatigue limit adjusted to staircase level (issue #102, decision d114)

- Adjusted the Janssens row tensile fatigue limit $\sigma_{-1}$ in `8_AISI316L.xlsx` via `set_summary`
  in `8_AISI316L.overlay.json` from $278.5\text{ MPa}$ ($241 \times 1.1558$) to $271.6\text{ MPa}$
  ($235 \times 1.1558$, the staircase load level of Janssens 2018 with PSI load-cell correction, as
  used in the PSI calibration of Tveit).
- 5 axial $R = -1$ failures previously sat below $\sigma_{-1}$ ($271.6$ and $277.4\text{ MPa}$); with
  $\sigma_{-1} = 271.6\text{ MPa}$, 3 failures sit above the threshold and 2 on it (one with
  $2.47 \times 10^6$ cycles outlasting the $2 \times 10^6$ reference life).
- `IMPLEMENTATION_REVISION` bumped to `2026-11-18`.
- 21 models with $\sigma_{-1}$ dependence improve their RMS load error on AISI316L; models without
  $\sigma_{-1}$ dependence and all other materials remain identical.
- Post-overlay dataset fingerprint updated to `5e2e99b9e3456f7d96a5dde8e66e4213e241cf88592a068d679506490824445d`.
- Regenerated AISI316L rows in `load_error_ranking.csv` and refreshed provenance manifests.

### Documentation consistency cleanup (no numerical change)

- MMCCM, MDC (guides, audit rows, sheets, MMCCM docstring): `multifat-pipeline` fills missing
  Basquin constants by identifiable fits on the R = -1 tests and, for MMCCM, the plastic terms by
  uniform-material-law estimates; otherwise the model is skipped (the text said nothing is fitted,
  or that runs stop with an error).
- References corrected to the Crossref records of `docs/research/doi-records.json`: MMCCM (IJF 83,
  135-149), MMMH (title; IJF 103, 294-308), MDC, RFB, OTT and COIN titles and authors, LIU2021
  authors; `references.bib` entries `wang2016`, `mertens2017`, `tveit2024`, `schmidt2023`.
- PCR docstring: PCr is maximised on the default 24 x 48 grid (PI keeps 48 x 192). SSF docstring:
  the published constants are in `Calibrations/42CrMo4_Reis/SSF_published.json` (also in the
  catalogue). DANGVAN docstrings cite the registered Dang Van (1973).
- MWHS guide and sheet state the shared identifiability rule; WANG2022 guide states the 1 %
  tolerance its test asserts; the MDC test docstring states the 4 % McDiarmid tolerance.

### EK: Eq. (18) mean-stress range enforced as the loading scope (maintainer decision, 17 November 2026)

- A case with sigma_m/sigma_bar outside -1.06 to 4.4 (Ellyin & Kujawski 1993, Eq. 18, p. 61) now
  raises `LoadingScopeError` and EK is skipped for that material. The guide stated the range, but
  the code evaluated such cases.
- No supplied workbook has a case outside the range. EK outputs on every workbook are identical
  before and after, so the revision and the ranking are unchanged.
- New audit row for EK, which had none; guide and LaTeX sheet updated; test of the bounds.

### Extension points: public registry accessors, Job and RowGroup dataclasses, FatPy mapping (#81)

- Replaced private cross-package imports (`_available_cp_methods`, `_available_direct_models`,
  `_PROBE_PARAMS`) with public registry functions `available_cp_methods()`,
  `available_direct_models()`, `base_models()`, `model_names()`, and `PROBE_PARAMS`.
- Replaced hard-coded model-name checks in `multifat-run` with capability queries (`load_error`
  capability in `ComparisonCapabilities`, with `TVEIT2025` override) and polymorphism
  (`effective_strengths` method on `BohmePapuga`/`CAIM`).
- Replaced 11-field job tuples in `multifat-pipeline` with a frozen `Job` dataclass.
- Replaced separate table row lists in `export_csv.save_fatigue_csv` with a frozen `RowGroup` dataclass.
- Added data correspondence mapping `src/multifat/data/fatpy_keys.json` and conversion utilities
  `to_fatpy`, `from_fatpy`, and `load_fatpy_mapping` in `src/multifat/data/fatpy.py`.
- Added tests in `tests/test_fatpy.py` enforcing bidirectional conversion, unit scaling (Young's
  modulus $E$ GPa/MPa), notation consistency with `docs/notation.md`, absence of private
  cross-package imports in `cli/`, `tools/`, and `examples/`, and absence of model-name branching
  in `cli/run.py`.

### MATAKE and MCD: circumscribed-circle shear amplitude by default (maintainer decision, 17 November 2026)

- MATAKE and MCD now use half the longest chord of the shear path out of phase
  (`np_shear_rule = "longest_chord"`), i.e. the minimum circumscribed circle for one harmonic. This
  is the amplitude with which Papadopoulos et al. (1997) restate both criteria and with which their
  Table 1 columns are reproduced (39/42 and 41/43 rows); RSS before. The other RSS methods (FIN,
  MDC, MGSE_YU, MGSE_ZHU, MSWT, GSE) keep RSS as a documented exception.
- Numerical change, out of phase only. Results on proportional datasets are unchanged (Santus,
  Berto, Meneghetti). Out-of-phase rows change on FABEST (23), Reis (7), AISI316L (12),
  Al7075_T6 (13), Al7075_T651_Wang (6) and the out-of-phase limit states of the fatigue-limit
  datasets. Measured with the pipeline on every workbook, before and after; no other model is
  touched.
- `IMPLEMENTATION_REVISION` 2026-11-17. Ranking, batch comparison and figures regenerated. Guides,
  LaTeX sheets, audit, `docs/shear-amplitude.md` and README updated.

### LaTeX method sheets for every registered model

- 46 new sheets in `docs/method-sheets/sheets/`; every registered model now has one (52 files,
  PDS_SWT/PDS_FS sharing one).
- `methods.tex` groups them by family. The pending list of `tests/test_method_sheets.py` is empty,
  so a new method cannot be registered without its sheet.
- Each sheet states, from the method's guide, docstring, audit row and tests:
  - the original equations and every symbol;
  - calibration and life law;
  - repository choices and verification.
- Worked-example numbers are those printed by the tested guides.
- `references.bib` holds the original papers (DOI-checked, duplicates merged).
- `docs/method-guides/SWTD.md`: Table 3 torsion rows match (2/sqrt3) SWTD = tau_a gamma_a/sqrt3,
  as the test asserts. The guide printed the inverse factor.

### New dataset: 7075-T651 of Wang, Shang, Sun & Chen (2018) (maintainer decision, 16 November 2026)

- New workbook `18_Al7075_T651_Wang.xlsx`: 52 constant-amplitude failures of thin-walled tubes.
  - Tests: uniaxial and torsion at R = -1 and with means; tension-torsion with means at 0-180
    deg.
  - Built by `tools/build_wang_al7075_t651.py` from Tables 1-5. The transcription is shared with
    the WANG2018 tests.
  - Sigm1 and Taum1 are the Table 2 S-N curves at 1e6 cycles; nu = E/(2G) - 1.
- Manifest, provenance and dataset tests (`tests/test_wang_al7075_dataset.py`).
- The dataset joins the homepage ranking: 45 models run on it.
  - Skipped for scope: WANG2022, GOLOS, KHALIJ, KAKUNO (no Sig0), RFB and SIH
    (t_-1/f_-1 = 0.573).
  - Skipped for missing constants: TVEIT2025, SSF.
  - WANG2018 ranks first on it; this is the dataset its authors used, not an independent check.
- `dataset_fingerprint` hashes the failure histories rounded to 1e-6 MPa when a failure has a
  phase that is not a multiple of 90 deg, as for the run-outs. Their cos/sin values differ in the
  last bit between platforms, which made the Windows fingerprint of this workbook differ. Every
  other recorded fingerprint is unchanged.

### Plane shear amplitude selectable for every method (maintainer decision, 16 November 2026)

- `--shear-amplitude-algorithm` of `multifat-run` and `multifat-calibrate` used to cover Findley
  only. It now applies to every model that exposes a plane shear amplitude: all critical-plane
  methods and SIH (`multifat.eval.shear_amplitude_option`).
- Choices:
  - `method_default`: each method's own amplitude; unchanged outputs;
  - `longest_chord` = `mcc_harmonic`: for one harmonic, half the longest chord equals the MCC
    radius;
  - `mce_harmonic` = `rss`.
- Models whose amplitude is fixed by their definition (BP, RFB, PAPADOPOULOS, LIUMAHADEVAN), or
  that use none, run unchanged with a warning.
- `multifat-calibrate` fits FS, FIN, LI, MGSE_ZHU and MSWT with the selection. An explicit
  selection records `amplitude_settings` for each of them.
- `multifat-pipeline --constants` refuses constants calibrated with an amplitude other than a
  method's own (generalising the Findley guard). The pipeline itself still evaluates
  `method_default`.

### SIH: shear stress intensity hypothesis, Zenner, Simbürger & Liu (2000) (#13, N15)

- New integral model `SIH` (`integral_methods/sih.py`): sigma_eq,a = {15/(8 pi) int [a tau_a^2
  (1 + m tau_m^2) + b sigma_a^2 (1 + n sigma_m)] dOmega}^(1/2) (Eq. 6). The constants a, b, m, n
  come from Table 1 and tau_Sch from Eq. (7). The model returns the index sigma_eq,a / sigma_W.
- Maintainer decisions:
  - shear amplitude on a plane: the longest-chord amplitude (`np_shear_rule`), because the
    available source does not define it for rotating shear vectors; `rss` is selectable;
  - mean stresses need `Sig0`.
- Material scope: 0.577 <= tau_W/sigma_W <= 0.866.
- Evidence reviewed by the maintainer, who confirmed the longest chord: the Section 3 statement that Eq. (5) does not depend on
  the phase holds exactly with `rss`, not with the longest chord (-7 % at 90 deg).
- Checked: the Table 1 factors (12/7, 4/21, 15/14) as sphere integrals; the four limits; Eqs.
  (5) and (5b).
- Guide, LaTeX sheet, catalogue N15, audit row; ranking and batch comparison regenerated.
- Runs on Reis, FABEST, CX18, 982FA and Delahay. Skipped on the others: mean stresses without
  Sig0, or the ratio outside the material scope.

### WANG2018: Wang, Shang, Sun & Chen (2018) weighted-plane criterion with mean shear stress (#13, N17)

- New case-level model `WANG2018` (`cp_methods/wang2018.py`):
  - the critical plane is the weighted mean of the instantaneous maximum-shear planes over one
    cycle (Eqs. 3-7);
  - of the candidates phi_hat and phi_hat + 90 deg, the larger tau_eq is kept;
  - tau_eq = sqrt([tau_a + k_tau |tau_m|]^2 + 1/3 [sigma_a + k_sigma sigma_m]^2) (Eq. 18), with
    k_tau and k_sigma from t_-1, f_0, t_0 (Eq. 15) or from the Eq. (19) estimates.
- Tension-torsion only, the paper's loading; other stress components are out of scope (maintainer
  decision). Constants that violate Eq. (17) are out of material scope.
- Life on the R = -1 torsion failures (`design_loading = "torsion"`), the paper's torsional S-N
  curve. The variable-amplitude procedure is not implemented.
- Paper data: Tables 3-5 transcribed in `tests/data/papers/wang2018_tables3_5.csv`. The Section 5.1
  estimated crack orientation (70 deg) is reproduced for every 30 and 90 deg out-of-phase test.
- Catalogue N17 is now implemented. Its family is corrected from integral and hybrid to critical
  plane.
- Also added: guide, LaTeX sheet, audit row. The ranking and the batch comparison were
  regenerated. WANG2018 runs on 6 life and 6 limit datasets. It is skipped on Ti6Al4V_Meneghetti
  (4 torsion tests) and on 25CrMo4 and 34Cr4 (biaxial normal stresses).

### LaTeX method sheets for the review article

- New `docs/method-sheets/`: one LaTeX document (`methods.tex`) with a shared notation section
  (frame, harmonic loading, invariants, plane quantities, material constants, master curve,
  load-equivalent error, plane search) and one sheet per method in `sheets/`, each with Source,
  Idea and scope, Inputs (every symbol and unit), Formulation, Constants and calibration, Life,
  Repository choices, Worked example and Verification; bibliography in `references.bib`.
- First sheets: KAKUNO, KHALIJ, MATSUBARA, PDS_SWT/PDS_FS. The other 46 models are listed as
  pending in `tests/test_method_sheets.py`, a list that may only shrink.
- `tests/test_method_sheets.py` checks coverage of the registry, `\input` targets, the block order,
  citation keys and balanced braces and environments.

### Homepage: cross-batch transfer figure (#71)

- New README section "Same steel, another batch" and figure `docs/images/batch_transfer.png`
  (`examples/plot_batch_transfer.py`, homepage style, part of `generate_gallery.py`): median and
  quartiles over the methods of the signed load error of each 42CrMo4 batch with its own constants
  and with those of the other two batches. A test ties the quoted percentages to the data.

### PDS_SWT, PDS_FS: Albinmousa & Al Hussain (2021) polar damage sum (#13, new N59/N60)

- New integral models `PDS_SWT` and `PDS_FS` (`integral_methods/polar_damage_sum.py`): the area
  1/2 int D(phi)^2 dphi of the polar curve of the SWT or Fatemi-Socie parameter over every plane
  perpendicular to the surface (Eqs. 8-9, 12-14), life via the master curve (Eqs. 17-19).
- Elastic, high-cycle reading of a strain-controlled LCF method (maintainer decision): linear-elastic
  strains replace the measured loops; documented in the module, guides and audit.
- New calibrated constant `k_PDS` (k_phi), fitted like `k_FS` on the R = -1 tension and torsion
  failures within [0, 10] (`cli/calibrate.py`, block `PDS_FS`) and set by the limit equalisation on
  fatigue-limit datasets (`eval/endurance.py`). Identified on FABEST (2.78), AISI316L (0.40),
  Ti6Al4V_Berto (9.06); not identifiable on Santus, Reis and Al7075_T6.
- New catalogue records N59 (PDS_SWT) and N60 (PDS_FS, variant), guides, audit rows; ranking and
  batch comparison regenerated.
- `docs/calibration-workflow.md`: the initial design-test tolerance is |R + 1| <= 1e-3 since #77
  (it still said 0.02).

### MATSUBARA: Matsubara, Hayashida & Kano (2018) criterion (#13, N09)

- New invariant model `MATSUBARA` (`invariant_methods/matsubara.py`): unified equivalent shear
  stress amplitude sqrt(J_2,amp,eq) = sigma_T sqrt(J_2,amp)/(sigma_T - S_max + sqrt3 sqrt(J_2,amp))
  (Eq. 15) with the Mamiya amplitude (Eq. 13) and the maximum S_max of the largest-range normal
  stress; life from the master curve on the R = -1 uniaxial tests (the paper's unified S-N curve).
- New optional constant `sigma_T` (true fracture strength, `docs/notation.md`); without it
  sigma_T = 1.84 Su for 400 <= Su <= 1000 MPa, else the material is out of scope (maintainer
  decision). Runs on 8 materials.
- Evidence (`tests/methods/test_matsubara.py`): Eq. (15) against Eq. (9) and the uniaxial line,
  the Mamiya amplitude in the time domain, S_max against a brute-force search. The paper prints no
  per-test value.
- Guide `docs/method-guides/MATSUBARA.md`, catalogue record N09, audit row; ranking and batch
  comparison regenerated.

### KAKUNO: Kakuno & Kawada (1979) criterion (#13, G04)

- New invariant model `KAKUNO` (`invariant_methods/kakuno_kawada.py`): Eq. (4), the octahedral
  shear stress amplitude against linear octahedral normal amplitude (B) and mean (C) terms, with
  A, B, C of Eqs. (7)-(9). Constants: `Sigm1`, `Taum1`; `Sig0` (pulsating limit) for mean normal
  stresses.
- Proportional loading only (principal stress amplitudes). Without `Sig0` a case with a mean
  normal stress is out of scope; the paper's k_1 = 0.83 for steels is not used (maintainer
  decision). Runs on 42CrMo4_Santus and Ti6Al4V_Berto.
- Evidence (`tests/methods/test_kakuno_kawada.py`): Eqs. (7)-(10) and (12) to machine precision;
  Fig. 5 tau intercept reproduced, sigma intercept a strict xfail. The paper prints no table.
- Guide `docs/method-guides/KAKUNO.md`, catalogue record G04, audit row; ranking and batch
  comparison regenerated.

### KHALIJ: Khalij, Pagnacco & Troian (2017) criterion (#13, N19)

- New invariant model `KHALIJ` (`invariant_methods/khalij.py`): Eq. (23) with the authors' OP
  setting alpha = 1/(3 r^2), evaluated in its frame-invariant form
  J_2,a/t_-1^2 + (1 - alpha) I_1,a/f_-1. Constants: `Sigm1`, `Taum1`.
- Table 6 (OP column) reproduced: 51 of 53 rows, with r from the Table 1 limits; two 3%C
  cast-iron rows are strict xfails (`tests/methods/test_khalij.py`).
- Zero-mean plane stress only (maintainer decision): the mean-stress terms of Eqs. (29)-(31) are
  not dimensionless as printed and no reading reproduces Table 7. Cases with a mean stress are out
  of scope, so the model runs on 42CrMo4_FABEST, 42CrMo4_Reis and the zero-mean limit datasets.
- Guide `docs/method-guides/KHALIJ.md`, catalogue record N19, audit row; ranking and batch
  comparison regenerated.

### Observation type of each Dataset row (decision d113, #18)

- New Dataset column `Observation type (failure/run-out/fatigue limit)`, read by the importer
  (`imp.schema.parse_observation_type`, stored as `meta["observation_type"]`); a type that
  contradicts Run out, or unknown text, stops the import.
- All 17 supplied workbooks carry it: run-outs of the S-N workbooks are `run-out` (126 rows), the
  rows of the 8 limit-only workbooks `fatigue limit` (147). `tools/add_observation_type.py` edits
  only the sheet and table XML; FABEST and Delahay are rebuilt by their scripts. Overlays are
  re-pinned to the new SHA-256 (`MANIFEST.json`, `PROVENANCE.md`).
- `imp.classification` takes declared types as given (blank: the rule before d113). A life
  workbook with declared fatigue limits is mixed: its limit states are also judged by the fatigue
  index (`FatigueLimit/`, own run report); `limit_state_runouts` selects them for the fatigue
  index path, so declared interrupted tests of a limit workbook are no longer taken as limits.
- No number changes: fingerprints, classification and fatigue-index outputs identical on all 17
  workbooks.

### Load error: shared scaled evaluations, no domain edge taken as a root (#78)

- `cli/run.py::load_errors` evaluates each scaled load once per case and reuses it across the
  variants and the master-curve and native-life searches (`_CachedModel`, keyed on the input
  tensors; skip errors are cached too). On 42CrMo4_Reis: SWT 4.2 s -> 1.1 s, FS 11.0 s -> 6.5 s.
- `eval/endurance.py::load_factor` returns NaN when the bisection closes on a sign change against
  +inf, i.e. on the edge of the criterion's domain (`scaled_damage` maps a LoadingScopeError to
  +inf): DP never reaches the target there. A finite jump across the target (a plane switch) is
  still a crossing and is kept. The native-life search treats a NaN life as a failed search.
- Numbers: `IMPLEMENTATION_REVISION` 2026-11-12. Only ZHANGYAO on Ti6Al4V_Berto changes: 9
  torsion tests (Torsion_Plain_12, 18-25) had a finite `Error_load` at the scope edge of the
  model (DP(s*) = 0.839 against the target 1.527 for Plain_18) and are now NaN. All other run and
  fatigue-index outputs are unchanged; ranking and batch comparison regenerated.

### One selector of the R = -1 uniaxial design tests (#77)

- `imp.experimental_data.is_design_test` (and `design_test_ratio`) is the single definition of a
  fully reversed uniaxial tension or pure torsion test: R = min/max of the loaded component over
  both load steps, whatever their order, within `DESIGN_RTOL_R` = 1e-3 of -1, every other
  component of both steps below 1e-3 of the peak. The master curve, the Basquin constants and the
  plane-constant calibration all use it; the calibration starts from the same 1e-3 (was 0.02)
  before its recorded relaxation (d107). Membership is unchanged on all 17 workbooks.

### Regression snapshot covers the pipeline (#73)

- `tools/regression_snapshot.py run` adds a pipeline stage: `multifat-pipeline --models all` on
  every workbook (calibrated FS/FIN/LI/MGSE/MSWT constants, the fatigue-index path of the
  fatigue-limit workbooks, sidecars and run reports) and one cross-batch transfer (Reis designed
  on Santus). The script-stage materials come from the workbooks' classification instead of a
  hand-written list. Volatile pipeline fields (dates, calibration ids, git state, timings) are
  left out of the hashes. `--skip-pipeline` runs the script stage only.

### Test: cross-batch transfer with a native-life model (#80)

- `test_transfer_native_life_model`: MOREL in the 42CrMo4_Santus -> 42CrMo4_Reis transfer keeps
  the design batch's own native lives on its design rows, and its 25 assessed rows have finite
  `Error_load`.

### Hand-entered constants accepted by `--constants` (#60)

- New constants kind `hand_entered_constants` (`cli/pipeline.py::load_constants`), as chosen by
  the TVEIT2025 owner: constants calibrated outside the pipeline, with `material_id` and the
  `provenance` of the entry, without the revision and fingerprint check of pipeline sidecars.
  The run sidecar labels them "hand-entered constants" and copies the provenance.
- The four `Calibrations/*/TVEIT2025_stage1.json` files use it; their constants, including
  `Sigm1_TVEIT`, `Taum1_TVEIT` and the numerical settings, and their provenance text are unchanged.
  The documented `multifat-pipeline ... --constants` command works again (it failed since #15).

### One visual system for the homepage

- Every README figure uses `examples/style.py`: `apply_homepage_style` (Inter, recessive axes),
  the validated colour-blind safe palette `SERIES`, one blue ramp for magnitudes (`blue_cmap`) and
  one title layout (`homepage_title`). Redrawn: the critical-plane sphere, the fatigue-limit
  envelopes (direct labels) and the load-error ranking.
- Removed: the polar damage landscape (repeated the sphere), the life-prediction correlation and
  the life-error benchmark (superseded by the load-equivalent ranking, #37), with
  `build_predictions_csv.py` and its CSV. Adding a method now refreshes the ranking only.

### Batch names after their author or project (decision d112)

- `42CrMo4_QT` is now `42CrMo4_Santus` (file `6_42CrMo4_Santus.xlsx`, `Calibrations/42CrMo4_Santus/`)
  and `Ti6Al4V` is now `Ti6Al4V_Meneghetti` (`1_Ti6Al4V_Meneghetti.xlsx`): every 42CrMo4 batch is
  +QT, so the heat treatment did not distinguish them. Workbook bytes and numbers are unchanged.
- The TVEIT2025 stage-1 constants move to `Calibrations/42CrMo4_Santus/` and
  `Calibrations/Ti6Al4V_Meneghetti/`; only their `workbook` path changes.
- The homepage ranking figure leaves out a material with fewer than 2 scored methods
  (Ti6Al4V_Meneghetti, only MOREL, #74) and names it in a note.

### Al7075_T6 material values (decision d111)

- An overlay sets the Pedranz 2025 fatigue limits at their own basis, 143.8 / 100.9 MPa at N = 5e6,
  R = -1 (Appendix B Table 5), instead of the unsourced 150 / 100 MPa, and adds the measured plastic
  strain-life terms of Benedetti 2021 Table 4 (eps'_f = 0.131, c = -0.674) in place of the steel
  uniform-material-law estimates. sigma'_f and b stay the Basquin fit of the HCF tests. Only
  Al7075_T6 outputs change; the homepage ranking is regenerated for it.

### Fatigue-limit datasets in the method ranking (decision d110, #37)

- The homepage ranking adds the fatigue-limit datasets (25CrMo4, 30NCD16, 34Cr4, 982FA, S65A,
  Ti6Al4V_Delahay) as a separate block of columns with its own mean rank. A limit state has the
  load-equivalent error $e_L = \ln(1 + FI_{load})$, exact for any damage parameter; the pure
  R = -1 axial limits are the reference and are left out. CX18 and GGG60 are left out until their
  data errors are resolved (#24).
- `load_error_ranking.csv` appends the column `Assessment` (fatigue_life / fatigue_limit); its
  life rows are unchanged. `dataset_fingerprint(block, runouts=True)` (opt-in) also hashes the
  run-outs (rounded to 1e-6 MPa, so that the last bits of cos/sin on Windows and Linux agree), so
  the manifest tracks the limit states of a limit dataset; the default is unchanged.

### Ti6Al4V_Delahay fatigue-limit dataset (decision d109)

- New workbook `17_Ti6Al4V_Delahay.xlsx`: the Ti-6Al-4V endurance limits of Delahay & Palin-Luc
  2006 (Tables 1-2), built by `tools/build_delahay_ti6al4v.py` from a transcription kept in
  `Experimental_Data/sources/Delahay2006/`. Sigm1 = 652 MPa (plane bending), Taum1 = 411 MPa;
  assessed points: plane bending + torsion at 0 and 90 deg. Tension (583) and rotating bending
  (602 MPa) are recorded as information. Existing outputs unchanged (regression snapshot).

### Fatigue-limit-only runs no longer crash in the comparison

- `multifat-pipeline` on fatigue-limit workbooks only (e.g. `2_CX18.xlsx`) ran every model and then
  stopped with `FileNotFoundError: No CSV files` in `multifat-compare`, which reads life results
  only. `compare` now reports that the folder holds only fatigue-limit results (their summary is
  already in `FatigueLimit/`) and exits cleanly; a folder without any result is still an error.
  Runs with life results are unaffected.

### Ti6Al4V material values (decision d108, #24)

- `1_Ti6Al4V.xlsx`: an overlay keeps only the Meneghetti 2018 properties (Su 980, Sy 915, Sigm1
  500, Taum1 422.5 MPa) instead of the average with the Delahay 2006 row (1035, 947.5, 576,
  416.75 MPa), a different batch. The workbook is unchanged.
- The same overlay removes the 4 Delahay rows (decision d109): Table 2 of Delahay & Palin-Luc 2006
  shows they are endurance limits of another batch, not run-outs, with plane bending entered as
  axial and the torsion limit as 441 instead of 411 MPa. New overlay operation `drop_case` clears
  a Dataset row after checking its authors, specimen and values. Failures and their errors are
  unchanged; GOLOS, no longer stopped by the 90° row, runs in the Berto -> Ti6Al4V transfer.
- Only Ti6Al4V outputs change (regression snapshot against main). MOREL, the one method with its
  own life law, now gives finite errors there (RMS of e_L 0.19 on 4 tests) and a Ti6Al4V cell in the
  homepage ranking; the other methods stay skipped (master curve not identifiable, #74).
- `plot_load_error_ranking`: a material on which fewer than 2 methods are scored is shown but left
  out of the ranks, so MOREL does not gain a free rank 1 on Ti6Al4V (its mean rank stays 26.3).

### Material validity ranges skip, not fail

- New typed skip error `MaterialScopeError` (`models/errors.py`): material constants outside a
  criterion's validity range. RFB raises it outside 0.577 < t_-1/s_-1 < 0.866 (Robert et al.,
  Appendix); the pipeline reports the model as skipped with eligibility `out_of_scope` instead of
  failed (AISI316L, t_-1/s_-1 = 0.517). It is not a `LoadingScopeError`, so load-scaling searches
  do not mistake it for an out-of-scope trial load. No number changes.

### Batched critical-plane evaluation (issue #69)

- Replaced `np.frompyfunc(lambda x: x**2, 1, 1)` in `multifat.models.utils.batch.pypow2`
  with an elementwise Python list comprehension over `np.asarray(x, dtype=float).ravel().tolist()`
  followed by `.reshape(arr.shape)`. This preserves bit-for-bit parity with CPython's
  `float ** 2` (avoiding the 1-ULP MSVC `pow(x, 2.0)` discrepancy) while eliminating NumPy
  `object`-dtype ufunc state implicated in intermittent `python312.dll` crashes (`0xc0000409`)
  during repeated pipeline runs on Windows Python 3.12 (`numpy==2.5.4`).
- Added direct unit and edge-case tests (`test_pypow2_matches_cpython_float_pow`,
  `test_pypow2_overflow_matches_cpython`) in `tests/test_batched_planes.py`.

### Calibration outside R = −1 flagged (issue #77)

- The plane-constant calibration may still relax its R = −1 tolerance or accept any stress ratio
  when a material has no fully reversed tests (accepted by the maintainer, decision d107), but the
  record now stores the tolerance used and whether any R was accepted (`selection.*_rtol_R`,
  `*_any_R`, `r_fallback`), and every result that adopted such a constant is reported as
  `ok (R fallback)` with a note. No current workbook needs the fallback: no output changes.

### Data

- `14_42CrMo4_Reis.xlsx`: R = −1 fatigue limits Sigm1 = 452 MPa and Taum1 = 372 MPa, the authors'
  trend-line values at 1e6 cycles (Anes et al., IJF 62 (2014), Table 5), added by an overlay
  (decision d100, #71). The 22 limit-based methods and WANG2022 now run on this batch.
- `42CrMo4_QT` (Santus 2018 + Benedetti 2022) is one batch (decision d101).

### Cross-batch comparison (#71)

- `multifat-pipeline --design-material A` predicts the failures of another batch from batch A: the
  Summary properties, calibrated constants, design set and master curve all come from A; outputs
  are `<material>_from_<A>` (decision d102). Extension point: `evaluate_model(..., design_data=...)`,
  `run_workbook(..., design_workbook=...)`, `run.load_transfer`. Default runs are unchanged.
- `examples/batch_comparison.py`: standalone pipeline of own-batch and cross-batch statistics for
  the material families of `examples/data/material_families.json` (bias with bootstrap intervals,
  Kruskal–Wallis with Benjamini–Hochberg, two-way decomposition into batch offsets and method ×
  batch interaction, transfer penalty, method-free S-N reference). Results in
  `examples/data/batch_comparison/`, figures `docs/images/batch_comparison_*.png`, guide
  `docs/batch-comparison.md`. `tests/test_batch_comparison.py` fails when a method or a family
  workbook changes until the results are regenerated.
- Statistical design revised (decision d104). New method-free part (`examples/batch_reference.py`):
  Basquin strength at 1e5 cycles per batch, path and series; the scatter between the same-steel
  series of the FABEST database as the yardstick of a batch difference; the relative strengths
  that predict the batch offsets and the transfer shifts without any method. The batch offsets
  get joint-bootstrap intervals (same specimens for all methods); the Kruskal–Wallis and
  Benjamini–Hochberg tests are removed in favour of effect sizes; run-outs enter a censored
  sensitivity check; a load factor covering 95 % of the tests compares own batch and transfer.
  Six new figures per family (`examples/batch_plots.py`).
- Material batch, test series and dataset are kept apart (decision d105): the FABEST database is one
  material batch (maintainer), so its between-series scatter is the test-series part, and
  `variance_components` estimates the material-batch part of the spread between datasets
  (`reference_components.csv`). Santus 2018 and Benedetti 2022 stay one batch with two series,
  decided on the data. Decision d104 is active.
- A transfer accepts `--constants` only for a published-constants record of the design batch
  (decision d103): SSF is assessed on FABEST and QT with the surface published for Reis.
- Second family Ti6Al4V (`Ti6Al4V`, `Ti6Al4V_Berto`). Batches are author groups, Santus 2018 and
  Benedetti 2022 being one (decision d101); `1_Ti6Al4V` mixes the Meneghetti and Delahay
  properties in its Summary (open, #24).
- `eval/load_error.py`: the secant search in ln s stops beyond |ln s| = 10 and falls back to the
  bisection, as documented, instead of overflowing (SSF extrapolated on FABEST). Regression
  snapshot: no existing output changes.


- Homepage data provenance manifests (issue #80): `examples/homepage_data.py` writes paired
  `.manifest.json` files alongside `examples/data/load_error_ranking.csv` and
  `examples/data/42CrMo4_QT_predictions.csv` recording their canonical LF SHA-256,
  `IMPLEMENTATION_REVISION`, registered models and post-overlay `dataset_fingerprint`s.
  `examples/build_load_error_ranking.py --check`, `examples/build_predictions_csv.py --check` and
  `tests/test_gallery_figures.py` fail when the committed CSVs go stale;
  `tests/test_dataset_overlays.py` now also covers invalid overlay operations, unsupported overlay
  schema versions and duplicate `base_sha256` aborts.
- Workbook integrity checks (issue #75): `Experimental_Data/MANIFEST.json` records the SHA-256 and
  post-overlay `dataset_fingerprint` of every supplied workbook (`tests/test_dataset_overlays.py`
  enforces that supplied bytes do not drift without updating the manifest), and loading a file
  whose hash does not match an overlay while its filename matches one now emits a `RuntimeWarning`
  and appends a warning to `ValidationReport.warnings`.
- New workbook `16_42CrMo4_FABEST.xlsx`, a third 42CrMo4 +QT batch, from the FABEST database
  (FAB_v20261014, doi:10.5281/zenodo.23359038, CC BY 4.0; #65). It holds the PBS series of plain hollow
  specimens: 79 failures and 10 run-outs, in axial, torsion, and in-phase and 90° out-of-phase
  tension–torsion loading, all at R = −1.
  - The original file is kept unchanged in `Experimental_Data/sources/FABEST/`. The workbook is
    generated by `tools/convert_fabest.py` and checked by `tests/test_fabest_conversion.py`.
  - Sigm1 and Taum1 are the authors' Basquin strengths at 1e6 cycles. They lie below the observed
    run-out levels; see `Experimental_Data/PROVENANCE.md`.

**Licence changed to MIT** (copyright Andrea Chiocca, Michele Sgamma and Tommaso Grossi), superseding
GPL-3.0-or-later. Incorporated third-party code is limited to permissive licences; the earlier rule that
GPL code is not ported still applies.

Branch `dev/tg-developments`, based on `main` at `e70a8af` (15 October 2026). Every method was
checked against its original publication; changes to numerical outputs are listed under
**Changed numerically** and were checked with the regression snapshot before each commit.

### Packaging and tooling

- Method documentation and docstrings (issue #32): added the remaining 16 method guides under
  `docs/method-guides/` (`LI`, `LIU2021`, `MGSE_YU`, `MGSE_ZHU`, `MSWT`, `MDC`, `MMCCM`,
  `SWTD`, `OTT`, `EK`, `RFB`, `VHN`, `PBP`, `COIN`, `MMMH`, `TVEIT2025`), bringing guide
  coverage to all 46 registered criteria; completed module docstrings (opening on line 1 with
  governing equations and `Repository choices:`) across all 19 remaining method modules and
  emptied `PENDING_DOCSTRINGS`; tightened `tests/test_method_guides.py::test_guides_exist` to
  require a guide for every registered criterion; updated `docs/method-guides/README.md`,
  `README.md`, and `docs/user-guide.md`.
- Implementation revision bump and policy (issue #76): `IMPLEMENTATION_REVISION` in
  `src/multifat/models/references.py` bumped to `"2026-11-11"` after merging PR #70
  (`origin/main` at `c09dc4c`); `CONTRIBUTING.md`, `docs/method-audit.md` and
  `.github/pull_request_template.md` state the bump rule for PRs that can alter numerical outputs.
  Pipeline sidecars calibrated under `"2026-10-07-audit1"` are rejected by
  `multifat.cli.records.check_record_compatibility`, and `load_constants` flags upstream
  `Calibrations/` records carrying an older `implementation_revision` with a provenance warning.
- Documentation and tool CLI checking (issue #79): `tools/check_docs.py` now parses fenced
  `tools/*.py` shell invocations against each script's `build_parser()`; `AGENTS.md`,
  `CONTRIBUTING.md`, `docs/user-guide.md` and `docs/adding-methods.md` are aligned on overlays,
  the `imp/schema.py` rule, `scripts/` compatibility shims, and `tools/regression_snapshot.py run <dir>`.
- The package is `multifat` in a `src/` layout with a `pyproject.toml` (hatchling) and a `uv.lock`;
  the requirements files are gone. Console scripts `multifat-calibrate`, `multifat-run`,
  `multifat-compare`, `multifat-pipeline`; the standalone scripts in `scripts/` are thin shims for them.
- Comprehensive User Guide (`docs/user-guide.md`) covering pipeline execution,
  single-step CLI tools, experimental workbook preparation, Python API recipes,
  and troubleshooting.
- All Python code is black-formatted (line length 100); the one-off reformat is listed in
  `.git-blame-ignore-revs`. CI runs ruff, black, the doc generators, the documentation check and the
  tests on Python 3.10-3.13 (and Windows on 3.12).
- `tools/regression_snapshot.py` runs every model on every workbook plus calibration and comparison
  and hashes all outputs, to prove a change output-identical.
- Batched plane evaluation and a per-calibration plane cache (bit-identical outputs): a full
  snapshot dropped from about 21 to about 4 minutes before the calibration moved to all planes.

### Master-curve identifiability and scatter (issue #74)

- **Changed numerically.** A master curve is used only if it is identifiable: at least 5 design
  tests, at least one decade of life, and a negative slope distinct from zero (|b| > 2 standard
  errors), the thresholds of the Basquin fits. Otherwise the model is skipped on that dataset with
  the reason (maintainer decision, d106). On the current data this skips the master-curve methods
  on `Ti6Al4V` (6 axial design tests over 0.39 decades; MWHS: 4 torsion tests).
- **Changed numerically.** The scatter of the master curve (`STD_ln_DP`, `STD`) is estimated with
  n − 2 degrees of freedom (two fitted parameters), so its square is unbiased. Median results
  (P_s = 0.5) and the calibrated constants are unchanged; curves with P_s > 0.5 shift slightly.
- New output `<material>_master_curves.csv` per pipeline run: n, decades, A, b with its standard
  error, R² in log space, scatter, identifiable and reason, for every fitted curve.
- Consequences: the homepage load-error ranking loses its Ti6Al4V column (no method keeps three
  assessed errors there); the cross-batch comparison of the Ti6Al4V family keeps only the
  transfer from Ti6Al4V_Berto. Three synthetic test fixtures that had four or fewer design tests
  were extended or their assertions updated.

### Survival probability (issue #16)

- `--survival-probability P_s` on `multifat-run` and `multifat-pipeline` (0.5 <= P_s < 1, default 0.5
  = median, unchanged results): the master curve is shifted by z = Phi^-1(P_s) standard deviations of
  the ln DP residuals (`multifat.eval.fitting.stdnum_for_survival`). Every result CSV ends with a
  `Survival_probability` column, the sidecar records it, and `multifat-compare` refuses to mix
  probabilities. Native life laws and fatigue-limit indices are not shifted.

### Load-equivalent error (issue #37)

- Result CSVs gain `Error_load` = ln(1/s*), the log of the factor on the whole load at which the
  method predicts the observed life (master curve or native law; failures only; not TVEIT2025).
  It does not depend on the S-N slope (Error_ln = k d e_L) nor on the degree d of the DP
  (Error_log10_dp = -d e_L / ln 10), so it is the recommended cross-method ranking metric.
  `multifat-compare` adds its mean, SD and the fractions within +-10 % and +-20 % of the load
  (`Within_10pct`, `Within_20pct` columns of the statistics tables). Existing columns unchanged.

### Calibration history (issue #15)

- Pipeline calibrations get a persistent `calibration_id` (UTC time + content signature) and are
  archived without overwriting (`--calibration-archive`, default `<results-dir>/calibrations`);
  equivalent reruns keep their id. `--calibration-id` reuses a record; records lacking provenance
  or made with another implementation revision are rejected, other data or procedure versions are
  flagged. The run report gains a `calibration_id` column; sidecars record
  `calibration_procedure` (2). Numerical results unchanged.

### Dataset workflow classification (issue #18)

- `multifat.imp.classification.classify_dataset` assigns each workbook a workflow
  (`fatigue_life`, `fatigue_limit`, `ambiguous`, `out_of_scope`) and each plain observation a role
  (`design`, `assessment`, `censored`, `limit_state`); the pipeline writes
  `<MATERIAL>_classification.json` and routes on it, with the previous routing preserved.
  Run-outs below the reference life of the Summary limits are reported (25CrMo4: all 11).

### README figures

- Header figure: Findley damage on the sphere of plane normals with the critical plane
  (`examples/plot_damage_sphere.py`); gallery item 5: all methods ranked by the RMS of the
  load-equivalent error across the six life datasets (`examples/plot_load_error_ranking.py`, data in
  `examples/data/load_error_ranking.csv`). Both are produced by `examples/generate_gallery.py`.

### FatPy alignment

multifat is expected to merge into FatPy (COST CA23109). This change documents the correspondence
of names and records one distinction FatPy makes; no multifat name was renamed.

- `docs/notation.md` has a *FatPy nomenclature* table. It lists the FatPy key (nomenclature v01)
  of each multifat name that has one, and the units that differ (FatPy gives E in GPa; multifat
  stays in MPa).
- A new optional Summary column, `Axial fatigue limit loading`, records for each source row how
  the axial R = -1 fatigue limit was obtained: `push-pull`, `rotating bending`, `plane bending`
  or `bending`.
  - The value is kept with the `Sigm1` source and reported in the validation summary.
  - It maps the limit onto `fat_lim_ten_m1` or `fat_lim_ben_m1`.
  - Limit values, and so numerical outputs, are unchanged.
- `docs/research/fatpy-nomenclature-v01-feedback.md` lists, for the FatPy coordinator:
  - keys that are missing (run-out flag, torsional ultimate strength, per-component loads,
    method constants);
  - ambiguities and unit inconsistencies;
  - spelling issues.

### Master-curve models (issue #38)

- `--sn-model {basquin,kohout_vechet,stromeyer}` on `multifat-run` and `multifat-pipeline`
  (`multifat.eval.sn_curves`): Basquin (default, unchanged), Kohout & Vechet (2001) Eq. (H2) and
  Stromeyer, fitted on ln DP with `scipy.optimize.least_squares` (SciPy is now a dependency,
  imported only for the nonlinear models). Unidentifiable fits are rejected with their reason;
  infinite lives below the fatigue limit are kept, counted in the run report and in the new
  `n_nonfinite` statistics column. CSVs end with `SN_model`. Plane-constant calibration and native
  life laws are unchanged. Known problems are documented in `docs/method-audit.md`.

### Methods

- New criterion `ZHANGYAO`, Zhang & Yao (2011), catalogue N06 (issue #13).
  - Formulation: a nonlinear combination of the shear amplitude $C_a$ and the Goodman-corrected
    normal stress on the maximum-$C_a$ plane.
  - Verification: 50 of the 66 fatigue limits of Tables 3–8 within 0.3 points.
  - Finding: the maximum-$C_a$ plane rule is ill-conditioned for this criterion, so results in
    out-of-phase loading depend on the tie rule.
  - Guide: `docs/method-guides/ZHANGYAO.md`.
- New criterion `WANG2022`, Wang et al. (2022), catalogue N05 (issue #13).
  - Formulation: the equivalent shear stress on the maximum-$C_a$ plane, with life from the
    torsional S–N curve; its master curve is fitted on the torsion failures.
  - Verification: the predicted lives of Tables 4–5 are reproduced within 1 %.
  - Finding: Eq. (18) as printed makes a tensile mean stress beneficial. By maintainer decision,
    loading with a mean normal stress on the critical plane is out of scope.
  - Guide: `docs/method-guides/WANG2022.md`.
- The native-life load error (`Error_load`) is NaN, instead of a failed model, when a trial load
  leaves the model's declared scope (`LoadingScopeError`), as the master-curve search already
  does.
- New criterion `LIUMAHADEVAN`, Liu & Mahadevan (2005), catalogue N16 (issue #13).
  - Formulation: normal, shear and hydrostatic amplitudes on the plane at a material angle from
    the maximum normal stress amplitude plane, with a normal mean-stress factor. Its native life
    law (Eq. 23) uses the tension and torsion Basquin curves.
  - Maintainer decisions: the worst plane on the $\alpha$-cone in 3D, the rss shear amplitude and
    the native life.
  - Verification: Eq. (19) and both reversed limits are exact (for $t_{-1}/f_{-1} \ge 0.5$). The
    published Figs. 2–5 are reproduced only in part (strict xfails, unresolved).
  - Guide: `docs/method-guides/LIUMAHADEVAN.md`.
- New criterion `SSF`, stress scale factor of Anes et al. (2014), catalogue N24 (issue #13).
  - Formulation: an equivalent shear stress with an ssf surface depending on the axial stress
    and the amplitude ratio, evaluated along the cycle; life on the pure-shear S–N curve.
  - By maintainer decision only the published 42CrMo4 surface is available
    (`Calibrations/42CrMo4_Reis/SSF_published.json`). On that workbook the results are in-sample.
  - Verification: Tables 5–8 reproduced from the Table 4 trend lines (two misprinted cells), and
    Eq. (31) within 0.06 of Table 9.
  - Guide: `docs/method-guides/SSF.md`.
- `--constants` accepts a published-constants record (`published_constants`, `material_id`,
  `reference`): constants printed in a paper, not a calibration of the data.
- New criterion `MMP`, Manson–McKnight–Papuga, Papuga & Fojtík (2017), catalogue G05 (issue #13).
  - Formulation: a von Mises-type equivalent amplitude and mean with weighted shear terms,
    combined by a Walker correction with the mean-stress factor $\beta_P$ (variant MMP_WP).
  - Verification: all 24 load cases of Table 11 within 0.1 point of the fatigue index error.
  - By maintainer decision, missing $f_0$ and $t_0$ are estimated (as for `PCR` and `MMMH`).
  - Guide: `docs/method-guides/MMP.md`.
- New criterion `KUJAWSKI2021`, Kujawski (2021) mSWT, catalogue N26 (issue #13).
  - Formulation: the SWT product with the maximum stress augmented by $|\sigma_m|/3$ when the
    mean stress is negative (Eq. 11b), applied on the planes of the SWT critical-plane extension
    (maintainer decision; the paper is uniaxial).
  - Verification: closed forms of Eqs. (11a)–(13b) and identity with `SWT` without a
    compressive normal mean; the paper validates on figures only.
  - Guide: `docs/method-guides/KUJAWSKI2021.md`.
- New criterion `GAM`, Gonçalves, Araújo & Mamiya (2005), catalogue N08 (issue #13). It is an
  endurance criterion for hard metals.
  - Formulation: the prismatic-hull deviatoric amplitude plus the maximum principal stress over the
    cycle, Eq. (11); the constants come from $f_{-1}$ and $t_{-1}$.
  - Verification: Tables 1–4 reproduce their 14 in-phase rows. The 27 out-of-phase rows of the
    tables used a phase-free peak principal stress, contrary to Eq. (11); they are strict xfails.
  - Guide: `docs/method-guides/GAM.md`.
- New criterion `NINIC`, the Ninić (2006) Damage Function, catalogue N11 (issue #13).
  - Formulation: a quadratic combination of $\tau_a/T_e$ and $\sigma_{n,\max}/S_e$, maximised
    over all planes in three dimensions.
  - Verification against the tables: 29 of 32 full-search rows of Tables 2–4 within 0.08 points,
    and Table 5 within 0.6 points.
  - Printed planes that are not the maximum of $D$, and load misprints, are strict xfails.
  - Finding: Eq. (4) for $k$ calibrates the uniaxial limit only for $S_e/T_e \ge \sqrt2$, and a
    warning is issued below that ratio.
  - Guide: `docs/method-guides/NINIC.md`.
- New criterion `GOLOS`, the Golos, Dębski & Dębski (2014) modified Gough limit curve, catalogue
  N13 (issue #13).
  - Scope: the curve is defined for in-phase, fully reversed bending–torsion only. A new
    `LoadingScopeError` makes the pipeline report such a model as `skipped` (eligibility
    `out_of_scope`) instead of `failed` when a case lies outside the scope.
  - Finding: the text's coefficients for $\kappa = 2$ contradict Eq. (30).
  - Guide: `docs/method-guides/GOLOS.md`.
- New criterion `MWHS`, the Wu, Li, Fang & Song (2017) modified Wu–Hu–Song strain parameter,
  catalogue N28 (issue #13).
  - Formulation: Eq. (19), with the longest-chord shear range.
  - Verification: Fig. 14 is reproduced, except the normal term at 90°.
  - Calibration: as in the source, the master curve is fitted on the fully reversed pure torsion
    failures. A model now declares this with `design_loading = "torsion"`, which is read through
    `extract_pure_torsion`. A workbook without such tests skips the model (`DesignSetError`).
  - Guide: `docs/method-guides/MWHS.md`.
- New criterion `LIU2025`, Liu, Xie & Hao (2025), a new catalogue record N58 from a
  maintainer-supplied paper (issue #13).
  - Formulation: a von Mises-like combination of the corrected plane amplitudes on the maximum
    shear plane, with Gerber and Serensen mean-stress corrections.
  - Verification: Fig. 3 is reproduced except two points.
  - Guide: `docs/method-guides/LIU2025.md`.
- Paper-reproduction tests for every registered method in `tests/methods/`, with the transcribed
  tables in `tests/data/papers/` (DOI, table, page). Confirmed paper misprints are strict `xfail`s
  with the reason (e.g. Papadopoulos 1997 Table 1 rows computed with sin δ, P1).
- Restored: SWTD, OTT, DANGVAN (Böhme et al. 2026 Eq. 3 form). ZHU_EDP is a kernel of Eqs. (19)-(23)
  and stays disabled (needs measured F_NP and ΔW0); CARSPA stays disabled.
- New registered criteria: CROSSLAND (N01), SINES (N02), PAPADOPOULOS 2001 (N03), MWCM (N04), VHN
  (N20), PBP (N37), PCR and PI (G01, G02), MDC (G06), MMMH (G07, endurance-limit Haigh line), COIN
  LiWI (G08), VSE (G15), EK (G16), MMCCM (G22), MOREL (G23), RFB (G41).
- Unregistered point kernels with their own tests: Flacelière-Morel (G24), Banvillet
  (G26, including the solid-bar tension-torsion closed form), Papadopoulos-Panoskaltsis (G29),
  Berto-Lazzarin SED (G30), Doudard-Hild-Calloch (N34). The catalogue records what each one needs
  before registration; G35 stays a specification.
- `optional_params` declares every constant a model reads with `params.get`; a missing required
  constant raises `MissingParamsError`. Critical-plane methods declare their out-of-phase shear rule
  (`np_shear_rule`) and can opt out of the `_ext` variant.
- Native life laws: MDC, MMCCM and MOREL export `Nf_native` and `Error_log10_native` next to the
  fitted universal curve.

### Experimental-data pipeline

- `multifat-pipeline <workbook> | --material NAME | --all [--jobs N]`: validates a workbook, calibrates
  the constants it lacks (sidecar JSON, never written into the workbook), runs every model and
  variant, compares, and reports skipped models with the reason. Outputs in `Results/current/`.
- Workbook contract in `src/multifat/imp/schema.py`, documented in the generated
  `docs/workbook-format.md`; contract violations raise `WorkbookError` with a validation report.
- Constants the pipeline supplies: Basquin curves fitted on the R = -1 tension and torsion tests
  (accepted only with at least 5 tests, one decade of life and a slope beyond 2 standard errors),
  MMCCM plastic terms from the uniform material law (flagged), k_MGSE from Zhu 2018 Eq. (12), the
  other plane constants by scatter minimisation on all planes.
- Run-outs are imported as censored observations instead of being dropped: evaluated with the
  same constants, exported with `Flag = 2` (recorded cycles as a lower bound, empty errors),
  drawn with arrows in the comparison plots (`--hide-runouts` to omit them), never used in
  calibration, curve fits or statistics. Run-out status also accepts true/false, 1/0 and
  runout/failure.
- The pipeline always writes `<MATERIAL>_calibration.json`, now with a provenance block (workbook
  SHA-256, failure and run-out counts, git commit, versions, numerical settings).
- Robustness: `--constants` takes exactly one workbook; `--material` accepts the full file stem
  and reports an unknown material without a traceback; Excel lock files (`~$*.xlsx`) are ignored;
  a one-workbook run keeps the eligibility matrix of an earlier `--all`; empty specimen or
  frequency cells are empty text; `tools/install_revised_dataset.py --apply` refuses to replace a
  workbook edited after installation unless `--replace-edited`.
- Result CSVs end with the observation identity columns (`Test_status`, `Specimen_id`, `Case_id`,
  `Observation_id`, `Dataset_fingerprint`); `multifat-compare --common-methods` compares methods on
  the failures they all assess; `scripts/3_review_subsets.py` forwards to it.
- `tools/regression_snapshot.py --data-folder DIR` (frozen inputs) and a compare on the base
  snapshot's columns; `multifat-calibrate --data-folder`.
- `Experimental_Data/PROVENANCE.md` traces every workbook row and constant to its source, with the
  open questions for the data owners.
- Notes written into outputs were reworded during the refactor (`constants_source.note` in
  `<MATERIAL>_calibration.json`, the run report, `<MATERIAL>_validation.txt` and the eligibility
  notes): the text changed, the values it describes did not (wording only, checked against the
  pre-refactor baseline).
- Calibration records one status per material, with its reason: `calibrated`, `not_applicable` (no
  plain failures: fatigue-limit or notched-only set), `infeasible` (no finite objective) or `error`.
  It is written to `<MATERIAL>_calibration.json` (`status`, `reason`), printed by
  `multifat-calibrate` as a summary table, and copied to `calibration_status` in the eligibility
  matrix and to `Comparison/calibration_status.csv`. The batch is isolated per material instead of
  stopping at the first failure; an infeasible fit names the block, the parameter and a diagnostic
  measured only when it fails, and the exit code is 1 when any material is `infeasible` or `error`.
  A material whose calibration failed keeps running on the constant's default, and the run report
  says `ok (default constant)` with the block and the default value; models skipped for that
  constant carry the same pointer in their reason. No fitted number, table or figure changes
  (`1_Ti6Al4V` keeps the constants it had); block isolation (`partial`) is deferred.

### Upstream sync (PRs #5-#9, up to ab0e828)

Merged on 19 October 2026 (PRs #5-#7, up to 71eaac3, merge commit a3c1d2e) and 20 October 2026
(PRs #8-#9, up to ab0e828, merge commit b8c84c9); upstream modules, originally developed for the legacy `fatigue/`
layout, were ported to `multifat`, each port naming its upstream commit.

- Data: the revised October 2026 collection (twelve numbered workbooks, `Al7075_T6.xlsx`, archives of
  the previous inputs) and the `Calibrations/` records, as supplied.
- Ported: revised workbook import, material IDs without ordering prefix, relative phases, exclusion of
  rows with unknown loading, identity and stress-ratio checks, optional E and nu (dfb514d, 904c619);
  observation identities, dataset fingerprint and common-subset comparison (`--common-methods`,
  904c619); the eligibility matrix (904c619, as pipeline output); run-out plot details (7a87256); the
  harmonic MCC/MCE shear amplitudes with their tests, unused by the criteria (88229e9, 2a17791); the
  dataset install and audit tools (904c619).
- PRs #8-#9: the isolated harmonic Findley damage evaluation `evaluate_findley` with MCC or MCE
  amplitude on the supplied planes (176ce53, 50dabaf; unregistered, nothing in production imports it),
  and the choice of the Findley shear amplitude in the k_FI calibration (36ff0a1, 0105b69):
  `multifat-calibrate --shear-amplitude-algorithm {method_default, mce_harmonic, rss, mcc_harmonic}`.
  The default and every production result are unchanged; MCC is the `longest_chord` plane rule through
  the established refined search, checked against `evaluate_findley`. The option is not exposed by
  `multifat-pipeline`, and `--constants` rejects a record whose Findley amplitude is not MCE (production
  Findley evaluates MCE). Docs: `docs/shear-amplitude.md`, `docs/calibration-workflow.md`.
- Used as inputs: the `Calibrations/` records through `multifat-pipeline --constants` (dc4fe6d); the
  failure-only `failure_fingerprint` in the calibration sidecar (our hash, not comparable with the
  records' `input_fingerprint`); a record made on a different number of failures is flagged.
- Not ported, by design: the versioned calibration store and record selection, with its method-scoped
  Findley compatibility and fingerprints (the pipeline calibrates in one run and records provenance), the `output_quantity` contract and `direct_life`
  (native lives are exported as `Nf_native`), `tools/migrate_calibrations.py` (the importer reads
  legacy coefficients directly) and `tools/run_tests.py` (pytest). His unittest suites are ported
  where the feature is.

### Fatigue-limit datasets

- A workbook whose plain rows are all limit states (run-outs at the endurance basis) with both R = -1
  fatigue limits is assessed by the fatigue index error of every model (Boehme & Papuga 2024,
  Eqs. 8-9), FI = DP / DP(axial R = -1 limit) - 1, with missing plane constants set so that the axial
  and torsional limits give the same DP (`eval/endurance.py`, outputs in `FatigueLimit/`). FI_load
  = 1/s - 1 (s: load factor at which DP reaches the axial reference) is the stress-scaled index
  for every damage parameter (d063).

### Changed numerically

- **TVEIT2025: end of the block integration** (issue #66, approved by the maintainer). The
  stopping rule compared only two consecutive blocks, so during a slow shakedown it stopped while
  the damage per block was still decaying to zero and extrapolated 9e7 to 2.4e9 cycles instead of
  an infinite life. Now a block without any damage increment is a shakedown (infinite life), and
  a cyclic steady state is accepted from block 400 on only when the back-stress has no net drift
  over 200 blocks and the mean damage per block is stable (`tveit2025.py` docstring). The rule
  reproduces 16 of 16 histories integrated to 50,000 blocks without a stopping rule
  (`tests/test_tveit2025_end_rule.py`). With the stage-1 constants the failures predicted
  infinite go from 3 to 18 (42CrMo4_QT), 2 to 14 (Al7075_T6) and 2 to 26 (Ti6Al4V_Berto);
  Ti6Al4V does not change. Nine Al7075_T6 tests loaded exactly at the model's fatigue limit stay
  finite but move from 4.5e8-1.9e9 to 4.6e9-8.0e9 cycles; all other lives move by at most 1 %.
  The load-equivalent error, computed from per-path curves outside the pipeline, changes by at
  most 0.006. New settings defaults: `min_blocks_TVEIT` 400, `max_blocks_TVEIT` 60000,
  `tol_alpha_TVEIT` 0.01 (net drift over 200 blocks, relative to the path length),
  `tol_damage_TVEIT` 0.001 (mean damage per block between two windows); the four
  `TVEIT2025_stage1.json` files carry them.
- **TVEIT2025: a creeping back-stress is no longer taken for a steady state** (follow-up to #67,
  found by the automated review). A back-stress that swings in every block can creep one way while
  its net drift over 200 blocks stays below 1 % of its path length; the peak of beta then sinks
  until a shakedown. The steady state is now refused while the trend of that peak reaches zero before
  the extrapolated failure and within the remaining blocks. Al7075_T6 stage-1 constants, torsion
  99.996 MPa with a constant out-of-plane shear of 0.895 MPa: 4.6e9 cycles before, infinite life
  (shakedown at block 2446) now. The 16 long histories are reproduced as before.
- **TVEIT2025: within-cycle failure time.** A failure inside a block is placed in the increment
  that crosses D = 1: (i - 1 + f)/(npts - 1) of the cycle, where the code counted one increment
  more ((i + f)/(npts - 1)). Lives shorten by at most 1/(npts - 1) = 0.001 cycle
  (`docs/audit-checkpoint.md`).
- `k_MGSE` no longer depends on the model selection: Zhu (2018) Eq. (12) on the Basquin fits is used
  whenever those fits are identifiable, also in a run that selects MGSE_ZHU without MDC or MMCCM
  (previously scatter minimisation, e.g. 1.239 instead of 0.988 on Al7075_T6). Full runs are unchanged.
- L_LI search interval widened from [0, 5] to [0, 20]: the minimum lay outside it on 42CrMo4_QT
  (L_LI 5.0 -> 7.07, life-scatter STD 1.153 -> 0.942) and 42CrMo4_Reis (5.0 -> 5.36, 0.368 -> 0.341).
  Interior L_LI values move by at most 3e-6 relative (grid refinement on the wider interval).
  A scatter minimum on an upper bound that only limits the search (all but a_MSWT = 1) is no longer
  used: it is recorded as `unidentified`, and models requiring it are skipped with that reason. On
  Ti6Al4V the wider L_LI interval makes the calibration feasible (L_LI = 11.5) but k_FI and k_FS end
  at 3, so FIN and FS stay skipped, now with that reason; LI, MGSE_ZHU and MSWT use the calibrated
  constants instead of their defaults. A model whose master curve gives no life predictions (slope
  `b >= 0` or `k = -1/b > 15`) is noted in the run report.
- Design flags of calibrated plane constants: the tension and torsion tests of the scatter
  objective behind a pipeline-calibrated `k_FI`, `k_FS`, `L_LI`, `k_MGSE` or `a_MSWT`, and the
  Basquin tests behind a `k_MGSE` from Zhu (2018) Eq. (12), are now `Flag = 0` for the models that
  adopt the constant; they were `Flag = 1` assessment rows, so the assessment statistics included
  the tests that fitted the constant. On the five life datasets 410 rows of FIN, FS, LI, MGSE_ZHU,
  MSWT and their `_ext` change flag; constants, master curves and predicted lives are unchanged
  (differences at most 6.7e-14 relative, from the design-row code path). `--constants` with a
  sidecar of this pipeline keeps its design tests when the `failure_fingerprint` matches.
- F1: result CSVs are written as UTF-8 (only text fields changed; SWT results now enter the
  comparison on Windows).
- F2: the critical-plane search refines the grid maximum (grid + analytical seeds + local
  refinement); the 522-plane grid under-estimated the continuous maximum by up to 1.2 %. F10:
  refinement also starts from tied maxima in different basins (workbook impact at most 2.2e-8).
- F3: the out-of-phase shear range follows each source: longest chord for FS, MKBM, LI and LIU2021,
  RSS for the other critical-plane methods.
- F5-F7, F9: SWTD, OTT and DANGVAN equations corrected before they were re-enabled.
- U2: MSWT uses the printed (1 - a)/2 shear weight of Jiang (2000); the workbook a_MSWT values were
  converted (a' = a/(2 - a)), so lives are unchanged and the damage parameter is rescaled.
- F4: calibration evaluates all planes and returns the smallest point of a flat minimum.
- Workbooks: 42CrMo4 `Plain_multiaxial_9` is proportional (phase 90° -> 0°, N1); AISI316L
  non-proportional rows get their 90° phase (N6), E = 193 GPa and ν = 0.3 from the literature (N7),
  and the 8 Janssens Table 4 points missing from the workbook (N8). Each changed cell carries an
  Excel comment with the source and the old value.
- New data (upstream 71eaac3): every result changes with the revised collection. 42CrMo4_QT and
  Ti6Al4V are the one-batch revised files; the Anes/Reis and Berto rows are the new workbooks
  `14_42CrMo4_Reis.xlsx` (with the printed strain-life constants) and `15_Ti6Al4V_Berto.xlsx` (d055).
  In the revised AISI316L, N6 and N8 are re-applied, Janssens E/nu restored (Braun's 219 GPa / 0.29 are
  XRD constants) and the Braun R = 0 mean stresses restored (d056); nu = 0.3 for four steels (d059).
  The plane constants of every material are calibrated, as the revised workbooks carry none.
- Revised-format import: an empty or `Null` stress excludes the row instead of reading 0; phases are
  relative to the axial phase (34Cr4 has nonzero axial phases); AISI316L's Sigm1, averaged over two
  stress ratios, is not used.
- Without E or nu, the critical-plane evaluator and the plane seeds skip the undefined strain
  tensors (finite inputs unchanged).

### Notation and documentation

- `docs/notation.md` is the symbol table behind the package vocabulary: paper symbol, Python name,
  workbook key and unit for every concept, plus the rule that a function implementing one published
  equation uses that equation's symbols. The code follows it end to end: `f_m1`/`t_m1` for the
  reversed limits, `sigma_z0`/`sigma_z1` for the load-step normal stress, `d_gam`/`d_tau` for the
  ranges, `DP`/`metric` for the plane scalars, `R_list`/`normals` for the plane frames,
  `case`/`refs`/`bounds` in the evaluation layer and `workbook`/`material`/`data`/`params` in the
  command line. The closing pass restored the coefficient sets of the 2026 Böhme-Papuga paper
  (`a_BP`..`d_BP`, `a_CAIM`..`d_CAIM`, `kappa_sq`), the MMMH surface stresses (`sigma_x`, `sigma_y`,
  `tau_xy`, `sigma_z`, `k_a`, `M`) and the evaluation-layer names (`block`, `case`, `key`, `target`,
  `df`, `material`, `amplitude`, `Nf`, `DP`, `eps_e`/`eps_p`). No algorithm, constant, public method
  name, registry key, CLI option or workbook key changed: `regression_snapshot.py` compares
  identical to the previous commit (1708 output files and 410 logs, no differences) and the
  pre-existing tests still pass.
- Python source, its docstrings and the tables generated from it spell the name `Boehme`;
  hand-written documentation keeps the author's spelling `Böhme`, including the method catalogue
  (`docs/catalogue/methods.json`) and everything generated from it, `docs/bibliography.md` and the
  README. The guard rejects `Böhme` in `src/**/*.py` and in `docs/workbook-format.md`.
- `tests/test_notation.py` guards the vocabulary: retired synonyms (read from `docs/notation.md`),
  `_old = _new` private aliases, module docstrings that do not open on their first line, string
  statements left outside the docstring slot, an override that renames the parameters of the method
  it overrides (a coefficient set may only take the subscript of its own criterion, `a_BP`..`d_BP`
  to `a_CAIM`..`d_CAIM`), the banned docstring phrases of `docs/adding-methods.md`, and every
  parameter of a public function or method in `eval/` and `models/` that appears neither in the
  table nor in `tests/data/notation_allow.txt`. That file holds exactly the names outside the
  table, one group per reason, so an entry whose name disappears fails like a missing one. The
  docstring rules themselves are in `docs/adding-methods.md`.
- The review found one logic difference, `pipeline.py` treating an empty constant mapping as absent
  and reloading it from the workbook: it is kept (`is None`) and pinned by
  `test_empty_data_mapping_is_used_as_given`.
- Method guides (issue #32), `docs/method-guides/`: a template, an index with the conventions
  shared by every method (load history, plane search and shear-amplitude rules, master curve,
  calibration of plane constants) and one guide each for FIN, FS, SWT, MATAKE, MCD, MWCM,
  CROSSLAND, SINES, DANGVAN and PAPADOPOULOS: original formulation in LaTeX, repository
  implementation and its departures, constants, original versus repository calibration, an
  executable worked example checked against a hand calculation, comparison capabilities and
  verification evidence. `tests/test_method_guides.py` runs every example and compares its printed
  output, checks the section order, the registry names on the command lines and the index, and
  ratchets the `Repository choices:` paragraph of the method-module docstrings (25 modules still
  pending). Module docstrings added or completed for FIN, FS, SWT, MATAKE, MCD, MWCM, CROSSLAND,
  SINES, DANGVAN and PAPADOPOULOS. Documentation only: `regression_snapshot.py` compares identical
  to `5cd5124`.
- Method guides for the strain and energy critical-plane criteria MKBM, VSE, GSE and GSA (issue
  #32). Each guide gives the torsion-to-tension limit ratio that the criterion implies with elastic
  strains, in closed form: about 0.68 for MKBM at $f/S_y = 0.75$, 0.88 or 0.62 for VSE (mode I or
  II), 0.63 for GSE and 0.55 for GSA. The guides also document two properties of the criteria:
  - VSE has no mean-stress effect;
  - GSA lives do not depend on $\sigma_f'$ when strains are elastic and $\tau_f'$ takes its
    default.

  The MKBM, VSE and GSA module docstrings now carry a `Repository choices:` paragraph, so 19
  modules remain pending. Documentation only, no numerical change.

### Repository

- Licensed under the GNU General Public License v3.0 or later (`LICENSE`, `pyproject.toml`);
  `CITATION.cff` gives the software citation.
- README rewritten for newcomers (quick start through the pipeline, coverage, how to add a method
  or a workbook); CONTRIBUTING covers the snapshot gate and workbook edits.
- The 88 pre-audit result CSVs moved to `Results/archive/pre-audit/`; `multifat-run` writes to
  `Results/current/` by default.
- Removed the 7 October audit checkpoint and the 15 October cleanup audit, superseded by this
  changelog; they remain in the git history. The remaining maintenance items they listed (live
  Issues, repository description and topics, branch protection) need repository-admin access.
