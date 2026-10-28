# Experimental inputs after revision, and the comparison subsets

Source: upstream PRs #5 and #6 (commits 904c619, 52b061b), ported to `multifat`. The import rules, eligibility outputs, observation identities and matched subsets follow the upstream design; the pipeline commands and the fatigue-limit assessment were built for the `multifat` package.

The October 2026 collection comprises twelve revised workbooks, installed without rewriting any measured cell, together with the retained `Al7075_T6.xlsx` and two batch workbooks introduced here (`14_42CrMo4_Reis.xlsx`, `15_Ti6Al4V_Berto.xlsx`). Byte-for-byte copies of the earlier workbooks remain under `Experimental_Data/archive/pre-dataset-migration/`; only the active folder is read by the importer. Every change applied here after the merge is recorded cell by cell in `Experimental_Data/PROVENANCE.md`.

```bash
uv run python tools/install_revised_dataset.py exp_data.zip
uv run python tools/install_revised_dataset.py exp_data.zip --apply
uv run python tools/audit_dataset_migration.py --old-dir Experimental_Data/archive/pre-dataset-migration
```

## Provenance and import

Header rules live in `multifat.imp.schema` (see [workbook format](workbook-format.md)). An ordering prefix such as `6_` is stripped from the file name, which yields the material ID `42CrMo4_QT`; two files carrying the same ID are rejected, and `--material 42CrMo4_QT` locates `6_42CrMo4_QT.xlsx`. Alloy, treatment and legacy axial/bending suffixes are left unmerged. When the source rows of a Summary name an alloy different from the file's, the case is reported as an identity issue, and the workbook itself is not rewritten.

Summary properties other than the R = -1 fatigue limits are taken from the named Average cell when that cell exists; otherwise they come from the row labelled Average. In a revised Summary with a `Stress ratio R` column, `Sigm1` and `Taum1` are arithmetic means of their respective positive, finite R = -1 source values whose material name matches the workbook. Each selected value retains the contributing Excel rows, references, authors, stress ratios and reference cycle counts in `fatigue_limit_sources`. A property with source values but no compatible R = -1 value is unavailable and reported. Legacy workbooks retain the Average-row rule. Properties that are absent remain absent. Strains are undefined without E or nu; stress observations are imported regardless, and only the models that require strains get skipped (`multifat.models.registry.stress_only_models()` finds the others by evaluation).

The source rows do not consistently identify metallurgical batches or whether a reported finite-life strength is an amplitude or a range. Same-ratio averaging is therefore a repository selection policy, not proof that sources describe one batch or a common endurance limit. In particular, Braun's R = 0 value for AISI316L is a stress range at 10^6 cycles; it is not assigned to `Sig0` or combined with Janssens' R = -1 amplitude. No mean-stress correction is applied. Review source compatibility before interpreting the resulting model comparisons.

Dataset columns are chosen by header, whatever their positions. Both `Amplitute` and `Amplitude` are accepted, and BENDING stress headers are read as the axial (longitudinal surface) stress. Forces in N and moments in Nmm stay in the workbook and are never substituted for stresses, so the stress-based reconstruction is unchanged, including the tensorial shear strain epsilon_xy = tau_xy / (2G).

In the revised format (a `Phase AXIAL load` column), radial and torsional phases are referenced to the axial phase, and a row is excluded, with its reason recorded in the validation report, when

- a stress cell is empty or `Null` (unknown loading, not an inactive component);
- active components have different known frequencies (the model is single-harmonic); frequency text such as `15-20` is retained as text and does not exclude a row;
- an internal pressure is given without a radial (hoop) stress.

Legacy workbooks keep their own conventions: a `Null` stress denotes an inactive component (0), and a missing phase is 0. Known failures and run-outs stay distinct; notched rows fall outside the scope.

## Input problems requiring the source

| Material | Constraint |
|---|---|
| CX18 | The Summary names 30NCD16 and holds its limits, yet rows 2-3 agree with FatLim's XC18; 3 of its 5 rows also occur in 30NCD16. Both R = -1 limits remain unavailable until the identity is resolved. |
| Ti6Al4V | One Summary source names Ti64Al4V (Delahay). Its values are not pooled with the matching Meneghetti row; the selected R = -1 limits are 500 and 422.5 MPa. The source identity still needs review. |
| AISI316L | The Average row mixes Janssens' R = -1 amplitude with Braun's R = 0 stress range. Source-row selection currently yields 241 and 144 MPa, but the newly uploaded Dataset contains corrected Janssens axial stresses while the Summary is unchanged. The axial limit and the distinct Braun batch require scientific review before using the resulting comparisons. |
| En3B, GJS600 | Every observation is notched, outside the plain-specimen scope. |
| Fatigue-limit sets | CX18, 30NCD16, 25CrMo4, 982FA, 34Cr4, GGG60 and S65A contain limit states only. In five sets their pure R = -1 rows contradict the Summary limits (see PROVENANCE.md), so their fatigue index errors await the source check. |

## Assessment and calibration

`multifat-pipeline` validates every workbook, calibrates only those constants a model needs but the workbook lacks, runs each model, and reports it as ok, skipped (with the missing inputs) or failed; a failure in one model never halts the others. `--constants FILE` draws constants from a `Calibrations/` record or an earlier sidecar instead of calibrating them (see [`Calibrations/README.md`](../Calibrations/README.md)).

When several workbooks are supplied (or `--all` is used), `Comparison/eligibility_matrix.csv` gives the method x material status, `method_availability.csv` the per-method counts, and `eligibility.json` the reasons:

| Status | Interpretation |
|---|---|
| available | the model ran and the material offers failures to compare |
| prediction_only | the model ran on run-outs alone |
| fatigue_index | the model ran on a fatigue-limit dataset (fatigue index error; no failures counted) |
| calibration_required | `--no-calibrate` was set, and only constants the pipeline can calibrate are missing |
| insufficient_inputs | a property or constant that the pipeline cannot supply is missing |
| out_of_scope | the workbook contains notched specimens only |
| failed | the model raised an error (reported, not hidden) |

Availability is not an accuracy ranking. If two selected workbooks share a material ID, the run stops before anything is written.

```bash
uv run multifat-pipeline --all --jobs 8 --results-dir Results/current
uv run multifat-pipeline --material 42CrMo4_QT --models BP,SWT --no-calibrate
```

## Fatigue-limit datasets

A workbook whose plain rows are all run-outs and whose Summary provides both R = -1 fatigue limits constitutes a set of limit states. Every model is then judged through its fatigue index error (Böhme & Papuga 2024, Eqs. 8-9): FI_axial = DP(case) / DP(axial R = -1 at Sigm1) - 1, and FI_torsion with the torsional limit. FI_load = 1/s - 1, where s is the factor on all loads at which DP(s * case) matches the axial reference, is the stress-scaled index for any damage parameter (it equals FI_axial for criteria linear in the loads, and sqrt(DP/DP_ref) - 1 for SWT-type parameters). A missing plane constant is set so that both limits yield the same DP; models that need life constants are skipped. Results are written to `FatigueLimit/<METHOD>_<MATERIAL>_FI.csv` and `<MATERIAL>_fatigue_index_summary.csv` (statistics over finite values; rows where a model is infinite by design on some loads, such as MWCM past its rho cut-off, are counted in `n_nonfinite`); the life comparison ignores that folder. See `multifat.eval.endurance`.

## Failure observations shared across methods

Result CSVs end with `Test_status`, `Specimen_id`, `Case_id`, `Observation_id` and `Dataset_fingerprint`. A case is identified by its specimen and stress history; an observation additionally by its cycle count, status and replicate occurrence. The dataset fingerprint hashes material properties and failure histories, so neither row order nor run-outs alter it.

`multifat-compare --common-methods A B ...` limits the comparison to those failure observations assessed by all selected methods, for the life (`--comparison-metric life`) or damage (`damage`) error, and writes `selected_common_observations.csv`, `subset_definition.txt`, `individual_availability.csv` and `common_subsets.json`. Duplicate method-observation results, or conflicting lives, raise an error instead of an inferred match. Run-outs stay in the plots as open markers with arrows, but never enter the statistics. `scripts/3_review_subsets.py` forwards to the same command.

```bash
uv run multifat-compare --results_dir Results/current --out_dir Results/current/Common \
    --common-methods BP SWT --comparison-metric life --flag_filter other
```

Neither a common subset nor a passing regression test counts as independent experimental validation.

## Checks

`uv run pytest` runs everything, including `tests/test_dataset_import.py` (revised headers, IDs, phases, exclusions, missing properties, batch workbooks), `tests/test_subsets.py`, `tests/test_dataset_tools.py` and `tests/test_endurance.py`.
