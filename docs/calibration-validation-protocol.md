# Calibration and independent-validation protocol

**Status:** Scientific review proposal for Issue #14. No production calibration procedure is changed by this document.

**Scope:** Multiaxial fatigue criteria implemented in \`multifat\`, initially for constant-amplitude loading and smooth specimens. The protocol applies per method, material identity and experimental campaign. It covers criterion-specific coefficients, fatigue-curve fitting and the provenance of independently assessed predictions.

## 1. Calibration principle

The **default** is the procedure prescribed by the original authors for the particular formulation implemented. This means reproducing their calibration equations, required experimental conditions, coefficient constraints and fitting objective on a new material, **not** copying numerical constants from a paper for an unrelated material.

An original procedure may be algebraic, based on specific fatigue limits, based on selected life curves, or based on a constrained optimisation. Where it prescribes no adjustable coefficient, no extra coefficient should be fitted merely to improve a benchmark. If the authors give multiple procedures, the implemented variant and the selected procedure must be identified explicitly.

The default cannot be marked \`source_verified\` unless the primary publication (or a documented, authoritative author clarification) actually supports its inputs, equations and objective. A repository adaptation must retain a different name and provenance. If the original procedure is unknown, inconsistent or unsupported by the available experiments, report \`source_review_required\` or \`insufficient_calibration_data\` rather than quietly using a generic optimiser.

## 2. Calibration inputs and fitted parameters

Keep separate:

1. **Measured material properties.** Examples include elastic modulus, Poisson's ratio, yield or tensile strength and experimentally established fatigue limits. These are supplied inputs with units, material and batch provenance; they are not automatically fit variables.
2. **Criterion-specific coefficients.** Examples in the current code include \`k_FS\`, \`k_FI\`, \`L_LI\`, \`k_MGSE\` and \`a_MSWT\`. Their calibration rules depend on the implemented criterion and its original source.
3. **Fatigue-curve parameters.** Basquin and strain-life coefficients, or a damage-parameter master curve, are a separate fitting task with its own observation-selection and validity requirements.
4. **Numerical controls.** Critical-plane grid resolution, search bounds, stopping tolerances and shear-amplitude rules must be recorded but are not material coefficients.

Do not estimate a missing physical property and label it as measured. Model-specific derived estimates, if scientifically permitted, must declare their equations and source. Calibrated criterion coefficients and fitted master curves must not be silently conflated.

## 3. Proposed user-selectable modes

| Mode | Selection and objective | Scientific interpretation |
| --- | --- | --- |
| \`author\` (default) | Use the original method's validated calibration selector, equations and objective. | Source-prescribed calibration on eligible material data. |
| \`all_eligible\` | Fit the explicitly designated free criterion coefficients using all *eligible experimental observations* from the permitted material/batch group. Declare the objective, weights and constraints. | User-forced optimisation, not automatically the original method's procedure. |
| \`subset\` | As above, but use a saved, explicit list of eligible observation IDs, loading conditions or named source groups. | User-forced subset optimisation with reproducible provenance. |
| \`supplied\` | Reuse explicitly provided constants after compatibility checks; do not fit them. | External or historically calibrated coefficients. |

The term "all available material properties" must not be interpreted as an instruction to optimise every property in the \`Summary\`. Users select *experimental observations* for fitting, and specify which criterion coefficients are free. The \`all_eligible\` mode does not mean fitting on every row regardless of loading mode, fatigue regime, material identity, missing tensors or censoring.

A user-forced mode may use all available proportional and non-proportional fatigue failures **only if** the model's objective is defined for them. A model whose source-defined calibration needs specific uniaxial/torsional properties cannot automatically exchange those for arbitrary multiaxial observations. An alternative multiaxial fitting objective must be explicitly identified as a repository choice.

The current \`multifat-pipeline --constants\` option supplies recorded coefficients; it is **not** an implementation of all four proposed modes. Software changes belong to subsequent, reviewed tasks.

## 4. Source, material and observation eligibility

Before fitting or assessing a criterion:

- Use one specified material identity, heat treatment and compatible experimental batch. Alternative source names are not evidence of equivalence. Keep unresolved identity/batch conflicts visible and prevent their use in pooled fitting unless explicitly accepted with a documented rationale. Coordinate the review with Issue #24.
- Record specimen geometry and exclusions. Notched specimens remain outside the present plain-specimen scope.
- Validate stress and strain components, units, harmonic load paths, load ratios and the definition of fatigue limits (amplitude versus range, reference cycles, and whether endurance is experimentally established).
- Preserve each observation's workbook, source row, reference, specimen/case ID, dataset fingerprint and failure/run-out status.
- Apply each original procedure's required loading-mode and stress-ratio selectors. Do not silently relax the allowed ratio or include a different loading mode. Any tolerance or allowed fallback needs explicit metadata and scientific approval.
- Treat run-outs as right-censored observations. They are **not** failures with life equal to the interruption cycle. Exclude them from ordinary failure-only regression and errors; include them in an estimation objective only if that objective explicitly handles censoring.
- Treat a fatigue-limit observation separately from a life observation. Survival up to a finite count does not establish an endurance limit.

Calibration subsets are selected deterministically and saved **before** evaluating prediction errors. If observations are added, corrected or removed, reconsider the fit and the calibration's compatibility fingerprint. Source/model revisions and optimisation settings must also be versioned.

## 5. Independent validation and fair comparisons

For each method and material, distinguish at least:

1. **Calibration fit**, evaluated on the observations used to estimate the constants or the master curve. Report these as in-sample, not independent predictive accuracy.
2. **Held-out evaluation**, on observations that were not used in *any* fitting step for the evaluated prediction. This applies both to criterion coefficients and to any downstream damage-parameter master curve.
3. **Source-prescribed versus user-optimised comparison**, reporting objective, observation count, fitted coefficients, eligible methods and prediction metrics on the *same* independent assessment set.

Use grouped validation by specimen batch, experimental campaign or publication source where trustworthy identifiers exist. Never split correlated replicates from one test campaign across training and validation merely to create more nominal folds. If independent material batches are not available, explicitly label within-batch validation and its limitations. For cross-material assessment, avoid fitting material-specific constants on observations of the held-out material unless the evaluation explicitly permits that transfer.

For fatigue-life methods, report native-life errors when the model predicts life itself, and master-curve-inverted life errors for damage-parameter methods, with the prediction basis clearly distinguished. For fatigue-limit methods, assess fatigue indices or limit-state measures; do not invent an unsupported S–N relationship. Mixed datasets require the separate assessment rules of Issue #18. Run-outs must not contribute ordinary life-error metrics.

**If \`all_eligible\` uses every failure observation available for one material, no independent within-material failure set remains.** Report training fit and, if possible, validation against a separate campaign, material batch or external dataset. Otherwise state that independent validation was unavailable.

## 6. Source-review inventory

The following entries describe the **current repository implementation or literature link**, not a completed verification of the original authors' calibration rules.

| Criterion / family | Repository behaviour or evidence | Primary-source verification needed |
| --- | --- | --- |
| Fatemi–Socie (\`FS\`) | The current calibrator fits \`k_FS\` through fatigue-scatter minimisation on selected reversed axial and torsional failures. | Fatemi & Socie (1988), DOI [10.1111/j.1460-2695.1988.tb01169.x](https://doi.org/10.1111/j.1460-2695.1988.tb01169.x). Verify the applicable variant, coefficient identification and curve-calibration sequence. The plane-selection convention may change the calibration procedure. |
| Findley (\`FIN\`) | \`k_FI\` is currently fitted with the same generic scatter objective; shear-amplitude definition is recorded. | Findley (1959), DOI [10.1115/1.4008327](https://doi.org/10.1115/1.4008327). Verify fatigue-strength inputs, normal/shear stress conventions and the calibration relation for this implementation. |
| Li (\`LI\`) | \`L_LI\` is fitted through the generic scatter objective. | Li et al. (2021), DOI [10.1016/j.euromechsol.2021.104353](https://doi.org/10.1016/j.euromechsol.2021.104353). Check the source's coefficient definition, constraints and required data. |
| MGSE and modified SWT (\`MGSE_ZHU\`, \`MSWT\`) | \`k_MGSE\` and \`a_MSWT\` can be fitted by scatter minimisation; the pipeline also provides a source-linked elastic \`k_MGSE\` estimate. | Zhu et al. (2018), DOI [10.1016/j.ijfatigue.2018.03.028](https://doi.org/10.1016/j.ijfatigue.2018.03.028). Establish whether the printed relation or an optimised variant is the source-defined default for each parameter. |
| Fatigue-limit criteria | Several criteria derive coefficients from endurance limits, without optimising life scatter. | Verify each formula and its stress-ratio, mean-stress and endurance-limit conventions. Do not impose a universal life-based calibration. |
| Methods with native fatigue life (\`MDC\`, \`MOREL\`, \`MMCCM\`, etc.) | Native-life outputs and fitted material curves coexist with the generic damage-parameter framework. | Verify the original life model's parameter identification independently of any repository master curve. |
| \`TVEIT2025\` in open PR #21 | The proposed implementation uses hand-entered \`A\`, \`C\`, \`K\`, \`L\`; stage-2 calibration is not yet implemented. | Tveit et al. (2025), DOI [10.1016/j.mechmat.2025.105452](https://doi.org/10.1016/j.mechmat.2025.105452), notably the identification equations and independent identifiability of its constants. Keep this entry provisional until PR #21 is integrated. |

For every implemented criterion, create a source-audit entry containing the primary DOI, equation/table/page, parameter names and units, required observations, calibration objective, numerical constraints, supported formulations, and independent reproduction check. Mark each \`source_verified\`, \`adaptation_documented\`, \`source_review_required\` or \`not_applicable\`. Do not promote a method to source-verified on the basis of a similar published implementation or a passing software test alone.

## 7. Proposed comparison experiment for Issue #14

Before changing the production calibration, run a focused, reproducible comparison on representative datasets:

- **AISI316L**, using the final Janssens-corrected workbook (39 failures and 13 run-outs after the documented overlay). Separate Janssens and Braun material/campaign inputs where needed. The latest Janssens \`R = -1\` reference limits are 278.5 MPa tensile and 144 MPa torsional.
- **42CrMo4_QT**, where axial and torsional fatigue data can support at least some of the existing calibration blocks, subject to source and batch checks.
- **Ti6Al4V / Ti6Al4V_Berto**, where the material-source separation must be made explicit before fitting or evaluating results.

For the methods whose original calibration is verified, evaluate (a) the authors' prescribed calibration, (b) \`all_eligible\` with a declared alternative objective, and (c) at least one predefined subset. Save selected observation IDs and grouped hold-out IDs for every fit. Compare fitted constants, fit residuals, independent prediction errors, sensitivity to the selection, eligibility failures and computation time. Never fabricate a validation result when all observations are used for fitting.

**This comparison remains to be executed.** No numerical comparison, scientific approval or new user-facing calibration mode is claimed in this first documentation phase.

## 8. Required calibration record

A future method-calibration record should capture at least:

- method ID, implemented formulation/version, primary reference and verification state;
- mode, optimisation objective, parameter bounds, constraints, weighting, algorithm, random seed (where applicable) and numerical settings;
- input workbook and immutable content identity, material/batch identity, measured/estimated property origins and units;
- exact calibration and independent-validation observation IDs, source references, selection predicates and exclusion reasons;
- separate criterion-coefficient and master-curve fitting stages, their data use, parameters, status, uncertainty/identifiability and reasons for failure;
- run-out/censoring treatment, fatigue-life or fatigue-limit metric definitions, software revision and compatible existing calibration IDs.

Do not overwrite historical records or silently reuse an incompatible material or formulation. Invalid source matching, missing required test types or a non-identifiable objective should produce an explicit diagnostic rather than a seemingly successful default fit.

## 9. Suggested subsequent software tasks

Following scientific review and the primary-source audit, prepare independently reviewable development work:

1. Add method-level calibration descriptors and source-verification status without altering the implemented equations.
2. Provide a deterministic observation/group selector and persistent training/validation split manifests.
3. Implement the source-prescribed calibrators method by method, with primary-source reproduction tests.
4. Add explicit \`author\`, \`all_eligible\`, \`subset\` and \`supplied\` controls with parameter/objective compatibility checks.
5. Keep model-constant fitting separate from master-curve fitting and ensure independent validation excludes both sets of design observations.
6. Produce matched comparison reports, calibrated constants, run-out diagnostics and auditable JSON records.
7. Coordinate classification and comparison semantics with Issue #18; coordinate batch identity checks with Issue #24.

**Review gate:** Issue #14 should remain open until author procedures are source-verified and the requested representative-material comparison has been performed. This document establishes the initial scientific protocol but does not yet satisfy those acceptance criteria.

## Related repository material

- [Current calibration workflow](calibration-workflow.md)
- [Method audit and implemented coefficient origins](method-audit.md)
- [Experimental-data and assessment rules](dataset-workflow.md)
- [Calibration record archive](../Calibrations/README.md)
- [Primary-reference catalogue](bibliography.md)
