# Numerical method audit — 7 October 2026

This audit covers all 20 previously registered base methods, the shared plane evaluator, calibration, and reporting. It corrects identifiable implementation errors and separates usable damage-parameter implementations from incomplete proxies. Passing analytical tests establishes the tested numerical properties; it does not establish predictive accuracy for every material or reproduce every author's validation dataset.

The source catalogue is [bibliography.md](bibliography.md), backed by `fatigue/models/references.py`. `get_models()` attaches a reference to each model, and new CSVs include `Method_reference` and `Implementation_revision`. No article PDFs are added to the repository.

## Retained methods

`tau_a`, `gamma_a`, and `epsilon_a` denote amplitudes; gamma is **engineering** shear strain (`gamma_ij = 2 epsilon_ij`). `sigma_max` is the plane-normal maximum. Except for FIN, existing tension-only normal-stress conventions remain as documented adaptations. The table distinguishes a checked numerical implementation from full independent source reproduction.

| Name | Implemented quantity / plane | Audit action and reference basis |
| --- | --- | --- |
| FS | `gamma_a (1 + k_FS sigma_max/Sy)` on maximum shear-strain plane | Retained. Tensor shear range already equals engineering shear amplitude, so no extra factor of two is applied. Original Fatemi–Socie plane rule checked against the publication abstract and equation reproduced by Yu (2017). Cyclic yield strength is approximated by workbook `Sy`. |
| FIN | `tau_a + k_FI sigma_max`, maximum damage plane | Fixed plane selection: previously maximized shear range only. Removed clipping of the signed normal maximum. Findley's maximum combined-parameter plane and stress form retained. |
| SWT | `sigma_max epsilon_a`, maximum normal-strain plane | Retained. Classical uniaxial SWT extended to candidate planes; compression clipped to zero. `_ext` selects maximum product instead. |
| MATAKE | `tau_a + (2 Taum1/Sigm1 - 1) sigma_max`, maximum shear plane | Retained fatigue-limit-calibrated form. Tension-only normal maximum is a repository adaptation; not claimed to reproduce the original paper for all mean stresses. |
| MCD | `tau_a + Taum1/(2 Su) sigma_max`, maximum shear plane | Retained common surface-crack form; torsional fatigue limit is used for the shear limit. Crack-type-specific calibration is not supplied by the workbooks. |
| MGSE_YU | `tau_max gamma_a + sigma_max epsilon_a`, maximum damage plane | Fixed engineering shear factor and actual shear-path maximum. Scalar equation checked against Yu (2017), Eq. (11). Numerical maximum-damage plane is the repository's explicit interpretation of the near-maximum-shear plane; it is not an exact reproduction of a restricted plane search. |
| MGSE_ZHU | `tau_max gamma_a + k_MGSE sigma_max epsilon_a`, maximum damage plane | Fixed engineering shear factor and actual shear-path maximum. Correct original article identified. The existing weighting and plane rule remain provisional: full original equation/plane-rule reproduction was not independently completed from the accessible article preview. |
| MSWT | `a sigma_max epsilon_a + (1-a) tau_a gamma_a`, maximum damage plane | Fixed engineering shear factor; reject invalid `a_MSWT` instead of silently clipping it. Associated with the modified SWT form reported in Zhu (2018), not presented as the complete incremental Jiang model. Original Jiang reference also listed in the bibliography. |
| MKBM | `gamma_a + (1 + sigma_max/Sy) epsilon_a`, maximum shear plane | Fixed engineering shear factor. Original Li (2011) article identified; analytical consistency checked against the equation recorded in the supplied implementation. Full original-paper reproduction remains unverified. |
| LIU2021 | Expanded modified shear-amplitude equation recorded in `liu2021.py`, maximum shear plane | Fixed engineering shear factor, Poisson-ratio bounds, and removable division by shear range. Preserve supplied `tau_f`. Original Liu et al. (2021) article identified. The interpretation of normal strain range and original Eq. (27) still require independent full-text verification. |
| LI | `(tau_a + norm(mean shear)) gamma_a + L sigma_max epsilon_a`, maximum shear-energy plane | Fixed mean shear: norm of the mean vector, not mean of endpoint norms. Fully reversed loading now has zero mean shear. Fixed engineering shear factor and explicit tie handling. Original Li (2021) article/plane description identified; exact coefficient formulation remains based on the supplied code, not independently reproduced full text. |
| GSE | Elastic/plastic shear and normal energy sum, maximum damage plane | Fixed engineering shear factor. Decompose real stress/strain harmonics before taking scalar ranges; synthetic extrema no longer create fictitious plastic strain. Scalar equation checked against Ince–Glinka form reproduced in Yu (2017), Eq. (3). |
| GSA | Stress-normalized elastic amplitudes plus plastic amplitudes, maximum damage plane | Same strain corrections as GSE. Preserve individually supplied fatigue-strength coefficients. Scalar equation checked against Ince–Glinka form reproduced in Yu (2017), Eq. (4). Missing axial coefficients still use the warned steel-based empirical estimate; missing shear coefficient uses the von Mises conversion. |
| BP | Sphere average of amplitude and mean-stress terms | Existing implementation retained. Original BP source is Böhme–Papuga (2024); implementation equations were supplied in Böhme, Papuga & Lange (2026). Four fatigue-limit cases and phase-origin invariance tested. |
| CAIM | Maximum amplitude term plus sphere-average mean term | Existing implementation retained. Coefficients/clipping checked against supplied `damage_criteria_CAIM.m`; independent dense angular reference and four fatigue-limit cases tested. This is distinct from BP. |

For FS, SWT, MATAKE, MCD, MKBM, LIU2021 and LI, `_ext` means maximum damage over the sampled planes. This is a repository variant, not automatically the cited author's plane-selection rule. FIN and other methods whose standard metric is already damage return identical standard/`_ext` values. Tied maximum metrics are resolved by the larger damage value. Zero harmonic shear correctly resets the encoded endpoints to their mean, even when the legacy component extrema would suggest a nonzero range.

## Disabled implementations

Both normal and `_ext` names are rejected with an explanatory error. Archived Python source remains for inspection and direct evaluation is guarded. Re-enabling requires a corrected formulation, source traceability, and numerical tests. Existing CSVs are historical outputs, not evidence that these formulations are correct.

| Name | Why disabled | What is required to restore it |
| --- | --- | --- |
| CARSPA | The scalar combination was evaluated on a maximum-shear plane or maximized over all planes. Neither implements Carpinteri–Spagnoli's material-dependent plane construction. | Implement the published plane-orientation procedure and validate against reference load cases. |
| DANGVAN | The evaluator reconstructs independent plane extrema, then the model computes hydrostatic stress from that synthetic tensor. Its trace can differ from the physical hydrostatic history. Adding a shear amplitude to independent hydrostatic extrema also loses simultaneity. | Implement the time-dependent mesoscopic/shear-hydrostatic criterion with its proper stress-path treatment. |
| OTT | A proportional cycle reduction is used as an arbitrary non-proportional model; the fatigue-ratio solver silently defaults to `m=2` when it cannot bracket a solution. | Implement/validate the full path model, or expose a clearly named proportional-only reduction with enforced parameter/loading limits. |
| SWTD | Stress and strain eigenvalues are independently sorted and multiplied by array index, which need not pair the same physical directions. The method ignores harmonic input. | Preserve physical directions and implement the published non-proportional applicability limits. |
| ZHU_EDP | Explicitly an “inspired” proxy: omits elastoplastic work `W0` and substitutes `abs(sin(phase))` for the published path factor. Basquin calibration cannot restore missing path physics. | Implement the actual EDP work/path formulation with required material data. Its parameter calibration has also been removed. |

The runner default is now `FS_ext`, because its former default `CARSPA_ext` is disabled.

## Shared numerical conventions and limitations

- Workbook strains are constructed using linear isotropic elasticity. The workflow does not simulate cyclic plasticity or non-proportional hardening. GSE/GSA therefore have zero plastic contribution for these imported histories; external total-strain histories can have a plastic part.
- Legacy CP methods retain the repository's ellipse convention: amplitude is the root-sum-square of the two shear semi-axes. BP/CAIM use the minimum-circumscribed-circle radius (major semi-axis). These are different non-proportional extensions and must not be described as identical original-paper conventions.
- Mean-shifted shear magnitude maxima use 720 time samples. GSE/GSA decompose the true harmonic tensors, not independently synthesized scalar extrema. Plane orientations remain a finite grid; this audit does not prove global continuous-plane maxima.
- BP/CAIM use 24 × 48 sphere quadrature; CAIM refines multiple grid maxima numerically. The MATLAB uses a different angular grid, so agreement is within discretization error, not bitwise equality.
- Missing R=0 strengths for BP/CAIM use the supplied 2026 article's estimates. They are maximum stresses, not amplitudes. These estimates and empirical strain-life properties are not measured values.
- Supplied material constants were calibrated with the old implementation. Recalibrate affected methods before drawing quantitative model rankings. Raw experimental workbooks and archived result CSVs are preserved.

## Damage-space error for every slope

`Error_log10_dp = log10(DP_e / DP)`, where `DP_e = DP_fit = A_surv * N_exp**b` and `DP = CP_value`. It is exported for calibration and assessment rows and included in statistics/plots for every slope. Positive values mean predicted allowable damage exceeds the evaluated point; negative values mean the point exceeds the curve. This sign convention is the reverse of the retained legacy `Error_ln_dp = ln(DP/DP_fit)`.

The default `STDNUM=0` uses the median fitted curve. For a nonzero shift, residual scatter is computed in log-damage space and the curve is shifted downward. This matches the previous life-space shift for negative slopes without dividing by a tiny slope. For flat or rising curves it is a damage-space convention, not an inferred life-survival probability. Nonpositive/nonfinite damage values yield NaN log ratios; they are not silently replaced by arbitrary positive values.

Life predictions are still suppressed when `b >= 0` or `k=-1/b > 15`. Exactly flat curves retain a finite constant `DP_e`. Constant-life calibration data cannot identify a slope and yield an invalid fit; the runner rejects it with an explanatory error. Parameter searches with no finite objective also fail explicitly instead of reporting an arbitrary bound as an optimum. Nonpositive evaluated DP values are retained with undefined log/life metrics rather than passed into life inversion. New columns do not change the meanings of `Error_ln` or `Error_ln_dp`.

## Result provenance

New CSVs carry revision `2026-10-07-audit1`. Comparison excludes disabled methods and pre-audit outputs by default. `--include-legacy` explicitly allows historical outputs from retained methods; it does not re-enable disabled methods. Do not mix corrected and old runs when ranking methods. See the README for regeneration commands.

## Verification

Run `python -m unittest discover -s tests -v`. Tests cover analytical plane values for all retained CP methods, mean-shear cancellation, Findley selection, plane ties, harmonic elastic/plastic splitting, shifted-ellipse maxima, disabled-name guards, bibliography coverage, flat/shallow fits, log base/sign, and export/statistics across both metric modes. BP/CAIM retain their separate limiting-case and angular-reference tests.

Full-text reproduction of the provisional methods listed above and independent validation against each publication's datasets remain outside the evidence established by these tests. They are explicitly retained with that limitation, rather than labelled scientifically validated.
