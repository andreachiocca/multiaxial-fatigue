# Method comparison capabilities

This register describes the **implemented case-level methods**, not every paper in the
scientific catalogue. It is the audit and user-facing counterpart of
`multifat.models.comparison.ComparisonCapabilities`. A method's declared capability
does not guarantee that a particular material has sufficient material properties,
independent calibration data, or a finite prediction.

## Access through the method API

Every registered model returned by `get_models` exposes the capability object
as `model.comparison_capabilities`. Methods without a native life law receive
the conventional default explicitly on registration. A new method that
implements `native_life` must declare its own capabilities and pass the
method-inventory regression test.

```python
from multifat.models.registry import get_models

method = get_models(["MOREL"])[0]
print(method.comparison_capabilities.life_prediction)  # native
print(method.comparison_capabilities.dp_comparison)    # True
print(method.comparison_capabilities.dp_reference)     # master_curve
```

## Output semantics

- **Life, fitted master curve**: the method computes a damage parameter `DP`, a
  data-dependent curve `DP = A N_f^b` is fitted, and the resulting
  `Nf_expected` is compared with the observed failure life. This is an
  *assessment-layer* life prediction, not necessarily an original-author life law.
- **Identifiable master curve (issue #74).** A fitted master curve is used only if its design
  set determines it: at least 5 tests, over at least one decade of life, with a negative slope
  distinct from zero ($|b| > 2$ standard errors), the thresholds of the Basquin fits. Otherwise
  the model is skipped on that dataset with the reason. Every fitted curve, used or not, is
  listed in `<material>_master_curves.csv` (n, decades, A, b and its standard error, R² in
  log space, scatter, identifiable, reason). The scatter is estimated with n − 2 degrees of
  freedom.
- **Life, native**: `Nf_native` is computed from the method's own fatigue-life law.
  The `Nf_expected` CSV column contains exactly the same native prediction for
  plot compatibility, with `Metric_mode=Nf_direct`. No supplementary generic
  Basquin inversion is used.
- **DP, fitted master curve**: `CP_value` (the method's DP at the applied load)
  is compared with `DP_fit = A N_{f,exp}^b`, calibrated on the selected design
  data. `DP_fit` is **not a directly measured experimental damage parameter**.
- **DP, model law**: the method itself provides the reference DP corresponding
  to the measured life. The reference is independent of the generic master curve.
- **DP, unsupported**: a diagnostic fatigue index must not be substituted for
  the physical damage parameter merely to fit the generic comparison interface.

A method supports a kind of plot only if it declares the appropriate scientific
capability **and the current result has the required values**. The `_ext` critical
plane variant inherits the corresponding method's capability.

## Audited native-life methods

| Method | Life route | DP comparison route | DP interpretation |
| --- | --- | --- | --- |
| `TVEIT2025` | Native continuum damage integration | **Not available** | The stage-1 auxiliary virgin endurance index ignores back-stress/damage evolution. |
| `MOREL` | Native Morel Eq. (38), fitted `p,q,r` | Empirical DP–life fit | `T_S/T_{S,lim}` is an endurance index. This separate empirical DP comparison does **not** assess the native Eq. (38) fit. |
| `MDC` | Native Eq. (8), torsional Basquin law | Model law, `tau'_f(2N_{f,exp})^{b_0}` | `tau_r`, equivalent shear stress in MPa. |
| `MMCCM` | Native inversion of the Manson–Coffin strain–life curve | Model definition, `1/N_{f,exp}` | `DP=1/N_{f,pred}`, reciprocal cycles. |

For `MOREL`, native-life fitting requires at least three suitable reversed
uniaxial design tests above the endurance threshold. If this fitting procedure
cannot be completed, its native prediction remains unavailable. A DP fit, if
available, is not silently substituted for it. For MDC and MMCCM, the original
model-law DP reference is independent of any fitted universal DP curve.
The published model parameters, their calibration populations, and held-out
validation populations must be recorded before independent error rankings.

## Other registered methods

All **31** other registered models currently use the implemented
`master_curve` life route and support calibrated `DP` comparison when a
usable dataset-level DP master curve exists. This includes every registered
method without a native-life implementation. They are listed below to make
the scope explicit and to support regression checks when new methods are added.

`BP`, `CAIM`, `COIN`, `CROSSLAND`, `DANGVAN`, `EK`, `FIN`, `FS`, `GSA`, `GSE`, `LI`, `LIU2021`, `MATAKE`, `MCD`, `MGSE_YU`, `MGSE_ZHU`, `MKBM`, `MMMH`, `MSWT`, `MWCM`, `OTT`, `PAPADOPOULOS`, `PBP`, `PCR`, `PI`, `RFB`, `SINES`, `SWT`, `SWTD`, `VHN`, `VSE`.

In total, the registry currently contains **35** methods (17 critical-plane
and 18 non-critical-plane). These are **implementation/comparison capabilities**,
not claims of complete reproduction of each original paper or of independent
experimental validation. See [the method index](METHOD_INDEX.md) and
[the scientific audit](method-audit.md) for those separate assessments.

## User commands

```bash
# Comparable life-based plots. Uses the native law where declared.
multifat-compare --results_dir Results/current --out_dir Results/Comparison --plot-kind life

# Only methods with an interpretable damage-parameter comparison.
multifat-compare --results_dir Results/current --out_dir Results/Comparison --plot-kind damage

# Both plot families (the default).
multifat-compare --results_dir Results/current --out_dir Results/Comparison --plot-kind both
```

## Load-equivalent error: the cross-method ranking metric

`Error_load` $= e_L = \ln(1/s^*)$, where $s^*$ scales the whole load history (means included;
strains follow elastically) so that the method predicts exactly the observed life: on the master
curve, $DP(s^*\boldsymbol\sigma) = DP_{fit}(N_{\exp})$; with a native law,
$N_{native}(s^*\boldsymbol\sigma) = N_{\exp}$. $e_L > 0$ means the method rates the test as more
damaging than it was (conservative). It is the life-data counterpart of `FI_load` on fatigue-limit
datasets ([issue #37](https://github.com/andreachiocca/multiaxial-fatigue/issues/37)).

Why it is the metric to rank methods with:

- **Independent of the S–N slope.** For a DP homogeneous of degree $d$ in the loads,
  $\mathrm{Error\_ln} = k\,d\,e_L$ with $k = -1/b$: the life error is the load error amplified by
  the flatness of the master curve. On AISI316L ($k \approx 22$) life errors are huge, or undefined
  (`DP_DIFF`), for every method, whatever its quality.
- **Independent of the DP's degree.** $\mathrm{Error\_log10\_dp} = -d\,e_L/\ln 10$: quadratic DPs
  (SWT, energy parameters) look twice as bad in DP space as linear ones of equal accuracy.
- **One scale for every route:** master-curve and native-life methods, including tests with an
  infinite native life prediction, which have no finite life error.

Homogeneous DPs use the closed form $e_L = \ln(DP/DP_{fit})/d$; others are solved by iteration on
the load factor (`multifat.eval.load_error`). Run-outs give bounds only and have no `Error_load`.
TVEIT2025 is not evaluated (each evaluation integrates many cycles). `multifat-compare` reports
mean, standard deviation and the fractions of tests within $\pm10\,\%$ and $\pm20\,\%$ of the
load (`Within_10pct`, `Within_20pct`) for every method. The complexity of a method (constants
fitted on fatigue data, required inputs) is to be reported alongside, and is best penalised by
held-out validation (#14).

## Cross-material ranking KPI and complexity axis (issue #57, decision d116)

Comparing fatigue criteria across heterogeneous materials requires handling differences in dataset difficulty, scale and sample size. `multifat-rank` (and `multifat-compare --rank`) implements a two-axis evaluation combining an accuracy KPI with explicit complexity counts and a reproducible selection rule.

### 1. Accuracy KPI: relative RMS load error

For method $m$ on material $k$, over the $n_{mk}$ assessment failures (`Flag = 1`, $n_{mk} \ge n_{\min} = 3$) with load-equivalent error $e_{L,i} = \ln(1/s^*_i)$:

$$R_{mk} = \sqrt{\frac{1}{n_{mk}}\sum_i e_{L,i}^2} = \sqrt{\bar e_L^{\,2} + \tfrac{n-1}{n}\,s_{e_L}^2}$$

Across materials, each method is scored relative to the median method on that material, preventing the hardest or largest dataset from dominating:

$$z_{mk} = \ln\frac{R_{mk}}{\operatorname{med}_{m'} R_{m'k}}, \qquad Z_m = \frac{1}{K_m}\sum_{k} z_{mk}$$

$e^{Z_m}$ reads directly as the typical relative accuracy: $Z = -0.39$ indicates a typical load error 32 % smaller than the median method on the same materials. The standard error $\mathrm{SE}(Z_m)$ is estimated by bootstrap resampling (1000 iterations) of tests within each material.

Additional diagnostics reported alongside $Z_m$:
- **Median bias:** $\operatorname{med}_k \bar e_L$;
- **Error bands:** fraction of tests within $\pm10\,\%$ ($|e^{e_L}-1| \le 0.10$) and $\pm20\,\%$;
- **Non-conservative share:** fraction of predictions where $e_L < -\ln(1.1)$ (unsafe by $> 10\,\%$);
- **Generalisation gap / Transferability:** $\Delta_m = \operatorname{med}_k \ln(R^{\mathrm{assessment}}_{mk} / R^{\mathrm{design}}_{mk})$.

### 2. Explicit complexity axis: $(p_{\mathrm{in}}, p_{\mathrm{fit}})$

Complexity is reported beside the accuracy KPI and never arbitrarily folded into a composite score:
- **$p_{\mathrm{fit}}$**: number of quantities fitted on fatigue data. Covers the master curve ($A, b$: 2) plus any calibrated plane constants (`k_FI`, `k_FS`, `L_LI`, `k_MGSE`, `a_MSWT`, `k_PDS`), or native-law parameters (MOREL $p,q,r$: 3; TVEIT2025: 4; MDC: 2; MMCCM: 4; LIUMAHADEVAN: 4).
- **$p_{\mathrm{in}}$**: number of experimental characterisation classes needed beyond the S–N data of the calibration paths:
  - axial fatigue limit (`Sigm1`);
  - torsional fatigue limit (`Taum1`);
  - mean-stress ($R = 0$) data (`Sig0`, `Tau0`);
  - static strength (`Su`, `Sy`, `Tu`);
  - S–N curves (`sigma_f`, `b_fat`, `tau_f`, `b0`);
  - strain-life/cyclic curve (`epsilon_f`, `c_fat`, `gamma_f`, `c0`, `K_prime`, `n_prime`);
  - elastic constants (any strain-based method: `E`, `nu`).

### 3. Family-wise 1-SE selection rule

To pick the best criteria within each family (Critical plane, Invariant, Integral, Energy):
1. Filter methods meeting minimum material coverage $K_m \ge K_{\min}$ (default 4).
2. Within each family, find the best accuracy $Z_{\mathrm{best}} = \min Z_m$ and its standard error $\mathrm{SE}_{\mathrm{best}}$.
3. Form the set of statistically indistinguishable candidates: $\{m : Z_m \le Z_{\mathrm{best}} + \mathrm{SE}_{\mathrm{best}}\}$.
4. Select the candidate with the lowest complexity: sort by $p_{\mathrm{in}}$, then $p_{\mathrm{fit}}$, then $Z_m$.
5. Test pick stability via leave-one-material-out replication.

## Error statistics and observation roles

For a valid DP comparison, `Error_log10_dp = log10(DP_fit / CP_value)`.
A positive value means that the reference DP at the observed life exceeds the
method's DP at the applied load. The retained `Error_ln_dp` column has the
opposite sign and a natural-log scale:
`Error_ln_dp = ln(CP_value / DP_fit) = -ln(10) Error_log10_dp`.
The life error is `Error_ln = ln(Nf_exp / Nf_expected)`. Undefined, nonpositive
or nonfinite quantities do not enter finite-error statistics. Infinite native
life predictions remain visible in the separate native-life summary.

`summary_error_stats.csv` and its styled table retain their historical
aggregate of eligible failed rows for compatibility. They may combine design
and other assessment observations and must not be labelled independent
validation. `summary_error_stats_by_role.csv` reports the same metrics
separately for `Flag=0` design observations and `Flag=1` other assessed
failures. Rows without a recognised flag are labelled `unclassified`;
`Flag=2` run-outs contribute to neither ordinary error table. The
`--common-methods` option restricts comparison to finite metrics on matched
failed observations assessed by every selected method. Tests behind a constant
fitted by the pipeline (Basquin constants and the plane constants `k_FI`, `k_FS`,
`L_LI`, `k_MGSE`, `a_MSWT`, `k_PDS`) are `Flag=0` for every model that adopts that constant.
A `Flag=1` row is still not automatically independent of all fits: the R = -1
fatigue limits, constants supplied in the workbook or by an upstream record, and
other experimental batches of the same source are not traced.
Establishing a source- and batch-aware held-out validation population remains
part of [Issue #14](https://github.com/andreachiocca/multiaxial-fatigue/issues/14).

The command writes `comparison_capabilities.csv` showing per-method and
per-result eligibility and DP reference source. Native-life models also write
`native_life_prediction_summary.csv`, which reports design observations,
assessment failures and run-outs separately, including the positive-infinite
native-life predictions in each group. Design observations are not independent
assessment failures. Infinite predictions are not silently
discarded as if the method made no prediction. Such counts are not an
independent validation statistic by themselves.

Comparison availability must not be confused with evidence of predictive
accuracy. Tests used to fit coefficients or DP master curves are *design data*,
not independent validation observations. For study-level comparisons, also
use matched source/batch-aware validation subsets.
