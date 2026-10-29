# Method comparison capabilities

This register describes the **implemented case-level methods**, not every paper in the
scientific catalogue. It is the audit and user-facing counterpart of
`multifat.models.comparison.ComparisonCapabilities`. A method's declared capability
does not guarantee that a particular material has sufficient material properties,
independent calibration data, or a finite prediction.

## Output semantics

- **Life, fitted master curve**: the method computes a damage parameter `DP`, a
  data-dependent curve `DP = A N_f^b` is fitted, and the resulting
  `Nf_expected` is compared with the observed failure life. This is an
  *assessment-layer* life prediction, not necessarily an original-author life law.
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

The command writes `comparison_capabilities.csv` showing per-method and
per-result eligibility and DP reference source. Native-life models also write
`native_life_prediction_summary.csv`, which counts experimentally failed
specimens for which the model predicts infinite life separately from run-outs
and from finite logarithmic life errors. Infinite predictions are not silently
discarded as if the method made no prediction. Such counts are not an
independent validation statistic by themselves.

Comparison availability must not be confused with evidence of predictive
accuracy. Tests used to fit coefficients or DP master curves are *design data*,
not independent validation observations. For study-level comparisons, also
use matched source/batch-aware validation subsets.
