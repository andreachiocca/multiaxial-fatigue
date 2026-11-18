# TVEIT2025 — Tveit, Ottosen, Erduran & Reyes (2025) Continuum Damage Model

**Registry name:** `TVEIT2025`  
**Family:** Incremental continuum damage, native life  
**Source file:** [`src/multifat/models/invariant_methods/tveit2025.py`](../../src/multifat/models/invariant_methods/tveit2025.py)

## Overview

`TVEIT2025` integrates an incremental continuum endurance-surface model cycle by cycle. A back-stress tensor $\boldsymbol{\alpha}$ and scalar damage $D$ evolve whenever the stress excursion leaves the moving yield-like endurance surface ($\beta > 0$, $d\beta > 0$), delivering native fatigue life when $D$ reaches $1$ or infinite life upon elastic shakedown.

## Original formulation

In Tveit, Ottosen, Erduran & Reyes (2025), Eqs. (1)–(12) and (24) (pp. 2–5), the moving endurance surface and damage evolution laws in deviatoric stress $\mathbf{s}$ and hydrostatic invariant $I_1$ are

```math
\beta = \frac{\sigma_e + A\,I_1 - f_{-1}}{f_{-1}},
\qquad
\sigma_e = \sqrt{3 J_2} - \lambda\!\left(1 - \frac{27}{4}\frac{J_3^2}{J_2^3}\right),
\qquad
\lambda = \sqrt{3}\,t_{-1} - f_{-1},
```

with invariants $J_2, J_3$ evaluated on the reduced deviator $\mathbf{s} - \boldsymbol{\alpha}$ and active evolution (when $\beta > 0$ and $d\beta > 0$):

```math
d\boldsymbol{\alpha} = C\,(\mathbf{s} - \boldsymbol{\alpha})\,d\beta,
\qquad
dD = K\,e^{L\beta}\,d\beta,
\qquad
d\beta = \frac{\mathbf{G} : d\mathbf{s} + A\,dI_1}{f_{-1} + C\,\mathbf{G} : (\mathbf{s} - \boldsymbol{\alpha})},
\qquad
\mathbf{G} = \frac{\partial\sigma_e}{\partial\mathbf{s}}.
```

Starting from $\boldsymbol{\alpha} = \mathbf{0}$ and $D = 0$, failure occurs at $D = 1$, while paths enclosed by the drifted endurance surface shakedown to $N_f = \infty$.

## Repository implementation

**Block integration and steady-state acceleration.** Each harmonic period is sampled at `npts_TVEIT = 1001` points and integrated explicitly block by block until one of four rules triggers:
1. *Failure:* $D \ge 1$ inside a block (linearly interpolated within the crossing increment).
2. *Shakedown:* a block (from block 2 onward) accumulates zero damage ($N_f = \infty$).
3. *Cyclic steady state:* after `min_blocks_TVEIT = 400` cycles, back-stress drift over the last 200 blocks is $\le 0.01$ of its path length and mean block damage agrees within $0.1\%$ of the preceding 200-block window—provided peak $\beta$ over the block is not trending down toward shakedown before $D = 1$. Life is then extrapolated as $N_f = n + (1 - D)/\overline{dD}$.
4. *Block cap:* after `max_blocks_TVEIT = 60000` cycles, life is extrapolated similarly and a warning is emitted.

**Repository choices and departures.**
- *Calibration by the authors' procedure:* without supplied constants, `multifat-pipeline` calibrates $A$, then $(C, K, L)$ by Eq. (25) (see [Calibration](#calibration)); a `--constants` file (e.g. [`Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json`](../../Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json)) takes precedence.
- *Constants-file limits and settings:* `Sigm1_TVEIT` and `Taum1_TVEIT` override workbook `Sigm1`/`Taum1` for `TVEIT2025` only so calibrated constants always match their reference limits. Optional `npts_TVEIT`, `min_blocks_TVEIT`, `max_blocks_TVEIT`, `tol_alpha_TVEIT`, and `tol_damage_TVEIT` may also be overridden in the file. With `1001` points per cycle, finite lives sit roughly $0.02$ decade shorter than converged high-resolution integration (Table 1 NiCrMo steel).
- *Symmetric `G` projection and DP:* `G` uses the deviatoric projection $\mathbf{s}\mathbf{s} - \frac{2}{3} J_2 \mathbf{I}$ of $\partial J_3 / \partial \mathbf{s}$ (Eq. 10 prints a different isotropic piece that cancels upon deviatoric contraction). `values[""]` holds the virgin-state maximum of $(\sigma_e + A I_1)/f_{-1}$ ($\beta + 1$ at $\boldsymbol{\alpha} = 0$), while `aux["Nf"]` returns native life.

## Material constants and inputs

| Symbol | Workbook / JSON key | Python name | Unit | Notes |
|---|---|---|---|---|
| $A, C, K, L$ | `A_TVEIT`, `C_TVEIT`, `K_TVEIT`, `L_TVEIT` | same | — | Required model constants: calibrated by the pipeline, or supplied via `--constants` |
| $f_{-1}, t_{-1}$ | `Sigm1_TVEIT`, `Taum1_TVEIT` | same | MPa | Optional limits replacing workbook `Sigm1`, `Taum1` for this model |
| solver controls | `npts_TVEIT`, `min_blocks_TVEIT`, `max_blocks_TVEIT`, `tol_alpha_TVEIT`, `tol_damage_TVEIT` | same | — | Optional integration settings (defaults `1001`, `400`, `60000`, `0.01`, `0.001`) |

## Calibration

`TVEIT2025` is not fitted by the scatter objective of `multifat-calibrate`. When no constants are
supplied, `multifat-pipeline` calibrates `A_TVEIT`, `C_TVEIT`, `K_TVEIT` and `L_TVEIT` by the
procedure of the authors (Tveit et al. 2025, Sect. 4; [`eval/tveit2025_calibration.py`](../../src/multifat/eval/tveit2025_calibration.py)),
on any fatigue-life workbook. The model provides it through the pipeline's model-calibration
extension point (`calibrated_params`, `calibrate_constants`); `--no-calibrate` turns it off.
Values entered by hand take precedence: a `--constants` file, or a Summary column named like a
declared parameter (`A_TVEIT`, `Sigm1_TVEIT`, `npts_TVEIT`, ...). The calibration uses none of
them, so it runs only when none is entered. With all four of `A_TVEIT`, `C_TVEIT`, `K_TVEIT` and
`L_TVEIT` entered, those are used; with only some of the parameters entered, the pipeline does
not calibrate and skips the model with the reason (status `partly_supplied`), because constants
fitted with other inputs than those of the evaluation would not hold. `multifat-run` skips models
without constants, and it does not compute the load error of `TVEIT2025` (its comparison
capabilities set `load_error = False`, `models/comparison.py`).

1. **Fatigue limits.** $f_{-1}$ and $t_{-1}$ are the workbook's own Summary values (R = −1 rows).
2. **Mean-stress constant $A$**, in the order of the authors:
   - a fatigue limit $(\sigma_{m,0}, \sigma_{a,0})$ at a tensile mean stress on the Haigh line
     $\sigma_a + A\,\sigma_m = f_{-1}$: $A = (f_{-1} - \sigma_{a,0})/\sigma_{m,0}$ (Eqs. 22–23), "the
     most reliable approach";
   - else a hand-entered $A$ (for example fitted to the life of a mean-stress test);
   - else the FKM estimate $A = a_M S_u/1000 + b_M$ (Eq. 26) for a declared material type with
     tabulated $a_M, b_M$ (steel: 0.35, −0.1);
   - else Goodman, $A = f_{-1}/S_u$, as the authors do where FKM gives no values.
3. **Three calibration points.** The axial R = −1 failures (shared selector `is_design_test`,
   #77) of the source of the fatigue limits (same DOI, else same authors). LOW and HIGH are the
   lowest and highest stress level above $f_{-1}$; MID is the level closest to
   $\sqrt{\text{LOW}\cdot\text{HIGH}}$; at each level, the replicate whose $\log N$ is closest to
   the level median. The papers prescribe three points of one S–N curve, not how they are chosen.
4. **Fit of $C, K, L$** by Eq. (25), $E = \sum_i (1 - \ln N_i/\ln N_{i,\exp})^2$, with the closed
   form of the periodic uniaxial cycle (Ottosen et al. 2008, Eqs. 23, 24, 28, 29). Uniaxial data do
   not identify $C$, so Eq. (25) is profiled over $C$ on a fixed grid (8 values per decade,
   $C \in [10^{-4}, 10^{3}]$), fitting $K \in [10^{-12}, 10^{-1}]$ and $L \in [0.05, 50]$ at each
   $C$. The authors used MATLAB `fmincon`; the grid gives the same answer on every run.
5. **Choice.** Among the grid points with $E \le 1.05\,E_{\min}$, the one of lowest $E$ whose
   closed form has a unique root at every point (G1) and whose time-integrated life agrees with the
   closed form within 0.10 decade (G2). The $C$ range of the equivalent sets is recorded.
6. **Record.** The constants, their sources, the points, the profile, the gates and the
   confidence go into the sidecar (`model_calibrations`) and an archived record; the calibration
   points are design rows (Flag 0) of `TVEIT2025`.

**Inputs file.** What a workbook cannot hold comes from the optional
`Calibrations/<material>/TVEIT2025_inputs.json`, for example:

```json
{
 "schema": "multifat/tveit2025-calibration-inputs/2",
 "material_id": "Al7075_T6",
 "mean_stress_fatigue_limit": {"sigma_a": 100.0, "R": 0.1, "source": "..."},
 "notes": "..."
}
```

Its other keys are `material_type` (for FKM), `A` (`{"value": ..., "source": ...}`),
`calibration_source` (authors or DOI of the calibration tests) and `numerical_settings`. Without
a file, $A$ is the Goodman estimate.

**Reuse.** A calibration takes from seconds to several minutes, because gate G2 integrates the
model over many cycles. A record in `Calibrations/<material>/` or in the run's calibration archive
with the same workbook bytes, overlay, inputs file, procedure version and implementation revision
is reused instead of repeated (the protocol allows identical deterministic reruns to reuse a
record).

**Confidence.** Each calibration lists what weakens it, each finding major or minor, with what to
try; the level is *high* without findings, *medium* with minor findings only and *low* with a major
one. A major finding means the calibration probably misrepresents the material (inconsistent data,
a degenerate fit); a minor one is a known limitation or an estimate the authors themselves use, so
that the level separates the calibrations with a real problem. The confidence is printed during
the run, stored in the record (judged again with the current rules when a record is reused) and
noted on the run-report row. The findings only describe the calibration inputs and the fit, never
the predictions:

| Finding | Severity | What to try |
|---|---|---|
| `failures_below_limit`: axial R = −1 failures below $f_{-1}$ before its reference life (the Summary "Number of Cycles"; unknown counts as before) | major | check $f_{-1}$ against the failures and run-outs of the source |
| `late_failures_below_limit`: failures below $f_{-1}$ after its reference life: they agree with a limit defined at that life, but $f_{-1}$ is a true threshold of the model | minor | for a material without a true fatigue limit, a limit at a longer reference life |
| `runouts_above_limit`: axial R = −1 run-outs above $f_{-1}$ | minor | the limit may be conservative |
| `limit_source_differs`: no test from the source of the limits | major | use the limits of the source of the tests |
| `limit_other_publication`: limits from another paper of the same authors | minor | confirm the same batch |
| `A_estimated`: $A$ from FKM or Goodman | minor | a fatigue limit at a tensile mean stress |
| `material_type_for_FKM`: Goodman because no material type is declared | minor | declare `material_type` |
| `A_supplied`: $A$ entered by hand | minor | a fatigue limit at a tensile mean stress |
| `point_not_R_minus_1`: a calibration point with a mean stress | major | use R = −1 failures |
| `narrow_stress_range`: HIGH/LOW below 1.10 | major | failures over a wider stress range |
| `low_point_near_limit`: LOW within 5 % of $f_{-1}$ | minor | a LOW level further above the limit |
| `single_specimen_level`: a chosen level with one specimen | minor | more replicates |
| `points_mixed_sources`: points from several sources | major | calibrate on one source |
| `poor_fit`: a point missed by more than 0.3 decade | minor | check the limit and the points |
| `constant_at_bound`: $K$ or $L$ on its search bound, or $C$ when the equivalent sets pin it | major | more informative tests |
| `C_not_identified`: $C$ over more than one decade among equivalent sets | minor | a mean-shear or 90° out-of-phase test |
| `point_at_block_limit`: the time integration of a calibration point reached `max_blocks_TVEIT`, so its life is extrapolated (gate G2 still checks it against the closed form) | minor | raise `max_blocks_TVEIT` in `numerical_settings` |

For example, Ti6Al4V_Meneghetti (its axial R = −1 failures lie between 523 and 551 MPa):

```text
TVEIT2025: calibration confidence low (2 major, 4 minor: constant_at_bound, narrow_stress_range,
  A_estimated, C_not_identified, low_point_near_limit, single_specimen_level)
```

**Calibration guide.** A step-by-step guide in LaTeX, with worked examples from the calibrations of
this repository and their confidence, is in
[`docs/calibration-guides/TVEIT2025/`](../calibration-guides/TVEIT2025/README.md).

## Worked example

Evaluate $90^\circ$ out-of-phase bending–torsion ($\sigma_a = 660\,\text{MPa}$, $\tau_a = 340\,\text{MPa}$) using the Stage 1 constants for 42CrMo4 Q&T:

```python
import numpy as np
from multifat.models.invariant_methods.tveit2025 import Tveit2025

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {
    "Sigm1_TVEIT": 433.0, "Taum1_TVEIT": 344.0,
    "A_TVEIT": 0.245625, "C_TVEIT": 0.0040672184319305,
    "K_TVEIT": 1.91047799401683e-05, "L_TVEIT": 10.1293188952565,
}
res = Tveit2025().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params=p,
    harmonics=dict(S_mean=Z, S_sin=bt(660.0, 0.0), S_cos=bt(0.0, 340.0)),
)
assert np.isclose(res.aux["Nf"], 263.0, atol=1.0)
print(f"index = {res.values['']:.4f}, Nf = {res.aux['Nf']:.0f}")
```

Output:

```text
index = 1.8986, Nf = 263
```

Batch CLI evaluation on 42CrMo4 Q&T with Stage 1 constants:

```bash
multifat-pipeline Experimental_Data/6_42CrMo4_Santus.xlsx --models TVEIT2025 --no-calibrate --constants Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json
```

## Comparison capabilities and outputs

`TVEIT2025` sets `ComparisonCapabilities(life_prediction="native", dp_comparison=False, dp_reference="none")`: only native life predictions are compared against experiments (`aux["Nf"]`), and no universal-curve DP comparison is performed because $(C, K, L)$ already encapsulate the damage-life law.

## Verification and limitations

- **Software verification (paper reproduction).** [`tests/methods/test_tveit2025.py`](../../tests/methods/test_tveit2025.py) reproduces all published $S$–$N$ curves in Figs. 5–6 across Table 1 materials within $1\%$ in stress, Eq. (24) torsional limit, infinite life inside the endurance surface and after mean-shear shakedown, and constants-file parameter overrides.
- **Calibration.** [`tests/methods/test_tveit2025_calibration.py`](../../tests/methods/test_tveit2025_calibration.py) checks Ottosen et al. (2008) Eq. (31) ($A = 0.225$ for SAE 4340), the closed form against an independent implementation (7 SAE 4340 cases) and against the time integration, the routes to $A$, the point rule, the source match and the confidence findings; with the stage-1 inputs, the profile reaches the $E_{\min}$ of the stage-1 calibrations. [`tests/test_tveit2025_calibration_pipeline.py`](../../tests/test_tveit2025_calibration_pipeline.py) runs it through the pipeline on a synthetic workbook with an exact fit.
- **Computational cost.** Finite-life evaluation integrates block by block until failure or steady state (up to 60,000 blocks unless bounded); a calibration repeats it at the calibration points for each equivalent set checked by G2, which is why identical calibrations are reused.
- **Identifiability.** Uniaxial data do not identify $C$: the equivalent sets span three to four decades of $C$ with the same uniaxial S–N curve. The chosen $C$ is set by the lowest-$E$ rule, not by the data; $C$ matters under mean shear stress and non-proportional loading.

## References

- Tveit, Ø., Ottosen, N. S., Erduran, E., & Reyes, A. (2025). Continuum mechanics fatigue model for arbitrary ratios of shear and normal stress fatigue strengths. *Mechanics of Materials*, 209, 105452. [doi:10.1016/j.mechmat.2025.105452](https://doi.org/10.1016/j.mechmat.2025.105452)
- Ottosen, N. S., Stenström, R., & Ristinmaa, M. (2008). Continuum approach to high-cycle fatigue modeling. *International Journal of Fatigue*, 30, 996–1006. [doi:10.1016/j.ijfatigue.2007.08.009](https://doi.org/10.1016/j.ijfatigue.2007.08.009)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
