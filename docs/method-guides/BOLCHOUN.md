# Substitute amplitude and out-of-phase damage (`BOLCHOUN`)

> **Registry name:** `BOLCHOUN` · **Catalogue entry:** [N53](../../METHODS.md#n53) ·
> **Source:** Bolchoun, Baumgartner & Kaufmann (2017), <https://doi.org/10.1016/j.ijfatigue.2017.01.012>
> (accepted manuscript) ·
> **Implementation:** [`integral_methods/bolchoun.py`](../../src/multifat/models/integral_methods/bolchoun.py)

## Overview

- **Type and family:** integral approach with a native life law. Every cycle is expressed by one
  substitute stress amplitude, the integral of the normal stress amplitudes over all planes. Its
  damage is read on a Wöhler curve interpolated between the axial and torsion curves, and an
  extra out-of-phase damage, calibrated on one experimental out-of-phase curve, is added.
- **Idea:** separate the in-phase damage, which depends on how torsional the loading is, from the
  effect of rotating shear vectors, measured by the area they sweep over all planes ({OP}$).
- **Proposed for:** variable-amplitude, in-phase and out-of-phase loading of thin-walled welded
  magnesium tubes (AZ31, AZ61), assessed with the notch stress approach ({ref}$ = 0.05 mm).
- **Status:** constant-amplitude branch, one harmonic, fully reversed loading.

## Original formulation

**Out-of-phase measure** (Eq. 9). $\boldsymbol\tau_{n,0}(t)$ is the shear stress vector on the plane
of normal $\mathbf n$ minus its mean. The integral runs over the unit hemisphere
(\mathbf n = \cos\psi\,d\theta\,d\psi$) and over the cycle:

```math
F_{OP} = \frac12\int_{t_1}^{t_2}\!\!\iint_{\mathbf n}
\left\|\boldsymbol\tau_{n,0}(t)\times\dot{\boldsymbol\tau}_{n,0}(t)\right\|\,d\mathbf n\,dt .
```

**Wöhler curve interpolation** (Eq. 11). With $\tau_{max,a}$ and $\sigma_{max,a}$ the largest shear
and normal stress amplitudes (Maximum Variance Method),

```math
f_I = \frac{\tau_{max,a}}{\sigma_{max,a}},\qquad \tfrac12 \le f_I \le 1,\qquad
X = 2\left[X_{ax}(1 - f_I) + X_{tors}\left(f_I - \tfrac12\right)\right],\quad X = k,\ N_k,\ L_k .
```

 = 1/2$ in pure tension and 1 in pure torsion; below 1/2 the axial curve is used.

**Substitute stress amplitude** (Eqs. 12–13):

```math
\sigma_{a,s} = \iint_{\mathbf n}\sigma_{n,a}\,d\mathbf n,\qquad \sigma_{n,a} = \sqrt{2\,\mathrm{Var}(\sigma_n(t))} .
```

**Out-of-phase damage** (Eq. 14), on the experimental out-of-phase curve (total damage) and on
the in-phase curve interpolated at its {I,OP}$:

```math
D_{OP}\big(F_{OP}(\sigma_{a,s})\big) = \sqrt{D_{total}^2(\sigma_{a,s}) - D_{IP}^2(\sigma_{a,s})} .
```

**Damage of a cycle** (Eq. 15), and the variant of Section 8 when out-of-phase loading is less
damaging than in-phase loading:

```math
D_{total,C} = \sqrt{D_{IP,C}^2 + \left(D_{OP,C}^{\,k_{I,C}/k_{I,OP}}\right)^2},\qquad
D_{total,C} = \sqrt{D_{IP,C}^2 - \left(D_{OP,C}^{\,k_{I,C}/k_{I,OP}}\right)^2} .
```

The damages are summed by Palmgren–Miner with Haibach's slope k - 2$ below the knee, and
 = L_S\,D_{al}/D_{total}$.

## Repository implementation

For one harmonic, $\sigma_{n,a} = \sqrt{(\mathbf n\cdot\mathbf S_{\sin}\mathbf n)^2 +
(\mathbf n\cdot\mathbf S_{\cos}\mathbf n)^2}$ and the shear vector traces an ellipse. The integral
of Eq. (9) over one period is the area of that ellipse, so

```math
F_{OP} = \iint_{\mathbf n}\pi\left\|\boldsymbol\tau_s(\mathbf n)\times\boldsymbol\tau_c(\mathbf n)\right\|d\mathbf n .
```

The hemisphere integrals use a 48 × 96 grid (\pi$ times the sphere average). $\tau_{max,a}$ is
the major semi-axis of the shear ellipse, the largest resolved shear amplitude the Maximum Variance
Method finds; both maxima are refined locally on the sphere. The life of a cycle is
 = 1/D_{total,C}$. The model reports $\sigma_{a,s}$ as its index.

**Repository choices**

- *Constant amplitude.* One cycle per cycle: no rainflow counting, no Miner sum, no Haibach slope,
  {al} = 1$.
- *Curves.* The axial and torsion curves are the pipeline's Basquin fits on the R = −1 failures,
  $\sigma_a = \sigma_f(2N)^{b}$ and $\tau_a = \tau_f(2N)^{b_0}$, converted to $\sigma_{a,s}$ by the
  constant factors of pure tension (\pi/3$) and pure torsion (/3$). There is no knee point:
  Eq. (11) interpolates  = -1/b$ and the amplitude $ at  = 10^6$ (maintainer decision,
  22 November 2026).
- *Out-of-phase curve.* The Basquin fit on the most populated fully reversed out-of-phase
  tension–torsion path of the dataset (`sigma_f_OP`, `b_OP`, axial amplitude), with its mean
  ratio `lambda_OP` $= \tau_a/\sigma_a$ and phase `delta_OP`. Its tests become design points
  (Flag 0); without an identifiable path (5 tests, 1 decade, slope distinct from 0) the model is
  skipped. {OP}(F)$ is evaluated exactly on that path at the amplitude with the same {OP}$
  (maintainer decision).
- *Sign of Eq. (15).* When the out-of-phase curve lies above the interpolated in-phase curve
  ({total} < D_{IP}$), the minus-sign variant of Section 8 is used, decided at the amplitude of
  each cycle (maintainer decision). A non-positive {total,C}^2$ is an infinite life.
- *Conditioning.* On the calibration path {I,C} = k_{I,OP}$ and Eq. (15) returns the
  out-of-phase curve exactly. With the minus sign and an out-of-phase curve far above the in-phase
  one, {IP,C}$ and {OP,C}$ nearly cancel, and a cycle whose $ differs slightly from
  {I,OP}$ can reach {total,C}^2 \le 0$. This is reported, not regularised.
- *Mean stresses.* Out of scope, as in the paper (Section 8 leaves them open): a case with a mean
  stress raises `LoadingScopeError` (maintainer decision).
- *Interpolation range.*  < 1/2$ uses the axial curve (paper);  > 1$, which the paper
  assumes not to exist, uses the torsion curve.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_f$, $ | `sigma_f`, `b_fat` | MPa, – | workbook or Basquin fit on R = −1 tension | yes |
| $\tau_f$, $ | `tau_f`, `b0` | MPa, – | workbook or Basquin fit on R = −1 torsion | yes |
| $\sigma_{f,OP}$, {OP}$ | `sigma_f_OP`, `b_OP` | MPa, – | Basquin fit on the out-of-phase path | yes |
| $\lambda$, $\delta$ | `lambda_OP`, `delta_OP` | –, degrees | that path | yes |
| $\boldsymbol\sigma(t)$ | `S_mean` (zero), `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** Axial, torsional and one out-of-phase Wöhler curve, all in
   $\sigma_{a,s}$, from experiments, literature or guidelines; {OP}(F_{OP})$ follows from
   Eq. (14) on the out-of-phase curve.
2. **Repository.** `multifat-pipeline` fills `sigma_f`, `b_fat`, `tau_f`, `b0` (workbook
   values first) and fits the out-of-phase curve (`eval/material_fits.out_of_phase_fill`). It runs
   on fully reversed datasets with an identifiable out-of-phase path, for example
   `multifat-pipeline --material 42CrMo4_Reis --models BOLCHOUN`.

## Worked example

An illustrative steel; the out-of-phase curve is fitted on a 90° path with
$\tau_a/\sigma_a = 0.577$ and lies above the in-phase curve.

```python
import numpy as np

from multifat.models.integral_methods.bolchoun import (
    WoehlerSet,
    bolchoun_quantities,
    tension_torsion_path,
)

# Illustrative steel: axial and torsion Basquin curves, and an out-of-phase curve fitted on a
# 90 deg path with tau_a / sigma_a = 0.577 (sigma_a = sigma_f_OP (2N)^b_OP, axial amplitude)
params = {
    "sigma_f": 1154.0, "b_fat": -0.061, "tau_f": 940.0, "b0": -0.065,
    "sigma_f_OP": 773.0, "b_OP": -0.036, "lambda_OP": 0.577, "delta_OP": 90.0,
}  # fmt: skip
curves = WoehlerSet(params)
zero = np.zeros((3, 3))

for lam, delta in ((0.0, 0.0), (0.577, 0.0), (0.577, 90.0)):
    S_sin, S_cos = tension_torsion_path(500.0, lam, delta)
    q = bolchoun_quantities(zero, S_sin, S_cos)
    life = curves.life(q["sigma_as"], q["F_OP"], q["f_I"])
    print(f"lambda = {lam:.3f}, delta = {delta:2.0f}: sigma_as = {q['sigma_as']:7.1f} MPa, "
          f"F_OP = {q['F_OP']:8.0f} MPa^2, f_I = {q['f_I']:.3f}, N = {life:9.0f}")  # fmt: skip

# torsion alone: the torsion Basquin curve
S_sin = np.zeros((3, 3))
S_sin[0, 1] = S_sin[1, 0] = 500.0
q = bolchoun_quantities(zero, S_sin, zero)
print(f"torsion tau_a = 500: N = {curves.life(q['sigma_as'], q['F_OP'], q['f_I']):.0f}, "
      f"Basquin {0.5 * (500.0 / 940.0) ** (1 / -0.065):.0f}")  # fmt: skip
# the out-of-phase path itself returns its own curve
print(f"out-of-phase curve at sigma_a = 500: N = {0.5 * (500.0 / 773.0) ** (1 / -0.036):.0f}")
```

Output:

```text
lambda = 0.000, delta =  0: sigma_as =  1047.2 MPa, F_OP =        0 MPa^2, f_I = 0.500, N =    450459
lambda = 0.577, delta =  0: sigma_as =  1244.5 MPa, F_OP =        0 MPa^2, f_I = 0.604, N =     25784
lambda = 0.577, delta = 90: sigma_as =  1403.2 MPa, F_OP =   356172 MPa^2, f_I = 0.577, N =     90113
torsion tau_a = 500: N = 8256, Basquin 8256
out-of-phase curve at sigma_a = 500: N = 90113
```

Pure tension has  = 1/2$ and no out-of-phase measure: it returns the axial curve, and pure
torsion returns the torsion curve. The in-phase path has  = r/(\sigma_a/2 + r)$ with
 = \sqrt{(\sigma_a/2)^2 + \tau_a^2}$. The out-of-phase path, at the same amplitudes, has a larger
$\sigma_{a,s}$ but a longer life: its curve lies above the in-phase one, so the minus-sign variant
applies, and on its own path the model returns exactly the fitted out-of-phase curve.

## Comparison capabilities and outputs

- **Life:** native (`Nf_native`, `Error_log10_native`); the out-of-phase calibration tests and the
  R = −1 tension and torsion failures are design points (Flag 0).
- **DP comparison:** not supported (`dp_comparison = False`); $\sigma_{a,s}$ is reported in
  `CP_value`, but the life depends also on {OP}$ and $.
- **Load histories:** native life law; refused by the generic accumulation of
  [load histories](../load-histories.md).
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_bolchoun.py`](../../tests/methods/test_bolchoun.py) checks:

- Eq. (12): the unit factors \pi/3$ (tension) and /3$ (torsion), and linearity;
- Eq. (9): the ellipse area equals the time integral on random planes, vanishes for proportional
  loading and scales with the square of the amplitude;
- $ = 1/2 and 1 at the end points, and its closed form for proportional loading;
- Eq. (11): the end points are the axial and torsion curves; clipping outside 1$;
- uniaxial tension and torsion return their Basquin curves;
- Eqs. (14)–(15): the out-of-phase path returns its own curve, with either sign;
- the out-of-phase fit (largest path, true mean ratio), the scope and missing-constant errors.

Not reproduced: the paper's results are figures for notched magnesium welds under variable
amplitude, with no tabulated per-test value for a reproducible stress state.

Limitations:

- constant amplitude only; no knee point;
- needs an out-of-phase path in the dataset, whose tests are then not predictions;
- the minus-sign variant is ill-conditioned when the out-of-phase curve lies far above the
  in-phase one;
- no mean-stress treatment.

## References

- A. Bolchoun, J. Baumgartner, H. Kaufmann (2017), A new method for fatigue life evaluation under
  out-of-phase variable amplitude loadings and its application to thin-walled magnesium welds,
  *International Journal of Fatigue* (accepted manuscript).
  <https://doi.org/10.1016/j.ijfatigue.2017.01.012>
