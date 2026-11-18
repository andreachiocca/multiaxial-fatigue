# Manson–McKnight–Papuga (`MMP`)

> **Registry name:** `MMP` · **Catalogue entry:** [G05](../../METHODS.md#g05) ·
> **Source:** Papuga & Fojtík (2017), <https://doi.org/10.1016/j.ijfatigue.2017.07.001> ·
> **Implementation:** [`invariant_methods/mmp.py`](../../src/multifat/models/invariant_methods/mmp.py)

## Overview

- **Type and family:** stress-based equivalent-stress criterion for fatigue limits (Manson–McKnight
  family).
- **Idea:** a von Mises-type equivalent stress is built separately from the component amplitudes
  and from the component means. The shear terms are weighted so that the torsion limits are met.
  The two are combined with a Walker mean-stress correction.
- **Mean stress:** Walker exponent $w$ from the repeated tension limit. A factor $\beta_P$ from the
  principal stress extremes gives the mean its sign and halves it in repeated torsion.
- **Proposed for:** a quick fatigue-limit estimate that "can be performed in an ordinary
  spreadsheet editor", validated on CSN 41 1523 structural steel (tension, bending, torsion,
  internal and external pressure).
- **Status:** implemented as MMP_WP, the paper's final variant. Eqs. (19)–(27) were checked on the
  paper. All 24 load cases of Table 11 are reproduced within 0.1 point.

## Original formulation

Papuga & Fojtík (2017, Eqs. 19–25):

```math
\sigma_a = \sqrt{\tfrac12\left[(\sigma_{x,a}-\sigma_{y,a})^2 + (\sigma_{y,a}-\sigma_{z,a})^2 + (\sigma_{z,a}-\sigma_{x,a})^2
  + 2\kappa^2\left(\tau_{xy,a}^2+\tau_{yz,a}^2+\tau_{zx,a}^2\right)\right]},
```

```math
\sigma_m = \sqrt{\tfrac12\left[(\sigma_{x,m}-\sigma_{y,m})^2 + (\sigma_{y,m}-\sigma_{z,m})^2 + (\sigma_{z,m}-\sigma_{x,m})^2
  + 2X_m^2\left(\tau_{xy,m}^2+\tau_{yz,m}^2+\tau_{zx,m}^2\right)\right]},
```

```math
D_p = \sigma_a^{\,w}\left(\sigma_a + \beta_P\,\sigma_m\right)^{1-w} \le f_{-1},
\qquad
\beta_P = \begin{cases}
\dfrac{\sigma_{1,\max}}{\sigma_{1,\max}-\sigma_{3,\min}} & |\sigma_{1,\max}| \ge |\sigma_{3,\min}| \\[2mm]
\dfrac{\sigma_{3,\min}}{\sigma_{1,\max}-\sigma_{3,\min}} & \text{otherwise}
\end{cases}
```

```math
\kappa = \frac{f_{-1}}{t_{-1}}, \qquad
w = \frac{\log(f_0/f_{-1})}{\log 2}, \qquad
X_m = 2\kappa\left[\left(\frac{2t_{-1}}{t_0}\right)^{\frac{1}{1-w}} - 1\right],
```

where:

- $f_{-1}$, $t_{-1}$ are the fully reversed axial and torsion limits;
- $f_0$, $t_0$ are the repeated ($R = 0$) limits, as the maximum stress of the cycle;
- $\sigma_{1,\max}$, $\sigma_{3,\min}$ are the largest and smallest principal stresses over the
  cycle.

**Calibration conditions.** The four limits used as input are met exactly:
- fully reversed tension and torsion, through $\kappa$;
- repeated tension ($\beta_P = 1$), through $w$;
- repeated torsion ($\beta_P = 1/2$), through $X_m$.

**Error index.** $\Delta FI = 100\,(D_p - f_{-1})/f_{-1}$ (Eq. 27).

**Phase.** The components are processed one by one, so out-of-phase and in-phase loading with the
same amplitudes give the same result. The paper names this as the method's weak point.

## Repository implementation

$\mathrm{DP} = D_p/f_{-1}$, the fatigue index (1 on the fatigue limit, $\Delta FI = 100(\mathrm{DP}-1)$).
$\sigma_{1,\max}$ and $\sigma_{3,\min}$ are searched over the harmonic cycle.

**Repository choices**

- *Repeated limits (maintainer decision).* No supplied workbook gives $f_0$ or $t_0$.
  - A missing `Sig0` is estimated from `Su` and $f_{-1}$ by Böhme et al. (2026) Eq. (1), as for
    `PCR` and `MMMH`.
  - A missing `Tau0` follows from $4t_{-1}/t_0 - 2f_{-1}/f_0 = 1$ (Zenner), quoted by the paper as
    its Eq. (23).
- *Signed amplitudes (maintainer decision).* With harmonics, each component amplitude is
  $\sqrt{s^2 + c^2}$ of its in-phase and 90° parts, with the sign of the in-phase part. This is
  exact for proportional loading and keeps the criterion phase-blind, as published.
- *Negative bracket.* If $\sigma_a + \beta_P\sigma_m \le 0$ (a large compressive mean), Eq. (20)
  is undefined; DP = 0.
- *Life.* The index is the damage parameter of the shared master curve $\mathrm{DP} = A N_f^b$, an
  extension of the fatigue-limit criterion
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $f_{-1}$ | `Sigm1` | MPa | workbook | yes |
| $t_{-1}$ | `Taum1` | MPa | workbook | yes |
| $f_0$ | `Sig0` | MPa | workbook, else estimated from `Su` | `Sig0` or `Su` |
| $t_0$ | `Tau0` | MPa | workbook, else Eq. (23) | no |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** Nothing is fitted. $\kappa$, $w$ and $X_m$ follow in closed form from
   the four limits $f_{-1}$, $t_{-1}$, $f_0$, $t_0$ (Table 10 of the paper). Plane-bending cases use
   the bending limits $b_{-1}$, $b_0$ in place of $f_{-1}$, $f_0$.
2. **Repository.** The same, with the workbook limits (and the estimates above). Only the
   master-curve constants are fitted, on the uniaxial $R = -1$ design tests.

## Worked example

CSN 41 1523 steel, constants and outer-surface stresses of the paper's Tables 7 and 10.

```python
import numpy as np

from multifat.models.invariant_methods.mmp import mmp_constants, mmp_index

# CSN 41 1523 steel, paper's Table 10: f_-1, t_-1, f_0, t_0 (MPa)
f_m1, t_m1, f_0, t_0 = 244.7, 169.9, 385.8, 330.8
kappa, w, X_m = mmp_constants(f_m1, t_m1, f_0, t_0)
print(f"kappa = {kappa:.3f}, w = {w:.3f} (Eq. 21), X_m = {X_m:.3f} (Eq. 25)")


def tensor(axial, tangential, radial, shear):
    """Stress tensor in cylindrical axes (a, t, r)."""
    return np.array([[axial, shear, 0.0], [shear, tangential, 0.0], [0.0, 0.0, radial]])


# Outer-surface stresses of the paper's Table 7: (mean, in-phase amplitude, 90 deg part)
cases = {
    "FF004 Ten+To 90 deg": (tensor(0, 0, 0, 0), tensor(194.0, 0, 0, 0), tensor(0, 0, 0, 140.6), 14.6),
    "FF005 Ten+To R = 0": (tensor(114.6, 0, 0, 107.5), tensor(114.6, 0, 0, 107.5), tensor(0, 0, 0, 0), -10.6),
    "FF059 Ten+To+RP": (tensor(227.5, -19.5, -40.0, 0), tensor(0, 0, 0, 172.7), tensor(0, 0, 0, 0), 23.7),
}
for label, (S_mean, S_sin, S_cos, printed) in cases.items():
    dFI = 100 * (mmp_index(S_mean, S_sin, S_cos, f_m1, t_m1, f_0, t_0) - 1)
    print(f"{label:20s}: Delta FI = {dFI:5.1f} % (paper {printed:5.1f} %)")
```

Output:

```text
kappa = 1.440, w = 0.657 (Eq. 21), X_m = 0.234 (Eq. 25)
FF004 Ten+To 90 deg : Delta FI =  14.6 % (paper  14.6 %)
FF005 Ten+To R = 0  : Delta FI = -10.6 % (paper -10.6 %)
FF059 Ten+To+RP     : Delta FI =  23.7 % (paper  23.7 %)
```

The three indices match Table 11. $X_m$ is 0.234 against the printed 0.233, within the rounding of
$t_{-1}$ and $t_0$.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_Santus --models MMP --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`), a repository extension.
- **DP comparison:** supported. The index is homogeneous of degree 1 in the load ($\beta_P$ is
  scale-free).
- **Fatigue-limit datasets:** $\mathrm{DP} - 1$ is the paper's $\Delta FI/100$.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_mmp.py`](../../tests/methods/test_mmp.py):

- **Table 11 (MMP_WP).** All 24 load cases are reproduced within 0.1 point of $\Delta FI$. Both
  surfaces are evaluated and the critical one is kept: FF061 is critical on the inner surface, as
  the paper marks.
- **Table 10.** $\kappa$ and $w$ to the printed digits; $X_m$ within 0.002, the effect of the
  0.1 MPa rounding of $t_0$.
- **Other checks:**
  - the four calibration limits;
  - $\beta_P$ (Table 9);
  - phase blindness;
  - the estimates of `Sig0` and `Tau0`.

Limitations:

- phase-blind for non-proportional loading (the paper's own limitation);
- validated by its authors on one steel only;
- estimated $f_0$, $t_0$ on the supplied workbooks;
- a fatigue-limit criterion, used for finite life only through the master curve.

## References

- J. Papuga, F. Fojtík (2017), Multiaxial fatigue strength of common structural steel and the
  response of some estimation methods, *International Journal of Fatigue* 104, 27–42.
  <https://doi.org/10.1016/j.ijfatigue.2017.07.001>
- H. Zenner, R. Simbürger, J. Liu (2000), On the fatigue limit of ductile metals under complex
  multiaxial loading, *International Journal of Fatigue* 22, 137–145.
  <https://doi.org/10.1016/S0142-1123(99)00107-3>
