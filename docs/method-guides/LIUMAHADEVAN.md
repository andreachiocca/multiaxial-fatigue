# Liu–Mahadevan (`LIUMAHADEVAN`)

> **Registry name:** `LIUMAHADEVAN` · **Catalogue entry:** [N16](../../METHODS.md#n16) ·
> **Source:** Liu & Mahadevan (2005), <https://doi.org/10.1016/j.ijfatigue.2005.01.003> ·
> **Implementation:** [`cp_methods/liu_mahadevan.py`](../../src/multifat/models/cp_methods/liu_mahadevan.py)

## Overview

- **Type and family:** stress-based critical-plane criterion with its own life model.
- **Idea:**
  - The fatigue fracture plane is the plane of maximum normal stress amplitude.
  - The critical plane lies at a material angle $\alpha$ from it, chosen so that the damage
    parameter meets both fully reversed limits while the hydrostatic term is as small as
    possible. The critical plane therefore depends on the material, not only on the stress state.
  - The damage combines the normal, shear and hydrostatic stress amplitudes.
- **Mean stress:** a normal mean on the critical plane multiplies the normal amplitude by
  $1 + \eta\,\sigma_{m,c}/f_{-1}$. Mean shear has no effect.
- **Proposed for:** high-cycle fatigue limits and lives of metals from ductile to extremely brittle
  ($t_{-1}/f_{-1}$ from 0.57 to 1.28), under in-phase and out-of-phase loading with mean stresses.
- **Status:** implemented with its native life law. Eqs. (6)–(25) were checked on the paper. The
  published fatigue-limit comparison (Figs. 2–5) is reproduced only in part; see
  [Verification](#verification-and-limitations).

## Original formulation

Liu & Mahadevan (2005, Eqs. 20–21 and Table 1), with $s = t_{-1}/f_{-1}$:

```math
\sqrt{\left[\frac{\sigma_{a,c}\left(1+\eta\,\sigma_{m,c}/f_{-1}\right)}{f_{-1}}\right]^2
+ \left(\frac{\tau_{a,c}}{t_{-1}}\right)^2 + k\left(\frac{\sigma^H_{a}}{f_{-1}}\right)^2} = \beta,
```

```math
s \le 1:\quad \cos 2\alpha = \frac{-2+\sqrt{4-4(1/s^2-3)(5-1/s^2-4s^2)}}{2(5-1/s^2-4s^2)},\quad
\beta = \sqrt{\cos^2(2\alpha)\,s^2+\sin^2(2\alpha)},\quad k = 0;
```

```math
s > 1:\quad \alpha = 0,\quad \beta = s,\quad k = 9(s^2-1);\qquad
\eta = \frac34 + \frac14\,\frac{\sqrt3 - f_{-1}/t_{-1}}{\sqrt3-1}\ (s\le1),\quad \eta = 1\ (s>1).
```

Here:

- $\sigma_{a,c}$, $\tau_{a,c}$ and $\sigma_{m,c}$ act on the critical plane;
- $\sigma^H_a$ is the hydrostatic stress amplitude.

**Life (Eqs. 23–24).** Replace the limits by the S–N strengths $f_{N}$ and $t_{N}$ at life $N_f$:

```math
\frac1\beta\sqrt{\left[\sigma_{a,c}\left(1+\eta_{N}\frac{\sigma_{m,c}}{f_{N}}\right)\right]^2
+ \left(\frac{f_N}{t_N}\right)^2\tau_{a,c}^2 + k\,(\sigma^H_a)^2} = f_N,
```

with $\alpha$, $\beta$, $k$ and $\eta$ evaluated at $s = t_N/f_N$, solved for $N_f$ by trial and
error. Any S–N form may be used.

**Error index (Eq. 25).** $I = (\text{LHS} - \text{RHS})/\text{LHS}$.

## Repository implementation

- `""` value: the fatigue-limit index $\sqrt{\cdots}/\beta$ of Eq. (20), 1 on the limit. It is a
  diagnostic only.
- Life: the native $N_f$ of Eq. (23), with $f_N = \sigma_f'(2N)^{b}$ and $t_N = \tau_f'(2N)^{b_0}$,
  solved by bisection in $\log N$. The comparison uses the native life
  (`dp_comparison = False`).

**Repository choices**

- *Plane rule in 3D (maintainer decision).* The paper rotates the plane in 2D. Here the critical
  plane is the worst plane on the cone of planes at angle $\alpha$ around the fracture-plane
  normal. This meets both reversed limits exactly for $s \ge 0.5$.
  - Below $s = 0.5$ (under the Tresca ratio) a tilted plane of the cone exceeds the torsion
    limit, so such materials raise `ValueError`.
  - Fracture planes of equal amplitude: the worst is kept.
- *Shear amplitude (maintainer decision).* $\tau_{a,c}$ is the rss of the harmonic shear
  components, chosen against Fig. 2: it matches hard steel at 0°, 60° and 90°, whereas the longest
  chord gives 4.6 % at 90°.
- *Eq. (19).* It is evaluated in the rationalised form
  $\cos2\alpha = -C/(1+\sqrt{1-AC})$, which removes the 0/0 at $s = 1/2$ and $s = 1$.
- *Mean stress.* $\eta$ from Eq. (21); experimental $\eta$ values are not used. The factor
  $1 + \eta\sigma_{m,c}/f$ is clipped at zero.
- *Life law (maintainer decision).* The native Eq. (23) with Basquin curves. The bisection is kept
  to lives where $t_N/f_N \ge 0.5$: a solution outside raises `LoadingScopeError`, and no solution
  below $10^{12}$ cycles gives an infinite life.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $f_{-1}$, $t_{-1}$ | `Sigm1`, `Taum1` | MPa | workbook | yes ($t_{-1}/f_{-1} \ge 0.5$) |
| $\sigma_f'$, $b$ | `sigma_f`, `b_fat` | MPa, – | workbook or fitted on the R = −1 tension failures | yes |
| $\tau_f'$, $b_0$ | `tau_f`, `b0` | MPa, – | workbook or fitted on the R = −1 torsion failures | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.**
   - $\alpha$, $\beta$ and $k$ follow in closed form from $t_{-1}/f_{-1}$.
   - $\eta$ is calibrated on uniaxial mean-stress tests, or taken from Eq. (21).
   - The life uses the material's uniaxial and torsional S–N curves.
2. **Repository.**
   - The same closed forms and Eq. (21).
   - The S–N curves are the Basquin constants of the workbook, or fitted by the pipeline on its
     $R = -1$ tension and torsion failures (`material_fits.basquin_fill`).
   - No master curve is fitted.

## Worked example

Hard steel of Papadopoulos et al. (1997), then the native life with illustrative curves.

```python
import numpy as np

from multifat.models.cp_methods.liu_mahadevan import (
    liu_mahadevan_index,
    liu_mahadevan_life,
    material_angle,
)

zero = np.zeros((3, 3))


def harmonic(sigma_a, tau_a, delta_deg, sigma_m=0.0):
    """sigma_x = sigma_m + sigma_a sin wt, tau_xy = tau_a sin(wt - delta)."""
    d = np.radians(delta_deg)
    S_mean = np.diag([sigma_m, 0.0, 0.0])
    S_sin = np.array([[sigma_a, tau_a * np.cos(d), 0], [tau_a * np.cos(d), 0, 0], [0, 0, 0]])
    S_cos = np.array([[0, -tau_a * np.sin(d), 0], [-tau_a * np.sin(d), 0, 0], [0, 0, 0]])
    return S_mean, S_sin, S_cos


# Hard steel of Papadopoulos et al. (1997), the paper's Table 2: f_-1 = 313.9, t_-1 = 196.2 MPa
f, t = 313.9, 196.2
alpha, beta, k = material_angle(t / f)
print(f"s = {t / f:.3f}: alpha = {np.degrees(alpha):.1f} deg, beta = {beta:.4f}, k = {k:.0f}")
for delta in (0, 90):  # Papadopoulos et al. (1997) Table 1, rows 1 and 4
    sigma_a, tau_a = (138.1, 167.1) if delta == 0 else (150.2, 181.7)
    index = liu_mahadevan_index(*harmonic(sigma_a, tau_a, delta), f, t)
    print(f"delta = {delta:2d} deg: index = {index:.3f}, I = {100 * (1 - 1 / index):5.1f} % (Eq. 25)")

# Native life, Eq. (23), with illustrative Basquin curves (2N convention)
sn = dict(sigma_f=1200.0, b_fat=-0.08, tau_f=700.0, b0=-0.06)
for delta in (0, 90):
    Nf = liu_mahadevan_life(*harmonic(400.0, 230.0, delta), **sn)
    print(f"sigma_a = 400, tau_a = 230 MPa, delta = {delta:2d} deg: N_f = {Nf:.3g}")
```

Output:

```text
s = 0.625: alpha = 39.2 deg, beta = 0.9875, k = 0
delta =  0 deg: index = 0.975, I =  -2.5 % (Eq. 25)
delta = 90 deg: index = 1.056, I =   5.3 % (Eq. 25)
sigma_a = 400, tau_a = 230 MPa, delta =  0 deg: N_f = 1.29e+04
sigma_a = 400, tau_a = 230 MPa, delta = 90 deg: N_f = 3.63e+04
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_Santus --models LIUMAHADEVAN --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** native (`Nf_expected` = Eq. 23). Load errors are found by bisection on the native life.
- **DP comparison:** not supported. The index is a fatigue-limit diagnostic, and the life-level
  equivalent stress depends on $N_f$.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_liu_mahadevan.py`](../../tests/methods/test_liu_mahadevan.py):

- **Eq. (19) / Table 1.** Values including the removable singular points; both reversed limits on
  the $\alpha$-cone for $0.5 \le s \le 1.28$.
- **Eq. (23).** It returns the input life on both S–N curves.
- **Figs. 2–5.** These are the mean absolute error per phase on the Papadopoulos et al. (1997)
  Table 1 data; the figures' Matake and McDiarmid curves equal that table's phase means, which
  confirms the data set.
  - Reproduced: hard steel at 0°, 60° and 90°, and 30NCD16 at 60°.
  - Not reproduced (strict xfails, unresolved): hard steel at 30°, 42CrMo4, 34Cr4 and 30NCD16.
    The balanced 90° out-of-phase rows come out about 10 % non-conservative here and small in
    the paper. No plane-rule or amplitude variant tried closes the gap.
- **Not reproducible.** Table 4, the fracture angles, needs Macha's loads. Figs. 10–13, the lives,
  need data sets not given in the paper.

Limitations:

- partial reproduction of the published comparison;
- 3D plane rule, shear amplitude and S–N form are repository choices;
- $t_{-1}/f_{-1} \ge 0.5$;
- slow load-error search (the life is re-solved at each trial load);
- a load error is undefined (NaN) when a trial load moves the life where $t_N/f_N < 0.5$: on
  Ti6Al4V_Berto (torsion slope $-0.035$ against $-0.11$ in tension) this happens for 17 of the 31
  assessed failures.

## References

- Y. Liu, S. Mahadevan (2005), Multiaxial high-cycle fatigue criterion and life prediction for
  metals, *International Journal of Fatigue* 27, 790–800.
  <https://doi.org/10.1016/j.ijfatigue.2005.01.003>
- I. V. Papadopoulos, P. Davoli, C. Gorla, M. Filippini, A. Bernasconi (1997), A comparative study
  of multiaxial high-cycle fatigue criteria for metals, *International Journal of Fatigue* 19(3),
  219–235. <https://doi.org/10.1016/S0142-1123(96)00064-3>
