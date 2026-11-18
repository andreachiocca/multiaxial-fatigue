# Zhang–Yao (`ZHANGYAO`)

> **Registry name:** `ZHANGYAO` (and `ZHANGYAO_ext`) · **Catalogue entry:** [N06](../../METHODS.md#n06) ·
> **Source:** Zhang & Yao (2011), <https://doi.org/10.1111/j.1460-2695.2010.01523.x> ·
> **Implementation:** [`cp_methods/zhang_yao.py`](../../src/multifat/models/cp_methods/zhang_yao.py)

## Overview

- **Type and family:** stress-based critical-plane criterion for fatigue limits.
- **Idea:** like McDiarmid, the critical plane is the plane of maximum shear stress amplitude
  $C_a$. Unlike McDiarmid's linear sum, the shear amplitude and the normal stress $N^*$ on that
  plane are combined **nonlinearly**. The combination is chosen so that both fully reversed
  limits are met exactly, which extends the criterion to $1.25 < f_{-1}/t_{-1} < 2$.
- **Mean stress:** the mean normal stress on the critical plane is included by a Goodman
  correction of the normal stress amplitude.
- **Proposed for:** high-cycle fatigue limits of ductile and hard metals under proportional and
  out-of-phase tension (or bending) and torsion, with mean stresses.
- **Status:** implemented for the declared repository scope. Eqs. (3), (6), (7) and (A1) were
  checked on the paper. 50 of the 66 fatigue limits of Tables 3–8 are reproduced within 0.3
  points; the others are explained below.

## Original formulation

Zhang & Yao (2011, Eqs. 3, 6, 7), on the plane $(\varphi_c, \theta_c)$ of maximum $C_a$:

```math
\sqrt{C_a^2 - C_a N^* + \frac{t_{-1}^2}{f_{-1}^2}\left(f_{-1} + 2C_a\right)N^*} \le t_{-1},
\qquad N^* = \frac{N_a}{1 - N_m/\sigma_u},
```

where:

- $C_a$ is the shear stress amplitude on the plane: the major semi-axis of the shear ellipse
  (Eq. A1, after Papadopoulos et al. 1997);
- $N_a$ and $N_m$ are the normal stress amplitude and mean on the plane;
- $\sigma_u$ is the ultimate tensile strength;
- $f_{-1}$ and $t_{-1}$ are the fully reversed bending (or tension) and torsion limits.

**Checks at the two limits.** Fully reversed torsion gives $C_a = t_{-1}$ and $N^* = 0$. Fully
reversed bending gives $C_a = N^* = f_{-1}/2$ on the 45° plane. Both make the left-hand side
equal to $t_{-1}$.

**Validity range.** The condition that a compressive mean stress be beneficial (Eqs. 8–9) leads to
$f_{-1}/t_{-1} < 2$. The lower bound 1.25 is taken from Papadopoulos.

**Error index.** $I = 100\,(\mathrm{LHS} - t_{-1})/t_{-1}$ (Eq. 10).

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\mathrm{DP} = \frac{1}{t_{-1}}\sqrt{C_a^2 - C_a N^* + \frac{t_{-1}^2}{f_{-1}^2}(f_{-1} + 2C_a)N^*},
```

which equals 1 on the fatigue limit, so $I = 100(\mathrm{DP} - 1)$.

- **Shear amplitude.** $C_a$ is half the longest chord of the shear path
  (`np_shear_rule = "longest_chord"`), which is the major semi-axis of Eq. (A1).
- **Plane selection.** `ZHANGYAO` evaluates DP on the plane of largest $C_a$. `ZHANGYAO_ext`
  takes the largest DP over all planes, a repository variant.

**Repository choices**

- *Ties of the maximum $C_a$.* They are resolved towards the larger DP. The paper gives no tie
  rule, and this criterion is sensitive to it (see [Verification](#verification-and-limitations)).
- *Ultimate strength $\sigma_u$.* `Su` is needed only when a plane carries a mean normal stress;
  without it such a case is skipped (`MissingParamsError`). Two cases are out of scope
  (`LoadingScopeError`): $N_m \ge \sigma_u$, or a negative radicand in Eq. (6).
- *Validity range.* A warning is issued outside $1.25 < f_{-1}/t_{-1} < 2$.
- *Finite life.* The index is used as the damage parameter of the shared master curve
  $\mathrm{DP} = A N_f^b$. This is an extension of the fatigue-limit criterion
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $f_{-1}$ | `Sigm1` | MPa | workbook | yes |
| $t_{-1}$ | `Taum1` | MPa | workbook | yes |
| $\sigma_u$ | `Su` | MPa | workbook | only with mean normal stresses |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** Nothing is fitted. Eq. (6) is built so that the two fully reversed
   limits $f_{-1}$ and $t_{-1}$ lie on it, and $\sigma_u$ enters only the mean-stress correction.
2. **Repository.** The same, with the workbook limits. Only the master-curve constants are fitted,
   on the uniaxial R = -1 design tests.

## Worked example

Hard steel and 42CrMo4 of the paper's Table 2.

```python
import numpy as np

from multifat.models.cp_methods.zhang_yao import zhang_yao_lhs
from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

zy = get_models(["ZHANGYAO"])[0]
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))


def index(params, sigma_a, tau_a, delta_deg=0.0, sigma_m=0.0, tau_m=0.0):
    """Index of sigma_x = sigma_m + sigma_a sin wt, tau_xy = tau_m + tau_a sin(wt - delta)."""
    d = np.radians(delta_deg)
    S_mean = np.array([[sigma_m, tau_m, 0.0], [tau_m, 0.0, 0.0], [0.0, 0.0, 0.0]])
    S_sin = np.array([[sigma_a, tau_a * np.cos(d), 0.0], [tau_a * np.cos(d), 0.0, 0.0], [0, 0, 0]])
    S_cos = np.array([[0.0, -tau_a * np.sin(d), 0.0], [-tau_a * np.sin(d), 0.0, 0.0], [0, 0, 0]])
    result = zy.evaluate_case(
        S0=S_mean + S_sin, S1=S_mean - S_sin, E0=zero, E1=zero, params=params, R_list=planes,
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos, "E_sin": zero, "E_cos": zero},
    )  # fmt: skip
    return result.values[""]


# Hard steel (paper's Table 2): f_-1 = 313.9, t_-1 = 196.2, sigma_u = 680 MPa
hard = {"Sigm1": 313.9, "Taum1": 196.2, "Su": 680.0}
# Eq. (6) puts both fully reversed limits on the threshold
assert abs(zhang_yao_lhs(196.2, 0.0, 313.9, 196.2) - 196.2) < 1e-9  # torsion
assert abs(zhang_yao_lhs(313.9 / 2, 313.9 / 2, 313.9, 196.2) - 196.2) < 1e-9  # bending, 45 deg

# 42CrMo4 (Table 2): f_-1 = 398, t_-1 = 260, sigma_u = 1025 MPa; Table 4 row 5 has tau_m = 128
cr = {"Sigm1": 398.0, "Taum1": 260.0, "Su": 1025.0}

print(f"Table 3 row 1 (in phase):     I = {100 * (index(hard, 138.1, 167.1) - 1):+.2f} % (paper  0.0)")
print(f"Table 3 row 4 (90 deg):       I = {100 * (index(hard, 150.2, 181.7, 90) - 1):+.2f} % (paper  8.5)")
print(f"Table 4 row 5 (tau_m = 128):  I = {100 * (index(cr, 266, 128, 0, 0, 128) - 1):+.2f} % (paper -9.9)")
```

Output:

```text
Table 3 row 1 (in phase):     I = -0.09 % (paper  0.0)
Table 3 row 4 (90 deg):       I = +8.68 % (paper  8.5)
Table 4 row 5 (tau_m = 128):  I = -9.85 % (paper -9.9)
```

All three rows lie within 0.2 points of the printed index.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models ZHANGYAO --all-variants --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`), a repository extension.
- **DP comparison:** supported. The index is homogeneous of degree 1 in the load only without mean
  stresses; use `Error_load` across methods.
- **Fatigue-limit datasets:** $\mathrm{DP} - 1$ is the paper's $I/100$.
- **`_ext`:** available.

## Verification and limitations

[`tests/methods/test_zhang_yao.py`](../../tests/methods/test_zhang_yao.py) checks Tables 3–8: 66
fatigue limits on hard steel, 42CrMo4, 34Cr4, 30NCD16, mild steel and LY12CZ.

- **Shear amplitude.** The longest-chord amplitude reproduces 50 of the 66 rows within 0.3
  points; `rss` reproduces only 6 within 0.06.
- **Mean-stress correction.** Eq. (7) as printed reproduces the 42CrMo4 mean-stress rows.
- **Plane rule.** The rule "maximum $C_a$" is ill-conditioned for this criterion. At the maximum,
  $C_a$ is stationary while $N^*$ changes to first order with the plane orientation.
  - A 1° plane grid already moves the in-phase indices by up to 0.3 points.
  - Over the planes within $10^{-4}$ of the maximum $C_a$, the index spans several points.
  - In 90° out-of-phase loading, where the shear amplitude can be equal on a continuum of planes,
    it spans tens of points.

  The printed values mostly sit at the top of the tied range, consistent with the
  larger-DP tie rule. This is the reason for the 0.3-point tolerance.
- **Strict expected failures:**
  - 90° rows whose maximum $C_a$ is shared by a continuum of planes (3-8, 5-3, 5-11, 6-6, 7-6,
    8-14, 8-15). The printed value lies inside the tied range.
  - Mean-stress rows of 34Cr4 and 30NCD16 (5-7, 5-8, 5-10, 6-3, 6-4, 6-7 to 6-10), not
    reproduced on any near-maximum plane with Eq. (7) and the $\sigma_u$ of Table 2. No single
    $\sigma_u$ fits them; unresolved.
- **Other checks:** Eq. (6) at both limits; the Goodman correction; the scope and missing-$\sigma_u$
  behaviour.

Limitations:

- a fatigue-limit criterion, used for finite life only through the master curve;
- results depend on the tie rule in out-of-phase loading;
- validity range $1.25 < f_{-1}/t_{-1} < 2$.

## References

- C.-C. Zhang, W.-X. Yao (2011), An improved multiaxial high-cycle fatigue criterion based on
  critical plane approach, *Fatigue & Fracture of Engineering Materials & Structures* 34(5),
  337–344. <https://doi.org/10.1111/j.1460-2695.2010.01523.x>
- I. V. Papadopoulos, P. Davoli, C. Gorla, M. Filippini, A. Bernasconi (1997), A comparative study
  of multiaxial high-cycle fatigue criteria for metals, *International Journal of Fatigue* 19(3),
  219–235. <https://doi.org/10.1016/S0142-1123(96)00064-3>
