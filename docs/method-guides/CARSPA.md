# Carpinteri–Spagnoli critical-plane criterion for hard metals (`CARSPA`)

> **Registry name:** `CARSPA` · **Catalogue entry:** [I_CARSPA](../../METHODS.md#i_carspa) ·
> **Source:** Carpinteri & Spagnoli (2001), <https://doi.org/10.1016/S0142-1123(00)00075-X> ·
> **Implementation:** [`cp_methods/car_spa.py`](../../src/multifat/models/cp_methods/car_spa.py)

## Overview

- **Type and family:** critical plane, stress based. The plane is *constructed*, not searched:
  its normal follows from the weighted mean principal stress directions.
- **Idea:** the fatigue fracture plane follows the averaged direction of the maximum principal
  stress; the critical plane, on which the assessment is made, is rotated from it by an angle
  $\delta$ that grows from 0 (very brittle metals) to 45° (the brittle–ductile threshold
  $\tau_{af}/\sigma_{af} = 1/\sqrt3$). On that plane the maximum normal stress and the shear
  amplitude are combined quadratically.
- **Proposed for:** high-cycle fatigue limits of hard metals, /\sqrt3 \le \tau_{af}/\sigma_{af}
  \le 1$, under in-phase and out-of-phase sinusoidal loading, with mean normal stresses.
- **Status:** implemented. Tables 2–4 of the 2001 paper are reproduced (see Verification).

## Original formulation

The principal directions 1, 2, 3 ($\sigma_1 \ge \sigma_2 \ge \sigma_3$) at each instant $ are
described by the principal Euler angles $\phi, \theta, \psi$ (rotations about $, '$ and the
3-axis; Carpinteri et al. 1999, matrix 2). In the notation of the 2000 and 2001 papers direction 1
is

```math
\hat{\mathbf 1} = (\cos\phi\sin\theta,\ \sin\phi\sin\theta,\ \cos\theta).
```

A two-stage procedure (1999, Fig. 3 and Table 1) reduces the angles to
/bin/bash \le \phi, \theta \le 90°$ and hmBc90° < \psi \le 90°$. The mean directions
$\hat 1, \hat 2, \hat 3$ follow from the weighted averages (2001, Eqs. 22–23)

```math
\hat\phi = \frac1W\sum_{k} W(t_k)\,\phi(t_k),\quad
\hat\theta = \frac1W\sum_{k} W(t_k)\,\theta(t_k),\quad
\hat\psi = \frac1W\sum_{k} W(t_k)\,\psi(t_k),\quad W = \sum_k W(t_k),
```

```math
W(t_k) = \begin{cases} 0 & \sigma_1(t_k) < c\,\sigma_{af} \[2pt]
\left(\dfrac{\sigma_1(t_k)}{\sigma_{af}}\right)^{m_\sigma} & \sigma_1(t_k) \ge c\,\sigma_{af}
\end{cases}\qquad 0 < c \le 1,
```

with \sigma = -1/m$, $ the slope of the R = −1 tension S–N curve. The critical-plane normal
$\mathbf w$ forms with $\hat 1$ the angle (Eq. 24)

```math
\delta = 45\cdot\frac32\left[1 - \left(\frac{\tau_{af}}{\sigma_{af}}\right)^2\right]\ \text{deg},
```

and the criterion is (Eqs. 25–27)

```math
\sigma_{eq} = \sqrt{N_{max}^2 + \left(\frac{\sigma_{af}}{\tau_{af}}\right)^2 C_a^2} \le \sigma_{af},
\qquad N_{max} = N_m + N_a,
```

where $, $ are the mean and amplitude of the normal stress on the critical plane and $
is the amplitude of the shear stress vector, the major semi-axis of the ellipse that its tip
describes (Eqs. 9–13). The mean shear stress does not enter. The paper excludes compressive means
 < -N_a$. For $\tau_{af}/\sigma_{af} = 1/\sqrt3$, Eq. (26) becomes {max}^2 + 3C_a^2 \le
\sigma_{af}^2$ (Eq. 28). The error index is  = 100\,(\sigma_{eq}-\sigma_{af})/\sigma_{af}$
(Eq. 30);  > 0$ is conservative.

## Repository implementation

- **Interface:** case-level model on the single-harmonic tensors `S_mean`, `S_sin`,
  `S_cos`; `sigma_eq` is the returned value (`""` variant). No plane search and no `_ext`
  variant.
- **Euler angles and reductions:** matrix (2) of 1999 with direction 1 in its third column, 2 in
  the first and 3 in the second (2000 notation). The first reduction is done by enumerating the
  column-sign changes with $\det = +1$ and the two equivalent Euler triples; the second applies
  Table 1 ($\phi, \theta \to |\phi|, |\theta|$, $\psi$ keeps its sign when $\phi\psi \ge 0$).
  This reproduces the explicit plane-stress cases (a1)–(c2) of the 2000 paper. At $\theta = 0$
  only $\phi + \psi$ is defined and is put in $\phi$.
- **Mirror and return (maintainer decision, 23 Nov 2026):** the second reduction folds directions
  into the first octant of the fixed frame, i.e. it mirrors them. Applied to the fixed frame alone,
  a load and its mirror image ($\tau \to -\tau$) would give different planes ( = +3.1\,\%$
  against hmBc1.7\,\%$ for hard steel test 3 with $\beta = 0°$ and 80°$). The fixed axes are
  therefore mirrored so that direction 1 at the instant of largest weight lies in the first
  octant, and the critical-plane normal is mirrored back. For the published tests the mirror is the
  identity.
- **Rotation sense:** $\mathbf w = \cos\delta\,\hat 1 + \sin\delta\,\hat 3$, in the plane
  $\hat 1\hat 3$. The paper gives the angle only; this sense reproduces Tables 2–4, the opposite
  one gives errors down to −15 %.
- **Sampling:** 720 instants per period. Where two principal stresses are equal the directions of
  the right limit in time are used. When $\sigma_1 < c\,\sigma_{af}$ at every instant (loads far
  below the limit, e.g. in a load-scaling search), the instants of largest $\sigma_1$ are used:
  the limit of Eq. (23) as the window shrinks to the peak.
- **Scope:** `MaterialScopeError` outside /\sqrt3 \le \tau_{af}/\sigma_{af} \le 1$;
  `LoadingScopeError` for  < -N_a$ on the critical plane.

## Material constants and inputs

| Symbol | Key | Notes |
| --- | --- | --- |
| $\sigma_{af}$ | `Sigm1` | R = −1 axial limit, required |
| $\tau_{af}$ | `Taum1` | R = −1 torsion limit, required |
| \sigma$ | `m_sigma` | optional; default hmBc1/b$ from `b_fat` |
| $ | `c_CARSPA` | optional; default 0.5 (Tables 2–4) |

\sigma$ only matters when the principal directions of the weighted instants differ (most
out-of-phase paths). If neither `m_sigma` nor `b_fat` is given, such a case raises
`MissingParamsError` and the model is skipped on that material; on proportional paths it runs
without them (maintainer decision).

## Calibration

None. The criterion uses the two fatigue limits; \sigma$ is a property of the tension S–N
curve. Finite-life comparisons use the shared master curve  = A N^b$ on $\sigma_{eq}$, a
repository extension.

## Worked example

Three Swedish hard steel limits of Table 2 (Nishihara–Kawamoto), with the constants of Table 1.

```python
import numpy as np

from multifat.models.cp_methods.car_spa import carspa_critical_plane
from multifat.models.registry import get_models

carspa = get_models(["CARSPA"])[0]
params = {"Sigm1": 313.9, "Taum1": 196.2, "m_sigma": 8.7}  # Swedish hard steel, Table 1


def bending_torsion(sigma_a, tau_a, beta):
    """Harmonics of sigma_l = sigma_a sin(wt) (axis y) and tau = tau_a sin(wt - beta)."""
    b = np.radians(beta)
    S_sin, S_cos = np.zeros((3, 3)), np.zeros((3, 3))
    S_sin[1, 1] = sigma_a
    S_sin[0, 1] = S_sin[1, 0] = tau_a * np.cos(b)
    S_cos[0, 1] = S_cos[1, 0] = -tau_a * np.sin(b)
    return {"S_mean": np.zeros((3, 3)), "S_sin": S_sin, "S_cos": S_cos}


# Table 2 rows 3, 8 and 12: (sigma_a, tau_a, beta; printed eta_cal and I)
for test, sigma_a, tau_a, beta, eta_p, I_p in (
    (3, 255.1, 127.5, 0, 23, 3),
    (8, 255.1, 127.5, 60, 16, -2),
    (12, 152.5, 184.2, 90, 38, 5),
):
    h = bending_torsion(sigma_a, tau_a, beta)
    res = carspa_critical_plane(h["S_mean"], h["S_sin"], h["S_cos"], params)
    eta = np.degrees(np.arccos(abs(res["axes"][1, 0])))
    s_eq = carspa.evaluate_case(S0=None, E0=None, S1=None, E1=None, params=params,
                                harmonics=h).values[""]  # fmt: skip
    print(f"test {test:2d}: eta_cal = {eta:4.1f} deg (Table 2: {eta_p}), N_max = "
          f"{res['N_m'] + res['N_a']:5.1f}, C_a = {res['C_a']:5.1f} MPa, "
          f"I = {100 * (s_eq - 313.9) / 313.9:+.1f} % (Table 2: {I_p:+d})")  # fmt: skip
print(f"delta = {np.degrees(res['delta']):.2f} deg")
```

Output:

```text
test  3: eta_cal = 22.5 deg (Table 2: 23), N_max = 151.8, C_a = 178.7 MPa, I = +3.1 % (Table 2: +3)
test  8: eta_cal = 17.0 deg (Table 2: 16), N_max = 184.4, C_a = 156.0 MPa, I = -1.1 % (Table 2: -2)
test 12: eta_cal = 37.5 deg (Table 2: 38), N_max = 153.7, C_a = 183.0 MPa, I = +5.3 % (Table 2: +5)
delta = 41.13 deg
```

The in-phase test 3 has a fixed principal frame, so $\hat 1$ is its direction. In tests 8 and 12
the frame rotates and the weight function decides $\hat 1$. Each error index lies within the
rounding of the printed one, given the printed $\eta_{cal}$.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_FABEST --models CARSPA --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on $\sigma_{eq}$ (`Nf_expected`), a repository extension.
- **DP comparison:** supported; $\sigma_{eq}$ is homogeneous in the load wherever the weighted
  window does not change, so `Error_load` and `Error_log10_dp` agree closely.
- **Fatigue-limit datasets:** $\sigma_{eq}/\sigma_{af} - 1$ is the paper's error index $.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_carspa.py`](../../tests/methods/test_carspa.py) checks:

- the criterion step: from the printed $\eta_{cal}$ of Tables 2–4, the 30 printed $ within
  their rounding (and that of $\eta_{cal}$);
- the averaging step: $\eta_{cal}$ within 1.5° for 27 of 30 rows, and the end-to-end $ for
  these rows. For the three $\beta = 90°$ rows with $\tau_a/\sigma_a = 0.2$ the tables print 0
  where the procedure gives 2.3–3.5°; no time step, threshold or exponent gives both these zeros
  and the printed 12° of hard steel test 11. They are strict xfails, unexplained;
- that the opposite rotation sense misses by about 15 points;
- the cases (a1)–(c2) of the 2000 paper, the reduction ranges on random frames, mirror
  invariance, the closed forms at $\tau_{af}/\sigma_{af} = 1/\sqrt3$ and 1, the scope errors and
  the `m_sigma` fallback.

Limitations:

- hard metals only; the later ductile-metal versions of the criterion are separate methods;
- single-harmonic constant-amplitude loading;
- a fatigue-limit criterion, so finite lives come only from the master curve.

## References

- A. Carpinteri, A. Spagnoli (2001), Multiaxial high-cycle fatigue criterion for hard metals,
  *International Journal of Fatigue* 23, 135–145. <https://doi.org/10.1016/S0142-1123(00)00075-X>
- A. Carpinteri, E. Macha, R. Brighenti, A. Spagnoli (1999), Expected principal stress directions
  under multiaxial random loading. Part I: theoretical aspects of the weight function method,
  *International Journal of Fatigue* 21, 83–88. <https://doi.org/10.1016/S0142-1123(98)00046-2>
- A. Carpinteri, R. Brighenti, A. Spagnoli (2000), A fracture plane approach in multiaxial
  high-cycle fatigue of metals, *Fatigue & Fracture of Engineering Materials & Structures* 23,
  355–364. <https://doi.org/10.1046/j.1460-2695.2000.00265.x>
