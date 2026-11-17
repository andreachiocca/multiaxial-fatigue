# Matake (`MATAKE`)

> **Registry name:** `MATAKE` (and `MATAKE_ext`) · **Catalogue entry:**
> [I_MATAKE](../../METHODS.md#i_matake) · **Source:** Matake (1977),
> <https://doi.org/10.1299/jsme1958.20.257> · **Implementation:**
> [`cp_methods/matake.py`](../../src/multifat/models/cp_methods/matake.py)

## Overview

- **Type and family:** stress-based critical-plane fatigue-limit criterion.
- **Idea:** like Findley, shear amplitude plus normal stress on a plane, but the critical plane is
  the plane of **maximum shear stress amplitude**, and the normal-stress coefficient follows
  directly from the two fatigue limits.
- **Proposed for:** the fatigue limit of metals under combined stresses.
- **Status:** fully implemented (tension-only repository adaptation). With the circumscribed-circle
  amplitude, the Matake column of Papadopoulos et al. (1997) Table 1 is reproduced on 39 of 42 rows.

## Original formulation

Matake (1977, Eq. 1, p. 5) writes, on the plane of maximum shear stress amplitude,

```math
\tau_a + \mu\,\sigma_{n,\max} = t_{-1}, \qquad \mu = \frac{2\,t_{-1}}{f_{-1}} - 1,
```

where $\tau_a$ is the shear stress amplitude and $\sigma_{n,\max}$ the maximum normal stress on that
plane. The coefficient (Appendix, p. 7) is the slope of the line through $(0, t_{-1})$ (torsion)
and $(f_{-1}/2, f_{-1}/2)$ (tension, 45° plane). The Appendix prints two different, inconsistent
expressions for it; both are misprints of $2t_{-1}/f_{-1} - 1$.

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\mathrm{DP} = \tau_a + \mu\,\sigma_{n,\max}, \qquad \mu = \frac{2\,t_{-1}}{f_{-1}} - 1,
\qquad \sigma_{n,\max} = \max(\sigma_{n,m} + \sigma_{n,a},\,0),
```

in MPa; $\mathrm{DP} = t_{-1}$ at both fatigue limits.

- `MATAKE`: DP on the plane of largest $\tau_a$ (Matake's rule); ties (e.g. the cone of 45° planes
  in tension) resolved towards the larger DP.
- `MATAKE_ext`: largest DP over all planes (repository variant).
- Out-of-phase shear amplitude: half the longest chord of the shear path, the circumscribed-circle
  radius for one harmonic (`longest_chord`) (maintainer decision, 17 November 2026; RSS before).

**Repository choices**

- *Compression clipped.* $\sigma_{n,\max} < 0$ is set to 0; Matake's linear relation is not claimed
  for compressive normal stresses.
- *Out-of-phase amplitude.* The circumscribed-circle amplitude (`longest_chord`) with which
  Papadopoulos et al. (1997) restate Matake and their Table 1 column is reproduced; Matake's paper
  does not define one for rotating shear. The RSS rule, used until 17 November 2026, misses the
  out-of-phase rows. Selectable with `--shear-amplitude-algorithm`.
- *Life.* Shared master curve on the DP
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $f_{-1}$ | `Sigm1` | MPa (amplitude) | workbook Summary, R = -1 | yes, $> 0$ |
| $t_{-1}$ | `Taum1` | MPa (amplitude) | workbook Summary, R = -1 | yes, $> 0$ |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $\mu$ from the two fatigue limits (closed form above).
2. **Repository.** The same closed form from the workbook limits. Nothing is fitted except the
   master-curve constants on the uniaxial R = -1 design tests. $\mu < 0$ when
   $t_{-1}/f_{-1} < 1/2$.

## Worked example

Synthetic demonstration material: $f_{-1} = 300$, $t_{-1} = 200$ MPa, so $\mu = 1/3$.

```python
import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

f_m1, t_m1 = 300.0, 200.0  # MPa, synthetic
mu = 2 * t_m1 / f_m1 - 1
matake = get_models(["MATAKE"])[0]
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))


def matake_dp(S_mean, S_sin):
    """MATAKE of S(t) = S_mean + S_sin sin(wt)."""
    result = matake.evaluate_case(
        S0=S_mean + S_sin,
        S1=S_mean - S_sin,
        E0=zero,
        E1=zero,
        params={"Sigm1": f_m1, "Taum1": t_m1},
        R_list=planes,
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": zero},
    )
    return result.values[""]


torsion = np.zeros((3, 3))
torsion[0, 1] = torsion[1, 0] = t_m1
assert abs(matake_dp(zero, torsion) - t_m1) < 1e-8 * t_m1
assert abs(matake_dp(zero, np.diag([f_m1, 0, 0])) - t_m1) < 1e-8 * t_m1

# Tension with a mean: on the 45 deg planes tau_a = s_a / 2, sigma_n,max = (s_m + s_a) / 2
s_a, s_m = 200.0, 100.0
with_mean = matake_dp(np.diag([s_m, 0, 0]), np.diag([s_a, 0, 0]))
assert abs(with_mean - (s_a / 2 + mu * (s_m + s_a) / 2)) < 1e-8 * with_mean

print(f"mu = {mu:.4f}; DP at both limits = t_-1 = {t_m1:.1f} MPa")
print(f"tension with mean: DP = {with_mean:.2f} MPa, fatigue index {with_mean / t_m1:.3f}")
```

Output:

```text
mu = 0.3333; DP at both limits = t_-1 = 200.0 MPa
tension with mean: DP = 150.00 MPa, fatigue index 0.750
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MATAKE --all-variants --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`).
- **DP comparison:** supported, against the master-curve `DP_fit`.
- **`_ext`:** available; maximum DP over planes, a repository variant.
- **Fatigue-limit datasets:** $\mathrm{DP}/t_{-1} - 1$ is Matake's error index.

## Verification and limitations

- [`tests/methods/test_matake.py`](../../tests/methods/test_matake.py): Matake (1977) Table 3 and
  the coefficient derivation (with the Appendix misprints documented).
- [`tests/methods/test_papadopoulos1997_matake_mcd.py`](../../tests/methods/test_papadopoulos1997_matake_mcd.py):
  Papadopoulos et al. (1997) Table 1 with the circumscribed-circle amplitude (39/42 rows).
- Limitations: tension-only normal stress; fatigue-limit criterion extended to finite life by the master curve.

## References

- T. Matake (1977), An explanation on fatigue limit under combined stress, *Bulletin of the JSME*. <https://doi.org/10.1299/jsme1958.20.257>
- I. V. Papadopoulos, P. Davoli, C. Gorla, M. Filippini, A. Bernasconi (1997), A comparative study
  of multiaxial high-cycle fatigue criteria for metals, *International Journal of Fatigue* 19(3),
  219–235. <https://doi.org/10.1016/S0142-1123(96)00064-3>
