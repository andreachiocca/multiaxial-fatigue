# McDiarmid (`MCD`)

> **Registry name:** `MCD` (and `MCD_ext`) · **Catalogue entry:** [I_MCD](../../METHODS.md#i_mcd) ·
> **Source:** McDiarmid (1991), <https://doi.org/10.1111/j.1460-2695.1991.tb00673.x> ·
> **Implementation:** [`cp_methods/mcdiarmid.py`](../../src/multifat/models/cp_methods/mcdiarmid.py)

## Overview

- **Type and family:** stress-based critical-plane criterion.
- **Idea:** shear amplitude on the plane of maximum shear amplitude, with a normal-stress term
  normalised by the ultimate tensile strength; the shear strength depends on whether cracks grow
  along the surface (case A) or into the depth (case B).
- **Proposed for:** high-cycle multiaxial fatigue of metals.
- **Status:** fully implemented for case A (surface-crack repository scope). With the
  circumscribed-circle amplitude, the McDiarmid column of Papadopoulos et al. (1997) Table 1 is
  reproduced on 41 of 43 rows.

## Original formulation

McDiarmid (1991, Eqs. 7–8, p. 437), on the plane of maximum shear stress amplitude:

```math
\frac{\tau_a}{t_{A,B}} + \frac{\sigma_{n,\max}}{2\,\sigma_u} = 1,
```

with $\sigma_{n,\max} = \sigma_{n,a} + \sigma_{n,m}$ the maximum normal stress on that plane,
$\sigma_u$ the ultimate tensile strength and $t_{A,B}$ the shear fatigue strength for case-A
(combined bending and twisting, cracks along the surface) or case-B cracks. The relation is stated
for $0 \le \sigma_{n,\max} \le \sigma_u$.

## Repository implementation

Multiplying by $t_A$ gives a stress-like DP (MPa) on every candidate plane:

```math
\mathrm{DP} = \tau_a + \frac{t_{-1}}{2\,S_u}\,\sigma_{n,\max}, \qquad
\sigma_{n,\max} = \max(\sigma_{n,m} + \sigma_{n,a},\,0),
```

which equals $t_{-1}$ on McDiarmid's limit.

- `MCD`: DP on the plane of largest $\tau_a$; ties resolved towards the larger DP.
- `MCD_ext`: largest DP over all planes (repository variant).
- Out-of-phase shear amplitude: half the longest chord of the shear path, the circumscribed-circle
  radius for one harmonic (`longest_chord`) (maintainer decision, 17 November 2026; RSS before).

**Repository choices**

- *Case A only:* $t_A = t_{-1}$ (fully reversed torsion limit); case-B strength is not supplied by
  the workbooks.
- *Validity range not enforced:* $\sigma_{n,\max} > \sigma_u$ is evaluated; $\sigma_{n,\max} < 0$ is
  clipped to 0.
- *Out-of-phase plane association:* McDiarmid's rule (plane of maximum principal-stress range) is
  not implemented; the shear amplitude is the circumscribed-circle amplitude with which Papadopoulos
  et al. (1997) restate McDiarmid and their Table 1 column is reproduced (RSS until 17 November
  2026).
- *Life:* shared master curve on the DP
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $t_A$ | `Taum1` | MPa (amplitude) | workbook Summary, R = -1 | yes, $> 0$ |
| $\sigma_u$ | `Su` | MPa | workbook Summary | yes, $> 0$ |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $t_A$ (or $t_B$) from torsion tests; $\sigma_u$ from a tensile test.
2. **Repository.** Workbook `Taum1` and `Su`; nothing else is fitted except the master-curve
   constants on the uniaxial R = -1 design tests. The tension fatigue limit is not an input: in
   reversed tension McDiarmid predicts $f_{-1} = t_A / (1/2 + t_A/(4\sigma_u))$, which a comparison
   with the measured `Sigm1` can test.

## Worked example

Synthetic demonstration material: $t_{-1} = 200$, $S_u = 800$ MPa.

```python
import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

t_m1, Su = 200.0, 800.0  # MPa, synthetic
mcd = get_models(["MCD"])[0]
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))


def mcd_dp(S_mean, S_sin):
    """MCD of S(t) = S_mean + S_sin sin(wt)."""
    result = mcd.evaluate_case(
        S0=S_mean + S_sin,
        S1=S_mean - S_sin,
        E0=zero,
        E1=zero,
        params={"Taum1": t_m1, "Su": Su},
        R_list=planes,
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": zero},
    )
    return result.values[""]


torsion = np.zeros((3, 3))
torsion[0, 1] = torsion[1, 0] = t_m1
assert abs(mcd_dp(zero, torsion) - t_m1) < 1e-8 * t_m1

# Reversed tension at the tension limit McDiarmid implies: DP must equal t_-1
f_implied = t_m1 / (0.5 + t_m1 / (4 * Su))
assert abs(mcd_dp(zero, np.diag([f_implied, 0, 0])) - t_m1) < 1e-8 * t_m1

# Tension with a mean: 45 deg planes, tau_a = s_a / 2, sigma_n,max = (s_m + s_a) / 2
s_a, s_m = 200.0, 100.0
with_mean = mcd_dp(np.diag([s_m, 0, 0]), np.diag([s_a, 0, 0]))
assert abs(with_mean - (s_a / 2 + t_m1 / (2 * Su) * (s_m + s_a) / 2)) < 1e-8 * with_mean

print(f"tension limit implied by McDiarmid: {f_implied:.1f} MPa")
print(f"tension with mean: DP = {with_mean:.2f} MPa, fatigue index {with_mean / t_m1:.3f}")
```

Output:

```text
tension limit implied by McDiarmid: 355.6 MPa
tension with mean: DP = 118.75 MPa, fatigue index 0.594
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MCD --all-variants --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`).
- **DP comparison:** supported, against the master-curve `DP_fit`.
- **`_ext`:** available; maximum DP over planes, a repository variant.

## Verification and limitations

- [`tests/methods/test_mcdiarmid.py`](../../tests/methods/test_mcdiarmid.py): Eq. (7)/(8) on bending,
  torsion and mean-stress cases.
- [`tests/methods/test_papadopoulos1997_matake_mcd.py`](../../tests/methods/test_papadopoulos1997_matake_mcd.py):
  Papadopoulos et al. (1997) Table 1 with the circumscribed-circle amplitude (41/43 rows).
- Limitations: case A only; validity range not enforced; the out-of-phase plane rule of the source
  is not implemented.

## References

- D. L. McDiarmid (1991), A general criterion for high cycle multiaxial fatigue failure, *Fatigue
  & Fracture of Engineering Materials & Structures* 14.
  <https://doi.org/10.1111/j.1460-2695.1991.tb00673.x>
- I. V. Papadopoulos, P. Davoli, C. Gorla, M. Filippini, A. Bernasconi (1997), A comparative study
  of multiaxial high-cycle fatigue criteria for metals, *International Journal of Fatigue* 19(3),
  219–235. <https://doi.org/10.1016/S0142-1123(96)00064-3>
