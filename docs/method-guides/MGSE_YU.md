# Yu et al. (2017) modified generalized strain energy (`MGSE_YU`)

> **Registry name:** `MGSE_YU` (and `MGSE_YU_ext`) · **Catalogue entry:** [I_MGSE_YU](../../METHODS.md#i_mgse_yu) ·
> **Source:** Yu, Zhu, Liu & Liu (2017), <https://doi.org/10.3390/ma10080923> ·
> **Implementation:** [`cp_methods/mgse_yu.py`](../../src/multifat/models/cp_methods/mgse_yu.py)

## Overview

- **Type and family:** energy-based critical-plane criterion with no material-specific parameters.
- **Idea:** combine the product of maximum shear stress and engineering shear strain amplitude with
  the product of maximum normal stress and normal strain amplitude on the same plane.
- **Mean stress:** enters through $\tau_{\max}$ and the tensile peak normal stress $\sigma_{n,\max}$.
- **Proposed for:** multiaxial fatigue life prediction of ductile and brittle metals under
  proportional and non-proportional loading without extra material constants.
- **Status:** implemented (`fully_implemented` with `basic_numerical_tests` in the catalogue).
  With linear-elastic workbook strains, `MGSE_YU` reduces to the same scalar expression as elastic
  [GSE](GSE.md).

## Original formulation

Yu et al. (2017, Eq. 11) define the modified generalized strain energy parameter on a critical
plane close to the maximum shear strain plane:

```math
\mathrm{MGSE} = \tau_{\max}\,\frac{\Delta\gamma}{2} + \sigma_{n,\max}\,\frac{\Delta\varepsilon_n}{2},
```

where $\tau_{\max}$ is the maximum shear stress magnitude on the plane over the cycle,
$\Delta\gamma/2$ is the engineering shear strain amplitude, $\sigma_{n,\max}$ is the maximum
normal stress, and $\Delta\varepsilon_n/2$ is the normal strain amplitude.

## Repository implementation

On every candidate plane (normal = local $z$ axis, shear components $xz, yz$):

```math
\mathrm{DP} = \tau_{\max}\,\gamma_a + \sigma_{n,\max}\,\varepsilon_{n,a},
\qquad \sigma_{n,\max} = \max(\sigma_{zz,0},\,\sigma_{zz,1},\,0),
```

with $\gamma_a = \Delta\gamma/2$ and $\varepsilon_{n,a} = \Delta\varepsilon_n/2$. DP has units of
MPa ($\mathrm{MJ/m^3}$). Both standard `MGSE_YU` and `MGSE_YU_ext` maximise DP directly, so the two
variants return identical values.

**Repository choices**

- *Maximum-damage critical plane.* The paper describes choosing a plane close to maximum shear; the
  repository explicitly maximises the full MGSE value over all candidate planes.
- *Shear path quantities.* Non-proportional shear range uses the default root-sum-square (`rss`)
  convention; $\tau_{\max}$ is evaluated over 720 time samples of the harmonic cycle (`prepare_planes`
  via `NEEDS_TAU_MAX`, falling back to endpoint maximum if `evaluate_on_plane` is called directly).
- *Compression clipped.* $\sigma_{n,\max} < 0$ is clipped to $0$.
- *Elastic strains and life curve.* Workbook strains are linear elastic (matching elastic `GSE`);
  life uses the shared master curve $\mathrm{DP} = A N_f^b$
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

`MGSE_YU` requires no criterion-specific material constants beyond the elastic constants used by the
importer to generate strains.

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook (for elastic strains) | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

`MGSE_YU` has no fitted criterion constant. Only the master-curve parameters $A$, $b$ are fitted on
the uniaxial $R = -1$ design tests.

With elastic strains in fully reversed tension ($\sigma_a$), on a plane at angle $\theta$ from the
load axis ($c = \cos^2\theta \ge \nu/(1+\nu)$):

```math
\mathrm{DP}(\theta) = \frac{\sigma_a^2}{E}\left[(2+\nu)\,c - (1+\nu)\,c^2\right],
\qquad \mathrm{DP}_{\max} = \frac{\sigma_a^2}{E}\,\frac{(2+\nu)^2}{4(1+\nu)},
```

and in fully reversed torsion $\mathrm{DP} = 2(1+\nu)\tau_a^2/E$. Equal DP at the two limits
implies $t_{-1}/f_{-1} = (2+\nu)/[2\sqrt{2}(1+\nu)] \approx 0.625$ for $\nu = 0.3$.

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$.

```python
import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

E, nu = 200000.0, 0.3  # MPa, - (synthetic)
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))
X = np.diag([1.0, 0.0, 0.0])  # unit uniaxial stress along x
T = np.zeros((3, 3))
T[0, 1] = T[1, 0] = 1.0  # unit shear stress tau_xy


def hooke(S):
    """Isotropic linear-elastic strain tensor (tensorial shear components)."""
    return ((1 + nu) * S - nu * np.trace(S) * np.eye(3)) / E


def evaluate(model, params, S_mean, S_sin, S_cos=zero):
    """Standard and _ext DP of S(t) = S_mean + S_sin sin(wt) + S_cos cos(wt), elastic strains."""
    result = model.evaluate_case(
        S0=S_mean + S_sin,
        S1=S_mean - S_sin,
        E0=hooke(S_mean + S_sin),
        E1=hooke(S_mean - S_sin),
        params=params,
        R_list=planes,
        harmonics={
            "S_mean": S_mean,
            "S_sin": S_sin,
            "S_cos": S_cos,
            "E_sin": hooke(S_sin),
            "E_cos": hooke(S_cos),
        },
    )
    return result.values[""], result.values["_ext"]


yu = get_models(["MGSE_YU"])[0]
s_a, t_a = 300.0, 200.0
dp_tension, ext_tension = evaluate(yu, {}, zero, s_a * X)
dp_torsion, _ = evaluate(yu, {}, zero, t_a * T)
dp_oop, ext_oop = evaluate(yu, {}, zero, 200.0 * X, 150.0 * T)

tension = s_a**2 / E * (2 + nu) ** 2 / (4 * (1 + nu))
torsion = 2 * (1 + nu) * t_a**2 / E
assert abs(dp_tension - tension) < 1e-7 * tension
assert abs(dp_tension - ext_tension) < 1e-12
assert abs(dp_torsion - torsion) < 1e-8 * torsion
assert abs(dp_oop - ext_oop) < 1e-12

print(f"MGSE_YU tension 300 MPa: {dp_tension:.6f} MPa (MGSE_YU_ext {ext_tension:.6f} MPa)")
print(f"MGSE_YU torsion 200 MPa: {dp_torsion:.6f} MPa")
print(f"MGSE_YU 90 deg OOP (200/150 MPa): {dp_oop:.6f} MPa (MGSE_YU_ext {ext_oop:.6f} MPa)")
```

Output:

```text
MGSE_YU tension 300 MPa: 0.457788 MPa (MGSE_YU_ext 0.457788 MPa)
MGSE_YU torsion 200 MPa: 0.520000 MPa
MGSE_YU 90 deg OOP (200/150 MPa): 0.492500 MPa (MGSE_YU_ext 0.492500 MPa)
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MGSE_YU --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`).
- **DP comparison:** supported against `DP_fit`.
- `MGSE_YU` and `MGSE_YU_ext` return the same value. On purely elastic workbook strains,
  `MGSE_YU` gives the same predictions as [GSE](GSE.md).

## Verification and limitations

- [`tests/test_method_audit.py`](../../tests/test_method_audit.py) verifies the plane expression.
- Maximum-damage selection is an explicit repository adaptation; Yu et al.'s restricted search near
  the maximum shear plane and their experimental datasets are not reproduced in `tests/methods/`.

## References

- Z. Y. Yu, S. P. Zhu, Q. Liu, Y. Liu (2017), Multiaxial fatigue damage parameter and life
  prediction without any additional material constants, *Materials* 10(8), 923.
  <https://doi.org/10.3390/ma10080923>
