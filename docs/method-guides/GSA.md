# Ince–Glinka generalized strain amplitude (`GSA`)

> **Registry name:** `GSA` · **Catalogue entry:** [I_GSA](../../METHODS.md#i_gsa) ·
> **Source:** Ince & Glinka (2014), <https://doi.org/10.1016/j.ijfatigue.2013.10.007> ·
> **Implementation:** [`cp_methods/gsa.py`](../../src/multifat/models/cp_methods/gsa.py)

## Overview

- **Type and family:** strain-based critical-plane criterion; the strain-amplitude form of the
  generalized strain energy [GSE](GSE.md).
- **Idea:** the GSE energy is divided by the fatigue strength coefficients. Each elastic strain
  amplitude is weighted by the ratio of the maximum stress to its fatigue strength coefficient,
  while the plastic strain amplitudes enter unweighted. The parameter is thus a strain amplitude
  that can be compared with a strain–life curve.
- **Proposed for:** the same scope as GSE: proportional and non-proportional loading, with mean
  stresses.
- **Status:** fully implemented for the declared repository scope. Eq. (5) is checked against the
  paper. No experimental dataset is reproduced.

## Original formulation

Ince & Glinka (2014, Eq. 5):

```math
\frac{\Delta\varepsilon^*}{2} = \left[\frac{\tau_{\max}}{\tau_f'}\frac{\Delta\gamma_e}{2}
+ \frac{\Delta\gamma_p}{2} + \frac{\sigma_{n,\max}}{\sigma_f'}\frac{\Delta\varepsilon_e}{2}
+ \frac{\Delta\varepsilon_p}{2}\right]_{\max},
```

with the symbols of [GSE](GSE.md#original-formulation). $\sigma_f'$ and $\tau_f'$ are the axial
and shear fatigue strength coefficients. The bracket is maximised over the planes, and the result
is related to life through a strain–life curve.

## Repository implementation

On every candidate plane, Eq. (5) is evaluated with the same strain decomposition, shear rule,
$\tau_{\max}$ and clipping of $\sigma_{n,\max}$ as [GSE](GSE.md#repository-implementation). DP is
dimensionless. The plane rule is the DP itself, so the `_ext` value is identical.

**Repository choices**

- *Fatigue strength coefficients.* If `sigma_f` is missing from the workbook, it is estimated with
  the uniform material law $\sigma_f' = 1.5\,S_u$ (Bäumel & Seeger, as summarised by Zhao et al.,
  2020). This estimate was developed for steels, and a warning is printed. If `tau_f` is missing,
  $\tau_f' = \sigma_f'/\sqrt3$ (von Mises); the paper gives no value.
- *Elastic strains on imported workbooks.* The plastic terms vanish, so

  ```math
  \mathrm{DP} = \frac{1}{\sigma_f'}\left(\sqrt3\,\tau_{\max}\gamma_a + \sigma_{n,\max}\varepsilon_{n,a}\right)
  ```

  when $\tau_f'$ takes its default. $\sigma_f'$ then only scales DP. The fitted master curve absorbs
  that scale, so the predicted lives do not depend on $\sigma_f'$. A separately supplied $\tau_f'$
  does change them.
- *Life.* The shared master curve replaces the strain–life curve
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_f'$ | `sigma_f` / Summary | MPa | workbook, else $1.5\,S_u$ | no (needs `Su` if missing) |
| $\tau_f'$ | `tau_f` / Summary | MPa | workbook, else $\sigma_f'/\sqrt3$ | no |
| $S_u$ | `Su` / Summary | MPa | workbook | only to estimate $\sigma_f'$ |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

GSA has no fitted criterion constant. $\sigma_f'$ and $\tau_f'$ are read or estimated, never fitted
by the pipeline. The master-curve constants $A$, $b$ are fitted on the uniaxial $R = -1$ design
tests.

With elastic strains and the default $\tau_f'$, let $a = 2\sqrt3(1+\nu) - \nu$ and
$c = \cos^2\theta$. Fully reversed tension gives

```math
\mathrm{DP}(\theta) = \frac{\sigma_a^2}{E\sigma_f'}\left[a\,c - (a-1)\,c^2\right],
\qquad \mathrm{DP}_{\max} = \frac{\sigma_a^2}{E\sigma_f'}\,\frac{a^2}{4(a-1)},
```

and fully reversed torsion gives $\mathrm{DP} = 2\sqrt3(1+\nu)\tau_a^2/(E\sigma_f')$. Equal damage at
the two limits implies $t_{-1}/f_{-1} \approx 0.55$ for $\nu = 0.3$. This is below the GSE value
(0.63), because the shear term carries the factor $\sigma_f'/\tau_f' = \sqrt3$.

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$,
$\sigma_f' = 1350$ MPa.

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


gsa = get_models(["GSA"])[0]
sigma_f = 1350.0  # MPa (synthetic); tau_f defaults to sigma_f / sqrt(3)
params = {"E": E, "nu": nu, "sigma_f": sigma_f}

s_a, t_a = 300.0, 200.0
dp_tension, _ = evaluate(gsa, params, zero, s_a * X)
dp_torsion, _ = evaluate(gsa, params, zero, t_a * T)

# Elastic strains, plane at theta from the axis (c = cos^2 theta):
# GSA = s^2/(E sigma_f) [a c - (a - 1) c^2], a = 2 sqrt(3) (1 + nu) - nu
a = 2 * np.sqrt(3) * (1 + nu) - nu
tension = s_a**2 / (E * sigma_f) * a**2 / (4 * (a - 1))
torsion = np.sqrt(3) * 2 * (1 + nu) * t_a**2 / (E * sigma_f)
assert abs(dp_tension - tension) < 1e-6 * tension
assert abs(dp_torsion - torsion) < 1e-9 * torsion

# With elastic strains and the default tau_f, sigma_f only scales the DP
doubled, _ = evaluate(gsa, {**params, "sigma_f": 2 * sigma_f}, zero, s_a * X)
assert abs(2 * doubled - dp_tension) < 1e-12

print(f"GSA tension 300 MPa: {dp_tension:.4e}")
print(f"GSA torsion 200 MPa: {dp_torsion:.4e}")
print(f"t/f at equal DP: {np.sqrt(a**2 / (4 * (a - 1)) / (2 * np.sqrt(3) * (1 + nu))):.3f}")
```

Output:

```text
GSA tension 300 MPa: 4.5963e-04
GSA torsion 200 MPa: 6.6716e-04
t/f at equal DP: 0.553
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models GSA --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`), not the strain–life curve.
- **DP comparison:** supported, against `DP_fit`. DP is quadratic in the load; use `Error_load` to
  compare across methods.
- **`_ext`:** identical to `GSA`.

## Verification and limitations

[`tests/methods/test_gse_gsa.py`](../../tests/methods/test_gse_gsa.py) checks:

- the four terms of Eq. (5) on a plane state with elastic and plastic strains;
- the default $\tau_f' = \sigma_f'/\sqrt3$.

Not verified: the experimental lives of Ince & Glinka (2014).

Limitations:

- elastic strains only on imported workbooks;
- the steel-based estimate of $\sigma_f'$ matters only through a supplied $\tau_f'$, or when plastic
  terms are present.

## References

- A. Ince, G. Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life
  prediction under proportional and non-proportional loadings, *International Journal of
  Fatigue* 62, 34–41. <https://doi.org/10.1016/j.ijfatigue.2013.10.007>
- Y. Zhao et al. (2020), Fatigue properties estimation and life prediction for steels under axial,
  torsional, and in-phase loading (source of the uniform-material-law estimate used for missing
  $\sigma_f'$).
