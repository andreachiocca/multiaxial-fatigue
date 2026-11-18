# Ince–Glinka generalized strain energy (`GSE`)

> **Registry name:** `GSE` · **Catalogue entry:** [I_GSE](../../METHODS.md#i_gse) ·
> **Source:** Ince & Glinka (2014), <https://doi.org/10.1016/j.ijfatigue.2013.10.007> ·
> **Implementation:** [`cp_methods/gse.py`](../../src/multifat/models/cp_methods/gse.py)

## Overview

- **Type and family:** energy-based critical-plane criterion (shear and normal energy terms).
- **Idea:** on each material plane, the damage is an energy density with four terms: elastic and
  plastic, shear and normal. The elastic terms are weighted by the **maximum** stresses, so that a
  mean stress increases damage. The plastic terms are weighted by the stress **amplitudes**. The
  critical plane is the plane on which this energy is largest.
- **Proposed for:** low- and high-cycle multiaxial fatigue under proportional and
  non-proportional loading, with mean stresses.
- **Status:** fully implemented for the declared repository scope. Eq. (1) is checked against the
  paper. No experimental dataset is reproduced. The strain-amplitude form of the same paper is the
  separate model [GSA](GSA.md).

## Original formulation

Ince & Glinka (2014, Eq. 1):

```math
W^*_{gen} = \left[\tau_{\max}\frac{\Delta\gamma_e}{2} + \frac{\Delta\tau}{2}\frac{\Delta\gamma_p}{2}
+ \sigma_{n,\max}\frac{\Delta\varepsilon_e}{2} + \frac{\Delta\sigma_n}{2}\frac{\Delta\varepsilon_p}{2}\right]_{\max},
```

where:

- $\Delta\gamma_e$ and $\Delta\gamma_p$ are the elastic and plastic **engineering** shear strain
  ranges;
- $\Delta\varepsilon_e$ and $\Delta\varepsilon_p$ are the elastic and plastic normal strain ranges;
- $\tau_{\max}$ and $\sigma_{n,\max}$ are the maximum shear and normal stresses on the plane;
- $\Delta\tau$ and $\Delta\sigma_n$ are the shear and normal stress ranges on the plane.

The bracket is maximised over the planes ("maximum damage plane"). The parameter is related to life
through an energy–life curve derived from the strain–life constants of the material.

## Repository implementation

On every candidate plane (normal = local $z$ axis), Eq. (1) is evaluated as written:

- the plastic strain is the total strain minus the isotropic elastic strain of the stresses,
  $\varepsilon_p = \varepsilon - \varepsilon_e(\boldsymbol\sigma; E, \nu)$;
- the shear ranges are taken on the shear path of the plane, by the default rule
  `np_shear_rule = "rss"`;
- $\tau_{\max}$ is the largest magnitude of the shear vector over the cycle, mean included,
  obtained from 720 time samples;
- $\sigma_{n,\max} = \max(\sigma_{n,m} + \sigma_{n,a},\,0)$.

DP has units of energy density, MPa (= MJ/m³). The plane rule and the DP coincide, so `GSE` is
already the maximum over planes. The `_ext` value is identical and is not documented separately.

**Repository choices**

- *Compression clipped.* $\sigma_{n,\max} < 0$ is set to 0.
- *Elastic strains on imported workbooks.* The workbook strains are elastic, so the two plastic
  terms vanish and
  $\mathrm{DP} = \tau_{\max}\,\gamma_{a} + \sigma_{n,\max}\,\varepsilon_{n,a}$.
  Plastic terms appear only with total-strain histories that carry a plastic part.
- *Harmonic decomposition.* The elastic/plastic split is applied to the harmonic tensors before the
  ranges are taken, so that synthetic extrema cannot create fictitious plastic strain.
- *Life.* The shared master curve replaces the energy–life curve
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook | yes (elastic/plastic split) |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

GSE has no criterion-specific constant. The master-curve constants $A$, $b$ are fitted on the
uniaxial $R = -1$ design tests.

With elastic strains the parameter fixes the ratio of the torsion and tension fatigue limits. In
fully reversed tension, on a plane whose normal makes an angle $\theta$ with the load
($c = \cos^2\theta$),

```math
\mathrm{DP}(\theta) = \frac{\sigma_a^2}{E}\left[(2+\nu)\,c - (1+\nu)\,c^2\right],
\qquad c^* = \frac{2+\nu}{2(1+\nu)},
\qquad \mathrm{DP}_{\max} = \frac{\sigma_a^2}{E}\,\frac{(2+\nu)^2}{4(1+\nu)} .
```

The critical plane is therefore inclined, at about 20° from the load axis for $\nu = 0.3$, not
normal to it. In fully reversed torsion $\mathrm{DP} = 2(1+\nu)\tau_a^2/E$. Equal damage at the two
limits gives

```math
\frac{t_{-1}}{f_{-1}} = \frac{2+\nu}{2\sqrt2\,(1+\nu)} \approx 0.63 \quad (\nu = 0.3).
```

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


gse = get_models(["GSE"])[0]
params = {"E": E, "nu": nu}

s_a, t_a = 300.0, 200.0
w_tension, _ = evaluate(gse, params, zero, s_a * X)
w_torsion, _ = evaluate(gse, params, zero, t_a * T)

# Tension, plane at theta from the axis (c = cos^2 theta), elastic strains:
# W = s^2/E [(2 + nu) c - (1 + nu) c^2], maximal at c* = (2 + nu) / (2 (1 + nu))
c_star = (2 + nu) / (2 * (1 + nu))
tension = s_a**2 / E * (2 + nu) ** 2 / (4 * (1 + nu))
theta_star = np.degrees(np.arccos(np.sqrt(c_star)))
assert abs(w_tension - tension) < 1e-6 * tension
# Torsion, plane of tau_xy: W = tau gamma_a = 2 (1 + nu) tau^2 / E
assert abs(w_torsion - 2 * (1 + nu) * t_a**2 / E) < 1e-9 * w_torsion

# Tensile mean stress: tau_max and sigma_n,max include the mean, the amplitudes do not
w_mean, _ = evaluate(gse, params, 100.0 * X, 250.0 * X)

print(f"GSE tension 300 MPa: {w_tension:.4f} MPa (plane at {theta_star:.1f} deg from the axis)")
print(f"GSE torsion 200 MPa: {w_torsion:.4f} MPa")
print(f"t/f at equal DP: {(2 + nu) / (2 * np.sqrt(2) * (1 + nu)):.3f}")
print(f"GSE tension 250 MPa with a 100 MPa mean: {w_mean:.4f} MPa")
```

Output:

```text
GSE tension 300 MPa: 0.4578 MPa (plane at 19.9 deg from the axis)
GSE torsion 200 MPa: 0.5200 MPa
t/f at equal DP: 0.626
GSE tension 250 MPa with a 100 MPa mean: 0.4451 MPa
```

The last line shows the mean-stress effect. A fully reversed amplitude of 250 MPa alone gives
$0.4578\,(250/300)^2 = 0.318$ MPa. The 100 MPa mean raises $\tau_{\max}$ and $\sigma_{n,\max}$, and
with them the DP, to 0.445 MPa.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models GSE --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`), not the energy–life curve.
- **DP comparison:** supported, against the master-curve `DP_fit`. DP is quadratic in the load;
  `Error_load` removes the factor 2.
- **`_ext`:** identical to `GSE`, because the plane rule is the DP itself.

## Verification and limitations

[`tests/methods/test_gse_gsa.py`](../../tests/methods/test_gse_gsa.py) checks the four terms of
Eq. (1) on a plane state with explicit elastic and plastic strain amplitudes.

Not verified: the experimental lives of Ince & Glinka (2014).

Limitations:

- with the workbooks' elastic strains, only the two elastic terms are active;
- the plastic terms need total-strain histories, which the workbook importer does not produce.

## References

- A. Ince, G. Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life
  prediction under proportional and non-proportional loadings, *International Journal of
  Fatigue* 62, 34–41. <https://doi.org/10.1016/j.ijfatigue.2013.10.007>
