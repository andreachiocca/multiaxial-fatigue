# Albinmousa–Al Hussain polar damage sum, SWT parameter (`PDS_SWT`)

> **Registry name:** `PDS_SWT` · **Catalogue entry:** [N59](../../METHODS.md#n59) ·
> **Source:** Albinmousa & Al Hussain (2021), <https://doi.org/10.1016/j.finmec.2021.100025> ·
> **Implementation:** [`integral_methods/polar_damage_sum.py`](../../src/multifat/models/integral_methods/polar_damage_sum.py)

## Overview

- **Type and family:** integral over planes (polar damage sum) of a critical-plane damage parameter.
- **Idea:** fatigue damage is not confined to the critical plane. Every plane perpendicular to the
  free surface contributes, and the damage is the **area enclosed by the polar curve** of the
  damage parameter (\varphi)$. The critical-plane value is the peak of that curve. The sibling
  model [`PDS_FS`](PDS_FS.md) uses the other damage parameter of the paper.
- **Proposed for:** strain-controlled low-cycle fatigue of six alloys (magnesium AZ31B and ZK60,
  steels S460N, 16MnR, AL6XN, aluminium 7075-T651), 241 tests under proportional and
  non-proportional tension–torsion paths, with measured stabilised hysteresis loops.
- **Status:** implemented as an **elastic, high-cycle reading** of the concept (maintainer
  decision): the strains are the linear-elastic strains of the imported stresses. The published
  procedure needs measured loops, so no per-test value of the paper can be reproduced.

## Original formulation

The plane at angle $\varphi$ is perpendicular to the free surface; its stresses and strains follow
from the in-plane transformations (Eqs. 7a–7d):

```math
\sigma_\varphi = \tfrac{\sigma_y}{2} + \tfrac{\sigma_y}{2}\cos2\varphi - \tau_{xy}\sin2\varphi, \qquad
\varepsilon_\varphi = \tfrac{\varepsilon_x+\varepsilon_y}{2} - \tfrac{\varepsilon_x-\varepsilon_y}{2}\cos2\varphi - \tfrac{\gamma_{xy}}{2}\sin2\varphi,
\qquad \tfrac{\gamma_\varphi}{2} = \tfrac{\gamma_{xy}}{2}\cos2\varphi - \tfrac{\varepsilon_x-\varepsilon_y}{2}\sin2\varphi,
```

with $\varepsilon_x = -(\nu^{el}\varepsilon_y^{el} + \nu^{pl}\varepsilon_y^{pl})$ (Eq. 7e). The damage parameter
on each plane is $D_{SWT}(\varphi) = \sigma_{max}(\varphi)\,\Delta\varepsilon_{max}(\varphi)/2$ (Eq. 8). The polar damage sum is the enclosed area (Eqs. 10–14)

```math
D^s = \frac12\int_0^{2\pi} D(\varphi)^2\,d\varphi ,
```

and the life follows from a power law ^s = \alpha N_f^{\beta}$ (Eqs. 17–19).

## Repository implementation

For one harmonic, on the plane of unit normal $\mathbf n = (\cos\psi, \sin\psi, 0)$ and tangent
$\mathbf t$, every plane quantity is  + a\sin\omega t + b\cos\omega t$; its maximum over the cycle is
 + \sqrt{a^2+b^2}$ and its half range $\sqrt{a^2+b^2}$. $\sigma(\psi) = \mathbf n\cdot\mathbf S\mathbf n$,
$\varepsilon(\psi) = \mathbf n\cdot\mathbf E\mathbf n$, $\gamma(\psi) = 2\,\mathbf t\cdot\mathbf E\mathbf n$, from the imported
harmonic stresses and their elastic strains. The integral uses 720 equally spaced planes.

**Repository choices**

- *Elastic version.* Linear-elastic strains replace the measured loops and the plastic Poisson
  ratio. This reads the concept in the high-cycle regime of the repository; it is not the
  published LCF procedure.
- *Negative plane damage.* (\varphi)$ is taken as $\max(D(\varphi), 0)$ before squaring, so a
  compressive $\sigma_{max}$ (SWT) or  + k_\varphi\sigma_{n,max}/\sigma_{yt} < 0$ (FS) adds no damage.
- *Life law.* The shared master curve, fitted on the R = −1 uniaxial tension failures.
- *Scope.* Plane stress on the x–y surface (radial/hoop, axial, torsion); a stress with a z
  component raises `LoadingScopeError`.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, - | workbook | yes (elastic strains) |
| $\boldsymbol\sigma(t)$, $\boldsymbol\varepsilon(t)$ | `S_mean`, `S_sin`, `S_cos`, `E_sin`, `E_cos` | MPa, - | Dataset sheet | yes |

## Calibration

Nothing is fitted in the damage parameter. The life law $D^s = \alpha_n N_f^{\beta_n}$ (Eq. 18) is the shared master curve, fitted on the R = −1 uniaxial tension failures (the paper fits the axial and shear tests together).

## Worked example

Three zero-mean paths with the same von Mises amplitude (400 MPa) on a steel; both polar sums are
printed, with the peak of the SWT polar curve (the critical-plane SWT value).

```python
import numpy as np

from multifat.models.integral_methods.polar_damage_sum import polar_damage
from multifat.models.registry import get_models

swt, fs = get_models(["PDS_SWT", "PDS_FS"])
E, nu = 206000.0, 0.3  # MPa, -; steel
params = {"Sy": 727.0, "k_PDS": 1.0}  # Sy of 42CrMo4_Santus, an illustrative k_phi
zero = np.zeros((3, 3))


def harmonics(S_sin, S_cos):
    """Harmonic stresses and their linear-elastic strains (zero mean)."""
    strain = lambda S: ((1 + nu) * S - nu * np.trace(S) * np.eye(3)) / E  # noqa: E731
    return {"S_mean": zero, "S_sin": S_sin, "S_cos": S_cos,
            "E_mean": zero, "E_sin": strain(S_sin), "E_cos": strain(S_cos)}  # fmt: skip


def plane(sy, txy):
    """Axial sigma_y and torsional tau_xy on the x-y surface."""
    return np.array([[0.0, txy, 0.0], [txy, sy, 0.0], [0.0, 0.0, 0.0]])


# three paths with the same von Mises amplitude, 400 MPa
for label, S_sin, S_cos in (
    ("tension 400", plane(400.0, 0.0), zero),
    ("torsion 231", plane(0.0, 400.0 / np.sqrt(3)), zero),
    ("tension 283 + torsion 163 at 90 deg", plane(283.0, 0.0), plane(0.0, 163.0)),
):
    h = harmonics(S_sin, S_cos)
    ds_swt = swt.evaluate_case(S0=None, S1=None, E0=None, E1=None, params={}, harmonics=h)
    ds_fs = fs.evaluate_case(S0=None, S1=None, E0=None, E1=None, params=params, harmonics=h)
    peak = polar_damage("SWT", zero, S_sin, S_cos, h["E_sin"], h["E_cos"]).max()
    print(f"{label:37s} D^s SWT = {ds_swt.values['']:.4f} MPa^2, peak D_SWT = {peak:.3f} MPa, "
          f"D^s FS = {ds_fs.values['']:.3e}")  # fmt: skip
```

Output:

```text
tension 400                           D^s SWT = 0.4778 MPa^2, peak D_SWT = 0.777 MPa, D^s FS = 1.646e-05
torsion 231                           D^s SWT = 0.1335 MPa^2, peak D_SWT = 0.337 MPa, D^s FS = 1.728e-05
tension 283 + torsion 163 at 90 deg   D^s SWT = 0.2124 MPa^2, peak D_SWT = 0.389 MPa, D^s FS = 1.862e-05
```

For tension the peak equals $\sigma_a^2/E = 0.777$ MPa, the critical-plane SWT value. At equal von
Mises amplitude the SWT polar sum is largest in tension and smallest in torsion, while the FS polar
sum, shear-driven, is larger in torsion than in tension and largest for the out-of-phase path.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_FABEST --models PDS_SWT --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on ^s$ (`Nf_expected`).
- **DP comparison:** supported. ^s$ is not homogeneous of degree one in the load (fourth power
  for SWT, quadratic times a stress term for FS), so use `Error_load`.
- **Fatigue-limit datasets:** evaluated like the other strain-based criteria (needs E and nu).
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_polar_damage_sum.py`](../../tests/methods/test_polar_damage_sum.py) checks:

- the plane transformations of Eqs. (7a), (7c), (7d);
- the closed form of a uniaxial elastic state, whose peak is the critical-plane SWT value (Fig. 2b);
- the convergence of the Riemann sum (Eq. 11) to the integral (Eq. 12) for an out-of-phase path;
- the effect of \varphi$, the models' interface and the plane-stress scope.

Not reproduced: Figs. 4 and 6–8 and Tables 3 and 5 rest on measured stabilised loops of LCF tests,
which the elastic version cannot reproduce; no per-test value is printed.

Limitations:

- elastic strains only (no plastic strain, no measured loops);
- plane stress on the surface only;
- the life law is fitted on the uniaxial tests, not on axial and shear together.

## References

- J. Albinmousa, M. Al Hussain (2021), Polar damage sum concept for constant amplitude
  proportional and nonproportional multiaxial fatigue analysis, *Forces in Mechanics* 4, 100025.
  <https://doi.org/10.1016/j.finmec.2021.100025>
- A. Fatemi, D. F. Socie (1988), A critical plane approach to multiaxial fatigue damage including
  out-of-phase loading, *Fatigue & Fracture of Engineering Materials & Structures* 11, 149–165.
- K. N. Smith, P. Watson, T. H. Topper (1970), A stress-strain function for the fatigue of metals,
  *Journal of Materials* 5, 767–778.
