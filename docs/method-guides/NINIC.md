# Ninić Damage Function (`NINIC`)

> **Registry name:** `NINIC` · **Catalogue entry:** [N11](../../METHODS.md#n11) ·
> **Source:** Ninić (2006), <https://doi.org/10.1016/j.ijfatigue.2005.04.014> ·
> **Implementation:** [`cp_methods/ninic.py`](../../src/multifat/models/cp_methods/ninic.py)

## Overview

- **Type and family:** stress-based critical-plane endurance criterion with a quadratic stress
  combination.
- **Idea:** on every plane, the shear stress amplitude and the maximum normal stress are each
  normalised by the corresponding endurance limit and combined as a vector sum. The plane of
  maximum damage is the critical plane, searched in three dimensions, so planes tilted into the
  surface are included.
- **Proposed for:** high-cycle fatigue limits of metals from brittle to ductile (mild and hard
  steel, cast iron, aluminium alloy) under in-phase and out-of-phase bending–torsion.
- **Status:** fully implemented for the declared repository scope. Eqs. (1)–(4) and Tables 1–5
  were checked on the paper.

## Original formulation

Ninić (2006, Eqs. 1–3) defines the Damage Function on the plane with normal direction cosines
$(l, m, n)$ and maximises it over all planes:

```math
D(l,m,n) = \left[\left(\frac{\tau_a(l,m,n)}{T_e}\right)^2 + k\left(\frac{\sigma_{n,\max}(l,m,n)}{S_e}\right)^2\right]^{1/2},
\qquad \max_{l,m,n} D = 1 \ \text{at the limit},
```

where:

- $T_e$ and $S_e$ are the fully reversed torsion and tension (bending) endurance limits;
- $\tau_a$ is the shear stress amplitude on the plane;
- $\sigma_{n,\max}$ is the maximum normal stress on the plane over the cycle.

The normal stress sensitivity factor is obtained by maximising $D$ for uniaxial loading (Eq. 4):

```math
k = \left(\frac{S_e}{T_e}\right)^2\left[1 - \left(\frac{S_e}{2T_e}\right)^2\right].
```

The error index is $I = 100\,(D - 1)$ (Eq. 14). For finite lives the paper suggests replacing
$S_e$ and $T_e$ by the fatigue strengths $S_N$ and $T_N$ at the life of interest.

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\mathrm{DP} = \sqrt{\left(\frac{\tau_a}{t_{-1}}\right)^2 + k\left(\frac{\max(\sigma_{n,\max}, 0)}{f_{-1}}\right)^2},
```

with $\tau_a$ the shear amplitude by the default `rss` rule. DP is both the plane metric and the
damage, so the critical plane is the plane of maximum DP. `NINIC_ext` equals `NINIC`.

**Repository choices**

- *Out-of-phase shear amplitude.* The paper does not define $\tau_a$ for non-proportional paths.
  The `rss` amplitude reproduces Tables 2–4; the longest chord misses rows 2-9 and 2-10 by
  1.7–1.9 points.
- *Compression clipped.* $\sigma_{n,\max} < 0$ is set to 0, because Eq. (1) squares it and a
  compressive normal stress would otherwise add damage. The paper's critical planes all carry
  tensile normal stresses.
- *Constant $k$.* $k$ follows Eq. (4) from `Sigm1` and `Taum1`, unless the Summary supplies `k_NIN`.
- *Finite life.* The index is used as the damage parameter of the shared master curve
  $\mathrm{DP} = A N_f^b$. This is an extension of the endurance criterion, not the paper's
  $S_N$, $T_N$ substitution ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $S_e$ | `Sigm1` / Summary R = -1 tension or bending limit | MPa | workbook | yes |
| $T_e$ | `Taum1` / Summary R = -1 torsion limit | MPa | workbook | yes |
| $k$ | `k_NIN` / Summary column `k_NIN` | – | Eq. (4) unless supplied | no |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $k$ from Eq. (4), so that the uniaxial limit gives $\max D = 1$. In
   torsion, the plane of $\tau_a = T_e$ carries no normal stress, so $D = 1$ independently of $k$.
2. **Repository.** The same, with the workbook limits. Nothing is fitted except the master-curve
   constants $A$ and $b$ on the uniaxial R = -1 design tests.

**Limit of Eq. (4).** For fully reversed uniaxial stress $S$ on a plane whose normal makes
$\theta$ with the load ($c = \cos^2\theta$):

```math
D^2 = \left(\frac{S}{T_e}\right)^2 c\,(1-c) + k\left(\frac{S}{S_e}\right)^2 c^2 .
```

With Eq. (4), the maximum lies at $c^* = 2\,(T_e/S_e)^2$, where $D(S_e) = 1$. This holds only if
$c^* \le 1$, i.e. $S_e/T_e \ge \sqrt2$. For more brittle behaviour ($S_e/T_e < \sqrt2$) the maximum
moves to the principal plane, $c = 1$, and the uniaxial limit gives $D = \sqrt{k} < 1$. This is the
paper's own cast-iron row 4-5, $I = -10.4$ %. The repository keeps Eq. (4) and warns when
$S_e/T_e < \sqrt2$. In that range, $k = 1$ would put both limits on $D = 1$, and it joins Eq. (4)
continuously at $S_e/T_e = \sqrt2$. This is a suggestion, not the source's rule; it can be
supplied as `k_NIN`.

## Worked example

Mild steel of the paper's Table 1 (Nishihara & Kawamoto): $S_e = 235.4$ MPa, $T_e = 137.3$ MPa.

```python
import numpy as np

from multifat.models.cp_methods.ninic import ninic_k
from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

ninic = get_models(["NINIC"])[0]
S_e, T_e = 235.4, 137.3  # MPa, mild steel of Nishihara & Kawamoto (paper's Table 1)
params = {"Sigm1": S_e, "Taum1": T_e}
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))


def damage(sigma_a, tau_a, phi_deg):
    """D for sigma_y = sigma_a sin wt, tau_xy = tau_a sin(wt - phi); y is the bar axis."""
    p = np.radians(phi_deg)
    S_sin = np.array([[0.0, tau_a * np.cos(p), 0.0], [tau_a * np.cos(p), sigma_a, 0.0], [0, 0, 0]])
    S_cos = np.array([[0.0, -tau_a * np.sin(p), 0.0], [-tau_a * np.sin(p), 0.0, 0.0], [0, 0, 0]])
    result = ninic.evaluate_case(
        S0=S_sin,
        S1=-S_sin,
        E0=zero,
        E1=zero,
        params=params,
        R_list=planes,
        harmonics={"S_mean": zero, "S_sin": S_sin, "S_cos": S_cos, "E_sin": zero, "E_cos": zero},
    )
    return result.values[""]


k = ninic_k(S_e, T_e)  # Eq. (4)
# S_e/T_e = 1.71 >= sqrt(2): both reference limits lie on D = 1
assert abs(damage(S_e, 0.0, 0) - 1.0) < 1e-6
assert abs(damage(0.0, T_e, 0) - 1.0) < 1e-9
# Uniaxial critical plane: cos^2(theta) = 2 (T_e/S_e)^2 from the bar axis
theta = np.degrees(np.arccos(np.sqrt(2 * (T_e / S_e) ** 2)))

print(f"k = {k:.4f}; uniaxial critical plane at {theta:.1f} deg from the axis")
print(f"Table 2 row 8 (phi = 90 deg): I = {100 * (damage(108.9, 131.5, 90) - 1):+.2f} %")
print(f"Table 2 row 9 (phi = 90 deg): I = {100 * (damage(201.1, 100.6, 90) - 1):+.2f} %")
```

Output:

```text
k = 0.7793; uniaxial critical plane at 34.4 deg from the axis
Table 2 row 8 (phi = 90 deg): I = +4.12 %
Table 2 row 9 (phi = 90 deg): I = +6.81 %
```

The uniaxial critical-plane normal lies at 34.4° from the bar axis, so the crack trace lies at
55.6°. This agrees with the paper's 56° for row 2-5. Rows 8 and 9 reproduce the printed error
indices. Row 9's critical plane is tilted into the surface, at $\theta_z = -23°$ in the paper.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models NINIC --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`), a repository extension.
- **DP comparison:** supported, against `DP_fit`. The index is homogeneous of degree 1 in the load.
- **Fatigue-limit datasets:** $\mathrm{DP} - 1$ is the paper's error index $I/100$.
- **`_ext`:** identical to `NINIC`.

## Verification and limitations

[`tests/methods/test_ninic.py`](../../tests/methods/test_ninic.py) checks:

- Eq. (4) against $k$ of Table 1, for the four materials;
- Tables 2–4: 29 of the 32 rows evaluated with the full plane search are reproduced within 0.08
  points;
- Table 5 (aluminium 6082-T6) within 0.6 points. Its printed plane stresses imply slightly
  different loads, e.g. $\sigma_{n,\max} = 16$ MPa for $\sigma_a = 15$ MPa.
- the re-evaluation restricted to planes normal to the surface, row 2-11: within 0.04 points of
  the text value $I = 0.0735$ % (the table prints 0.735);
- that the printed plane stresses give the printed index on every row;
- the uniaxial limit of Eq. (4) on both sides of $S_e/T_e = \sqrt2$, the torsion limit, and the
  clipping of compression.

Strict expected failures:

- **Rows 2-7 and 3-12:** the printed plane is not the maximum of $D$; the search finds 0.24 and
  1.8 points more.
- **Row 3-13:** the surface re-evaluation of the row 3-12 load is 1.8 points off, like row 3-12.
- **Row 3-4:** the printed $\tau_{xy,a} = 21.9$ MPa contradicts its printed ratio 0.2. Gonçalves et
  al. (2005) list 62.8 MPa for this test.
- **Row 5-6:** the printed $\sigma_{n,\max}$ exceeds every principal stress of the printed load.

Limitations:

- an endurance criterion;
- Eq. (4) is not consistent for $S_e/T_e < \sqrt2$;
- the out-of-phase shear amplitude is a repository rule, not a definition of the source.

## References

- D. Ninić (2006), A stress-based multiaxial high-cycle fatigue damage criterion, *International
  Journal of Fatigue* 28(2), 103–113. <https://doi.org/10.1016/j.ijfatigue.2005.04.014>
- T. Nishihara, M. Kawamoto (1945), The strength of metals under combined alternating bending and
  torsion with phase difference, *Memoirs of the College of Engineering, Kyoto Imperial
  University* 11, 85–112.
