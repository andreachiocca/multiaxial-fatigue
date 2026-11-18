# Modified Kandil–Brown–Miller (`MKBM`)

> **Registry name:** `MKBM` (and `MKBM_ext`) · **Catalogue entry:** [I_MKBM](../../METHODS.md#i_mkbm) ·
> **Source:** Li, Zhang, Sun & Li (2011), <https://doi.org/10.1016/j.ijfatigue.2010.07.003> ·
> **Implementation:** [`cp_methods/mkbm.py`](../../src/multifat/models/cp_methods/mkbm.py)

## Overview

- **Type and family:** strain-based critical-plane criterion (shear strain plus normal strain,
  with a normal-stress correction).
- **Idea:** the Kandil–Brown–Miller parameter combines the maximum shear strain and the normal
  strain on the same plane. Li et al. multiply the normal-strain term by
  $1 + \sigma_{n,\max}/\sigma_y$, so that the additional hardening and the crack opening produced by
  non-proportional loading also enter the parameter.
- **Proposed for:** low-cycle multiaxial fatigue of several metals (steels, aluminium and titanium
  alloys) under proportional and out-of-phase tension–torsion.
- **Status:** implemented for the declared repository scope ("requires review" in the catalogue).
  The plane strains of the thin-walled tube (Eqs. 7–8, 13, 15–16) and Eq. (27) are reproduced. No
  experimental dataset of the paper is reproduced.

## Original formulation

Li et al. (2011, Eq. 27) define an equivalent strain amplitude on the plane of maximum shear
strain range:

```math
\frac{\Delta\varepsilon^*_{eq}}{2} = \frac{\Delta\gamma_{\max}}{2}
+ \left(1 + \frac{\sigma_{n,\max}}{\sigma_y}\right)\frac{\Delta\varepsilon_n}{2},
```

where:

- $\Delta\gamma_{\max}$ is the maximum **engineering** shear strain range;
- $\Delta\varepsilon_n$ is the normal strain range on that plane;
- $\sigma_{n,\max}$ is the maximum normal stress on that plane;
- $\sigma_y$ is the yield strength.

For general paths the shear range is the longest chord of the shear-strain path on the plane
(Eq. 19). The equivalent strain amplitude is related to life through the uniaxial strain–life
curve of the material.

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\mathrm{DP} = \gamma_a + \left(1 + \frac{\sigma_{n,\max}}{S_y}\right)\varepsilon_{n,a},
\qquad \sigma_{n,\max} = \max(\sigma_{n,m} + \sigma_{n,a},\,0),
```

where:

- $\gamma_a = \Delta\gamma/2$ is the engineering shear strain amplitude (twice the norm of the
  tensorial shear range, halved);
- $\varepsilon_{n,a} = \Delta\varepsilon_n/2$.

DP is dimensionless.

- `MKBM`: DP on the plane of largest $\Delta\gamma$ (the source's plane rule). Among tied planes,
  the one with the larger DP is used.
- `MKBM_ext`: the largest DP over all planes. This is a repository variant.
- Out-of-phase loading: $\Delta\gamma$ is the longest chord of the shear-strain path
  (`np_shear_rule = "longest_chord"`), as in Eq. (19).

**Repository choices**

- *Compression clipped.* $\sigma_{n,\max} < 0$ is set to 0, so a compressive normal stress never
  lowers the factor below 1.
- *Monotonic yield strength.* The workbook `Sy` is used for $\sigma_y$.
- *Elastic strains.* Strains follow from the stresses through isotropic linear elasticity; the
  cyclic plasticity and additional hardening that motivate the original correction are not
  simulated.
- *Life.* The shared master curve $\mathrm{DP} = A N_f^b$ replaces the strain–life curve
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_y$ | `Sy` / Summary yield strength | MPa | workbook | yes, $> 0$ |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook | yes (strains) |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

MKBM has no criterion-specific constant. The only fitted quantities are the master-curve constants
$A$, $b$ from the uniaxial $R = -1$ design tests. In the original formulation, the strain–life
constants $\sigma_f'$, $b$, $\varepsilon_f'$, $c$ play that role.

With elastic strains and no mean stress, the parameter fixes the ratio of the torsion and tension
fatigue limits:

- tension, 45° plane: $\mathrm{DP} = (1+\nu)\frac{f}{E} + \left(1 + \frac{f}{2S_y}\right)\frac{(1-\nu)f}{2E}$;
- torsion: $\mathrm{DP} = 2(1+\nu)t/E$.

The ratio therefore depends on $f/S_y$, which is about 0.68 in the example below. Torsion tests of
materials with a different ratio are expected to lie off the tension master curve.

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$,
$S_y = 400$ MPa.

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


Sy = 400.0  # MPa (synthetic)
mkbm = get_models(["MKBM"])[0]
params = {"Sy": Sy, "E": E, "nu": nu}

s_a, t_a = 300.0, 200.0
dp_tension, ext_tension = evaluate(mkbm, params, zero, s_a * X)
dp_torsion, _ = evaluate(mkbm, params, zero, t_a * T)

# Tension, 45 deg plane: gamma_a = (1+nu) s/E, eps_n,a = (1-nu) s/(2E), sigma_n,max = s/2
tension = (1 + nu) * s_a / E + (1 + s_a / 2 / Sy) * (1 - nu) * s_a / (2 * E)
# Torsion, plane of tau_xy: gamma_a = tau/G, no normal strain or stress
torsion = 2 * (1 + nu) * t_a / E
assert abs(dp_tension - tension) < 1e-9 * tension
assert abs(dp_torsion - torsion) < 1e-9 * torsion

# Fatigue-limit ratio t/f predicted by MKBM for f = 300 MPa: equal DP in tension and torsion
f = s_a
ratio = (tension * E / (2 * (1 + nu))) / f

# 90 deg out-of-phase tension-torsion: sigma_a = 200, tau_a = 100 MPa
dp_oop, ext_oop = evaluate(mkbm, params, zero, 200.0 * X, 100.0 * T)

print(f"MKBM tension 300 MPa: {dp_tension:.6f} (MKBM_ext {ext_tension:.6f})")
print(f"MKBM torsion 200 MPa: {dp_torsion:.6f}")
print(f"t/f at equal DP for f = 300 MPa: {ratio:.3f}")
print(f"MKBM 90 deg out of phase: {dp_oop:.6f} (MKBM_ext {ext_oop:.6f})")
```

Output:

```text
MKBM tension 300 MPa: 0.002672 (MKBM_ext 0.003372)
MKBM torsion 200 MPa: 0.002600
t/f at equal DP for f = 300 MPa: 0.685
MKBM 90 deg out of phase: 0.002800 (MKBM_ext 0.002800)
```

**Tension.** `MKBM_ext` is about 26 % higher than `MKBM`. The largest DP lies on a plane closer to
the load axis, where the normal strain and stress are larger and the shear strain is slightly lower.

**Out-of-phase tension–torsion** ($\sigma_a = 2\tau_a$, 90°). Every plane containing the
specimen's radial direction has the same shear range. The tie is resolved towards the largest DP,
which is the plane normal to the axis: there $\varepsilon_{n,a} = \sigma_a/E$ and
$\sigma_{n,\max} = \sigma_a$. As a result, `MKBM` and `MKBM_ext` coincide.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MKBM --all-variants --results-dir Results/current/guides
```

MKBM is skipped on a workbook without `Sy`, `E` or `nu`; the reason is written to
`<MATERIAL>_run_report.tsv`.

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`), not the strain–life curve.
- **DP comparison:** supported, against the master-curve `DP_fit`. DP is a sum of strain terms with
  a stress-dependent factor, so it is not homogeneous in the load. The load-equivalent error
  `Error_load` is the scale-free metric to compare with other methods.
- **`_ext`:** available; maximum DP over planes, a repository variant.

## Verification and limitations

[`tests/methods/test_mkbm.py`](../../tests/methods/test_mkbm.py) checks:

- the thin-walled tube plane strains, Eqs. (7)–(8), at every instant of the cycle;
- the harmonic evaluator against Eqs. (15)–(16);
- the critical-plane angle of Eq. (13);
- the normal-stress factor of Eq. (27), including clipping of compression;
- the longest-chord shear range of Eq. (19).

The paper swaps the mean-strain terms of Eqs. (9) and (10); this misprint is recorded in the test
module.

Not verified: reproduction of the experimental lives of Li et al. (2011).

Limitations:

- elastic strains only, so the non-proportional hardening captured by the original parameter
  through the measured stresses is absent;
- monotonic instead of cyclic yield strength.

## References

- J. Li, Z. Zhang, Q. Sun, C. Li (2011), Multiaxial fatigue life prediction for various metallic
  materials based on the critical plane approach, *International Journal of Fatigue* 33(2),
  90–101. <https://doi.org/10.1016/j.ijfatigue.2010.07.003>
- M. W. Brown, K. J. Miller (1973), A theory for fatigue failure under multiaxial stress–strain
  conditions, *Proceedings of the Institution of Mechanical Engineers* 187, 745–755.
- F. A. Kandil, M. W. Brown, K. J. Miller (1982), Biaxial low-cycle fatigue fracture of 316
  stainless steel at elevated temperatures, *Mechanical Behaviour and Nuclear Applications of
  Stainless Steel at Elevated Temperatures*, The Metals Society, London, 203–210.
