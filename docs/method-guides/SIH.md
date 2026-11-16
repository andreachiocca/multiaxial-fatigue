# Shear stress intensity hypothesis (`SIH`)

> **Registry name:** `SIH` · **Catalogue entry:** [N15](../../METHODS.md#n15) ·
> **Source:** Zenner, Simbürger & Liu (2000), <https://doi.org/10.1016/S0142-1123(99)00107-3>
> (the hypothesis of Zenner, Heidenreich & Richter) ·
> **Implementation:** [`integral_methods/sih.py`](../../src/multifat/models/integral_methods/sih.py)

## Overview

- **Type and family:** integral approach. The squared shear and normal stress amplitudes, with
  their mean-stress corrections, are averaged over all planes through the point.
- **Idea:** Novozhilov showed that the mean square of the shear stress over all planes is
  proportional to $J_2'$. Integrating amplitudes instead of instantaneous stresses lets every plane
  contribute. A phase shift that rotates the principal directions then spreads the shear amplitude
  over more planes, instead of lowering its maximum.
- **Proposed for:** fatigue limits of ductile steels ($R_m$ = 400–1600 MPa) with
  $0.577 \le \tau_W/\sigma_W \le 0.866$. The source covers alternating and static normal and shear
  stresses, phase shifts, different frequencies and non-sinusoidal waveforms. It was compared with
  179 test series (Table 2).
- **Status:** single-harmonic loading. The shear amplitude of a rotating shear vector is not
  defined in the available source (see the repository choices).

## Original formulation

**Equivalent stress amplitude** (Eq. 6). $\gamma$ and $\varphi$ are the polar angles of the plane
normal. $\tau_{\gamma\varphi,a}$, $\tau_{\gamma\varphi,m}$, $\sigma_{\gamma\varphi,a}$ and
$\sigma_{\gamma\varphi,m}$ are the shear and normal stress amplitudes and means on that plane.

```math
\sigma_{eq,a} = \left\{\frac{15}{8\pi}\int_0^{\pi}\!\!\int_0^{2\pi}
\left[a\,\tau_{\gamma\varphi,a}^2\left(1 + m\,\tau_{\gamma\varphi,m}^2\right)
 + b\,\sigma_{\gamma\varphi,a}^2\left(1 + n\,\sigma_{\gamma\varphi,m}\right)\right]
\sin\gamma\,d\gamma\,d\varphi\right\}^{1/2},
\qquad \sigma_{eq,a} = \sigma_W \text{ at the fatigue limit.}
```

**Constants** (Table 1). They come from the reversed axial and torsion limits $\sigma_W$,
$\tau_W$ and the pulsating (R = 0) limits $\sigma_{Sch}$, $\tau_{Sch}$, written as maximum
stresses:

```math
a = \tfrac15\left[3\left(\frac{\sigma_W}{\tau_W}\right)^2 - 4\right],\qquad
b = \tfrac15\left[6 - 2\left(\frac{\sigma_W}{\tau_W}\right)^2\right],
```

```math
a\,m = \frac{\sigma_W^2 - (\sigma_W/\tau_W)^2(\tau_{Sch}/2)^2}{\tfrac{12}{7}(\tau_{Sch}/2)^4},\qquad
b\,n = \frac{\sigma_W^2 - (\sigma_{Sch}/2)^2 - \tfrac{4}{21}\,a\,m\,(\sigma_{Sch}/2)^4}{\tfrac{15}{14}(\sigma_{Sch}/2)^3}.
```

**Pulsating torsion limit** (Eq. 7), when it is not measured:

```math
\tau_{Sch} = \frac{4\,\tau_W}{2\sigma_W/\sigma_{Sch} + 1}.
```

With $a = 1$ and $b = 0$, Eq. (6) without means reduces to Eq. (5), the shear-only integral. This
gives $\tau_W/\sigma_W = 1/\sqrt3$. For constant principal directions it equals
$\tfrac{3}{\sqrt2}\tau_{oct,a}$ (Eq. 5b).

## Repository implementation

On each plane of normal $\mathbf n$, the normal stress of one harmonic has amplitude
$\sqrt{(\mathbf n\cdot\mathbf S_{\sin}\mathbf n)^2 + (\mathbf n\cdot\mathbf S_{\cos}\mathbf n)^2}$ and
mean $\mathbf n\cdot\mathbf S_m\mathbf n$. The shear vector traces an ellipse. The integral is
evaluated on a 32 × 64 Gauss–Legendre × azimuth grid. The factor $15/(8\pi)$ times the sphere area
is 7.5. The model returns the index $\sigma_{eq,a}/\sigma_W$, which equals 1 on the fatigue limit.

**Repository choices**

- *Shear amplitude.* When the shear vector rotates, the available source does not define
  $\tau_{\gamma\varphi,a}$; it cites Liu (1991) and Liu & Zenner (1993). By maintainer decision the
  model uses half the longest chord of the ellipse, i.e. its minimum circumscribed circle radius
  (`np_shear_rule = "longest_chord"`). The quadratic amplitude $\sqrt{s^2 + c^2}$ (`"rss"`) is
  selectable. **Evidence:** Section 3 states that the integral of Eq. (5) does not depend on the
  phase between $\sigma_x$ and $\tau_{xy}$ (Fig. 2). That holds exactly with `"rss"`. With the
  longest chord, Eq. (5) drops by 7 % at 90° for $\tau_{xy,a} = \sigma_{x,a}/2$. The maintainer confirmed the
  longest chord after reviewing this evidence (16 November 2026).
- *Mean shear.* $\tau_{\gamma\varphi,m}$ is the magnitude of the mean shear vector.
- *Mean stresses.* They need $\sigma_{Sch}$ (workbook `Sig0`); without it, a case with a mean
  stress raises `LoadingScopeError`. $\tau_{Sch}$ is `Tau0`, otherwise Eq. (7).
- *Material scope.* $\tau_W/\sigma_W$ outside $[0.577, 0.866]$ raises `MaterialScopeError`
  (Section 4).
- *Clipping.* A negative integral (large compressive mean) is set to zero.
- *Life.* The source is a fatigue-limit criterion. Finite lives use the shared master curve on the
  index. Section 6 says lives need damage sums per plane instead; that is not implemented.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_W$ | `Sigm1` | MPa | workbook | yes |
| $\tau_W$ | `Taum1` | MPa | workbook | yes |
| $\sigma_{Sch}$ | `Sig0` (maximum stress of the R = 0 limit) | MPa | workbook | for mean stresses |
| $\tau_{Sch}$ | `Tau0` (maximum stress of the R = 0 limit) | MPa | workbook | optional (Eq. 7) |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $a$ and $b$ make the reversed axial and torsion limits exact; $m$ and
   $n$ make the pulsating limits exact. The source says that $\sigma_{Sch}$ can be estimated from
   the mean stress sensitivity $M$, but gives no formula.
2. **Repository.** The same, from the workbook limits. Nothing is fitted. The master curve for
   finite lives is fitted on the R = −1 uniaxial failures.

## Worked example

Hard steel of Nishihara & Kawamoto. These are the fully reversed limits and two combined limits of
Papadopoulos et al. (1997), Table 1.

```python
import numpy as np

from multifat.models.integral_methods.sih import sih_constants, sih_index

# Hard steel of Nishihara & Kawamoto (Papadopoulos et al. 1997, Table 1): fully reversed limits
f_W, t_W = 313.9, 196.2
k = sih_constants(f_W, t_W)
print(f"tau_W/sigma_W = {t_W / f_W:.3f}: a = {k['a']:.4f}, b = {k['b']:.4f}")


def bending_torsion(s_a, t_a, delta_deg):
    """Harmonic parts of sigma_yy = s_a sin wt and tau_xy = t_a sin(wt + delta)."""
    d = np.radians(delta_deg)
    S_sin, S_cos = np.zeros((3, 3)), np.zeros((3, 3))
    S_sin[1, 1] = s_a
    S_sin[0, 1] = S_sin[1, 0] = t_a * np.cos(d)
    S_cos[0, 1] = S_cos[1, 0] = t_a * np.sin(d)
    return np.zeros((3, 3)), S_sin, S_cos


for s_a, t_a, delta in ((313.9, 0, 0), (0, 196.2, 0), (138.1, 167.1, 0), (150.2, 181.7, 90)):
    parts = bending_torsion(s_a, t_a, delta)
    lc = sih_index(*parts, f_W, t_W)  # longest chord (default)
    rss = sih_index(*parts, f_W, t_W, shear_rule="rss")
    print(f"sigma_a = {s_a:5.1f}, tau_a = {t_a:5.1f}, delta = {delta:2d}: "
          f"index {lc:.3f} (longest chord), {rss:.3f} (rss)")  # fmt: skip
```

Output:

```text
tau_W/sigma_W = 0.625: a = 0.7358, b = 0.1761
sigma_a = 313.9, tau_a =   0.0, delta =  0: index 1.000 (longest chord), 1.000 (rss)
sigma_a =   0.0, tau_a = 196.2, delta =  0: index 1.000 (longest chord), 1.000 (rss)
sigma_a = 138.1, tau_a = 167.1, delta =  0: index 0.959 (longest chord), 0.959 (rss)
sigma_a = 150.2, tau_a = 181.7, delta = 90: index 0.998 (longest chord), 1.042 (rss)
```

The two reversed limits are exact by construction. In phase, both amplitudes agree, since the shear
path is a segment. At 90° out of phase, the measured limit gives 0.998 with the longest chord and
1.042 with the quadratic amplitude.

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`); a repository extension.
- **DP comparison:** supported without mean stresses, where the index is homogeneous in the load.
  With means, use `Error_load`.
- **Fatigue-limit datasets:** the index is the fatigue index directly.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_sih.py`](../../tests/methods/test_sih.py) checks:

- that the Table 1 factors 12/7, 4/21 and 15/14 are the sphere integrals of the uniaxial and
  torsional plane stresses;
- that the four calibration limits give index 1 exactly, with given or estimated $\tau_{Sch}$, for
  both amplitude rules;
- Eq. (5): $a = 1$, $b = 0$ at $\tau_W/\sigma_W = 1/\sqrt3$;
- Eq. (5b): the octahedral form for proportional loading;
- the Section 3 phase independence, exact with `"rss"`, and the 7 % drop with the longest chord;
- the scope errors, and that the model follows its `np_shear_rule`.

Not reproduced: the paper prints only statistics of 179 test series (Table 2), not the series.

Limitations:

- the original shear amplitude definition for rotating shear vectors is not available;
- mean stresses need `Sig0`;
- material scope $0.577 \le \tau_W/\sigma_W \le 0.866$.

## References

- H. Zenner, A. Simbürger, J. Liu (2000), On the fatigue limit of ductile metals under complex
  multiaxial loading, *International Journal of Fatigue* 22, 137–145.
  <https://doi.org/10.1016/S0142-1123(99)00107-3>
- R. Heidenreich, H. Zenner, I. Richter (1983), Schubspannungsintensitätshypothese —
  Dauerschwingfestigkeit bei mehrachsiger Beanspruchung, *Forschungshefte FKM* 105 (not consulted).
- J. Liu, H. Zenner (1993), Berechnung der Dauerschwingfestigkeit bei mehrachsiger Beanspruchung,
  *Mat.-wiss. u. Werkstofftech.* 24, 240–249 (not consulted).
