# Matsubara–Hayashida–Kano unified equivalent shear stress (`MATSUBARA`)

> **Registry name:** `MATSUBARA` · **Catalogue entry:** [N09](../../METHODS.md#n09) ·
> **Source:** Matsubara, Hayashida & Kano (2018), <https://doi.org/10.1016/j.ijfatigue.2017.12.001> ·
> **Implementation:** [`invariant_methods/matsubara.py`](../../src/multifat/models/invariant_methods/matsubara.py)

## Overview

- **Type and family:** stress invariants. An equivalent shear stress amplitude (Mamiya's
  prismatic hull) is combined with a Stage II term, the maximum of the normal stress whose range is
  the largest.
- **Idea:** every stress state is mapped, at the same fatigue life, onto the fully reversed
  uniaxial line $\sqrt{J_{2,amp}}/S_{max} = 1/\sqrt3$. The mapped value, the **unified equivalent
  shear stress amplitude**, is then read on the reversed axial S–N curve with ordinate
  $\sigma_a/\sqrt3$. Only axial data and the true fracture strength are needed: no torsion test.
- **Proposed for:** high-cycle fatigue limits and lives of steels and aluminium alloys under
  combined axial, torsion and biaxial loading, in and out of phase, with mean stresses and with
  sinusoidal, triangular and trapezoidal waveforms.
- **Status:** implemented for sinusoidal loading. The paper prints no per-test value, so the
  evidence is the closed forms (see below).

## Original formulation

**Equivalent shear stress amplitude** (Mamiya, Eqs. 7, 8, 13): with the deviatoric stress $s_{ij}$,

```math
S_1 = \tfrac{\sqrt3}{2}s_{xx},\; S_2 = \tfrac12(s_{yy}-s_{zz}),\; S_3 = s_{xy},\; S_4 = s_{yz},\; S_5 = s_{zx},
\qquad
\sqrt{J_{2,amp}} = \sqrt{\sum_{i=1}^{5}\left(\frac{S_{i,\max} - S_{i,\min}}{2}\right)^2} .
```

**Stage II term** (Section 2.2.3.2): $S_{max}$ is the maximum over the cycle of the normal stress
that has the largest range among all planes of the evaluated area.

**Fatigue limit** (Eq. 9), with the reversed axial limit $\sigma_w$ and the true fracture strength
$\sigma_T$:

```math
\sqrt{J_{2,amp}} + \frac{\sigma_w}{\sqrt3}\,\frac{1}{\sigma_T-\sigma_w}\,S_{max}
 = \frac{\sigma_w}{\sqrt3}\,\frac{\sigma_T}{\sigma_T-\sigma_w} .
```

**Unified equivalent shear stress amplitude** (Eq. 15): the line of constant life through the state
$(S_{max}, \sqrt{J_{2,amp}})$ passes through $(\sigma_T, 0)$ (Fig. 3) and meets the uniaxial line at

```math
\sqrt{J_{2,amp,eq}} = \frac{\sigma_T\sqrt{J_{2,amp}}}{\sigma_T - S_{max} + \sqrt3\,\sqrt{J_{2,amp}}} .
```

**Life** (Section 2.2.3.4.1): the reversed axial S–N curve, with $\sqrt{J_{2,amp,eq}} = \sigma_a/\sqrt3$
on its ordinate, gives $N_f$. Without fatigue data the paper also estimates this curve from
$\sigma_u$ and $\sigma_w$ (Eqs. 16–17). When $\sigma_T$ is unknown and
$400 \le \sigma_u \le 1000$ MPa, it uses $\sigma_T = 1.84\,\sigma_u$.

## Repository implementation

For one harmonic the half-ranges of the $S_i$ are the amplitudes of the components, so
$J_{2,amp} = \tfrac12(\lVert\operatorname{dev}\mathbf S_{\sin}\rVert^2 +
\lVert\operatorname{dev}\mathbf S_{\cos}\rVert^2)$, the Mamiya and Li values coinciding (Section
2.2.3.1). On the plane of unit normal $\mathbf n$ the normal stress has amplitude
$\sqrt{(\mathbf n\cdot\mathbf S_{\sin}\mathbf n)^2 + (\mathbf n\cdot\mathbf S_{\cos}\mathbf n)^2}$ and
mean $\mathbf n\cdot\mathbf S_m\mathbf n$. $S_{max}$ is their sum on the plane of the largest
amplitude, found by a search over the unit sphere (a quasi-uniform grid, the principal directions
of the harmonic parts, and local refinement). The model returns $\sqrt{J_{2,amp,eq}}$ (MPa).

**Repository choices**

- *True fracture strength.* The optional constant `sigma_T`. Without it,
  $\sigma_T = 1.84\,S_u$, only for $400 \le S_u \le 1000$ MPa as the paper states; outside that range
  the material is out of scope (`MaterialScopeError`). Maintainer decision.
- *Life.* The shared master curve $\mathrm{DP} = A N_f^b$ with $\mathrm{DP} = \sqrt{J_{2,amp,eq}}$,
  fitted on the R = −1 uniaxial tension failures. This is the paper's unified equivalent S–N
  curve from the reversed axial S–N curve, because those tests have
  $\sqrt{J_{2,amp,eq}} = \sigma_a/\sqrt3$. The estimate from $S_u$ and $\sigma_w$ (Eqs. 16–17) is
  not implemented.
- *Ties.* When several planes share the largest normal stress range, the largest $S_{max}$ is taken
  (conservative; the paper does not discuss ties). Example: torsion with a mean shear stress.
- *Undefined states.* $S_{max} \ge \sigma_T + \sqrt3\sqrt{J_{2,amp}}$ leaves Eq. (15) undefined and
  raises `LoadingScopeError`.
- *Waveforms.* The harmonic input represents sinusoidal loading only.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_T$ | `sigma_T` / Summary column `sigma_T` | MPa | workbook | optional |
| $\sigma_u$ | `Su` / Summary tensile strength | MPa | workbook | if `sigma_T` is absent |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

$\sigma_w$ does not enter the damage parameter: the fatigue limit is the endurance of the master
curve, and in the fatigue-limit path the reference is the workbook's `Sigm1`.

## Calibration

1. **Original procedure.** $\sigma_T$ from tensile data (or $1.84\,\sigma_u$); the unified S–N curve
   from the reversed axial S–N curve. Nothing is fitted on multiaxial or torsion tests.
2. **Repository.** The same: the master curve is fitted on the R = −1 uniaxial failures, i.e. on
   the paper's own calibration set.

## Worked example

Hard steel of Nishihara & Kawamoto (Table 2: $\sigma_u = 680$, $\sigma_w = 314$, $\tau_w = 196$ MPa).

```python
import numpy as np

from multifat.models.invariant_methods.matsubara import (
    equivalent_shear_amplitude,
    max_range_normal_stress,
    true_fracture_strength,
)
from multifat.models.registry import get_models

matsubara = get_models(["MATSUBARA"])[0]
params = {"Su": 680.0}  # hard steel of Nishihara & Kawamoto (Table 2): sigma_w = 314, tau_w = 196
sigma_T = true_fracture_strength(params)
zero = np.zeros((3, 3))
print(f"sigma_T = 1.84 Su = {sigma_T:.1f} MPa; fatigue limit sigma_w/sqrt3 = {314 / np.sqrt(3):.1f} MPa")


def bt(s, t):
    return np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])


for label, S_sin, S_cos in (
    ("bending 314 MPa", bt(314.0, 0.0), zero),
    ("torsion 196 MPa", bt(0.0, 196.0), zero),
    ("bending 150.2 + torsion 181.7 at 90 deg", bt(150.2, 0.0), bt(0.0, 181.7)),
):
    out = matsubara.evaluate_case(
        S0=None, S1=None, E0=None, E1=None, params=params,
        harmonics={"S_mean": zero, "S_sin": S_sin, "S_cos": S_cos},
    )  # fmt: skip
    print(f"{label:40s} sqrt(J2,amp) = {equivalent_shear_amplitude(S_sin, S_cos):6.1f}, "
          f"S_max = {max_range_normal_stress(zero, S_sin, S_cos):6.1f}, "
          f"sqrt(J2,amp,eq) = {out.values['']:6.1f} MPa")  # fmt: skip
```

Output:

```text
sigma_T = 1.84 Su = 1251.2 MPa; fatigue limit sigma_w/sqrt3 = 181.3 MPa
bending 314 MPa                          sqrt(J2,amp) =  181.3, S_max =  314.0, sqrt(J2,amp,eq) =  181.3 MPa
torsion 196 MPa                          sqrt(J2,amp) =  196.0, S_max =  196.0, sqrt(J2,amp,eq) =  175.8 MPa
bending 150.2 + torsion 181.7 at 90 deg  sqrt(J2,amp) =  201.3, S_max =  199.5, sqrt(J2,amp,eq) =  179.9 MPa
```

The bending limit maps onto itself. The torsion limit gives 175.8 MPa, 3 % below the fatigue
limit $\sigma_w/\sqrt3$, and the 90° out-of-phase fatigue limit 1 % below: the criterion predicts
both from axial data alone.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_Santus --models MATSUBARA --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on $\sqrt{J_{2,amp,eq}}$ (`Nf_expected`), the paper's own route.
- **DP comparison:** supported. $\sqrt{J_{2,amp,eq}}$ is not homogeneous in the load (the
  denominator holds $\sigma_T$), so use `Error_load`.
- **Fatigue-limit datasets:** the index uses the workbook `Sigm1` as reference, consistent with
  Eq. (9).
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_matsubara.py`](../../tests/methods/test_matsubara.py) checks:

- that R = −1 uniaxial states map onto $\sigma_a/\sqrt3$ (Fig. 4) for any $\sigma_T$;
- that states on the Eq. (9) limit line, with and without mean stresses, give $\sigma_w/\sqrt3$;
- Eq. (13) against a time-domain evaluation of the Mamiya components;
- $S_{max}$ against a brute-force search over 200 000 planes, and the tie rule;
- the $\sigma_T$ rules and the undefined-state error.

Not reproduced: the paper reports its verification as error histograms (Figs. 6–7) and life scatter
plots (Figs. 8–12); no per-test value is printed.

Limitations:

- sinusoidal loading only;
- $\sigma_T$ estimated only within $400 \le S_u \le 1000$ MPa;
- the S–N estimate without fatigue data (Eqs. 16–17) is not available.

## References

- G. Matsubara, A. Hayashida, D. Kano (2018), Predicting the multiaxial fatigue limit and the
  multiaxial high-cycle fatigue life based on the unified equivalent shear stress from axial
  strength characteristics with various waveforms, *International Journal of Fatigue* (accepted
  manuscript). <https://doi.org/10.1016/j.ijfatigue.2017.12.001>
- E. N. Mamiya, J. A. Araújo (2002), Fatigue limit under multiaxial loadings: on the definition of
  the equivalent shear stress, *Mechanics Research Communications* 29, 141–151.
