# Khalij–Pagnacco–Troian improved Gough / Nishihara–Kawamoto criterion (`KHALIJ`)

> **Registry name:** `KHALIJ` · **Catalogue entry:** [N19](../../METHODS.md#n19) ·
> **Source:** Khalij, Pagnacco & Troian (2017), <https://doi.org/10.1016/j.ijfatigue.2017.02.016> ·
> **Implementation:** [`invariant_methods/khalij.py`](../../src/multifat/models/invariant_methods/khalij.py)

## Overview

- **Type and family:** stress invariants. An equivalent shear stress amplitude, measured by the
  prismatic hull in the principal coordinate system (PHPS), is combined with the amplitude of the
  first stress invariant.
- **Idea:** the Gough and Nishihara–Kawamoto bending–torsion criteria are generalised to plane
  stress. The weight between the shear and hydrostatic terms depends on the limit ratio
  $r = t_{-1}/s_{-1}$, so one formula moves from a von Mises ellipse (mild metals) towards a
  Rankine-like envelope (cast irons).
- **Proposed for:** high-cycle fatigue limits of metals from mild steels to brittle cast irons,
  $0.5 \lesssim r \le 1$, under in-phase and out-of-phase sinusoidal plane stress.
- **Status:** implemented for zero-mean loading (the authors' "OP" setting). Table 6 is reproduced;
  the mean-stress extension is not implemented (see below).

## Original formulation

For sinusoidal plane stress with amplitudes $s_{xx,a}$, $s_{yy,a}$, $s_{xy,a}$ and a phase shift
$\phi_{yy}$ between the two normal stresses, the criterion is (Eq. 23)

```math
\mathfrak D = \alpha(r)\left[\left(\frac{s_{xx,a}}{s_{-1}}\right)^2 + \left(\frac{s_{yy,a}}{s_{-1}}\right)^2
 - \frac{s_{xx,a}}{s_{-1}}\frac{s_{yy,a}}{s_{-1}}\cos\phi_{yy}\right]
 + \left(\frac{s_{xy,a}}{t_{-1}}\right)^2 + \bigl(1-\alpha(r)\bigr)\,\frac{I_{1,a}}{s_{-1}} \le 1,
```

with the amplitude of the first invariant (Eq. 22)

```math
I_{1,a} = \sqrt{s_{xx,a}^2 + s_{yy,a}^2 + 2\,s_{xx,a}s_{yy,a}\cos\phi_{yy}} .
```

The phase of the shear stress does not enter: the PHPS amplitude of Eq. (12) depends only on
$\phi_{yy}$. The authors' setting ("OP") makes the criterion independent of the coordinate system
for pure shear (Eqs. 24–26):

```math
\alpha(r) = \frac{1}{3r^2}, \qquad r = \frac{t_{-1}}{s_{-1}} .
```

The Gough ($\alpha = 1/r - 1$) and Nishihara–Kawamoto ($\alpha = \tfrac12(1/r^2 - 1)$) settings are
used in the paper only for comparison.

With $\alpha = 1/(3r^2)$ the first two groups equal $J_{2,a}/t_{-1}^2$, where
$J_{2,a} = \tfrac13(s_{xx,a}^2 + s_{yy,a}^2 - s_{xx,a}s_{yy,a}\cos\phi_{yy}) + s_{xy,a}^2$ is the
PHPS amplitude (Eq. 12). The criterion is therefore

```math
\mathfrak D = \frac{J_{2,a}}{t_{-1}^2} + \left(1 - \frac{1}{3r^2}\right)\frac{I_{1,a}}{s_{-1}} \le 1 .
```

## Repository implementation

The model evaluates the invariant form above from the harmonic stress tensors:
$J_{2,a} = \tfrac12\left(\lVert\operatorname{dev}\mathbf S_{\sin}\rVert^2 +
\lVert\operatorname{dev}\mathbf S_{\cos}\rVert^2\right)$ and
$I_{1,a} = \sqrt{(\operatorname{tr}\mathbf S_{\sin})^2 + (\operatorname{tr}\mathbf S_{\cos})^2}$.
In plane stress both equal the paper's expressions (a test checks Eq. 23 term by term). It returns
$\mathfrak D$, which is 1 on the limit.

**Loading scope.** A case with a nonzero mean stress, or with an out-of-plane ($z$) stress
component, raises a `LoadingScopeError`. The model is then **skipped** for that material, with
the reason in the run report, and the eligibility matrix records `out_of_scope`.

**Repository choices**

- *Setting.* Only the OP setting, $\alpha = 1/(3r^2)$, the authors' proposal.
- *No mean stresses.* The paper extends the criterion to mean stresses with Liu's correction
  (Eqs. 29–31), using the pulsating limits $s_0$ and $\tau_0$. That part is not implemented, by
  maintainer decision:
  - the printed term $\gamma_2(r)\,s_{ms}/s_{-1}^2$ is not dimensionless, because $s_{ms}$ is
    linear in the mean stresses;
  - Eq. (30) for $\gamma_2$ reproduces the pulsating tension limit only for a purely quadratic base
    criterion, not for Eq. (23);
  - none of four readings (as printed, $/s_{-1}$, $\times s_{xx,a}$, $\times$ amplitude norm; means
    fixed or scaled with the amplitudes) reproduces Table 7 (at best 2 of 18 OP rows).
- *Ratio.* $r = t_{-1}/s_{-1}$ of the material. Table 6 was computed with this unrounded ratio, not
  with the two-digit $r$ printed beside it.
- *Finite life.* $\mathfrak D$ is used as the damage parameter of the shared master curve
  $\mathrm{DP} = A N_f^b$, a repository extension. Because of the linear $I_{1,a}$ term,
  $\mathfrak D$ is not homogeneous in the load, so the load-equivalent error is the appropriate
  cross-method metric.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $s_{-1}$ | `Sigm1` / Summary R = -1 tension or bending limit | MPa | workbook | yes |
| $t_{-1}$ | `Taum1` / Summary R = -1 torsion limit | MPa | workbook | yes |
| $\boldsymbol\sigma(t)$ | `S_mean` (zero), `S_sin`, `S_cos` (plane stress) | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $\alpha$ follows from $r = t_{-1}/s_{-1}$; nothing is fitted.
2. **Repository.** The same, with the workbook limits. Only the master-curve constants are fitted.

For $r < 1/\sqrt3$ (mild metals), $\alpha > 1$ and the weight of $I_{1,a}$ is negative: a hydrostatic
amplitude then lowers $\mathfrak D$. This is the published setting and is not altered.

## Worked example

Three in-phase fatigue limits of Table 6, with the limits of Table 1.

```python
import numpy as np

from multifat.models.invariant_methods.khalij import khalij_alpha
from multifat.models.registry import get_models

khalij = get_models(["KHALIJ"])[0]
zero = np.zeros((3, 3))


def damage(params, sxx, syy, sxy):
    """KHALIJ damage of in-phase, zero-mean plane stress amplitudes (MPa)."""
    S = np.array([[sxx, sxy, 0.0], [sxy, syy, 0.0], [0.0, 0.0, 0.0]])
    out = khalij.evaluate_case(
        S0=S, S1=-S, E0=None, E1=None, params=params,
        harmonics={"S_mean": zero, "S_sin": S, "S_cos": zero},
    )  # fmt: skip
    return out.values[""]


# Table 1 limits and three Table 6 rows (s_xx,a/s_-1, s_yy,a/s_-1, s_xy,a/t_-1; printed OP Damage)
for name, s_m1, t_m1, (x, y, t, printed) in (
    ("0.1%C steel", 235.4, 137.3, (0.43, 0.0, 0.89, 0.98)),
    ("CrV-V", 530.0, 309.0, (0.24, -0.84, 0.0, 0.95)),
    ("Silal cast iron", 205.5, 177.9, (0.66, 0.0, 0.67, 1.01)),
):
    params = {"Sigm1": s_m1, "Taum1": t_m1}
    D = damage(params, x * s_m1, y * s_m1, t * t_m1)
    print(f"{name:16s} r = {t_m1 / s_m1:.3f}, alpha = {khalij_alpha(s_m1, t_m1):.3f}: "
          f"D = {D:.3f} (Table 6: {printed})")  # fmt: skip
```

Output:

```text
0.1%C steel      r = 0.583, alpha = 0.980: D = 0.982 (Table 6: 0.98)
CrV-V            r = 0.583, alpha = 0.981: D = 0.958 (Table 6: 0.95)
Silal cast iron  r = 0.866, alpha = 0.445: D = 1.009 (Table 6: 1.01)
```

The three limits give $\mathfrak D$ between 0.96 and 1.01, each within one unit of the last digit
printed in Table 6.
For the cast iron, $\alpha = 0.445$: the first-invariant term carries more than half of the
damage, which is what lets the criterion follow brittle metals.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_FABEST --models KHALIJ --results-dir Results/current/guides
```

FABEST holds zero-mean tests only. On a workbook with any mean stress the model is skipped as out
of scope, with the reason in the run report.

## Comparison capabilities and outputs

- **Life:** master curve on $\mathfrak D$ (`Nf_expected`), a repository extension.
- **DP comparison:** supported. $\mathfrak D$ is not homogeneous in the load; use `Error_load`.
- **Fatigue-limit datasets:** $\mathfrak D - 1$ is the paper's own measure (Table 6).
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_khalij.py`](../../tests/methods/test_khalij.py) checks:

- Table 6, OP column: 51 of 53 rows to one unit of the last printed digit. Two 3%C cast-iron rows
  give 1.035 and 1.063 against the printed 1.06 and 1.01, unexplained; they are strict xfails;
- that the printed two-digit $r$ would miss 7 rows, so the unrounded ratio is the one used;
- that the invariant form equals Eq. (23) for several $\phi_{yy}$, with a phase-shifted shear that
  must not matter;
- both limits, pure shear in its principal frame, and invariance under an in-plane rotation;
- the rejection of mean stresses and of out-of-plane components.

Limitations:

- zero-mean loading only;
- plane stress only;
- a fatigue-limit criterion, so finite lives are available only through the master curve.

## References

- L. Khalij, E. Pagnacco, R. Troian (2017), Fatigue criterion improvement of Gough and Nishihara &
  Kawamoto to predict the fatigue damage of a wide range of metallic materials, *International
  Journal of Fatigue* 99, 137–150. <https://doi.org/10.1016/j.ijfatigue.2017.02.016>
- L. Khalij, E. Pagnacco, S. Lambert (2010), A measure of the equivalent shear stress amplitude
  from a prismatic hull in the principal coordinate system, *International Journal of Fatigue* 32,
  1977–1984.
