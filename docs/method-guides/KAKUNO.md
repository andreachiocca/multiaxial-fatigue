# Kakuno–Kawada octahedral criterion (`KAKUNO`)

> **Registry name:** `KAKUNO` · **Catalogue entry:** [G04](../../METHODS.md#g04) ·
> **Source:** Kakuno & Kawada (1979), <https://doi.org/10.1111/j.1460-2695.1979.tb01358.x> ·
> **Implementation:** [`invariant_methods/kakuno_kawada.py`](../../src/multifat/models/invariant_methods/kakuno_kawada.py)

## Overview

- **Type and family:** stress invariants (octahedral stresses), a correction of Sines' criterion.
- **Idea:** the octahedral shear stress amplitude at the fatigue limit decreases linearly with the
  octahedral normal stress, but with **different slopes for its amplitude and its mean**. Sines
  keeps only the mean term, which forces $t_{-1}/f_{-1} = 1/\sqrt3$ without mean stresses.
- **Proposed for:** fatigue limits of round bars under combined repeated bending and torsion with
  superimposed static bending and torsion; plain and moderately notched members.
- **Status:** implemented for proportional loading. The paper validates the criterion only on
  figures, so the evidence is the closed forms and Fig. 5 (see below).

## Original formulation

With principal stress amplitudes $\sigma_{a_i}$ and principal mean stresses $\sigma_{m_i}$
(Eq. 4):

```math
\tfrac13\sqrt{(\sigma_{a_1}-\sigma_{a_2})^2 + (\sigma_{a_2}-\sigma_{a_3})^2 + (\sigma_{a_3}-\sigma_{a_1})^2}
 = A - B\,(\sigma_{a_1}+\sigma_{a_2}+\sigma_{a_3}) - C\,(\sigma_{m_1}+\sigma_{m_2}+\sigma_{m_3}) .
```

The constants follow from the fully reversed torsion limit $\tau_w$, the fully reversed bending limit
$\sigma_{wb}$ and the pulsating bending limit $\sigma_{up}$ (Eqs. 7–9):

```math
A = \sqrt{\tfrac23}\,\tau_w, \qquad
B = \frac{\sqrt2}{3}\left(\sqrt3\,\frac{\tau_w}{\sigma_{wb}} - 1\right), \qquad
C = \sqrt{\tfrac23}\,\frac{1-k_1}{k_1}\,\frac{\tau_w}{\sigma_{wb}}, \qquad k_1 = \frac{\sigma_{up}}{2\sigma_{wb}} .
```

For bending $\sigma_a$, $\sigma_m$ and torsion $\tau_a$, $\tau_m$ this becomes Eq. (10):

```math
\left(\frac{\tau_a}{\tau_w}\right)^2\frac{1}{p^2}
+ \left(\frac{\sigma_a}{\sigma_{wb}}\right)^2\frac{\frac{2}{\sqrt3}\frac{\sigma_{wb}}{\tau_w}-1}{p^2}
+ \frac{\sigma_a}{\sigma_{wb}}\,\frac{2-\frac{2}{\sqrt3}\frac{\sigma_{wb}}{\tau_w}}{p} = 1,
\qquad p = 1 - \frac{\sigma_m}{\sigma_{wb}}\,\frac{1-k_1}{k_1} .
```

The torsional mean stress does not enter, as observed for plain specimens. Without mean stresses
$p = 1$ (Eq. 11). The paper recommends $k_1 = 0.83$ as the average of structural carbon and alloy
steels.

## Repository implementation

The model returns the fatigue index

```math
\frac{\tau_{oct,a} + B\,\lvert\operatorname{tr}\mathbf S_a\rvert + C\operatorname{tr}\mathbf S_m}{A},
\qquad \tau_{oct,a} = \sqrt{\tfrac23}\sqrt{J_2(\mathbf S_a)},
```

which is 1 on the limit; $\mathbf S_a$ is the amplitude tensor and $\mathbf S_m$ the mean tensor.

**Loading scope.** The principal stress amplitudes of Eq. (4) exist only for proportional loading.
A case whose harmonic parts $\mathbf S_{\sin}$ and $\mathbf S_{\cos}$ are not parallel raises a
`LoadingScopeError`, so the model is **skipped** for that material (`out_of_scope`). The mean
stress may have any direction.

**Repository choices**

- *Amplitude of the first invariant:* $\lvert\operatorname{tr}\mathbf S_a\rvert$, positive like any
  amplitude.
- *Limits:* $\sigma_{wb}$, $\tau_w$ and $\sigma_{up}$ are the workbook values `Sigm1`, `Taum1` and
  `Sig0` (the maximum stress of the R = 0 cycle). The paper uses plane-bending limits; the workbook
  records whether its axial limit is tension or bending.
- *Without `Sig0`:* $C$ is undefined. A case with a mean normal stress
  ($\operatorname{tr}\mathbf S_m \ne 0$) is out of scope; cases without one are evaluated. The
  paper's $k_1 = 0.83$, stated for steels, is not used as a default (maintainer decision).
- *Finite life:* the index is used as the damage parameter of the shared master curve
  $\mathrm{DP} = A N_f^b$, a repository extension. The index is not homogeneous in the load (linear
  terms), so the load-equivalent error is the appropriate cross-method metric.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_{wb}$ | `Sigm1` / Summary R = -1 bending or tension limit | MPa | workbook | yes |
| $\tau_w$ | `Taum1` / Summary R = -1 torsion limit | MPa | workbook | yes |
| $\sigma_{up}$ | `Sig0` / Summary R = 0 limit (maximum stress) | MPa | workbook | for mean normal stresses |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` (proportional) | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $A$, $B$, $C$ from the three limits (Eqs. 7–9); nothing is fitted.
2. **Repository.** The same, with the workbook limits. Only the master-curve constants are fitted.

$B < 0$ when $\tau_w/\sigma_{wb} < 1/\sqrt3$, and $C < 0$ when $k_1 > 1$: the published constants are
used as they come.

## Worked example

The Ni–Cr–Mo steel of Fig. 5, with stresses normalised by $\sigma_{wb}$.

```python
import numpy as np

from multifat.models.invariant_methods.kakuno_kawada import kakuno_constants
from multifat.models.registry import get_models

kakuno = get_models(["KAKUNO"])[0]
# Ni-Cr-Mo steel of Fig. 5, normalised by s_wb: t_w/s_wb = 0.617, k_1 = s_up/(2 s_wb) = 0.83
s_wb = 1.0
params = {"Sigm1": s_wb, "Taum1": 0.617, "Sig0": 2 * 0.83 * s_wb}
A, B, C = kakuno_constants(params["Sigm1"], params["Taum1"], params["Sig0"])
print(f"A = {A:.4f}, B = {B:.4f}, C = {C:.4f}  (Eqs. 7-9, stresses / s_wb)")


def index(sigma_a, tau_a, sigma_m, tau_m=0.0):
    """KAKUNO index of in-phase bending-torsion with static bending and torsion."""
    bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])  # noqa: E731
    out = kakuno.evaluate_case(
        S0=None, S1=None, E0=None, E1=None, params=params,
        harmonics={"S_mean": bt(sigma_m, tau_m), "S_sin": bt(sigma_a, tau_a),
                   "S_cos": np.zeros((3, 3))},
    )  # fmt: skip
    return out.values[""]


# Fig. 5 (sigma_m/s_wb = 0.46): the Eq. (10) curve meets the tau axis at tau_a = p t_w
p = 1 - 0.46 * (1 - 0.83) / 0.83
print(f"tau_a = {p * 0.617:.3f}, sigma_m = 0.46: index {index(0.0, p * 0.617, 0.46):.4f}")
print(f"same with tau_m = 0.46 t_w: index {index(0.0, p * 0.617, 0.46, 0.46 * 0.617):.4f}")
print(f"sigma_a = 0.67, tau_a = 0.44 (Gough point): index {index(0.67, 0.44, 0.46):.3f}")
```

Output:

```text
A = 0.5038, B = 0.0324, C = 0.1032  (Eqs. 7-9, stresses / s_wb)
tau_a = 0.559, sigma_m = 0.46: index 1.0000
same with tau_m = 0.46 t_w: index 1.0000
sigma_a = 0.67, tau_a = 0.44 (Gough point): index 1.087
```

The point where the Eq. (10) curve meets the $\tau$ axis has index 1, with or without a static
torsion. The Gough test point read from Fig. 5 lies 9 % outside the curve, as drawn.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_Santus --models KAKUNO --results-dir Results/current/guides
```

42CrMo4_Santus holds proportional tests only; its torsion tests with a mean shear stress need no
`Sig0`, because a torsional mean does not enter. On a workbook with out-of-phase tests, such as
42CrMo4_FABEST, the model is skipped as out of scope.

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`), a repository extension.
- **DP comparison:** supported. The index is not homogeneous in the load; use `Error_load`.
- **Fatigue-limit datasets:** index − 1 measures the distance from the limit along the index.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_kakuno_kawada.py`](../../tests/methods/test_kakuno_kawada.py) checks:

- that the three reference limits give index 1 (Eqs. 7–9);
- that points on the Eq. (10) curve give index 1, with and without a bending mean stress;
- that a torsional mean stress does not change the index;
- Eq. (12): $B = C$ forces $k_1 = 1$ for $\sigma_{wb}/\tau_w = \sqrt3$, as the paper states;
- Fig. 5: the $\tau$ intercept, 0.559 against 0.557 read from the figure. The $\sigma$ intercept,
  0.906 computed against about 0.93 drawn, is a strict xfail: no single set of the caption values
  matches both intercepts;
- the rejection of non-proportional loading and of a mean normal stress without `Sig0`.

Limitations:

- proportional loading only;
- validated in the source on figures only;
- a fatigue-limit criterion, so finite lives are available only through the master curve.

## References

- H. Kakuno, Y. Kawada (1979), A new criterion of fatigue strength of a round bar subjected to
  combined static and repeated bending and torsion, *Fatigue of Engineering Materials and
  Structures* 2, 229–236. <https://doi.org/10.1111/j.1460-2695.1979.tb01358.x>
- G. Sines (1955), Failure of materials under combined repeated stresses with superimposed static
  stresses, NACA Technical Note 3495.
