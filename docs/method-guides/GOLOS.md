# Golos–Dębski–Dębski modified Gough curve (`GOLOS`)

> **Registry name:** `GOLOS` · **Catalogue entry:** [N13](../../METHODS.md#n13) ·
> **Source:** Golos, Dębski & Dębski (2014), <https://doi.org/10.1016/j.tafmec.2014.07.005> ·
> **Implementation:** [`invariant_methods/golos.py`](../../src/multifat/models/invariant_methods/golos.py)

## Overview

- **Type and family:** stress-based fatigue-limit curve for combined bending and torsion.
- **Idea:** two classical limit curves are interpolated by one corrected ellipse, with
  coefficients linear in the limit ratio $\kappa = f_{-1}/t_{-1}$:
  - the Gough ellipse, which corresponds to the von Mises (HMH) equivalent stress, $\kappa = \sqrt3$;
  - the curve of Dębski's equivalent stress, $\kappa = 4/3$.
- **Proposed for:** high-cycle fatigue limits under **in-phase, fully reversed** bending and
  torsion. Six steels with $1.53 \le \kappa \le 1.75$ are used in the paper.
- **Status:** implemented for the source's loading scope. Cases outside it are not extended, so
  the model is skipped on them (see below). The paper's validation is graphical, so no table is
  reproduced.

## Original formulation

Golos et al. (2014) generalise Gough's ellipse (Eq. 29):

```math
A\left(\frac{\sigma}{f_{-1}}\right)^2 + B\,\frac{\sigma}{f_{-1}} + C\left(\frac{\tau}{t_{-1}}\right)^2 = 1,
```

where $\sigma$ and $\tau$ are the bending and torsion stress amplitudes. The coefficients are fixed
by the two reference curves (Eq. 30):

```math
A = \frac{2}{3\sqrt3 - 4}\,\kappa + \frac{\sqrt3 - 4}{3\sqrt3 - 4}, \qquad
B = -\frac{2}{3\sqrt3 - 4}\,\kappa + \frac{2\sqrt3}{3\sqrt3 - 4}, \qquad C = 1 .
```

Their engineering approximation gives the final criterion, Eqs. (31)–(32):

```math
\left(\frac{\sigma}{f_{-1}}\right)^2 (1.67\kappa - 1.90) + \frac{\sigma}{f_{-1}}(-1.67\kappa + 2.90)
+ \left(\frac{\tau}{t_{-1}}\right)^2 = 1 .
```

Properties of these coefficients:

- **Bending limit.** $A + B = 1$ in both forms, so the bending limit lies on the curve for every
  $\kappa$.
- **Reference curves.** $\kappa = \sqrt3$ gives $A = 1$, $B = 0$, the Gough ellipse (Eq. 26).
  $\kappa = 4/3$ gives $A = 1/3$, $B = 2/3$, i.e. the curve of
  $\sigma_{eq} = \sigma/3 + \tfrac23\sqrt{\sigma^2 + 4\tau^2}$ (Eqs. 27–28).

## Repository implementation

The model returns the left-hand side of Eq. (32), which equals 1 on the limit curve.

**Loading scope.** A fully reversed, in-phase plane stress whose principal amplitudes satisfy
$\sigma_1 \ge 0 \ge \sigma_2$ is exactly a bending–torsion state in a rotated frame, with

```math
\sigma = \sigma_1 + \sigma_2, \qquad \tau = \sqrt{-\sigma_1\sigma_2}.
```

The repository uses this frame-invariant reduction. Any other case raises a `LoadingScopeError`:
a mean stress, a phase shift, a triaxial state, or principal amplitudes of equal sign. The model
is then **skipped** for that material, with the reason in `<MATERIAL>_run_report.tsv`, and the
eligibility matrix records `out_of_scope`. It is never silently extended.

**Repository choices**

- *Coefficients.* Eq. (32), i.e. the approximated coefficients of Eq. (31), as the paper's final
  criterion. Eq. (30) is available as `golos_coefficients(..., exact=True)`. The two differ by less
  than 0.01 for $1.5 \le \kappa \le 1.75$.
- *Frame-invariant reduction.* The reduction above, without extending the source's loading scope.
- *Finite life.* The index is used as the damage parameter of the shared master curve
  $\mathrm{DP} = A N_f^b$. This is an extension of the limit curve. Because of the linear term
  $B\,\sigma/f_{-1}$, the index is not homogeneous in the load, so the load-equivalent error is the
  appropriate cross-method metric
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $f_{-1}$ | `Sigm1` / Summary R = -1 bending or tension limit | MPa | workbook | yes |
| $t_{-1}$ | `Taum1` / Summary R = -1 torsion limit | MPa | workbook | yes |
| $\boldsymbol\sigma(t)$ | `S_mean` (zero), `S_sin`, `S_cos` (in phase) | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $A$ and $B$ follow from $\kappa = f_{-1}/t_{-1}$; $C = 1$. Nothing is
   fitted.
2. **Repository.** The same, with the workbook limits. Only the master-curve constants are fitted.

**Finding.** The text also states $A = 1$, $B = 0$ for $\kappa = 2$. Eq. (30) is the straight line
through the $\kappa = 4/3$ and $\kappa = \sqrt3$ points and gives $A = 1.45$, $B = -0.45$ at
$\kappa = 2$, so the statement is not consistent with the equations.

For $\kappa > 1.737$, Eq. (31) gives $B < 0$. A small bending stress then lowers the index below the
torsion-only value; the source's six steels reach $\kappa = 1.75$.

## Worked example

Hard steel of Nishihara & Kawamoto ($f_{-1} = 313.9$ MPa, $t_{-1} = 196.2$ MPa), with the three
in-phase fatigue limits of that dataset as tabulated by Gonçalves et al. (2005).

```python
import numpy as np

from multifat.models.invariant_methods.golos import golos_coefficients
from multifat.models.registry import get_models

golos = get_models(["GOLOS"])[0]
f_m1, t_m1 = 313.9, 196.2  # MPa, hard steel of Nishihara & Kawamoto (paper's Table 1)
params = {"Sigm1": f_m1, "Taum1": t_m1}
zero = np.zeros((3, 3))
A, B = golos_coefficients(f_m1, t_m1)  # Eq. (31), kappa = 1.60


def index(sigma_a, tau_a):
    """GOLOS index of in-phase, fully reversed bending sigma_a and torsion tau_a."""
    S = np.array([[sigma_a, tau_a, 0.0], [tau_a, 0.0, 0.0], [0.0, 0.0, 0.0]])
    out = golos.evaluate_case(
        S0=S, S1=-S, E0=None, E1=None, params=params,
        harmonics={"S_mean": zero, "S_sin": S, "S_cos": zero},
    )  # fmt: skip
    return out.values[""]


print(f"A = {A:.3f}, B = {B:.3f}, C = 1")
for sigma_a, tau_a in ((138.1, 167.1), (245.3, 122.6), (299.1, 62.8)):
    gough = (sigma_a / f_m1) ** 2 + (tau_a / t_m1) ** 2  # Gough ellipse, Eq. (2)
    print(f"sigma_a = {sigma_a:5.1f}, tau_a = {tau_a:5.1f}: GOLOS {index(sigma_a, tau_a):.3f}, "
          f"Gough {gough:.3f}")  # fmt: skip
```

Output:

```text
A = 0.772, B = 0.228, C = 1
sigma_a = 138.1, tau_a = 167.1: GOLOS 0.975, Gough 0.919
sigma_a = 245.3, tau_a = 122.6: GOLOS 1.040, Gough 1.001
sigma_a = 299.1, tau_a =  62.8: GOLOS 1.021, Gough 1.010
```

All three fatigue limits give an index within 4 % of 1. Compared with the Gough ellipse, the
correction raises the index of the torsion-dominated test, where it is most visible.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models GOLOS --results-dir Results/current/guides
```

On a workbook with mean stresses or out-of-phase tests, which includes every supplied workbook, the
model is skipped as out of scope, with the reason in the run report.

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`), a repository extension.
- **DP comparison:** supported. The index is not homogeneous in the load; use `Error_load`.
- **Fatigue-limit datasets:** $\mathrm{DP} - 1$ measures the distance from the curve along the
  index, not the normal distance $r_i$ the paper uses in its scatter index (Fig. 8).
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_golos.py`](../../tests/methods/test_golos.py) checks:

- $A + B = 1$ for Eqs. (30) and (31);
- the special cases $\kappa = \sqrt3$ and $\kappa = 4/3$ (Eq. 28 coincides with Eq. 27);
- the accuracy of Eq. (31) on the paper's range;
- the inconsistency at $\kappa = 2$;
- the frame-invariant reduction of a rotated bending–torsion state;
- the rejection of every out-of-scope case.

[`tests/test_pipeline.py`](../../tests/test_pipeline.py) checks that an out-of-phase workbook
skips the model, with `mode = "out_of_scope"`, and that an in-phase workbook runs it.

Not reproduced: the scatter index of Table 2. The experimental points of Figs. 2–7 are not
tabulated, and $r_i$ is the normal distance to the curve in the normalised plane.

Limitations:

- in-phase, fully reversed bending–torsion only;
- no mean-stress or phase effect;
- a fatigue-limit curve, so finite lives are available only through the master curve.

## References

- K. M. Golos, D. K. Dębski, M. A. Dębski (2014), A stress-based fatigue criterion to assess
  high-cycle fatigue under in-phase multiaxial loading conditions, *Theoretical and Applied
  Fracture Mechanics* 73, 3–8. <https://doi.org/10.1016/j.tafmec.2014.07.005>
- H. J. Gough, H. V. Pollard (1935), The strength of metals under combined alternating stresses,
  *Proceedings of the Institution of Mechanical Engineers* 131.
