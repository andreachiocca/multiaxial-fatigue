# Stress scale factor (`SSF`)

> **Registry name:** `SSF` · **Catalogue entry:** [N24](../../METHODS.md#n24) ·
> **Source:** Anes, Reis, Li, Fonte & de Freitas (2014), <https://doi.org/10.1016/j.ijfatigue.2013.05.004> ·
> **Implementation:** [`invariant_methods/ssf.py`](../../src/multifat/models/invariant_methods/ssf.py)

## Overview

- **Type and family:** stress-based equivalent shear stress with a material-specific,
  load-dependent scale factor between axial and shear stress.
- **Idea:** criteria such as Crossland, Matake or Dang Van weight the axial stress against the
  shear stress with a constant derived from the fatigue limits. Anes et al. measure this "stress
  scale factor" (ssf) on proportional tension–torsion tests and find that it depends on the
  axial stress level and on the amplitude ratio. They fit an ssf surface. Along a non-proportional
  path, the instantaneous ssf turns each instant of the axial stress into an equivalent shear
  stress.
- **Mean stress:** not addressed (all tests at $R = -1$).
- **Proposed for:** finite-life tension–torsion fatigue of high-strength 42CrMo4.
- **Status:** implemented with the published surface only (maintainer decision). Eqs. (29)–(31)
  and Tables 4–10 were checked on the paper.

## Original formulation

Anes et al. (2014, Eqs. 29–31):

```math
\tau_{eqv} = \max_t\left[\tau(t) + \mathrm{ssf}(\lambda,\sigma)\,\sigma(t)\right],
\qquad \lambda = \tan^{-1}\frac{\tau}{\sigma},
```

```math
\mathrm{ssf}(\sigma,\lambda) = a + b\sigma + c\sigma^2 + d\sigma^3 + f\lambda^2 + g\lambda^3 + h\lambda^4 + i\lambda^5,
```

with $a = 2.69$, $b = -9.90\times10^{-3}$, $c = 1.69\times10^{-5}$, $d = -9.52\times10^{-9}$,
$f = -5.99$, $g = 11.72$, $h = -8.04$, $i = 1.63$ ($\sigma$ in MPa, $\lambda$ in rad; there is no
$e$ term).

**Calibration of the surface (Tables 4–9).**
1. For each proportional case, axial ($AD$) and shear ($AB$) S–N trend lines are fitted, together
   with the pure-shear reference curve ($AC$, case 2: $\tau_a = 864.78\,N_f^{-0.061}$).
2. At each life, $\mathrm{ssf} = (AC - AB)/AD$, the shear increment needed to reach the
   pure-shear curve per unit axial stress.
3. Eq. (31) is regressed on these points ($R^2 = 0.98$).

**Life.** $\tau_{eqv}$ is compared with the pure-shear S–N curve.

**Non-proportional paths (Section 4.5).** The instantaneous stresses enter Eq. (31): $\lambda(t)$
and $\mathrm{ssf}(t)$ vary along the cycle (Figs. 12–14).

## Repository implementation

DP $= \tau_{eqv}$ (MPa), with the cycle sampled at 4096 points:

```math
\tau_{eqv} = \max_t\left[\,|\tau(t)| + \mathrm{ssf}\!\left(|\sigma(t)|, \tan^{-1}\tfrac{|\tau(t)|}{|\sigma(t)|}\right)|\sigma(t)|\,\right].
```

**Repository choices**

- *Published surface only (maintainer decision).*
  - The eight coefficients are required constants (`a_SSF` … `i_SSF`), so without them the model
    is skipped.
  - The paper's values are in
    [`Calibrations/42CrMo4_Reis/SSF_published.json`](../../Calibrations/42CrMo4_Reis/SSF_published.json).
    It is a published-constants record, read by `--constants`.
  - The per-material recalibration of the surface is not implemented.
- *Coefficients as printed (maintainer decision).* No refit. They miss the paper's Table 9 by up
  to 0.06.
- *Magnitudes.* $\lambda \in [0, \pi/2]$ and the sum uses $|\tau|$ and $|\sigma|$. The result
  therefore does not depend on axis sign conventions. It equals the signed form of Fig. 14
  wherever $\tau$ and $\sigma$ share the sign.
- *Scope.* Only fully reversed loading with one normal and one shear stress acting in the same
  plane (tension–torsion); anything else raises `LoadingScopeError`.
- *Range warning.* A warning is issued when $|\sigma|$ exceeds 737 MPa, the largest axial
  amplitude behind the fit.
- *Life.* The master curve is fitted on the pure torsion $R = -1$ failures
  (`design_loading = "torsion"`), the paper's reference curve. A power law is used, as in the
  paper.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $a, b, c, d, f, g, h, i$ | `a_SSF` … `i_SSF` | MPa$^{-k}$, – | `--constants` record or workbook Summary | yes |
| $\boldsymbol\sigma(t)$ | `S_mean` (zero), `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** The ssf surface needs proportional tension–torsion S–N data at several
   amplitude ratios plus pure tension and pure torsion (cases 1–5 of the paper).
2. **Repository.** The published surface for the paper's 42CrMo4.
   - The workbook `14_42CrMo4_Reis` holds the tests the surface was fitted on. Its SSF results are
     therefore **in-sample**: their errors measure the fit, not prediction (the paper's case 6 is
     the only test not used in the fit).
   - The master curve is fitted on the workbook's pure torsion failures.
3. **Cross-batch use (decision d103, issue #71).** The surface may also be applied to another
   batch of the same steel through a transfer with 42CrMo4_Reis as the design batch, which is the
   only context in which the record is used outside its own workbook:

       multifat-pipeline --material 42CrMo4_FABEST --design-material 42CrMo4_Reis --models SSF \
           --constants Calibrations/42CrMo4_Reis/SSF_published.json

   The master curve is then Reis's pure torsion curve, and the results are **out of sample**.
   - FABEST reaches axial amplitudes of 800 MPa, beyond the 737 MPa of the fit, so Eq. (31) is
     extrapolated there; the method warns.
   - The extrapolated cubic is not monotonic in the load, so the load-error search falls back to
     bisection (`eval/load_error.py`).
   - On 42CrMo4_Santus, SSF is skipped: its torsion tests with a mean shear are outside the R = −1
     scope.
   - Results: [cross-batch comparison](../batch-comparison.md).

## Worked example

Cases 4 and 6 of the paper with the published coefficients.

```python
import numpy as np

from multifat.models.invariant_methods.ssf import PUBLISHED, ssf_equivalent_shear, stress_scale_factor

zero = np.zeros((3, 3))


def harmonic(sigma_a, tau_a, delta_deg):
    """sigma_y = sigma_a sin wt, tau_xy = tau_a sin(wt + delta): (S_mean, S_sin, S_cos)."""
    d = np.radians(delta_deg)
    S_sin = np.array([[0, tau_a * np.cos(d), 0], [tau_a * np.cos(d), sigma_a, 0], [0, 0, 0]])
    S_cos = np.array([[0, tau_a * np.sin(d), 0], [tau_a * np.sin(d), 0, 0], [0, 0, 0]])
    return zero, S_sin, S_cos


# Proportional case 4 of the paper at 1e5 cycles (Table 7): sigma_a = 480, tau_a = 277 MPa
lam = np.arctan(277 / 480)
print(f"case 4: lambda = {lam:.2f} rad, ssf (Eq. 31) = {stress_scale_factor(480, lam, PUBLISHED):.3f} "
      f"(Table 7: 0.32), tau_eqv = {ssf_equivalent_shear(*harmonic(480, 277, 0), PUBLISHED):.1f} MPa "
      f"(pure-shear reference 428 MPa)")
# Non-proportional case 6 (Section 4.5): sigma_a = 510, tau_a = 294 MPa, 90 deg
print(f"case 6: tau_eqv = {ssf_equivalent_shear(*harmonic(510, 294, 90), PUBLISHED):.1f} MPa "
      f"(text p.31: 434 MPa)")
```

Output:

```text
case 4: lambda = 0.52 rad, ssf (Eq. 31) = 0.279 (Table 7: 0.32), tau_eqv = 411.0 MPa (pure-shear reference 428 MPa)
case 6: tau_eqv = 415.3 MPa (text p.31: 434 MPa)
```

The printed coefficients give 0.28 instead of Table 7's 0.32 for case 4. For case 6 they give
415 MPa instead of the 434 MPa of the text. Both gaps come from the rounding of the coefficients.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material 42CrMo4_Reis --models SSF --constants Calibrations/42CrMo4_Reis/SSF_published.json --results-dir Results/current/guides
```

The torsion S–N curve of this workbook is flat ($k \approx 16 > 15$, the paper's exponent is
$-0.061$), so only damage-parameter and load errors are reported.

## Comparison capabilities and outputs

- **Life:** master curve on $\tau_{eqv}$ fitted on the torsion failures.
- **DP comparison:** supported. $\tau_{eqv}$ is not homogeneous in the load (the ssf depends on
  $\sigma$), so use `Error_load` across methods.
- **Design set:** torsion instead of tension, as for `WANG2022`.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_ssf.py`](../../tests/methods/test_ssf.py):

- **Tables 5–8.** They are reproduced from the Table 4 trend lines, i.e. the calibration
  procedure. Two printed cells are misprints and are strict xfails:
  - Table 5 prints $BC = 0$ for case 1, where it should be $AC - AB = AC$;
  - Table 8 prints $AC = 508$ at $10^3$ cycles, where it should be 567.

  In both cases the table's own ssf column uses the correct value.
- **Table 9.** Eq. (31) as printed matches it within 0.06, with an RMS gap below 0.035.
- **Figs. 13–14 (case 6).**
  - The ssf history spans 0.37–0.79, as in Fig. 13.
  - The peak $\tau_{eqv}$ is 415 MPa, within the coefficient gap of the printed 434 MPa. Matching
    within 1 % is a strict xfail.
- **Other checks:**
  - proportional peaks;
  - pure torsion as the reference;
  - the scope errors;
  - the published-constants record.

Limitations:

- one material, with in-sample results on its workbook;
- fully reversed tension–torsion only;
- printed coefficients rounded to three digits;
- a polynomial outside its fitted range is unreliable (warning above 737 MPa).

## References

- V. Anes, L. Reis, B. Li, M. Fonte, M. de Freitas (2014), New approach for analysis of complex
  multiaxial loading paths, *International Journal of Fatigue* 62, 21–33.
  <https://doi.org/10.1016/j.ijfatigue.2013.05.004>
