# Gonçalves–Araújo–Mamiya (`GAM`)

> **Registry name:** `GAM` · **Catalogue entry:** [N08](../../METHODS.md#n08) ·
> **Source:** Gonçalves, Araújo & Mamiya (2005), <https://doi.org/10.1016/j.ijfatigue.2004.05.006> ·
> **Implementation:** [`invariant_methods/goncalves_araujo_mamiya.py`](../../src/multifat/models/invariant_methods/goncalves_araujo_mamiya.py)

## Overview

- **Type and family:** stress-based endurance criterion, of the effective-stress (invariant-like)
  type.
- **Idea:** two measures of the stress history are combined.
  - Shear solicitation: the size of the rectangular prism that encloses the deviatoric stress path.
    It is larger for a non-proportional path than for a proportional one with the same $\sqrt{J_{2,a}}$.
  - Normal solicitation: the **maximum principal stress** over the cycle, in place of the
    hydrostatic stress. The worst-oriented embryo crack, normal to the largest principal stress,
    governs the opening effect.
- **Proposed for:** high-cycle fatigue limits of hard metals, $1.3 \le f_{-1}/t_{-1} < \sqrt3$,
  under in-phase and out-of-phase sinusoidal loading of the same frequency, with mean stresses.
- **Status:** fully implemented for the declared repository scope. Eqs. (11)–(13) and Appendix A
  were checked on the paper. Tables 1–4 are reproduced for the in-phase rows; the out-of-phase rows
  of the tables use a different peak stress (see [Verification](#verification-and-limitations)).

## Original formulation

Gonçalves et al. (2005, Eq. 12) state the endurance condition

```math
\sqrt{\sum_{i=1}^{5} a_i^2} + \kappa\,\sigma_{p\max} \le \lambda .
```

The two terms are defined as follows.

- **Deviatoric amplitudes $a_i$ (Eq. A.17).** The deviatoric stress
  $\mathbf S(t) = \boldsymbol\sigma(t) - \tfrac13 \mathrm{tr}\,\boldsymbol\sigma(t)\,\mathbf I$ is
  written in an orthonormal basis of the deviatoric space as $\mathbf S(t) = \sum_i s_i(t)\mathbf N_i$
  (Eqs. A.1–A.3). Then
  $a_i = \tfrac12\left(\max_t s_i(t) - \min_t s_i(t)\right)$, the half-sides of the enclosing
  rectangular prism.
- **Peak principal stress (Eq. 11).**
  $\sigma_{p\max} = \max_t \max_{i} \sigma_{p,i}(t)$, the largest principal stress over the cycle.

Proposition 1 (Eq. A.5) shows that, for a path enclosed by an ellipsoid with semi-axes $\lambda_i$,
$\sum_i a_i^2 = \sum_i \lambda_i^2$ in **any** orthonormal basis. This identity avoids the
construction of the minimum enclosing ellipsoid (Eqs. 8–9).

The constants follow from the fully reversed bending and torsion limits (Eq. 13):

```math
\kappa = \frac{\sqrt2}{f_{-1} - t_{-1}}\left(t_{-1} - \frac{f_{-1}}{\sqrt3}\right), \qquad
\lambda = \sqrt2\,\frac{t_{-1} f_{-1}}{f_{-1} - t_{-1}}\left(1 - \frac{1}{\sqrt3}\right).
```

The error index $I = 100\,(\mathrm{LHS} - \lambda)/\lambda$ (Eq. 14) is positive on the
conservative side.

## Repository implementation

For a single-harmonic history $\boldsymbol\sigma(t) = \boldsymbol\sigma_m + \boldsymbol\sigma_s\sin\omega t + \boldsymbol\sigma_c\cos\omega t$,
each component is $s_i(t) = m_i + p_i\sin\omega t + q_i\cos\omega t$. Hence $a_i^2 = p_i^2 + q_i^2$,
and summing over an orthonormal basis gives

```math
\sqrt{\sum_{i=1}^{5} a_i^2} = \sqrt{\|\mathrm{dev}\,\boldsymbol\sigma_s\|_F^2 + \|\mathrm{dev}\,\boldsymbol\sigma_c\|_F^2},
```

with $\|\cdot\|_F$ the Frobenius norm. This is exact and basis-independent. $\sigma_{p\max}$ is the
largest eigenvalue of $\boldsymbol\sigma(t)$, maximised over the cycle: 720 time samples followed by
a golden-section refinement, the routine shared with [MMMH](../../src/multifat/models/invariant_methods/mmmh.py).

The model returns the fatigue index

```math
\mathrm{DP} = \frac{\sqrt{\sum a_i^2} + \kappa\,\sigma_{p\max}}{\lambda},
```

which equals 1 on the fatigue limit, so $I = 100(\mathrm{DP} - 1)$.

**Repository choices**

- *Eq. (11) as written.* The peak principal stress is taken over the cycle. The paper's tables used
  a phase-free peak instead (see [Verification](#verification-and-limitations)); the repository
  follows the equation.
- *Range of validity.* $f_{-1} \le t_{-1}$ is rejected, because $\kappa$ and $\lambda$ are then
  undefined. Outside the hard-metal range $[1.3, \sqrt3)$ the criterion is still evaluated, with a
  warning. For $f_{-1}/t_{-1} \ge \sqrt3$, $\kappa \le 0$, so a tensile principal stress would no
  longer be detrimental.
- *Finite life.* The index is used as the damage parameter of the shared master curve
  $\mathrm{DP} = A N_f^b$. This extends an endurance criterion beyond its source and is not part
  of it ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $f_{-1}$ | `Sigm1` / Summary R = -1 tension or bending limit | MPa | workbook | yes, $> t_{-1}$ |
| $t_{-1}$ | `Taum1` / Summary R = -1 torsion limit | MPa | workbook | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** $\kappa$ and $\lambda$ follow in closed form from $f_{-1}$ and $t_{-1}$
   (Eq. 13). They place the fully reversed bending limit and the fully reversed torsion limit
   exactly on the threshold:
   - bending: $\sqrt{\sum a_i^2} = \sqrt{2/3}\,f_{-1}$ and $\sigma_{p\max} = f_{-1}$;
   - torsion: $\sqrt{\sum a_i^2} = \sqrt2\,t_{-1}$ and $\sigma_{p\max} = t_{-1}$.
2. **Repository.** The same closed form, with the R = -1 limits of the workbook Summary. Nothing is
   fitted except the master-curve constants $A$ and $b$, on the uniaxial R = -1 design tests.
3. A material whose limits are missing or flagged by the importer is skipped. The reason is written
   to `<MATERIAL>_run_report.tsv`.

## Worked example

Hard steel of the paper's Table 1 (Nishihara & Kawamoto): $f_{-1} = 313.9$ MPa,
$t_{-1} = 196.2$ MPa.

```python
import numpy as np

from multifat.models.invariant_methods.goncalves_araujo_mamiya import gam_constants
from multifat.models.registry import get_models

gam = get_models(["GAM"])[0]
f_m1, t_m1 = 313.9, 196.2  # MPa, hard steel of Nishihara & Kawamoto (paper's Table 1)
params = {"Sigm1": f_m1, "Taum1": t_m1}
weight, limit = gam_constants(f_m1, t_m1)  # kappa and lambda of Eq. (13)


def bending_torsion(sigma, tau):
    return np.array([[sigma, tau, 0.0], [tau, 0.0, 0.0], [0.0, 0.0, 0.0]])


def index(sigma_a, tau_a, beta_deg, sigma_m=0.0, tau_m=0.0):
    """GAM index of sigma = sigma_m + sigma_a sin wt, tau = tau_m + tau_a sin(wt - beta)."""
    b = np.radians(beta_deg)
    S_mean = bending_torsion(sigma_m, tau_m)
    S_sin = bending_torsion(sigma_a, tau_a * np.cos(b))
    S_cos = bending_torsion(0.0, -tau_a * np.sin(b))
    result = gam.evaluate_case(
        S0=S_mean + S_sin,
        S1=S_mean - S_sin,
        E0=None,
        E1=None,
        params=params,
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": S_cos},
    )
    return result.values[""]


# Eq. (13) puts both reference limits on the threshold
assert abs(index(f_m1, 0.0, 0) - 1.0) < 1e-12
assert abs(index(0.0, t_m1, 0) - 1.0) < 1e-12

# Table 1, rows 1-1 (in phase) and 1-4 (90 deg out of phase)
in_phase = index(138.1, 167.1, 0)
out_of_phase = index(150.2, 181.7, 90)

# Row 1-4 by hand: sqrt(sum a_i^2) = sqrt(2/3 sigma_a^2 + 2 tau_a^2); sigma_p,max is the peak of
# sigma_1(t) = sigma/2 + sqrt(sigma^2/4 + tau^2) over the cycle, found here on a fine time grid
wt = np.linspace(0.0, 2 * np.pi, 200001)
sigma, tau = 150.2 * np.sin(wt), -181.7 * np.cos(wt)
peak = np.max(sigma / 2 + np.sqrt(sigma**2 / 4 + tau**2))
amp = np.sqrt(2 / 3 * 150.2**2 + 2 * 181.7**2)
assert abs(out_of_phase - (amp + weight * peak) / limit) < 1e-6
# The paper's Table 1 combines the peaks as if in phase: sigma_p,max = 271.7 MPa
table_peak = 150.2 / 2 + np.hypot(150.2 / 2, 181.7)
table_index = (amp + weight * table_peak) / limit

print(f"kappa = {weight:.4f}, lambda = {limit:.2f} MPa")
print(f"row 1-1, in phase:     I = {100 * (in_phase - 1):+.2f} %")
print(f"row 1-4, Eq. (11):     I = {100 * (out_of_phase - 1):+.2f} % (sigma_p,max = {peak:.1f} MPa)")
print(f"row 1-4, paper's peak: I = {100 * (table_index - 1):+.2f} % (sigma_p,max = {table_peak:.1f} MPa)")
```

Output:

```text
kappa = 0.1799, lambda = 312.76 MPa
row 1-1, in phase:     I = -1.91 %
row 1-4, Eq. (11):     I = +2.51 % (sigma_p,max = 199.5 MPa)
row 1-4, paper's peak: I = +6.66 % (sigma_p,max = 271.7 MPa)
```

The two indices of row 1-4 differ only in $\sigma_{p\max}$:

- **Eq. (11).** At 90° the bending and torsion peaks do not coincide, and the largest principal
  stress over the cycle is 199.5 MPa.
- **The paper's Table 1.** The peaks are combined as if in phase, giving 271.7 MPa, and with it the
  printed $I = 6.66$ %.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models GAM --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve on the index (`Nf_expected`), a repository extension.
- **DP comparison:** supported, against the master-curve `DP_fit`. The index is homogeneous of
  degree 1 in the load.
- **Fatigue-limit datasets:** $\mathrm{DP} - 1$ is directly the paper's error index $I/100$.
- No `_ext` variant.

## Verification and limitations

[`tests/methods/test_goncalves_araujo_mamiya.py`](../../tests/methods/test_goncalves_araujo_mamiya.py)
checks Tables 1–4, column $I^d$ (41 fatigue limits, four steels):

- all 14 in-phase rows are reproduced to 0.02 points by Eq. (12) as written;
- all 41 rows are reproduced to 0.02 points if $\sigma_{p\max}$ is replaced by the phase-free peak
  $s/2 + \sqrt{s^2/4 + x^2}$, with $s = \sigma_a + \sigma_m$ and $x = \tau_a + \tau_m$;
- with Eq. (11), the 27 out-of-phase rows are 0.3 to 13.4 points lower, so less conservative. The
  largest difference is row 3-7 at $\beta = 180°$, where the bending and torsion peaks never
  coincide.

These rows are kept as strict expected failures, recording a calculation inconsistency of the
tables, not of the criterion. The paper's remark that the model is conservative in the presence of
mean stresses partly reflects this phase-free peak.

The same module also tests:

- that Eq. (13) puts both reference limits on the threshold;
- that the closed-form $\sqrt{\sum a_i^2}$ equals the half-ranges of Eq. (A.17) computed in the
  basis of Eq. (A.3) for a random out-of-phase path;
- the rejection of $f_{-1} \le t_{-1}$ and the hard-metal warning.

Limitations:

- an endurance criterion for hard metals;
- restricted by the source to iso-frequency sinusoidal loading, as the repository data are;
- finite-life use only through the master curve.

## References

- C. A. Gonçalves, J. A. Araújo, E. N. Mamiya (2005), Multiaxial fatigue: a stress based criterion
  for hard metals, *International Journal of Fatigue* 27(2), 177–187.
  <https://doi.org/10.1016/j.ijfatigue.2004.05.006>
- E. N. Mamiya, J. A. Araújo (2002), Fatigue limit under multiaxial loadings: on the definition of
  the equivalent shear stress, *Mechanics Research Communications* 29, 141–151.
