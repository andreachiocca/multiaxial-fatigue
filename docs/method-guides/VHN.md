# VHN — Vu, Halm & Nadot (2010) Stress Invariant Criterion

**Registry name:** `VHN`  
**Family:** Invariant, stress, master-curve life  
**Source file:** [`src/multifat/models/invariant_methods/vu_halm_nadot.py`](../../src/multifat/models/invariant_methods/vu_halm_nadot.py)

## Overview

`VHN` is a stress-invariant endurance criterion designed to capture out-of-phase loading sensitivity without spherical integration. It blends the peak alternating octahedral shear stress $\max J'_2(t)$ with the cycle-mean octahedral shear stress $J_{2,\text{mean}}$ and a branched hydrostatic term $I_f(I_{1,a}, I_{1,m})$ selected according to ultimate tensile strength.

## Original formulation

In Vu, Halm & Nadot (2010), Eqs. (6)–(21), the fatigue criterion requires $f_{\max} = \max_t f(t) \le \beta$ with $\beta = t_{-1}$ (Eq. 12) and

```math
f(t) = \sqrt{\gamma_1\,J'_2(t)^2 + \gamma_2\,J_{2,\mathrm{mean}}^2 + \gamma_3\,I_f(I_{1,a}, I_{1,m})},
\qquad
\gamma_3 = \frac{t_{-1}^2 - f_{-1}^2/3}{f_{-1}},
```

where $J'_2(t) = \sqrt{\tfrac{1}{2}\,\mathbf{s}_a(t) : \mathbf{s}_a(t)}$ is evaluated on the alternating deviatoric stress $\boldsymbol{\sigma}(t) - \boldsymbol{\sigma}_m$, $J_{2,\mathrm{mean}}$ is its time average over one period (Eq. 7), and $I_{1,a}$, $I_{1,m}$ are the amplitude and mean of $I_1$. The coefficients satisfy $\gamma_1 + 4\gamma_2/\pi^2 = 1$ (Eq. 13) with two UTS regimes (Eqs. 17–20):
- $R_m < 750\,\text{MPa}$: $\gamma_1 = 0.65$, $\gamma_2 = 0.8636$, $I_f = I_{1,a} + I_{1,m}$;
- $R_m > 750\,\text{MPa}$: $\gamma_1 = 0.3$, $\gamma_2 = 1.7272$, $I_f = I_{1,a} + \frac{f_{-1}}{t_{-1}}\,I_{1,m}$.

## Repository implementation

**Analytical ellipse quadrature.** For a single harmonic, the alternating deviatoric stress traces an ellipse in 5-D deviatoric space with $J_2$ semi-axes $p \ge q$. The repository computes $\max_t J'_2(t) = p$ and the exact period average $J_{2,\mathrm{mean}} = \frac{2}{\pi}\,p\,E(1 - q^2/p^2)$ via an arithmetic–geometric-mean complete elliptic integral (`_elliptic_second_kind`), eliminating time discretization error.

**Repository choices and departures.**
- *Boundary at `Su = 750 MPa`:* Values with `Su < 750.0` take the low-strength weights and `Su >= 750.0` take the high-strength branch.
- *Clamping negative radicands:* Strong compressive mean hydrostatic stress that turns the radicand negative clamps it to $0$ before taking the square root.
- *Single-harmonic scope:* Tables 13–14 of Vu et al. (2010) involve two-frequency non-periodic signals outside the single-harmonic workbook contract and are excluded.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $f_{-1}$ | `Sigm1` | `Sigm1` | MPa | Required reversed axial/bending endurance limit ($> 0$) |
| $t_{-1}$ | `Taum1` | `Taum1` | MPa | Required reversed torsion endurance limit ($> 0$) |
| $R_m$ | `Su` | `Su` | MPa | Required ultimate tensile strength selecting the $\gamma_1,\gamma_2$ branch |

## Calibration

No additional optional parameters are defined (`optional_params = None`).

## Worked example

Evaluate an out-of-phase cycle ($\sigma_{xx,m} = 100\,\text{MPa}$, $\sigma_{xx,a} = 260\,\text{MPa}$, $\tau_{xy,a} = 150\,\text{MPa}$ at $90^\circ$ shift) on high-strength steel ($f_{-1} = 398\,\text{MPa}$, $t_{-1} = 260\,\text{MPa}$, $R_m = 1025\,\text{MPa}$):

```python
import numpy as np
from multifat.models.invariant_methods.vu_halm_nadot import VuHalmNadot

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"Sigm1": 398.0, "Taum1": 260.0, "Su": 1025.0}
h = dict(S_mean=bt(100.0, 0.0), S_sin=bt(260.0, 0.0), S_cos=bt(0.0, 150.0))
val = VuHalmNadot().evaluate_case(S0=Z, E0=Z, S1=Z, E1=Z, params=p, harmonics=h).values[""]
assert np.isclose(val, 0.9500, atol=1e-3)
print(f"VHN = {val:.4f}")
```

Output:

```text
VHN = 0.9500
```

Batch CLI evaluation:

```bash
multifat-run --model VHN
```

## Comparison capabilities and outputs

`VHN` sets `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and outputs the dimensionless fatigue-limit index $f_{\max}/t_{-1}$ under variant `""`. The paper's percent error index (Eq. 21) is $100\,(\text{index} - 1)$.

## Verification and limitations

- **Software verification (paper reproduction).** In [`tests/methods/test_vu_halm_nadot.py`](../../tests/methods/test_vu_halm_nadot.py), 117 of 119 rows across Tables 1–12 of Vu et al. (2010) match printed $I^e$ within `0.06` percentage points. Two confirmed printing errors are strict `xfail`s (Table 2 row 9 prints $\delta_{xy} = 0$ where all printed indices match $90^\circ$; Table 11 row 1 prints $\Sigma_{xx,a} = 52.90$ where $I^a, I^e$ match $552.90$). Calibration identities and 40,000-step time-domain $J_{2,\mathrm{mean}}$ quadrature agree to `1e-10`.
- **Frequency limitation.** Requires `Su` in the workbook and applies only to synchronous single-harmonic cycles.

## References

- Vu, Q. H., Halm, D., & Nadot, Y. (2010). Multiaxial fatigue criterion for complex loading based on stress invariants. *International Journal of Fatigue*, 32(7), 1004–1014. [doi:10.1016/j.ijfatigue.2009.11.006](https://doi.org/10.1016/j.ijfatigue.2009.11.006)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
