# RFB — Robert, Fogue & Bahuaud (1994) Integral Criterion

**Registry name:** `RFB`  
**Family:** Integral, stress, master-curve life  
**Source file:** [`src/multifat/models/integral_methods/robert_fogue_bahuaud.py`](../../src/multifat/models/integral_methods/robert_fogue_bahuaud.py)

## Overview

`RFB` integrates a plane-wise linear combination of maximum projected shear-stress amplitude, normal stress amplitude, and normal mean stress over the unit sphere in root-mean-square form. Three weighting parameters $(a, b, d)$ calibrate the integral to uniaxial reversed, torsional reversed, and pulsating tension limits.

## Original formulation

In Robert, Fogue & Bahuaud (1994), Eqs. (8), (15)–(25), and the Appendix (reduced with $c = 0$, Eq. 25), the plane damage measure for normal $\mathbf{h}$ is

```math
E_h = a \max_{\psi} \tau_{h\psi,a} + b\,\sigma_{hh,a} + d\,\sigma_{hh,m},
```

and the multiaxial failure condition averages $E_h^2$ uniformly over the sphere $S$ ($|S| = 4\pi$):

```math
E = \frac{1}{s_{-1}} \sqrt{\frac{1}{4\pi}\int_S E_h^2\,dS} = 1.
```

Setting $k = s_{-1}^2 / t_{-1}^2$, Eqs. (20), (21), and (24) fix $(a, b, d)$ for $0.577 < t_{-1}/s_{-1} < 0.866$:

```math
D = 225\left(1 - \frac{8(k - 3)^2}{25}\right), \qquad
b^2 = \frac{15 - \sqrt{D}}{2}, \qquad
a^2 = \frac{12 k - 21 + b^2}{2},
```

```math
d = \frac{-(3b + 2a) + \sqrt{(3b + 2a)^2 + 3\,(60 s_{-1}^2 / s_0^2 - 15)}}{3}.
```

## Repository implementation

**Spherical quadrature and shear projection.** The unit hemisphere is integrated on a fixed $48 \times 192$ Gauss–Legendre grid (`GRID = (48, 192)`). For single-harmonic stress cycles, $\max_\psi \tau_{h\psi,a}$ equals the minimum circumscribed circle radius `tau_a` on plane $\mathbf{h}$.

**Repository choices and departures.**
- *Published misprint corrections:* Source Eq. (16) prints $1/\pi$ instead of $1/(4\pi)$ (violating Eqs. 8 and 18–24), the Appendix omits the factor $1/2$ under the square root for $a$ (violating Eqs. 20 and 24), and Eq. (18) prints $E$ where $E^2$ is intended. The repository enforces the mathematically consistent $1/(4\pi)$ normalization and the $1/2$ factor in $a$.
- *Strength inputs:* Fixed-life strengths at target $N$ in the paper are mapped to workbook endurance limits (`Sigm1`, `Taum1`, and `Sig0` or the Böhme et al. 2026 Eq. 1 estimate of $s_0$ from `Su`), returning the scalar validity index $E$ (`I` in Eq. 27).
- *Ratio domain:* Ratios outside $0.577 < t_{-1}/s_{-1} < 0.866$ raise `MaterialScopeError`: the model is skipped on that material (eligibility `out_of_scope`), not reported as failed (e.g. AISI316L, $t_{-1}/s_{-1} = 0.517$).

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $s_{-1}$ | `Sigm1` | `Sigm1` | MPa | Required reversed tension/bending endurance limit |
| $t_{-1}$ | `Taum1` | `Taum1` | MPa | Required reversed torsion endurance limit ($0.577 < t_{-1}/s_{-1} < 0.866$) |
| $s_0$ | `Sig0` | `Sig0` | MPa | Optional $R = 0$ endurance limit; estimated via Böhme et al. (2026) from `Su` when absent |
| $\sigma_u$ | `Su` | `Su` | MPa | Optional UTS used only to estimate `Sig0` |

## Calibration

No model-specific parameters are fitted statistically; $(a, b, d)$ follow analytically from the three fatigue limits.

## Worked example

Evaluate combined reversed bending and torsion ($\sigma_a = 250\,\text{MPa}$, $\tau_a = 120\,\text{MPa}$) on the steel parameters of Table 2 ($s_{-1} = 361\,\text{MPa}$, $t_{-1} = 228\,\text{MPa}$, $s_0 = 600\,\text{MPa}$):

```python
import numpy as np
from multifat.models.integral_methods.robert_fogue_bahuaud import RobertFogueBahuaud

S = np.array([[250.0, 120.0, 0.0], [120.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"Sigm1": 361.0, "Taum1": 228.0, "Sig0": 600.0}
val = RobertFogueBahuaud().evaluate_case(S0=S, S1=-S, E0=Z, E1=Z, params=p).values[""]
assert np.isclose(val, 0.8698, atol=1e-3)
print(f"RFB = {val:.4f}")
```

Output:

```text
RFB = 0.8698
```

Batch CLI evaluation:

```bash
multifat-run --model RFB
```

## Comparison capabilities and outputs

`RFB` sets `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and outputs a single dimensionless index `""` calibrated to $1$ at the endurance limit.

## Verification and limitations

- **Software verification (paper reproduction).** In [`tests/methods/test_robert_fogue_bahuaud.py`](../../tests/methods/test_robert_fogue_bahuaud.py), Eqs. (19)–(24) hold to `4e-4`, and 52 of 59 rows in Tables 1–3 of Robert et al. (1994) reproduce $I$ within rounding (`0.006`). Six biaxial tension–pressure rows in Table 1 (rows 1, 2, 5, 6, 7, 12) differ by `0.006`–`0.015` and Table 3 row 11 ($\delta = 45^\circ$, printed `1.04` vs `1.073`, breaking its own phase trend) remain strict `xfail`s.
- **Strength ratio limitation.** Milder or hyper-ductile metals with $t_{-1}/s_{-1} \notin (0.577,\,0.866)$ fall outside the real-coefficient domain of Eq. (20).

## References

- Robert, J.-L., Fogue, M., & Bahuaud, J. (1994). Fatigue life prediction under periodical or random multiaxial stress states. In *Automation in Fatigue and Fracture: Testing and Analysis*, ASTM STP 1231 (pp. 369–387). [doi:10.1520/STP13959S](https://doi.org/10.1520/STP13959S)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
