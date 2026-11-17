# OTT — Cycle-Based Reduced Ottosen Criterion

**Registry name:** `OTT`  
**Family:** Invariant, master-curve life  
**Source file:** [`src/multifat/models/invariant_methods/ottosen.py`](../../src/multifat/models/invariant_methods/ottosen.py)

## Overview

`OTT` implements the cycle-based reduced form of the continuum Ottosen endurance surface formulated by Tveit, Reyes & Erduran (2024). It pairs a non-quadratic Hershey–Hosford deviatoric norm (controlled by $t_{-1}/b_{-1}$) with the mean first stress invariant $I_{1,m}$ calibrated against zero-to-maximum ($R = 0$) bending strength.

## Original formulation

In Tveit, Reyes & Erduran (2024), Eqs. (11)–(17) (p. 6), fatigue endurance requires

```math
\sigma_{e,a} + A\,I_{1,m} \le S_0,
```

where the mean hydrostatic stress is $I_{1,m} = \frac{1}{2}(\min_P I_1 + \max_P I_1)$ (Eq. 12), the proportional Hershey–Hosford effective stress amplitude of exponent $m$ is (Eq. 13)

```math
\sigma_{e,a} = \left(\frac{|s_{1,a} - s_{2,a}|^m + |s_{2,a} - s_{3,a}|^m + |s_{3,a} - s_{1,a}|^m}{2}\right)^{\!1/m},
```

and the three material parameters satisfy Eqs. (14)–(17):

```math
S_0 = b_{-1}, \qquad A = \frac{2 b_{-1}}{b_0} - 1, \qquad
\kappa = \frac{t_{-1}}{b_{-1}} = \left(1 + 2^{m-1}\right)^{-1/m}.
```

## Repository implementation

The entire left-hand side $\sigma_{e,a} + A\,I_{1,m}$ is reported in MPa as the damage parameter correlated with fatigue life. For proportional cycles `ott_parameter` recovers the closed-form Eq. (13) value identically.

**Repository choices and departures.**
- *Mean-stress strength `b_0` fallback:* $b_0$ is read via explicit `A_OTT` when provided; otherwise from `Sig0` (the $R = 0$ strength), or from the Böhme et al. (2026) Eq. (1) estimate built from `Su` and `Sigm1`.
- *Exponent `m` inversion and cap:* When `m_OTT` is omitted, $m$ is solved by bisecting $\kappa = \tau_{-1}/\sigma_{-1}$ over $[1.0, 2.7670]$; values with $\kappa \ge 0.58524$ are capped at `M_CAP = 2.7670` (the maximum of Eq. 17) and values with $\kappa \le 0.5$ clamp to $m = 1.0$.
- *Non-proportional bounding ball:* Tveit et al. (2024) request the smallest hypersphere enclosing the deviatoric path without specifying an algorithm. Because the Hershey–Hosford effective stress is a centrally symmetric norm and single-harmonic paths are symmetric about `S_mean`, `ott_parameter` places the ball centre at `S_mean` and maximizes $\sigma_e(S(t) - S_{\text{mean}})$ over 720 phase angles across one period.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $b_{-1}$ | `Sigm1` | `Sigm1` | MPa | Required fully reversed bending/axial endurance limit |
| $t_{-1}$ | `Taum1` | `Taum1` | MPa | Required fully reversed torsional endurance limit |
| $m$ | `m_OTT` | `m_OTT` | — | Optional Hershey–Hosford exponent (capped at `2.767`) |
| $A$ | `A_OTT` | `A_OTT` | — | Optional mean-stress sensitivity factor |
| $b_0$ | `Sig0` | `Sig0` | MPa | Optional $R = 0$ normal endurance limit |
| $\sigma_u$ | `Su` | `Su` | MPa | Optional UTS used only for Böhme $b_0$ estimate |

## Calibration

No additional optional parameters are fitted statistically; `m_OTT` and `A_OTT` are derived analytically from standard workbook constants.

## Worked example

Compare in-phase and $90^\circ$ out-of-phase combined bending–torsion with mean axial tension on steel ($b_{-1} = 433\,\text{MPa}$, $t_{-1} = 250\,\text{MPa}$, $b_0 = 700\,\text{MPa}$):

```python
import numpy as np
from multifat.models.invariant_methods.ottosen import OttosenReduced, m_of_kappa

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"Sigm1": 433.0, "Taum1": 250.0, "Sig0": 700.0}
m_val = m_of_kappa(p["Taum1"] / p["Sigm1"])
v_ip = OttosenReduced().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params=p,
    harmonics=dict(S_mean=bt(100.0, 0.0), S_sin=bt(300.0, 150.0), S_cos=Z),
).values[""]
v_np = OttosenReduced().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params=p,
    harmonics=dict(S_mean=bt(100.0, 0.0), S_sin=bt(300.0, 0.0), S_cos=bt(0.0, 150.0)),
).values[""]
assert np.isclose(m_val, 2.0, atol=1e-3) and v_ip > v_np
print(f"m = {m_val:.4f}, OTT (in-phase) = {v_ip:.2f} MPa, OTT (90-deg) = {v_np:.2f} MPa")
```

Output:

```text
m = 2.0007, OTT (in-phase) = 420.57 MPa, OTT (90-deg) = 323.71 MPa
```

Batch CLI evaluation:

```bash
multifat-run --model OTT
```

## Comparison capabilities and outputs

`OTT` sets `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and returns a single unlabelled metric `""` in MPa.

## Verification and limitations

- **Software verification (paper reproduction).** In [`tests/methods/test_ottosen.py`](../../tests/methods/test_ottosen.py), every in-phase row of Tveit et al. (2024) Table A.1 matches the printed error index `Err` for both modified ($m$) and original ($m = 2$) variants within rounding tolerance. Out-of-phase Table A.1 rows agree within $0.8$ percentage points against the paper's incremental Runge–Kutta steady state (except 30NCD16 at $\phi = 90^\circ$, where Table A.1 prints `0.00` while the reduced closed form returns $-26\%$ to $-34\%$, marked strict `xfail`).
- **Exponent range restriction.** Because Eq. (17) reaches its global maximum $\kappa_{\max} \approx 0.58524$ at $m = 2.7670$, steels with $t_{-1}/b_{-1} > 0.585$ clamp to `M_CAP = 2.7670`.

## References

- Tveit, Ø., Reyes, A., & Erduran, E. (2024). A continuum approach to multiaxial high-cycle fatigue modeling for ductile metallic materials. *Results in Engineering*, 22, 102171. [doi:10.1016/j.rineng.2024.102171](https://doi.org/10.1016/j.rineng.2024.102171)
- Ottosen, N. S., Stenström, R., & Ristinmaa, M. (2008). Continuum approach to high cycle fatigue modeling. *International Journal of Fatigue*, 30(6), 996–1006.
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
