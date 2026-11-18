# MMMH — Modified Mohr–Mises Hypothesis

**Registry name:** `MMMH`  
**Family:** Invariant (signed Mohr–Mises), stress, master-curve life  
**Source file:** [`src/multifat/models/invariant_methods/mmmh.py`](../../src/multifat/models/invariant_methods/mmmh.py)

## Overview

`MMMH` evaluates the signed Modified Mohr–Mises equivalent stress history of Mertens & Liebich (2017) for smooth components, pairing the maximum of the quadratic Mertens form $\sqrt{I_1^2 - k_a^2 I_2}$ with an equivalent mean normal stress $\sigma_{\mathrm{eq},m}$ defined via principal stresses on an endurance-limit Haigh line.

## Original formulation

In Mertens & Liebich (2017), Eqs. (1), (4), (9), and (22)–(25) for smooth surfaces ($m_a = n_{\mathrm{el}} = m_{o1} = m_{a1} = 1$), with $k_a = s_{-1}/t_{-1}$, Mohr circle centre $M_M = \frac{1}{2}(\sigma_1 + \sigma_2) - \sigma_3$, and signed radius $|R_M| = \frac{1}{2}|\sigma_1 - \sigma_2|$, the alternating equivalent history is

```math
\sigma_{\mathrm{eq,MMM}}(\eta) = V(\eta)\,\sqrt{(4 - k_a^2)\,M_M^2 + k_a^2\,R_M^2},
\qquad
V(\eta) = \operatorname{sgn}\!\left(\sqrt{4 - k_a^2}\,M_M + k_a\,R_M\right),
```

with alternating amplitude $\sigma_{\mathrm{eq},a} = \frac{1}{2}(\max_\eta \sigma_{\mathrm{eq,MMM}} - \min_\eta \sigma_{\mathrm{eq,MMM}})$ and equivalent mean stress $\sigma_{\mathrm{eq},m} = \sigma_{o1} - \sigma_{a1}$ (the cycle peak of maximum principal total stress minus the half-range of signed alternating $\sigma_1$).

## Repository implementation

**Single-harmonic identities.** Because a single-harmonic alternating tensor flips sign after half a period, both $M_M$ and signed $R_M$ reverse sign identically; hence $\sigma_{\mathrm{eq},a} = \max_t \sqrt{I_1^2 - k_a^2 I_2}$ coincides identically with `COIN` at $k_\phi = 0.5$ (`mmmh_amplitude`), and $\sigma_{a1}$ equals the peak alternating principal stress over one period (`mmmh_mean` via `peak_principal_stress`).

**Repository choices and departures.**
- *Endurance-limit Haigh line instead of unobtainable FKM $M_N$:* The paper's cycle-dependent Haigh diagram references ($M_N$ from refs. [16, 23]) cannot be obtained, while the paper explicitly allows guideline endurance values. The registered model returns $(\sigma_{\mathrm{eq},a} + M\,\sigma_{\mathrm{eq},m})/s_{-1}$, reading optional `M_MMMH` or setting $M = 2 s_{-1}/s_0 - 1$ (Dubbel 2004, E1) from `Sig0` or the Böhme et al. (2026) Eq. (1) estimate built from `Su` and `Sigm1`.
- *Domain constraint and text misprint:* Requires $0 < s_{-1}/t_{-1} \le 2$ (real $\sqrt{4 - k_a^2}$). For $k_a^2 = 0$, Eq. (8) gives $|s_{1a} + s_{2a}|$ rather than the paper's text misprint $|s_{1a} - s_{2a}|$.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $s_{-1}$ | `Sigm1` | `Sigm1` | MPa | Required reversed axial/bending endurance limit ($> 0$) |
| $t_{-1}$ | `Taum1` | `Taum1` | MPa | Required reversed torsion endurance limit ($0 < s_{-1}/t_{-1} \le 2$) |
| $M$ | `M_MMMH` | `M` | — | Optional Haigh slope (`2 Sigm1/Sig0 - 1` when omitted) |
| $s_0$ | `Sig0` | `Sig0` | MPa | Optional $R = 0$ endurance limit (Böhme estimate from `Su` when absent) |
| $\sigma_u$ | `Su` | `Su` | MPa | Optional UTS used only to estimate `Sig0` |

## Calibration

No statistical calibration is performed; `M_MMMH` defaults analytically from `Sig0` or `Su` via `Sigm1`.

## Worked example

Evaluate $90^\circ$ out-of-phase bending–torsion with axial mean stress ($\sigma_{xx,m} = 80\,\text{MPa}$, $\sigma_{xx,a} = 200\,\text{MPa}$, $\tau_{xy,a} = 120\,\text{MPa}$) on steel ($s_{-1} = 313.9\,\text{MPa}$, $t_{-1} = 196.2\,\text{MPa}$, $s_0 = 500\,\text{MPa}$):

```python
import numpy as np
from multifat.models.invariant_methods.mmmh import MMMH

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"Sigm1": 313.9, "Taum1": 196.2, "Sig0": 500.0}
h = dict(S_mean=bt(80.0, 0.0), S_sin=bt(200.0, 0.0), S_cos=bt(0.0, 120.0))
val = MMMH().evaluate_case(S0=Z, E0=Z, S1=Z, E1=Z, params=p, harmonics=h).values[""]
assert np.isclose(val, 0.7023, atol=1e-3)
print(f"MMMH = {val:.4f}")
```

Output:

```text
MMMH = 0.7023
```

Batch CLI evaluation:

```bash
multifat-run --model MMMH
```

## Comparison capabilities and outputs

`MMMH` sets `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and outputs the dimensionless endurance-limit Haigh index `""`.

## Verification and limitations

- **Software verification.** [`tests/methods/test_mmmh.py`](../../tests/methods/test_mmmh.py) verifies single-harmonic reduction against 20,000-step time-domain `mmmh_equivalent_history` and `COIN` ($k_\phi = 0.5$) to `1e-6`, proportional principal-stress limits (Eq. 8), literal Eqs. (20)–(23) non-proportional mean-stress checks, and the three Haigh calibration points ($s_{-1}$, $t_{-1}$, $s_0$).
- **Smooth-surface and Haigh-line scope.** Stress gradient and notch support factors ($m_a, n_{\mathrm{el}}, m_{o1}, m_{a1}$) are fixed at $1$, and Mertens & Liebich (2017) publish no specimen tables for direct reproduction.

## References

- Mertens, H., & Liebich, R. (2017). A time dependent equivalent stress function for proportional and non-proportional transient loaded and notched metallic components. *International Journal of Fatigue*, 103, 294–308. [doi:10.1016/j.ijfatigue.2017.06.003](https://doi.org/10.1016/j.ijfatigue.2017.06.003)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
