# COIN — Complex Invariants of the Modified Mohr–Mises Criterion (LiWI)

**Registry name:** `COIN`  
**Family:** Invariant (complex frequency domain), stress, master-curve life  
**Source file:** [`src/multifat/models/invariant_methods/coin_liwi.py`](../../src/multifat/models/invariant_methods/coin_liwi.py)

## Overview

`COIN` formulates the quadratic Modified Mohr–Mises invariant expression $I_1^2 - k_a^2 I_2$ directly in complex stress-amplitude space for single-frequency multiaxial histories, blending the Hermitian magnitude $H$ (invariant RMS) and complex square $|C|$ through a phase-sensitivity weight $k_\phi$.

## Original formulation

In Schmidt & Kraft (2023), Eqs. (21)–(32) (pp. 5–6), for a single-frequency stress state $\boldsymbol{\sigma}(t) = \operatorname{Re}(\mathbf{Z}\,e^{j\omega t})$ with complex amplitude tensor $\mathbf{Z}$ and strength ratio $k_a = s_{-1}/t_{-1}$ (Eq. 22), the quadratic Mertens form $\sigma_{\mathrm{eq},ka}^2 = I_1^2 - k_a^2 I_2$ (Eq. 21) defines a Hermitian magnitude $H$ and complex square $C$:

```math
H = |\operatorname{tr}\mathbf{Z}|^2 - \frac{k_a^2}{2}\!\left(|\operatorname{tr}\mathbf{Z}|^2 - \operatorname{tr}(\mathbf{Z}\,\overline{\mathbf{Z}})\right),
\qquad
C = (\operatorname{tr}\mathbf{Z})^2 - \frac{k_a^2}{2}\!\left((\operatorname{tr}\mathbf{Z})^2 - \operatorname{tr}(\mathbf{Z}^2)\right),
```

combined into the equivalent amplitude

```math
\sigma_{\mathrm{COIN}} = \sqrt{k_\phi\,H + (1 - k_\phi)\,|C|}.
```

Writing $Q(t) = I_1(t)^2 - k_a^2 I_2(t)$, these quantities satisfy $H = 2\,\overline{Q}$ and $|C| = 2(\max_t Q - \overline{Q})$: $k_\phi = 1$ recovers the invariant RMS, $k_\phi = 0.5$ equals the time-domain maximum $\sqrt{\max_t Q}$ (matching `MMMH` alternating amplitude), and $k_\phi = 0$ cancels $90^\circ$ phase-shifted invariants.

## Repository implementation

**Complex phasor construction and 3-D invariants.** Harmonic components `(S_sin, S_cos)` map to $\mathbf{Z} = \mathbf{S}_{\cos} - j\,\mathbf{S}_{\sin}$. Although derived for plane stress in Eqs. (24)–(30), the invariant traces carry over to general 3-D symmetric stress tensors without modification, and the output is normalized by $s_{-1}$ (`sigma_COIN / Sigm1`).

**Repository choices and departures.**
- *Default `k_phi = 1.0`:* Schmidt & Kraft (2023) fit $k_\phi$ empirically per material (Table 3 reports $0.85$–$1.10$) without giving an analytical rule. When optional `k_phi` is absent from the workbook, the repository defaults to `1.0` (RMS formulation).
- *No amplitude mean-stress term:* As in the paper (where mean stress is injected separately via PRAM/PRAJ damage laws), `COIN` includes no mean-stress contribution.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $s_{-1}$ | `Sigm1` | `Sigm1` | MPa | Required reversed axial/bending endurance limit ($> 0$) |
| $t_{-1}$ | `Taum1` | `Taum1` | MPa | Required reversed torsion endurance limit ($> 0$) |
| $k_\phi$ | `k_phi` | `k_phi` | — | Optional phase-sensitivity blend weight (default `1.0`) |

## Calibration

No statistical calibration is implemented for `k_phi` (`optional_params` documents the `1.0` default); `k_phi` is read from the workbook if present.

## Worked example

Compare the default RMS weight (`k_phi = 1.0`) against the time-domain-maximum weight (`k_phi = 0.5`) under $90^\circ$ out-of-phase bending–torsion ($\sigma_a = 200\,\text{MPa}$, $\tau_a = 120\,\text{MPa}$) on steel ($s_{-1} = 313.9\,\text{MPa}$, $t_{-1} = 196.2\,\text{MPa}$):

```python
import numpy as np
from multifat.models.invariant_methods.coin_liwi import CoinLiwi

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"Sigm1": 313.9, "Taum1": 196.2}
h = dict(S_mean=Z, S_sin=bt(200.0, 0.0), S_cos=bt(0.0, 120.0))
c1 = CoinLiwi().evaluate_case(S0=Z, E0=Z, S1=Z, E1=Z, params=p, harmonics=h).values[""]
c05 = CoinLiwi().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params={**p, "k_phi": 0.5}, harmonics=h
).values[""]
assert np.isclose(c1, 0.8832, atol=1e-3) and np.isclose(c05, 0.6371, atol=1e-3)
print(f"COIN (k_phi=1.0) = {c1:.4f}, COIN (k_phi=0.5) = {c05:.4f}")
```

Output:

```text
COIN (k_phi=1.0) = 0.8832, COIN (k_phi=0.5) = 0.6371
```

Batch CLI evaluation:

```bash
multifat-run --model COIN
```

## Comparison capabilities and outputs

`COIN` sets `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and returns the normalized index $\sigma_{\mathrm{COIN}} / s_{-1}$ under key `""`.

## Verification and limitations

- **Software verification.** [`tests/methods/test_coin_liwi.py`](../../tests/methods/test_coin_liwi.py) verifies closed-form plane-stress Eqs. (27) and (30) across arbitrary component phases, equivalence against 20,000-step time-domain $Q(t)$ at $k_\phi = 1.0$ (`sqrt(2 mean Q)`) and $k_\phi = 0.5$ (`sqrt(max Q)`), exact cancellation at $k_\phi = 0$ for $\sigma_x = k_a \tau_{xy}$ shifted by $90^\circ$, and uniaxial/torsional unit limits.
- **No printed specimen tables.** Schmidt & Kraft (2023) report summary scatter-band life statistics (Table 4) on external datasets without specimen tables, and the pure-index output carries zero mean-stress sensitivity unless handled upstream.

## References

- Schmidt, A. T., & Kraft, J. (2023). A new equivalent stress approach based on complex invariants: The COIN LiWI method. *International Journal of Fatigue*, 168, 107474. [doi:10.1016/j.ijfatigue.2022.107474](https://doi.org/10.1016/j.ijfatigue.2022.107474)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
