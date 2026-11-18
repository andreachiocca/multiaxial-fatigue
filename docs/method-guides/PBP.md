# PBP — Projection-by-Projection Criterion

**Registry name:** `PBP`  
**Family:** Invariant (deviatoric projection), stress, master-curve life  
**Source file:** [`src/multifat/models/invariant_methods/pbp.py`](../../src/multifat/models/invariant_methods/pbp.py)

## Overview

`PBP` resolves the 5-dimensional deviatoric stress vector path onto its principal covariance directions, extracts each projection's cycle amplitude, and recombines them quadratically into an equivalent deviatoric amplitude $\Gamma$ compared against a linear function of the hydrostatic–deviatoric stress ratio $\rho_{\text{ref}}$.

## Original formulation

In Cristofori, Susmel & Tovo (2008), Eqs. (4)–(13) and Appendix B, the deviatoric stress tensor is mapped to a 5-vector $\mathbf{s}(t)$ with $\|\mathbf{s}\|^2 = J_2$ (Eq. A.3) and projected onto the eigenvectors of its covariance matrix (Eqs. 4–5). Combining the cycle-counted projection amplitudes yields (Eqs. 12–13)

```math
\Gamma = \sqrt{\sum_{i=1}^5 \left(\sqrt{J_{2a,i}}\right)^{\!2}},
```

with hydrostatic stress ratio (Eq. 6) and reference fatigue strength (Eq. 8)

```math
\rho_{\mathrm{ref}} = \frac{\sqrt{3}\,\sigma_{H,\mathrm{ref}}}{\Gamma}, \qquad
\sqrt{J_{2A,\rho}} = t_{-1} + \rho_{\mathrm{ref}}\!\left(\frac{f_{-1}}{\sqrt{3}} - t_{-1}\right),
```

calibrated to yield $\rho_{\mathrm{ref}} = 1$ at reversed axial loading ($\sqrt{J_{2A,1}} = f_{-1}/\sqrt{3}$) and $\rho_{\mathrm{ref}} = 0$ at reversed torsion ($\sqrt{J_{2A,0}} = t_{-1}$).

## Repository implementation

**Closed-form single-harmonic reduction.** For a single harmonic $\mathbf{s}(t) = \mathbf{a}\sin\omega t + \mathbf{b}\cos\omega t$, every orthonormal projection $i$ contains a single harmonic cycle of amplitude $\sqrt{(\mathbf{e}_i\cdot\mathbf{a})^2 + (\mathbf{e}_i\cdot\mathbf{b})^2}$. Summing the squared amplitudes over any complete orthonormal basis yields $\Gamma^2 = J_2(A) + J_2(B)$ identically without needing covariance eigendecomposition (`pbp_amplitude`), while the standalone covariance routine `projection_amplitudes` is also retained. Mean and alternating hydrostatic stress sum to $\sigma_{H,\mathrm{ref}} = \sigma_{H,a} + \sigma_{H,m}$.

**Repository choices and departures.**
- *Retained $\sqrt{3}$ in Eq. (6):* Eq. (B.5) in the paper's appendix drops the $\sqrt{3}$ factor of Eq. (6); the repository retains $\sqrt{3}$ so pure reversed axial tension produces $\rho_{\mathrm{ref}} = 1$.
- *Endurance-limit index only:* The repository exposes the normalized endurance-limit index $\Gamma / \sqrt{J_{2A,\rho}}$ (`0.0` if $\Gamma = 0$, `inf` if $\sqrt{J_{2A,\rho}} \le 0$). The paper's finite-life slope Eqs. (9)–(11) with $k_\rho$ is not implemented; life prediction uses the shared master curve.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $f_{-1}$ | `Sigm1` | `Sigm1` | MPa | Required reversed axial/bending endurance limit ($> 0$) |
| $t_{-1}$ | `Taum1` | `Taum1` | MPa | Required reversed torsion endurance limit ($> 0$) |

## Calibration

No additional optional parameters exist (`optional_params = None`).

## Worked example

Verify that in-phase and $90^\circ$ out-of-phase bending–torsion ($\sigma_a = 200\,\text{MPa}$, $\tau_a = 120\,\text{MPa}$) yield identical `PBP` single-harmonic indices on medium-carbon steel ($f_{-1} = 313.9\,\text{MPa}$, $t_{-1} = 196.2\,\text{MPa}$):

```python
import numpy as np
from multifat.models.invariant_methods.pbp import ProjectionByProjection

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"Sigm1": 313.9, "Taum1": 196.2}
v_ip = ProjectionByProjection().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params=p, harmonics=dict(S_mean=Z, S_sin=bt(200.0, 120.0), S_cos=Z)
).values[""]
v_np = ProjectionByProjection().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params=p,
    harmonics=dict(S_mean=Z, S_sin=bt(200.0, 0.0), S_cos=bt(0.0, 120.0)),
).values[""]
assert np.isclose(v_ip, v_np, rtol=1e-12) and np.isclose(v_ip, 0.8962, atol=1e-3)
print(f"PBP (in-phase) = {v_ip:.4f}, PBP (90-deg) = {v_np:.4f}")
```

Output:

```text
PBP (in-phase) = 0.8962, PBP (90-deg) = 0.8962
```

Batch CLI evaluation:

```bash
multifat-run --model PBP
```

## Comparison capabilities and outputs

`PBP` uses `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and outputs a dimensionless endurance-limit index `""`.

## Verification and limitations

- **Software verification.** [`tests/methods/test_pbp.py`](../../tests/methods/test_pbp.py) verifies the 5-vector $J_2$ norm identity, equivalence between 8,000-point sampled covariance projection and closed-form $\sqrt{J_2(A)+J_2(B)}$, exact Eq. (6)–(8) axial and torsional calibration ($1.0$), phase-shift invariance, and mean-stress $\rho_{\mathrm{ref}}$ dependence. Cristofori et al. (2008) print no per-specimen numerical validation table beyond Wöhler constants and plots.
- **Phase insensitivity on single harmonics.** In bending–torsion single-harmonic loading, the closed form (Appendix Eq. B.9) removes phase-angle sensitivity between normal and shear channels.

## References

- Cristofori, A., Susmel, L., & Tovo, R. (2008). A stress invariant based criterion to estimate fatigue damage under multiaxial loading. *International Journal of Fatigue*, 30(9), 1646–1658. [doi:10.1016/j.ijfatigue.2007.11.006](https://doi.org/10.1016/j.ijfatigue.2007.11.006)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
