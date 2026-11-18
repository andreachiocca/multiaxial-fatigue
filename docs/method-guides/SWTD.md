# SWTD — Deviatoric Smith–Watson–Topper Energy Parameter

**Registry name:** `SWTD`  
**Family:** Energy-based, deviatoric, master-curve life  
**Source file:** [`src/multifat/models/energy_based_methods/swtd.py`](../../src/multifat/models/energy_based_methods/swtd.py)

## Overview

`SWTD` replaces the normal-plane stress and strain of the classic Smith–Watson–Topper parameter with deviatoric stress and deviatoric strain components in order to apply an energy-type parameter even under compressive mean stresses where classical $\sigma_{n,\max}$ drops to zero.

## Original formulation

In Kujawski (2014), Eq. (25) (p. 7), the deviatoric Smith–Watson–Topper parameter is

```math
\mathrm{SWTD} = \max_{i \in \{1,2,3\}} \max\!\left(S_{\max,i}\,|e_{a,i}|,\; S_{\min,i}\,|e_{a,i}|\right),
```

where $S_{\max}, S_{\min}$ are the deviatoric stress tensors at the loading extrema and $e_{a,i}$ are the principal deviatoric strain amplitudes. Following Eqs. (23)–(24) of the paper, $S_{\max,i}$, $S_{\min,i}$, and $e_{a,i}$ are taken along the same principal axes.

## Repository implementation

For two-step histories `(S0, E0, S1, E1)`, the implementation diagonalizes the deviatoric strain amplitude $e_a = \frac{1}{2}\operatorname{dev}(E_0 - E_1)$ once, rotates $\operatorname{dev}(S_0)$ and $\operatorname{dev}(S_1)$ into that eigenvector basis (`shared eigenframe`), takes the diagonal entries $S_{\max,i}, S_{\min,i}$, and clamps negative products to $0$.

**Repository choices and departures.**
- *Shared eigenframe vs independent eigenvalue sorting:* Projecting $S_{\max}$ and $S_{\min}$ onto the principal axes of $e_a$ matches Eqs. (23)–(24) and reproduces Table 2 uniaxial compression–compression (independent eigenvalue sorting would double compressive transverse terms).
- *Non-proportional single harmonics:* Because Kujawski (2014) gives no non-proportional procedure, single-harmonic histories scan 360 phase angles $\theta \in [0,\pi)$ and evaluate `swtd_eq25` on each opposing pair $(\theta, \theta+\pi)$, returning the maximum. Under proportional or in-phase cycles this reduces identically to the two-step evaluation.
- *Pure torsion Eq. (25) vs printed Table 3:* In pure reversed torsion Eq. (25) yields $\tau_a \gamma_a / 2$, whereas Kujawski (2014) plots and tabulates $\tau_a \gamma_a / \sqrt{3}$ in Table 3 (path A); the repository preserves the published Eq. (25) formula without an empirical $2/\sqrt{3}$ rescaling.

## Material constants and inputs

No material constants are required (`required_params = set()`, `optional_params = None`). Only the full stress and strain extrema `(S0, E0, S1, E1)` or single-harmonic tensors (`S_mean`, `S_sin`, `S_cos`, `E_mean`, `E_sin`, `E_cos`) are read from the workbook. Note that under elastic–plastic compression–compression tests where Kujawski (2014) Table 2 evaluates $\nu_{\text{eff}} = 0.5 - (0.5 - \nu)\,\Delta\sigma_{\text{eq}} / (E\,\Delta\varepsilon_{\text{eq}})$, the workbook elastic strains reflect $\nu$ unless full strains are imported.

## Calibration

`SWTD` has no model-specific parameters to calibrate. Master-curve $S$–$N$ parameters are fitted from the damage values across tests when `multifat-run` predicts life.

## Worked example

Compare a two-step proportional cycle and a $90^\circ$ out-of-phase cycle in Al 2024-T351 ($E = 71.7\,\text{GPa}$, $\nu = 0.3$):

```python
import numpy as np
from multifat.models.energy_based_methods.swtd import SWTD

E, nu = 71700.0, 0.3
strain = lambda S: ((1.0 + nu) * S - nu * np.trace(S) * np.eye(3)) / E
bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))

S0, S1 = bt(300.0, 120.0), bt(-100.0, -120.0)
sw_ip = SWTD().evaluate_case(S0=S0, E0=strain(S0), S1=S1, E1=strain(S1), params={}).values[""]

h_np = dict(
    S_mean=Z, S_sin=bt(250.0, 0.0), S_cos=bt(0.0, 150.0),
    E_mean=Z, E_sin=strain(bt(250.0, 0.0)), E_cos=strain(bt(0.0, 150.0)),
)
sw_np = SWTD().evaluate_case(S0=Z, E0=Z, S1=Z, E1=Z, params={}, harmonics=h_np).values[""]
assert sw_ip > 0.0 and sw_np > 0.0
print(f"SWTD (two-step) = {sw_ip:.4f} MJ/m^3, SWTD (90-deg NP) = {sw_np:.4f} MJ/m^3")
```

Output:

```text
SWTD (two-step) = 0.8186 MJ/m^3, SWTD (90-deg NP) = 0.5110 MJ/m^3
```

Batch CLI evaluation:

```bash
multifat-run --model SWTD
```

## Comparison capabilities and outputs

`SWTD` uses `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and produces a single unlabelled output variant `""` representing energy density ($\text{MPa}\cdot\text{strain} = \text{MJ/m}^3$).

## Verification and limitations

- **Software verification (paper reproduction).** In [`tests/methods/test_swtd.py`](../../tests/methods/test_swtd.py), 10 of 11 rows of Table 2 (uniaxial compression–compression using $\nu_{\text{eff}}$) match within $0.7\%$ (row 11 is a strict `xfail` confirmed misprint: `0.335` printed vs `0.324` from Eq. 25), Table 3 torsion rows match $\frac{2}{\sqrt{3}}\,\text{SWTD} = \tau_a\gamma_a/\sqrt3$ for 9 of 10 rows (row 10 is an inconsistent strict `xfail`: `0.58` printed vs `0.61`), and rotation plus phase-origin invariance hold.
- **Paths B–D limitation.** Table 3 paths B, C, D in Kujawski (2014) do not reconcile with a single tensor convention and the paper omits any out-of-phase algorithm, so experimental non-proportional validation relies on the repository's 360-chord single-harmonic extension.

## References

- Kujawski, D. (2014). A deviatoric version of the SWT parameter. *International Journal of Fatigue*, 67, 95–102. [doi:10.1016/j.ijfatigue.2013.12.002](https://doi.org/10.1016/j.ijfatigue.2013.12.002)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
