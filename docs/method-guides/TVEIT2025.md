# TVEIT2025 — Tveit, Ottosen, Erduran & Reyes (2025) Continuum Damage Model

**Registry name:** `TVEIT2025`  
**Family:** Incremental continuum damage, native life  
**Source file:** [`src/multifat/models/invariant_methods/tveit2025.py`](../../src/multifat/models/invariant_methods/tveit2025.py)

## Overview

`TVEIT2025` integrates an incremental continuum endurance-surface model cycle by cycle. A back-stress tensor $\boldsymbol{\alpha}$ and scalar damage $D$ evolve whenever the stress excursion leaves the moving yield-like endurance surface ($\beta > 0$, $d\beta > 0$), delivering native fatigue life when $D$ reaches $1$ or infinite life upon elastic shakedown.

## Original formulation

In Tveit, Ottosen, Erduran & Reyes (2025), Eqs. (1)–(12) and (24) (pp. 2–5), the moving endurance surface and damage evolution laws in deviatoric stress $\mathbf{s}$ and hydrostatic invariant $I_1$ are

```math
\beta = \frac{\sigma_e + A\,I_1 - f_{-1}}{f_{-1}},
\qquad
\sigma_e = \sqrt{3 J_2} - \lambda\!\left(1 - \frac{27}{4}\frac{J_3^2}{J_2^3}\right),
\qquad
\lambda = \sqrt{3}\,t_{-1} - f_{-1},
```

with invariants $J_2, J_3$ evaluated on the reduced deviator $\mathbf{s} - \boldsymbol{\alpha}$ and active evolution (when $\beta > 0$ and $d\beta > 0$):

```math
d\boldsymbol{\alpha} = C\,(\mathbf{s} - \boldsymbol{\alpha})\,d\beta,
\qquad
dD = K\,e^{L\beta}\,d\beta,
\qquad
d\beta = \frac{\mathbf{G} : d\mathbf{s} + A\,dI_1}{f_{-1} + C\,\mathbf{G} : (\mathbf{s} - \boldsymbol{\alpha})},
\qquad
\mathbf{G} = \frac{\partial\sigma_e}{\partial\mathbf{s}}.
```

Starting from $\boldsymbol{\alpha} = \mathbf{0}$ and $D = 0$, failure occurs at $D = 1$, while paths enclosed by the drifted endurance surface shakedown to $N_f = \infty$.

## Repository implementation

**Block integration and steady-state acceleration.** Each harmonic period is sampled at `npts_TVEIT = 1001` points and integrated explicitly block by block until one of four rules triggers:
1. *Failure:* $D \ge 1$ inside a block (linearly interpolated within the crossing increment).
2. *Shakedown:* a block (from block 2 onward) accumulates zero damage ($N_f = \infty$).
3. *Cyclic steady state:* after `min_blocks_TVEIT = 400` cycles, back-stress drift over the last 200 blocks is $\le 0.01$ of its path length and mean block damage agrees within $0.1\%$ of the preceding 200-block window—provided peak $\beta$ over the block is not trending down toward shakedown before $D = 1$. Life is then extrapolated as $N_f = n + (1 - D)/\overline{dD}$.
4. *Block cap:* after `max_blocks_TVEIT = 60000` cycles, life is extrapolated similarly and a warning is emitted.

**Repository choices and departures.**
- *No internal `C, K, L` optimizer:* Calibration of $(C, K, L)$ (Eq. 25, after Ottosen et al. 2008) is not implemented inside `multifat-calibrate`; constants must be supplied in a JSON `--constants` file (e.g., [`Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json`](../../Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json)).
- *Constants-file limits and settings:* `Sigm1_TVEIT` and `Taum1_TVEIT` override workbook `Sigm1`/`Taum1` for `TVEIT2025` only so calibrated constants always match their reference limits. Optional `npts_TVEIT`, `min_blocks_TVEIT`, `max_blocks_TVEIT`, `tol_alpha_TVEIT`, and `tol_damage_TVEIT` may also be overridden in the file. With `1001` points per cycle, finite lives sit roughly $0.02$ decade shorter than converged high-resolution integration (Table 1 NiCrMo steel).
- *Symmetric `G` projection and DP:* `G` uses the deviatoric projection $\mathbf{s}\mathbf{s} - \frac{2}{3} J_2 \mathbf{I}$ of $\partial J_3 / \partial \mathbf{s}$ (Eq. 10 prints a different isotropic piece that cancels upon deviatoric contraction). `values[""]` holds the virgin-state maximum of $(\sigma_e + A I_1)/f_{-1}$ ($\beta + 1$ at $\boldsymbol{\alpha} = 0$), while `aux["Nf"]` returns native life.

## Material constants and inputs

| Symbol | Workbook / JSON key | Python name | Unit | Notes |
|---|---|---|---|---|
| $A, C, K, L$ | `A_TVEIT`, `C_TVEIT`, `K_TVEIT`, `L_TVEIT` | same | — | Required model constants supplied via `--constants` |
| $f_{-1}, t_{-1}$ | `Sigm1_TVEIT`, `Taum1_TVEIT` | same | MPa | Optional limits replacing workbook `Sigm1`, `Taum1` for this model |
| solver controls | `npts_TVEIT`, `min_blocks_TVEIT`, `max_blocks_TVEIT`, `tol_alpha_TVEIT`, `tol_damage_TVEIT` | same | — | Optional integration settings (defaults `1001`, `400`, `60000`, `0.01`, `0.001`) |

## Calibration

`TVEIT2025` is not auto-fitted by `multifat-calibrate`; `A_TVEIT`, `C_TVEIT`, `K_TVEIT`, and `L_TVEIT` are marked as supplied constants (`--constants`). Without them, `multifat-run` automatically skips datasets lacking `TVEIT2025` constants (`LOAD_ERROR_SKIP = {"TVEIT2025"}`).

## Worked example

Evaluate $90^\circ$ out-of-phase bending–torsion ($\sigma_a = 660\,\text{MPa}$, $\tau_a = 340\,\text{MPa}$) using the Stage 1 constants for 42CrMo4 Q&T:

```python
import numpy as np
from multifat.models.invariant_methods.tveit2025 import Tveit2025

bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {
    "Sigm1_TVEIT": 433.0, "Taum1_TVEIT": 344.0,
    "A_TVEIT": 0.245625, "C_TVEIT": 0.0040672184319305,
    "K_TVEIT": 1.91047799401683e-05, "L_TVEIT": 10.1293188952565,
}
res = Tveit2025().evaluate_case(
    S0=Z, E0=Z, S1=Z, E1=Z, params=p,
    harmonics=dict(S_mean=Z, S_sin=bt(660.0, 0.0), S_cos=bt(0.0, 340.0)),
)
assert np.isclose(res.aux["Nf"], 263.0, atol=1.0)
print(f"index = {res.values['']:.4f}, Nf = {res.aux['Nf']:.0f}")
```

Output:

```text
index = 1.8986, Nf = 263
```

Batch CLI evaluation on 42CrMo4 Q&T with Stage 1 constants:

```bash
multifat-pipeline Experimental_Data/6_42CrMo4_Santus.xlsx --models TVEIT2025 --no-calibrate --constants Calibrations/42CrMo4_Santus/TVEIT2025_stage1.json
```

## Comparison capabilities and outputs

`TVEIT2025` sets `ComparisonCapabilities(life_prediction="native", dp_comparison=False, dp_reference="none")`: only native life predictions are compared against experiments (`aux["Nf"]`), and no universal-curve DP comparison is performed because $(C, K, L)$ already encapsulate the damage-life law.

## Verification and limitations

- **Software verification (paper reproduction).** [`tests/methods/test_tveit2025.py`](../../tests/methods/test_tveit2025.py) reproduces all published $S$–$N$ curves in Figs. 5–6 across Table 1 materials within $1\%$ in stress, Eq. (24) torsional limit, infinite life inside the endurance surface and after mean-shear shakedown, and constants-file parameter overrides.
- **Computational cost & external calibration.** Finite-life evaluation integrates block by block until failure or steady state (up to 60,000 blocks unless bounded), and `A, C, K, L` calibration remains external to the automated CLI pipeline.

## References

- Tveit, Ø., Ottosen, N. S., Erduran, E., & Reyes, A. (2025). Continuum mechanics fatigue model for arbitrary ratios of shear and normal stress fatigue strengths. *Mechanics of Materials*, 209, 105452. [doi:10.1016/j.mechmat.2025.105452](https://doi.org/10.1016/j.mechmat.2025.105452)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
