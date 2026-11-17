# MMCCM — Modified Manson–Coffin Curve Method

**Registry name:** `MMCCM`  
**Family:** Critical plane, strain–stress, native life  
**Source file:** [`src/multifat/models/cp_methods/mmccm.py`](../../src/multifat/models/cp_methods/mmccm.py)

## Overview

`MMCCM` extends Susmel and Lazzarin's Modified Wöhler Curve Method to the low- and medium-cycle strain–life regime by interpolating the Manson–Coffin axial and torsional strain–life constants using the plane stress ratio $\rho = \sigma_{n,\max} / \tau_a$.

## Original formulation

In Wang & Susmel (2016), Eqs. (1)–(8), the critical plane is the plane experiencing the largest engineering shear-strain amplitude $\gamma_a$. On that plane, with stress ratio $\rho = \sigma_{n,\max} / \tau_a$ (Eq. 1), the life $N_f$ solves the interpolated Coffin–Manson equation (Eq. 2)

```math
\gamma_a = \frac{\tau'_f(\rho)}{G}\,(2 N_f)^{b(\rho)} + \gamma'_f(\rho)\,(2 N_f)^{c(\rho)},
```

where the coefficients vary linearly with $\rho$ (Eqs. 5–6)

```math
\frac{\tau'_f(\rho)}{G} = \rho\,(1 + \nu_e)\,\frac{\sigma'_f}{E} + (1 - \rho)\,\frac{\tau'_f}{G}, \qquad
\gamma'_f(\rho) = \rho\,(1 + \nu_p)\,\varepsilon'_f + (1 - \rho)\,\gamma'_f,
```

and the exponents follow weighted harmonic means (Eqs. 7–8)

```math
b(\rho) = \frac{b\,b_0}{(b_0 - b)\,\rho + b}, \qquad
c(\rho) = \frac{c\,c_0}{(c_0 - c)\,\rho + c}.
```

At $\rho = 1$ (reversed uniaxial loading) and $\rho = 0$ (reversed torsion), Eqs. (5)–(8) recover the axial (Eq. 3) and torsional (Eq. 4) Manson–Coffin curves. Ties among planes sharing the maximum $\gamma_a$ are broken in favour of the plane with the largest $\rho$ (p. 5).

## Repository implementation

**Critical plane and tie breaking.** `MMCCM` overrides the out-of-phase shear rule to `longest_chord`, which matches the maximum-variance half-range on single-harmonic paths. Ties on $\gamma_a$ are broken by selecting the largest per-plane damage $1 / N_f$; because larger $\rho$ lowers the fatigue curve at fixed $\gamma_a$, this automatically chooses the plane with the largest $\rho$.

**Repository choices and departures.**
- *Strains and elastic modulus:* The repository evaluates workbook (elastic) strain tensors with $\nu_e = \nu$ and requires all eight uniaxial and torsional strain–life constants (no Universal Material Law estimation fallback).
- *Plastic Poisson's ratio:* $\nu_p$ is not specified in Wang & Susmel (2016); the repository defaults to $\nu_p = 0.5$ unless overridden by optional `nu_p`.
- *Life inversion and singularities:* Eq. (2) is inverted by bisection in $\log(2 N_f)$ over $[1,\,2\times 10^{12}]$, returning damage $0$ beyond $10^{12}$ cycles. No cut-off on $\rho$ is applied; if high mean stress triggers a pole in Eqs. (7)–(8) (non-positive coefficient or non-negative exponent), damage is set to $+\infty$. Because high-$\rho$ poles already maximize damage, no `_ext` variant is exposed (`ext_variant = False`).

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $E, \nu$ | `E`, `nu` | `E`, `nu` | MPa, — | Elastic modulus and Poisson's ratio ($G = E / (2(1+\nu))$) |
| $\sigma'_f, \varepsilon'_f, b, c$ | `sigma_f`, `epsilon_f`, `b_fat`, `c_fat` | `sigma_f`, `epsilon_f`, `b_fat`, `c_fat` | MPa, — | Axial Manson–Coffin constants |
| $\tau'_f, \gamma'_f, b_0, c_0$ | `tau_f`, `gamma_f`, `b0`, `c0` | `tau_f`, `gamma_f`, `b0`, `c0` | MPa, — | Torsional Manson–Coffin constants |
| $\nu_p$ | `nu_p` | `nu_p` | — | Optional plastic Poisson's ratio (default `0.5`) |

## Calibration

Direct runs with `multifat-run` require all eight strain–life constants to be supplied in the workbook (`CALIBRATED_OPTIONAL` exposes no `MMCCM` plane parameters). However, `multifat-pipeline` (and `multifat-calibrate --basquin`) fills missing constants: the four Basquin parameters ($\sigma'_f, b, \tau'_f, b_0$) from fits on the workbook's $R=-1$ tension and torsion tests (only when identifiable; the torsion tests become design rows), and the four plastic terms ($\varepsilon'_f, c, \gamma'_f, c_0$) via Uniform Material Law estimates (steel form), because elastic high-cycle data cannot identify them.

## Worked example

Evaluate combined reversed bending and torsion ($\sigma_a = 280\,\text{MPa}$, $\tau_a = 140\,\text{MPa}$) on a steel parameter set:

```python
import numpy as np
from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.mmccm import MMCCM
from multifat.planes.sampling import build_R_list

p = {
    "E": 200000.0, "nu": 0.3,
    "sigma_f": 900.0, "epsilon_f": 0.4, "b_fat": -0.09, "c_fat": -0.6,
    "tau_f": 520.0, "gamma_f": 0.69, "b0": -0.08, "c0": -0.55,
}
S = np.array([[280.0, 140.0, 0.0], [140.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
Eps = ((1.0 + p["nu"]) * S - p["nu"] * np.trace(S) * np.eye(3)) / p["E"]
res = best_cp_for_case(
    S0=S, S1=-S, E0=Eps, E1=-Eps, R_list=build_R_list(), method=MMCCM(), params=p
)
nf = 1.0 / res.by_metric
assert np.isclose(nf, 110221.0, rtol=1e-3)
print(f"damage = {res.by_metric:.6e}, Nf = {nf:.0f}")
```

Output:

```text
damage = 9.072641e-06, Nf = 110221
```

Batch CLI evaluation:

```bash
multifat-run --model MMCCM
```

## Comparison capabilities and outputs

`MMCCM` sets `ComparisonCapabilities(life_prediction="native", dp_comparison=True, dp_reference="model_law")`. The criterion outputs dimensionless damage $D = 1 / N_f$ directly from Eq. (2), and predicted life is $1 / D$.

## Verification and limitations

- **Software verification.** Exact reduction to Eqs. (3)–(4) at $\rho \in \{0, 1\}$, harmonic-mean interpolation (Eqs. 7–8), critical-plane tie breaking on largest $\rho$, run-out truncation at $10^{12}$ cycles, and singular $\rho$ pole handling are verified in [`tests/methods/test_mmccm.py`](../../tests/methods/test_mmccm.py).
- **Elastic workbook strain scope.** Workbook imports supply elastic strains only; no cyclic plasticity solver or strain estimation fallback is included for unlabelled low-cycle plastic strain inputs. Also note there is no `MMCCM_ext` variant.

## References

- Wang, Y., & Susmel, L. (2016). The Modified Manson–Coffin Curve Method to estimate fatigue lifetime under complex constant and variable amplitude multiaxial fatigue loading. *International Journal of Fatigue*, 83, 135–149. [doi:10.1016/j.ijfatigue.2015.10.005](https://doi.org/10.1016/j.ijfatigue.2015.10.005)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
