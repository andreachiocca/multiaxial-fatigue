# Li et al. (2021) improved strain energy density (`LI`)

> **Registry name:** `LI` (and `LI_ext`) · **Catalogue entry:** [I_LI](../../METHODS.md#i_li) ·
> **Source:** Li et al. (2021), <https://doi.org/10.1016/j.euromechsol.2021.104353> ·
> **Implementation:** [`cp_methods/li.py`](../../src/multifat/models/cp_methods/li.py)

## Overview

- **Type and family:** energy-based critical-plane criterion combining shear strain-energy density
  (with mean shear) and a weighted normal strain-energy term.
- **Idea:** add the magnitude of the mean shear stress $\lvert\tau_m\rvert$ to the alternating
  shear stress $\tau_a$ in the shear energy term, and weight the normal term $\sigma_{n,\max}\varepsilon_{n,a}$
  by a material constant $L$. The critical plane is the plane of maximum shear strain-energy density.
- **Mean stress:** both positive and negative mean shear stresses increase damage via
  $\lvert\tau_m\rvert$; normal mean stress enters through $\sigma_{n,\max}$ (tensile only).
- **Proposed for:** multiaxial high- and low-cycle fatigue of metals under proportional and
  non-proportional loading with mean stresses.
- **Status:** implemented for the declared repository scope ("requires review" in the catalogue).
  The mean-shear norm cancellation on reversed loading and plane-value assertions are tested;
  full-text equation-level reproduction of the paper's dataset tables remains outside current tests.

## Original formulation

Li et al. (2021, Eq. 3):

```math
\left(\frac{\Delta\tau_t}{2} + \lvert\tau_m\rvert\right)\frac{\Delta\gamma_{\max}}{2}
+ L\,\sigma_{n,\max}\,\frac{\Delta\varepsilon_n}{2} = f(N_f),
```

where on the candidate plane:

- $\Delta\tau_t/2$ and $\Delta\gamma_{\max}/2$ are the shear stress and **engineering** shear
  strain amplitudes;
- $\lvert\tau_m\rvert$ is the magnitude of the mean shear stress vector;
- $\sigma_{n,\max}$ and $\Delta\varepsilon_n/2$ are the maximum normal stress and normal strain
  amplitude;
- $L$ is a material constant calibrated from axial and torsional fatigue data (Section 2.4,
  Eqs. 31–32).

The critical plane is defined as the plane maximising the shear strain-energy term
$(\Delta\tau_t/2 + \lvert\tau_m\rvert)\Delta\gamma_{\max}/2$, resolving ties toward the larger
normal term.

## Repository implementation

On every candidate plane (normal = local $z$ axis, shear components $xz, yz$):

```math
\mathrm{DP} = (\tau_a + \lVert\boldsymbol\tau_m\rVert)\,\gamma_a
+ L\,\sigma_{n,\max}\,\varepsilon_{n,a},
```

with $\tau_a = \Delta\tau/2$, $\lVert\boldsymbol\tau_m\rVert = \sqrt{\tau_{xz,m}^2 + \tau_{yz,m}^2}$,
$\gamma_a = \Delta\gamma/2$, $\varepsilon_{n,a} = \Delta\varepsilon_n/2$, and
$\sigma_{n,\max} = \max(\sigma_{zz,0},\,\sigma_{zz,1},\,0)$. DP has units of MPa ($\mathrm{MJ/m^3}$).

- `LI`: DP on the plane maximising the shear energy term $(\tau_a + \lVert\boldsymbol\tau_m\rVert)\gamma_a$,
  with ties resolved toward the larger total DP.
- `LI_ext`: largest DP over all candidate planes (repository variant).
- Out-of-phase loading: shear range uses the longest chord of the shear path on the plane
  (`np_shear_rule = "longest_chord"`), matching the paper's cycle-range definition.

**Repository choices**

- *Mean shear vector norm.* $\lVert\boldsymbol\tau_m\rVert$ is the Euclidean norm of the mean shear
  vector $(\tau_{xz,m}, \tau_{yz,m})$, rather than the mean of endpoint norms, so fully reversed
  loading gives $\lVert\boldsymbol\tau_m\rVert = 0$.
- *Compression clipped.* $\sigma_{n,\max} < 0$ is clipped to $0$.
- *Optional constant $L$.* Read from `L_LI`, falling back to `1.0` with a `RuntimeWarning` when
  omitted in single-model evaluation (`multifat-run`); calibrated by `multifat-pipeline` when absent.
- *Elastic strains and life curve.* Workbook strains are isotropic linear elastic; the shared
  master curve $\mathrm{DP} = A N_f^b$ replaces the paper's torsional life curve
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $L$ | `L_LI` / Summary column `L LI` | – | workbook, `--constants`, pipeline calibration ($[0, 20]$), or default `1.0` | no (default `1.0`) |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook (for elastic strains) | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure (Li et al., 2021, Section 2.4, Eqs. 31–32).** First fit the four torsional
   strain–life constants ($\tau_f', \gamma_f', b_0, c_0$) on fully reversed torsion tests; then
   evaluate $L \in [0, 10]$ in steps of $0.01$ to minimise the axial fitting error against the
   torsional life curve.
2. **Repository (`multifat-pipeline`).** When `L_LI` is absent from the workbook and `--constants`,
   the pipeline searches $L_\mathrm{LI} \in [0, 20]$ to minimise the $\ln N_f$ scatter of the
   $R = -1$ tension and torsion failures around one master curve. (The search stops as
   `unidentified` if it hits the search bound.)
3. **Fatigue-limit datasets** (run-outs only): $L_\mathrm{LI} \in [0, 20]$ is chosen so that the
   axial ($f_{-1}$) and torsional ($t_{-1}$) limits yield equal DP. With elastic strains on the
   45° tension plane, $L_\mathrm{LI} = \frac{4(1+\nu)}{1-\nu}\left(\frac{t_{-1}^2}{f_{-1}^2} - \frac{1}{2}\right)$.

All tests contributing to the fitted $L_\mathrm{LI}$ are exported as design rows (`Flag = 0`) in
the CSV.

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$, $L = 1.5$.

```python
import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

E, nu = 200000.0, 0.3  # MPa, - (synthetic)
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))
X = np.diag([1.0, 0.0, 0.0])  # unit uniaxial stress along x
T = np.zeros((3, 3))
T[0, 1] = T[1, 0] = 1.0  # unit shear stress tau_xy


def hooke(S):
    """Isotropic linear-elastic strain tensor (tensorial shear components)."""
    return ((1 + nu) * S - nu * np.trace(S) * np.eye(3)) / E


def evaluate(model, params, S_mean, S_sin, S_cos=zero):
    """Standard and _ext DP of S(t) = S_mean + S_sin sin(wt) + S_cos cos(wt), elastic strains."""
    result = model.evaluate_case(
        S0=S_mean + S_sin,
        S1=S_mean - S_sin,
        E0=hooke(S_mean + S_sin),
        E1=hooke(S_mean - S_sin),
        params=params,
        R_list=planes,
        harmonics={
            "S_mean": S_mean,
            "S_sin": S_sin,
            "S_cos": S_cos,
            "E_sin": hooke(S_sin),
            "E_cos": hooke(S_cos),
        },
    )
    return result.values[""], result.values["_ext"]


li = get_models(["LI"])[0]
L = 1.5
params = {"L_LI": L}

s_a, t_a, t_m = 300.0, 200.0, 50.0
dp_tension, ext_tension = evaluate(li, params, zero, s_a * X)
dp_torsion, _ = evaluate(li, params, t_m * T, t_a * T)
dp_torsion_neg, _ = evaluate(li, params, -t_m * T, t_a * T)

# Tension on 45 deg max-shear-energy plane: tau_a = s_a/2, gamma_a = (1+nu)s_a/E,
# sigma_n,max = s_a/2, eps_n,a = (1-nu)s_a/(2E)
tension = (1 + nu) * s_a**2 / (2 * E) + L * (1 - nu) * s_a**2 / (4 * E)
# LI_ext reaches the normal plane (theta = 0 deg) for L = 1.5: DP = L * s_a^2 / E
ext_expected = L * s_a**2 / E
# Torsion with mean shear: DP = (tau_a + |tau_m|) * 2(1+nu)tau_a/E
torsion = (t_a + abs(t_m)) * 2 * (1 + nu) * t_a / E

assert abs(dp_tension - tension) < 1e-8 * tension
assert abs(ext_tension - ext_expected) < 1e-8 * ext_expected
assert abs(dp_torsion - torsion) < 1e-8 * torsion
assert abs(dp_torsion - dp_torsion_neg) < 1e-12

print(f"LI tension 300 MPa (L=1.5): {dp_tension:.6f} MPa (LI_ext {ext_tension:.6f} MPa)")
print(f"LI torsion 200 +/- 50 MPa: {dp_torsion:.6f} MPa (sign-reversed {dp_torsion_neg:.6f} MPa)")
```

Output:

```text
LI tension 300 MPa (L=1.5): 0.410625 MPa (LI_ext 0.675000 MPa)
LI torsion 200 +/- 50 MPa: 0.650000 MPa (sign-reversed 0.650000 MPa)
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models LI --all-variants --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`), fitted on the $R = -1$ uniaxial design set.
- **DP comparison:** supported against `DP_fit`. Because `LI` is homogeneous of degree 2 in the
  stress amplitude on elastic paths, $e_L = \tfrac{1}{2}\ln(\mathrm{DP}_{\text{fit}}/\mathrm{DP})$
  under proportional reversed loading (bisection when mean stress is present).
- Both `LI` and `LI_ext` are registered; `LI` selects maximum shear energy, `LI_ext` maximises
  the total DP.

## Verification and limitations

- [`tests/test_method_audit.py`](../../tests/test_method_audit.py) checks the analytical plane
  expression and confirms that reversed loading produces zero mean shear regardless of sign.
- Exact published parameter tables from Li et al. (2021) are not reproduced in `tests/methods/`;
  implementation status in the catalogue remains `requires_review`.
- Workbook strains are purely elastic; cyclic plasticity is not modelled.

## References

- P. Li, Y. Li, Y. Yang, B. Li, F. Wang (2021), Multiaxial fatigue life prediction for metals by
  means of an improved strain energy density-based critical plane criterion, *European Journal of
  Mechanics - A/Solids* 90, 104353. <https://doi.org/10.1016/j.euromechsol.2021.104353>
