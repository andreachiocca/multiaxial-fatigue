# Zhu et al. (2018) weighted generalized strain energy (`MGSE_ZHU`)

> **Registry name:** `MGSE_ZHU` (and `MGSE_ZHU_ext`) · **Catalogue entry:** [I_MGSE_ZHU](../../METHODS.md#i_mgse_zhu) ·
> **Source:** Zhu et al. (2018), <https://doi.org/10.1016/j.ijfatigue.2018.03.028> ·
> **Implementation:** [`cp_methods/mgse_zhu.py`](../../src/multifat/models/cp_methods/mgse_zhu.py)

## Overview

- **Type and family:** energy-based critical-plane criterion introducing a material weighting factor
  $k_\mathrm{MGSE}$ on the normal strain-energy term.
- **Idea:** generalise [MGSE_YU](MGSE_YU.md) across ductile and brittle metals by weighting the
  normal energy contribution $\sigma_{n,\max}\varepsilon_{n,a}$ with $k_\mathrm{MGSE}$.
- **Mean stress:** enters through $\tau_{\max}$ and the tensile peak normal stress $\sigma_{n,\max}$.
- **Proposed for:** multiaxial fatigue life prediction of both ductile and brittle metals.
- **Status:** implemented for the declared repository scope ("requires review" in the catalogue).
  Plane values, `zhu_k_mgse`, and pipeline calibration are tested; full-source reproduction of the
  paper's life benchmarks remains outside current tests.

## Original formulation

Zhu et al. (2018, Eqs. 7–12, pp. 7–8):

```math
\mathrm{MGSE} = \tau_{\max}\,\frac{\Delta\gamma}{2} + k_\mathrm{MGSE}\,\sigma_{n,\max}\,\frac{\Delta\varepsilon_n}{2},
```

where $k_\mathrm{MGSE}$ is a material-dependent weighting constant. In Eq. (12), $k_\mathrm{MGSE}$
is expressed as a function of life $N_f$ from the axial and torsional strain–life constants, and
Section 4.1 (pp. 12–13) averages it over a chosen finite-life range ($5000$–$50\,000$ cycles for
most metals in the paper), recommending $k_\mathrm{MGSE} = 1$ when fatigue properties are unavailable.

## Repository implementation

On every candidate plane (normal = local $z$ axis, shear components $xz, yz$):

```math
\mathrm{DP} = \tau_{\max}\,\gamma_a + k_\mathrm{MGSE}\,\sigma_{n,\max}\,\varepsilon_{n,a},
\qquad \sigma_{n,\max} = \max(\sigma_{zz,0},\,\sigma_{zz,1},\,0),
```

with $\gamma_a = \Delta\gamma/2$ and $\varepsilon_{n,a} = \Delta\varepsilon_n/2$. DP has units of
MPa ($\mathrm{MJ/m^3}$). Because the selection metric is DP itself, `MGSE_ZHU` and `MGSE_ZHU_ext`
return identical values.

**Repository choices**

- *Critical plane.* Maximises the full weighted MGSE value over all candidate planes; the
  plane rule and weighting remain provisional pending full-source review.
- *Shear path quantities.* Uses the default `rss` non-proportional shear range and the 720-sample
  harmonic maximum $\tau_{\max}$ via `NEEDS_TAU_MAX` (falling back to endpoint max outside
  `prepare_planes`).
- *Compression clipped.* $\sigma_{n,\max} < 0$ is clipped to $0$.
- *Default `k_MGSE`.* Optional; in single-model runs (`multifat-run`) without a workbook `k MGSE`,
  it defaults to `1.0` with a `RuntimeWarning`.
- *Elastic strains and life curve.* Workbook strains are linear elastic; life uses the shared
  master curve $\mathrm{DP} = A N_f^b$
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $k_\mathrm{MGSE}$ | `k_MGSE` / Summary column `k MGSE` | – | workbook, `--constants`, Zhu Eq. (12), scatter calibration ($[0, 5]$), or default `1.0` | no (default `1.0`) |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook (for elastic strains) | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure (Zhu et al., 2018, Eq. 12, Section 4.1).** Compute $k_\mathrm{MGSE}(N)$ from
   the axial and torsional strain–life constants and average it over a regime-dependent life
   interval (5,000–50,000 cycles for typical metals), or take $k_\mathrm{MGSE} = 1$ for ductile
   materials when data are missing.
2. **Repository (`multifat-pipeline`).** When `k_MGSE` is absent from the workbook and `--constants`:
   - First, if identifiable Basquin curves (`sigma_f`, `b_fat`, `tau_f`, `b0`) can be fitted on the
     workbook's $R = -1$ tension and torsion failures, the pipeline computes `k_MGSE` with
     `multifat.eval.material_fits.zhu_k_mgse` using the **elastic terms only** of Eq. (12) averaged
     over $N \in [10^4, 2\times 10^6]$ (independent of which other models are selected).
   - Otherwise, it falls back to minimising the $\ln N_f$ master-curve scatter over
     $k_\mathrm{MGSE} \in [0, 5]$ on the $R = -1$ tension and torsion failures.
3. **Fatigue-limit datasets** (run-outs only): $k_\mathrm{MGSE} \in [0, 5]$ equalises the axial
   and torsional limit DP.

The observations behind `k_MGSE` (Basquin or scatter tests) are exported as design rows (`Flag = 0`).

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$,
$k_\mathrm{MGSE} = 0.8$, plus `zhu_k_mgse` on sample Basquin constants.

```python
import numpy as np

from multifat.eval.material_fits import zhu_k_mgse
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


zhu = get_models(["MGSE_ZHU"])[0]
k = 0.8
params = {"k_MGSE": k}

s_a, t_a = 300.0, 200.0
dp_tension, ext_tension = evaluate(zhu, params, zero, s_a * X)
dp_torsion, _ = evaluate(zhu, params, zero, t_a * T)

# Elastic tension on plane c = cos^2 theta: s_a^2/E * [ (2(1+nu)-k*nu)*c - (2-k)*(1+nu)*c^2 ]
A_z = 2 * (1 + nu) - k * nu
B_z = (2 - k) * (1 + nu)
c_star = min(1.0, A_z / (2 * B_z))
tension = s_a**2 / E * (A_z * c_star - B_z * c_star**2)
torsion = 2 * (1 + nu) * t_a**2 / E

assert abs(dp_tension - tension) < 1e-7 * tension
assert abs(dp_tension - ext_tension) < 1e-12
assert abs(dp_torsion - torsion) < 1e-8 * torsion

k_eq12 = zhu_k_mgse(1200.0, -0.08, 700.0, -0.075, E, nu)
assert 0.0 < k_eq12 < 5.0

print(f"MGSE_ZHU tension 300 MPa (k=0.8): {dp_tension:.6f} MPa (MGSE_ZHU_ext {ext_tension:.6f} MPa)")
print(f"MGSE_ZHU torsion 200 MPa (k=0.8): {dp_torsion:.6f} MPa")
print(f"Zhu Eq. (12) k_MGSE (sample Basquin): {k_eq12:.4f}")
```

Output:

```text
MGSE_ZHU tension 300 MPa (k=0.8): 0.401654 MPa (MGSE_ZHU_ext 0.401654 MPa)
MGSE_ZHU torsion 200 MPa (k=0.8): 0.520000 MPa
Zhu Eq. (12) k_MGSE (sample Basquin): 0.7063
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MGSE_ZHU --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`).
- **DP comparison:** supported against `DP_fit`.
- `MGSE_ZHU` and `MGSE_ZHU_ext` are identical because the selection metric equals the DP.

## Verification and limitations

- [`tests/test_method_audit.py`](../../tests/test_method_audit.py) checks the closed-form plane
  value; [`tests/test_pipeline.py`](../../tests/test_pipeline.py) checks that `k_MGSE` calibration
  does not depend on whether other models (`MDC`) are also requested.
- Elastic-only `zhu_k_mgse` averaged over $10^4$–$2\times 10^6$ cycles is a repository convention
  differing from Zhu et al.'s $5000$–$50\,000$ cycle range; experimental benchmarks from the paper
  are not reproduced.

## References

- S.-P. Zhu, Z.-Y. Yu, J. Correia, A. De Jesus, F. Berto (2018), Evaluation and comparison of
  critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials,
  *International Journal of Fatigue* 113, 213–229.
  <https://doi.org/10.1016/j.ijfatigue.2018.03.028>
