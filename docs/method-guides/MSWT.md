# Modified Smith–Watson–Topper (`MSWT`)

> **Registry name:** `MSWT` (and `MSWT_ext`) · **Catalogue entry:** [I_MSWT](../../METHODS.md#i_mswt) ·
> **Source:** Zhu et al. (2018), <https://doi.org/10.1016/j.ijfatigue.2018.03.028> ·
> **Implementation:** [`cp_methods/mswt.py`](../../src/multifat/models/cp_methods/mswt.py)

## Overview

- **Type and family:** energy-based critical-plane criterion combining normal and shear strain-energy
  ranges with a weighting parameter $a$.
- **Idea:** evaluate the scalar range form of Jiang's plastic-work increment
  (Jiang, 2000, Eq. 6; Jiang et al., 2007, Eq. 4) as written in Zhu et al. (2018, Eq. 13), so the
  parameter weights tensile versus shear cracking modes through $a \in [0, 1]$.
- **Mean stress:** enters through the tensile peak normal stress $\langle\sigma_{n,\max}\rangle$.
- **Proposed for:** multiaxial fatigue life prediction across tensile-, mixed-, and shear-cracking
  metals.
- **Status:** implemented for the declared repository scope ("requires review" in the catalogue).
  Implements the scalar range form of Zhu et al. (2018), not Jiang's incremental cyclic plasticity
  criterion.

## Original formulation

Zhu et al. (2018, Eq. 13, p. 8) write the modified SWT parameter on the maximum-damage plane as:

```math
\mathrm{DP}_{\mathrm{MSWT}} = 2\,a\,\Delta\varepsilon\,\langle\sigma_{\max}\rangle
+ \frac{1 - a}{2}\,\Delta\gamma\,\Delta\tau,
```

where $\langle x \rangle = \max(x, 0)$, $\Delta\varepsilon$ is the normal strain range,
$\sigma_{\max}$ is the maximum normal stress, and $\Delta\gamma$ and $\Delta\tau$ are the
engineering shear strain and shear stress ranges on the candidate plane. Under plastic
incompressibility ($\nu_p = 0.5$), $a \le 0.375$ favours shear planes, $0.375 < a < 0.5$ gives
mixed planes, and $a \ge 0.5$ favours tensile planes.

## Repository implementation

On every candidate plane (normal = local $z$ axis, shear components $xz, yz$):

```math
\mathrm{DP} = 2\,a\,\Delta\varepsilon_n\,\max(\sigma_{zz,0},\,\sigma_{zz,1},\,0)
+ \frac{1 - a}{2}\,\Delta\gamma\,\Delta\tau
= 4\,a\,\sigma_{n,\max}\,\varepsilon_{n,a} + 2(1 - a)\,\tau_a\,\gamma_a,
```

with $\Delta\gamma = 2\lVert(\Delta\varepsilon_{xz},\,\Delta\varepsilon_{yz})\rVert$ and
$\Delta\tau = \lVert(\Delta\tau_{xz},\,\Delta\tau_{yz})\rVert$. DP has units of MPa ($\mathrm{MJ/m^3}$).
Because the selection metric equals DP itself, `MSWT` and `MSWT_ext` return identical values.

**Repository choices**

- *Scalar range form.* Evaluates Zhu et al.'s range equation on elastic workbook strains rather than
  integrating incremental cyclic plasticity.
- *Shear rule.* Non-proportional shear range uses the default `rss` convention.
- *Constant validation and default.* `a_MSWT` is validated finite in $[0, 1]$ (raising `ValueError`
  outside); when omitted in `multifat-run`, it defaults to `0.5` with a `RuntimeWarning`.
- *Weighting correction (2026-10-18, finding U2).* Earlier revisions halved the shear coefficient,
  giving half its intended weight; when loading pre-fix upstream calibration records, `load_constants`
  automatically converts old values via $a' = a / (2 - a)$.
- *Life curve.* Uses the shared power-law master curve $\mathrm{DP} = A N_f^b$
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $a$ | `a_MSWT` / Summary column `a MSWT` | – | workbook, `--constants`, pipeline calibration ($[0, 1]$), or default `0.5` | no ($0 \le a \le 1$, default `0.5`) |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook (for elastic strains) | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure (Zhu et al., 2018, Section 4.1, p. 14).** Identify $a$ by bringing the
   axial and torsional damage-parameter–versus–life curves as close together as possible.
2. **Repository (`multifat-pipeline`).** When `a_MSWT` is absent from the workbook and `--constants`,
   the pipeline searches $a_\mathrm{MSWT} \in [0, 1]$ to minimise the $\ln N_f$ scatter of the
   $R = -1$ tension and torsion failures around one master curve. Unlike arbitrary search bounds,
   $a_\mathrm{MSWT} = 1.0$ is accepted as a valid physical upper limit (`PHYSICAL_UPPER`).
3. **Fatigue-limit datasets** (run-outs only): $a_\mathrm{MSWT} \in [0, 1]$ equalises the axial and
   torsional limit DP.

All tests contributing to the fitted $a_\mathrm{MSWT}$ are exported as design rows (`Flag = 0`).

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$. Notice that
in pure torsion ($\tau_a = 200$ MPa), the 0° shear plane gives $4(1-a)(1+\nu)\tau_a^2/E$ while the
45° principal-stress plane gives $4a(1+\nu)\tau_a^2/E$: for $a < 0.5$ the shear plane controls
and for $a > 0.5$ the normal plane controls, giving identical DP for $a = 0.4$ and $a = 0.6$.

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


mswt = get_models(["MSWT"])[0]
s_a, t_a = 300.0, 200.0

dp_tension, ext_tension = evaluate(mswt, {"a_MSWT": 0.4}, zero, s_a * X)
dp_tor_04, _ = evaluate(mswt, {"a_MSWT": 0.4}, zero, t_a * T)
dp_tor_06, _ = evaluate(mswt, {"a_MSWT": 0.6}, zero, t_a * T)

# Elastic tension over c = cos^2 theta in [nu/(1+nu), 1]
cs = np.linspace(nu / (1 + nu), 1.0, 10001)
tension = float(
    np.max(4 * s_a**2 / E * (0.4 * cs * ((1 + nu) * cs - nu) + 0.6 * (1 + nu) * cs * (1 - cs)))
)
torsion_max = 4 * 0.6 * (1 + nu) * t_a**2 / E

assert abs(dp_tension - tension) < 1e-7 * tension
assert abs(dp_tension - ext_tension) < 1e-12
assert abs(dp_tor_04 - torsion_max) < 1e-8 * torsion_max
assert abs(dp_tor_06 - torsion_max) < 1e-8 * torsion_max

print(f"MSWT tension 300 MPa (a=0.4): {dp_tension:.6f} MPa (MSWT_ext {ext_tension:.6f} MPa)")
print(f"MSWT torsion 200 MPa (a=0.4): {dp_tor_04:.6f} MPa; (a=0.6): {dp_tor_06:.6f} MPa")
```

Output:

```text
MSWT tension 300 MPa (a=0.4): 0.720000 MPa (MSWT_ext 0.720000 MPa)
MSWT torsion 200 MPa (a=0.4): 0.624000 MPa; (a=0.6): 0.624000 MPa
```

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MSWT --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`).
- **DP comparison:** supported against `DP_fit`.
- `MSWT` and `MSWT_ext` return identical values.

## Verification and limitations

- [`tests/test_method_audit.py`](../../tests/test_method_audit.py) verifies the Eq. (13) scalar
  value; [`tests/test_pipeline.py`](../../tests/test_pipeline.py) verifies the $a/(2-a)$ migration
  on legacy calibration records.
- Attribution/scope review and experimental benchmark reproduction from Zhu et al. (2018) remain
  outside current tests.

## References

- S.-P. Zhu, Z.-Y. Yu, J. Correia, A. De Jesus, F. Berto (2018), Evaluation and comparison of
  critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials,
  *International Journal of Fatigue* 113, 213–229.
  <https://doi.org/10.1016/j.ijfatigue.2018.03.028>
- Y. Jiang (2000), A fatigue criterion for general multiaxial loading, *Fatigue & Fracture of
  Engineering Materials & Structures* 23(1), 19–32.
  <https://doi.org/10.1046/j.1460-2695.2000.00247.x>
