# Liu virtual strain energy (`VSE`)

> **Registry name:** `VSE` (and `VSE_ext`) · **Catalogue entry:** [G15](../../METHODS.md#g15) ·
> **Source:** Liu (1993), ASTM STP 1191, <https://doi.org/10.1520/STP24796S> ·
> **Implementation:** [`cp_methods/liu_vse.py`](../../src/multifat/models/cp_methods/liu_vse.py)

## Overview

- **Type and family:** energy-based critical-plane criterion.
- **Idea:** on a material plane, the "virtual" strain energy is the product of the stress ranges
  and the strain ranges acting on that plane. It is split into:
  - a tensile part, Mode I: normal stress range times normal strain range;
  - a shear part, Mode II: shear stress range times shear strain range.

  Failure is governed by the plane on which the energy of the dominant cracking mode is largest.
- **Proposed for:** low- and high-cycle multiaxial fatigue of metals whose cracking mode (tensile
  or shear) depends on material and temperature.
- **Status:** fully implemented for the declared repository scope. Eqs. (7)–(11) checked on the
  PDF. The paper reports its results as figures only, so no table is reproduced.

## Original formulation

Liu (1993, Eq. 7) defines on each plane the virtual strain energy

```math
\Delta W = \Delta W_n + \Delta W_s = \Delta\sigma_n\,\Delta\varepsilon_n + \Delta\tau\,\Delta\frac{\gamma}{2},
```

where:

- $\Delta\sigma_n$ and $\Delta\varepsilon_n$ are the normal stress and normal strain ranges on the
  plane;
- $\Delta\tau$ and $\Delta(\gamma/2)$ are the shear stress range and the tensorial shear strain
  range on the plane.

Two failure modes are distinguished:

- **Mode I (Eq. 8):** the plane on which $\Delta W_n$ is largest; the damage is the total $\Delta W$
  on that plane.
- **Mode II (Eq. 9):** the plane on which $\Delta W_s$ is largest, with the same total.

Under uniaxial loading, $\Delta W$ reduces to $\Delta\sigma\,\Delta\varepsilon$ (Eq. 11), and in
torsion the two modes coincide (p. 73). Life follows from an energy–life curve with elastic and
plastic power-law terms (Eq. 6), fitted to the uniaxial data. The paper gives no rule to choose
the mode a priori: the mode depends on the material and the temperature.

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\Delta W_n = \left|\Delta\sigma_{zz}\,\Delta\varepsilon_{zz}\right|, \qquad
\Delta W_s = \left|\Delta\sigma_{xz}\,\Delta\varepsilon_{xz} + \Delta\sigma_{yz}\,\Delta\varepsilon_{yz}\right|,
\qquad \mathrm{DP} = \Delta W_n + \Delta W_s .
```

The shear term is the dot product of the shear stress range vector and the tensorial shear strain
range vector. DP has units of energy density, MPa (= MJ/m³).

- `VSE`: DP on the plane maximising $\Delta W_n$ (mode I, default) or $\Delta W_s$ (mode II,
  `vse_mode = 2`).
- `VSE_ext`: the largest total $\Delta W$ over all planes, i.e. Eq. (7) read literally.
  In general it differs from Eq. (8).

**Repository choices**

- *Ranges at the cycle endpoints.* All ranges are differences between the two endpoints of the
  repository cycle. For out-of-phase loading these endpoints are those of the longest chord of the
  shear path (`np_shear_rule = "longest_chord"`). Liu's rules for out-of-phase loading
  (pp. 80–81) are stated only in words; this rule is a documented adaptation.
- *Mode choice by the user.* One model is registered, mode I by default. Mode II is selected per
  workbook with the Summary column `vse_mode` set to 2.
- *Elastic strains.* With the workbook's elastic strains, the plastic term of Eq. (6) is absent.
- *Life.* The shared master curve replaces Eq. (6)
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| mode | `vse_mode` / Summary column | – | workbook | no; 1 (default) or 2 |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook | yes (strains) |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

VSE has no criterion-specific constant. The master-curve constants $A$, $b$ are fitted on the
uniaxial $R = -1$ design tests, the role played by Eq. (6) in the paper.

Consequences to keep in mind when reading results:

- **Torsion-to-tension limit ratio.** With elastic strains:
  - tension: mode I gives $\Delta W = 4f^2/E$ and mode II gives $2f^2/E$;
  - torsion: both modes give $4(1+\nu)t^2/E$.

  Equal damage at the two limits therefore implies $t/f = 1/\sqrt{1+\nu} \approx 0.88$ in mode I
  and $1/\sqrt{2(1+\nu)} \approx 0.62$ in mode II, for $\nu = 0.3$. Ductile metals, whose ratio is
  close to $1/\sqrt3$, are better described by mode II.
- **No mean-stress effect.** Only ranges enter the parameter, so a mean stress does not change
  $\Delta W$. Tests with a mean stress are predicted with the life of the corresponding fully
  reversed load.

## Worked example

Synthetic demonstration material (not experimental data): $E = 200$ GPa, $\nu = 0.3$.

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


vse = get_models(["VSE"])[0]
mode_1 = {"E": E, "nu": nu}  # vse_mode absent: mode I, Eq. (8)
mode_2 = {"E": E, "nu": nu, "vse_mode": 2}  # mode II, Eq. (9)

s_a, t_a = 300.0, 200.0
w1_tension, ext_tension = evaluate(vse, mode_1, zero, s_a * X)
w2_tension, _ = evaluate(vse, mode_2, zero, s_a * X)
w1_torsion, _ = evaluate(vse, mode_1, zero, t_a * T)
w2_torsion, _ = evaluate(vse, mode_2, zero, t_a * T)

# Eq. (11): plane normal to the load, Delta W = (2 s)(2 s / E)
assert abs(w1_tension - 4 * s_a**2 / E) < 1e-9 * w1_tension
# Mode II in tension: 45 deg plane, half of the mode I value (elastic strains)
assert abs(w2_tension - 2 * s_a**2 / E) < 1e-9 * w2_tension
# Torsion: both modes give Delta tau Delta(gamma/2) = 4 tau^2 (1+nu)/E
assert abs(w1_torsion - 4 * (1 + nu) * t_a**2 / E) < 1e-9 * w1_torsion
assert abs(w2_torsion - w1_torsion) < 1e-9 * w1_torsion

# Ranges only: a tensile mean stress leaves Delta W unchanged
w_mean, _ = evaluate(vse, mode_1, 150.0 * X, s_a * X)
assert abs(w_mean - w1_tension) < 1e-9 * w1_tension

print(f"VSE mode I, tension 300 MPa: {w1_tension:.3f} MPa (VSE_ext {ext_tension:.3f})")
print(f"VSE mode II, tension 300 MPa: {w2_tension:.3f} MPa")
print(f"VSE torsion 200 MPa, modes I and II: {w1_torsion:.3f}, {w2_torsion:.3f} MPa")
print(f"mode I t/f at equal DP: {1 / (1 + nu) ** 0.5:.3f}; mode II: {(0.5 / (1 + nu)) ** 0.5:.3f}")
print(f"VSE mode I with a 150 MPa mean: {w_mean:.3f} MPa")
```

Output:

```text
VSE mode I, tension 300 MPa: 1.800 MPa (VSE_ext 1.800)
VSE mode II, tension 300 MPa: 0.900 MPa
VSE torsion 200 MPa, modes I and II: 1.040, 1.040 MPa
mode I t/f at equal DP: 0.877; mode II: 0.620
VSE mode I with a 150 MPa mean: 1.800 MPa
```

**On a supplied workbook** (mode I unless the Summary sets `vse_mode`):

```bash
uv run multifat-pipeline --material Al7075_T6 --models VSE --all-variants --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`), not Eq. (6).
- **DP comparison:** supported, against the master-curve `DP_fit`. DP is quadratic in the load, so
  a DP error is about twice the corresponding stress error. `Error_load` removes this factor.
- **`_ext`:** available; Eq. (7) maximised over all planes.
- **Mixed results:** both modes write to `VSE_<MATERIAL>.csv`, so keep mode I and mode II runs in
  separate results folders.

## Verification and limitations

[`tests/methods/test_liu_vse.py`](../../tests/methods/test_liu_vse.py) checks:

- Eq. (11) in tension, for both modes;
- the coincidence of the two modes in torsion;
- that `VSE_ext` (Eq. 7) bounds the mode I value (Eq. 8).

The single fitted curve printed in the paper (Fig. 9) has coefficients 34 435 and 1 126, which do
not match its plot; 344.35 and 11.26 do.

Not verified: the experimental lives of Liu (1993).

Limitations:

- elastic strains only;
- no mean-stress effect;
- the mode must be chosen by the user;
- the out-of-phase ranges follow a repository rule, not an equation of the paper.

## References

- K. C. Liu (1993), A method based on virtual strain-energy parameters for multiaxial fatigue life
  prediction, in *Advances in Multiaxial Fatigue*, ASTM STP 1191, 67–84.
  <https://doi.org/10.1520/STP24796S>
