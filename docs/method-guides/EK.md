# EK — Ellyin & Kujawski (1993) Distortional Strain Energy

**Registry name:** `EK`  
**Family:** Energy-based, distortional, master-curve life  
**Source file:** [`src/multifat/models/energy_based_methods/ellyin_kujawski.py`](../../src/multifat/models/energy_based_methods/ellyin_kujawski.py)

## Overview

`EK` scores fatigue severity using the cyclic distortional strain energy range $\Delta W_d$ scaled by a hydrostatic mean-stress multiplier and divided by a multiaxial triaxiality constraint factor $\bar{\rho}^n$ constructed from in-plane free-surface principal strains.

## Original formulation

In Ellyin & Kujawski (1993), Eqs. (5b) and (11)–(18), the damage parameter $\psi$ is

```math
\psi = \frac{\Delta W_d}{\bar{\rho}^{\,n}}
       \left(1 + \alpha\,\frac{\sigma_m}{\bar{\sigma}}\right),
\qquad
E^* = \frac{3 E}{2(1 + \nu)},
```

where $\Delta W_d$ is the distortional strain-energy range per reversal (elastic contribution $\Delta W_d^e = \Delta\bar{\sigma}^2 / (2 E^*)$ with $\Delta\bar{\sigma} = 2\bar{\sigma}$), $\sigma_m = \sigma_{1,m} + \sigma_{2,m} + \sigma_{3,m}$ is the sum of the three principal mean stresses (Eq. 5b, without a $1/3$ factor), and the triaxiality constraint factor is

```math
\bar{\rho} = (1 + \bar{\nu})\,\frac{\varepsilon_{\max}}{\gamma_{\max}},
```

with effective Poisson's ratio $\bar{\nu} = \nu_e$ in the high-cycle elastic limit ($\bar{\rho} = 1$ uniaxial, $1+\nu$ torsion, $1-\nu$ equibiaxial). Default recommended constants on p. 61 (Eq. 18) are $\alpha = 0.5$ and $n = 1.0$.

## Repository implementation

**Repository choices and departures.**
- *High-cycle elastic default with optional plasticity:* Because the cyclic stress–strain curve (`K_prime`, `n_prime`) is rarely present in HCF workbooks and is never fitted, `EK` defaults to the elastic distortional energy $\Delta W_d^e$. When both `K_prime` and `n_prime` are supplied in the workbook, a Masing plastic term $\Delta W_d^p = \frac{1 - n'}{2(1+n')}\,\Delta\bar{\sigma}\,\Delta\varepsilon_p$ with $\varepsilon_{p,a} = (\bar{\sigma}/K')^{1/n'}$, $\nu_p = 0.5$, and strain-weighted $\bar{\nu}$ (matching Eq. 13 in uniaxial tension) is added under assumed stress control.
- *Free-surface convention:* The repository assumes the $x$–$y$ plane lies on the free surface ($z$ normal) and derives out-of-plane elastic normal strain via $\varepsilon_{zz} = -\frac{\nu}{1-\nu}(\varepsilon_{xx} + \varepsilon_{yy})$ inside `constraint_factor`.
- *Non-proportional single harmonics:* Ellyin & Kujawski (1993) leave non-proportional cycles open. For single-harmonic histories, the repository defines $\bar{\sigma} = \sqrt{3}$ times the longest deviatoric stress half-chord and evaluates $\varepsilon_{\max}$ and $\gamma_{\max}$ as maxima over 360 phase samples around one cycle.
- *Mean-stress validity bounds:* Eq. (18) holds for $-1.06 \le \sigma_m/\sigma_a \le 4.4$ (p. 61). The repository reads the ratio as $\sigma_m/\bar{\sigma}$, the same ratio as the mean-stress factor, and a case outside raises `LoadingScopeError`: the model is skipped for that material (maintainer decision, 17 November 2026; such cases were evaluated before). No supplied workbook has such a case, so no output changed. A non-positive factor $(1 + \alpha\,\sigma_m/\bar{\sigma})$, possible only with a user `alpha_EK`, still clamps $\psi$ to $0$.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $E, \nu$ | `E`, `nu` | `E`, `nu` | MPa, — | Required elastic constants |
| $\alpha$ | `alpha_EK` | `alpha` | — | Optional mean-stress weight (default `0.5`) |
| $n$ | `n_EK` | `n` | — | Optional constraint exponent (default `1.0`) |
| $K', n'$ | `K_prime`, `n_prime` | `K_prime`, `n_prime` | MPa, — | Optional cyclic Ramberg–Osgood pair activating plasticity |

## Calibration

No optional parameters of `EK` are calibrated statistically; only the shared master Basquin curve is fitted across experimental tests.

## Worked example

Verify the constraint factor limits ($\bar{\rho} = 1.00$ uniaxial, $1.30$ torsion) and distortional energy parameter $\psi$ in steel ($E = 206\,\text{GPa}$, $\nu = 0.3$):

```python
import numpy as np
from multifat.models.energy_based_methods.ellyin_kujawski import EllyinKujawski, constraint_factor

E, nu = 206000.0, 0.3
strain = lambda S: ((1.0 + nu) * S - nu * np.trace(S) * np.eye(3)) / E
bt = lambda s, t: np.array([[s, t, 0.0], [t, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
p = {"E": E, "nu": nu}

for label, S in [("uniaxial", bt(300.0, 0.0)), ("torsion", bt(0.0, 200.0))]:
    Eps = strain(S)
    p_bar = constraint_factor(Eps, Z, nu)
    psi = EllyinKujawski().evaluate_case(
        S0=Z, E0=Z, S1=Z, E1=Z, params=p,
        harmonics=dict(S_mean=Z, S_sin=S, S_cos=Z, E_sin=Eps, E_cos=Z),
    ).values[""]
    print(f"{label}: p_bar = {p_bar:.2f}, psi = {psi:.5f} MJ/m^3")
```

Output:

```text
uniaxial: p_bar = 1.00, psi = 0.75728 MJ/m^3
torsion: p_bar = 1.30, psi = 0.77670 MJ/m^3
```

Batch CLI evaluation:

```bash
multifat-run --model EK
```

## Comparison capabilities and outputs

`EK` sets `ComparisonCapabilities(life_prediction="master_curve", dp_comparison=True, dp_reference="master_curve")` and outputs a single scalar damage parameter `""` in $\text{MJ/m}^3$ ($\text{MPa}\cdot\text{strain}$).

## Verification and limitations

- **Software verification (analytical paper limits).** In [`tests/methods/test_ellyin_kujawski.py`](../../tests/methods/test_ellyin_kujawski.py), $\bar{\rho} \in \{1, 1+\nu, 1-\nu\}$ holds across uniaxial, torsional, and equibiaxial stress states; Eq. (18) uniaxial mean-stress equivalence holds from $\sigma_m/\sigma_a = -0.8$ to $+4.0$; torsion/tension ratios match $\sqrt{(1+\nu)/3}$ ($n=1$) and $1/\sqrt{3}$ ($n=0$); mean shear stress has no effect; and Masing plasticity recovers Eq. (13) and uniaxial closed-form energy.
- **No printed per-test tables.** Ellyin & Kujawski (1993) publish analytical forms and figures without tabulated per-specimen $\psi$ values, and the free-surface plane assumption presupposes plane-stress surface tensors in the $x$–$y$ frame.

## References

- Ellyin, F., & Kujawski, D. (1993). A multiaxial fatigue criterion including mean-stress effect. In *Advances in Multiaxial Fatigue*, ASTM STP 1191 (pp. 55–66). [doi:10.1520/STP24795S](https://doi.org/10.1520/STP24795S)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
