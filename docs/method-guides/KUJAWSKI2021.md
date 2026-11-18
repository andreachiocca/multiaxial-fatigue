# Kujawski 2021 mSWT (`KUJAWSKI2021`)

> **Registry name:** `KUJAWSKI2021` (and `KUJAWSKI2021_ext`) · **Catalogue entry:** [N26](../../METHODS.md#n26) ·
> **Source:** Kujawski (2021), <https://doi.org/10.1111/ffe.13488> ·
> **Implementation:** [`cp_methods/kujawski2021.py`](../../src/multifat/models/cp_methods/kujawski2021.py)

## Overview

- **Type and family:** energy-type (stress × strain) mean-stress function, used here in the
  critical-plane extension of the repository's [SWT](SWT.md).
- **Idea:** the SWT product $\sigma_{\max}\varepsilon_a$ is read as the total damaging energy
  density $W_D$ in the first quadrant of shifted "damaging" axes. Under a negative mean stress the
  hysteresis loop keeps about one third of the mean after relaxation, so the SWT maximum stress is
  augmented by $|\sigma_m|/3$.
- **Mean stress:** identical to SWT for $\sigma_m \ge 0$; more damaging than SWT for
  $\sigma_m < 0$, where SWT is nonconservative ($R < -1$). Compression–compression cycles can be
  damaging.
- **Proposed for:** uniaxial R-ratio effects in ductile iron and 7075-T651 aluminium alloy.
- **Status:** implemented for the declared repository scope. Eqs. (9)–(13b) were checked on the
  paper. The paper validates the function on figures only, so the tests are closed forms.

## Original formulation

Kujawski (2021, Eqs. 11a–11b):

```math
W_D = \sigma_{\max}\,\varepsilon_a \quad (\sigma_m \ge 0),
\qquad
W_D = \left(\sigma_{\max} + \frac{|\sigma_m|}{3}\right)\varepsilon_a \quad (\sigma_m < 0),
```

where:

- $\sigma_{\max}$, $\sigma_m$ are the maximum and mean stress of the cycle;
- $\varepsilon_a$ is the strain amplitude;
- $|\sigma_m|/3$ is the stabilised part of a negative mean stress (Fig. 7, strain ratio
  $R_\varepsilon = -2$).

**Stress-based forms (Eqs. 12a–13b).** The equivalent fully reversed stress amplitude is
$\sigma_{ar} = \sqrt{\sigma_{\max}\sigma_a}$ for $\sigma_m \ge 0$ and
$\sigma_{ar} = \sqrt{(\sigma_{\max} + |\sigma_m|/3)\,\sigma_a}$ for $\sigma_m < 0$; with the strain
amplitude, $\sigma_{ar} = \sqrt{E\varepsilon_a(\sigma_{\max} + |\sigma_m|/3)}$. Life follows from
the fully reversed curve, $\sigma_{ar}$ or $W_D$ against $2N_f$.

The paper is uniaxial. It only says that the interpretation could show how the function "can be
extended to more complex and multiaxial loading situations".

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\mathrm{DP} = W_D = \max\!\left(\sigma_{n,\max} + \frac{|\sigma_{n,m}|}{3}\,[\sigma_{n,m} < 0],\ 0\right)\varepsilon_{n,a},
\qquad \varepsilon_{n,a} = \tfrac12|\Delta\varepsilon_{zz}|.
```

DP is in MPa (stress × strain).

- `KUJAWSKI2021`: DP on the plane of largest normal strain range, the plane rule of `SWT`; ties
  resolved towards the larger DP.
- `KUJAWSKI2021_ext`: largest DP over all planes (repository variant).

**Repository choices**

- *Multiaxial extension (maintainer decision).* The source is uniaxial. The function is applied
  on the planes of the repository's SWT critical-plane extension, with the normal stress and
  strain of each plane. Without a compressive normal mean on the plane, the result is identical
  to `SWT`.
- *Clip.* The bracket $\sigma_{\max} + |\sigma_m|/3 = \sigma_a - 2|\sigma_m|/3$ is clipped at
  zero, as SWT clips a fully compressive cycle; DP = 0 when $|\sigma_m| \ge 1.5\,\sigma_a$.
- *Mean stress.* The applied mean is used. The factor 1/3 already stands for its stabilised
  relaxation, and multifat strains are linear-elastic.
- *Life.* DP is $W_D$, not $\sigma_{ar}$ (a monotone transform); the shared master curve
  $\mathrm{DP} = A N_f^b$ replaces the reference curve
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)).

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook | yes (strains) |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

1. **Original procedure.** Nothing is fitted in the function. The life curve is the fully
   reversed ($R = -1$) $\sigma_{ar}$– or $W_D$–life curve of the material.
2. **Repository.** The same: only the master-curve constants are fitted, on the uniaxial
   $R = -1$ design tests, where $W_D$ equals the SWT value.

## Worked example

7075-T651 elastic constants; uniaxial cycles with positive, zero and negative mean stress.

```python
import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

mswt, swt = get_models(["KUJAWSKI2021", "SWT"])
planes = build_R_list(n1=20, n2=40)
E, nu = 71700.0, 0.33  # 7075-T651 (the paper's Figs. 8, 10, 11)


def eps(s):
    """Linear-elastic strain tensor."""
    return ((1 + nu) * s - nu * np.trace(s) * np.eye(3)) / E


def dp(model, sigma_a, sigma_m, tau_a=0.0):
    """DP for sigma_x = sigma_m + sigma_a sin wt, tau_xy = tau_a sin wt."""
    S_mean = np.diag([sigma_m, 0.0, 0.0])
    S_sin = np.array([[sigma_a, tau_a, 0.0], [tau_a, 0.0, 0.0], [0.0, 0.0, 0.0]])
    zero = np.zeros((3, 3))
    result = model.evaluate_case(
        S0=S_mean + S_sin, S1=S_mean - S_sin, E0=eps(S_mean + S_sin), E1=eps(S_mean - S_sin),
        params={"E": E, "nu": nu}, R_list=planes,
        harmonics={"S_mean": S_mean, "S_sin": S_sin, "S_cos": zero,
                   "E_sin": eps(S_sin), "E_cos": zero},
    )  # fmt: skip
    return result.values[""]


for sigma_a, sigma_m in [(200.0, 100.0), (200.0, 0.0), (200.0, -90.0), (150.0, -180.0)]:
    w_k, w_s = dp(mswt, sigma_a, sigma_m), dp(swt, sigma_a, sigma_m)
    print(
        f"sigma_a = {sigma_a:5.1f}, sigma_m = {sigma_m:6.1f}: W_D = {w_k:.4f} MPa "
        f"(SWT {w_s:.4f}), sigma_ar = sqrt(E W_D) = {np.sqrt(E * w_k):5.1f} MPa"
    )
print(f"in-phase, sigma_m = 100: mSWT {dp(mswt, 150, 100, 90):.4f} = SWT {dp(swt, 150, 100, 90):.4f}")
```

Output:

```text
sigma_a = 200.0, sigma_m =  100.0: W_D = 0.8368 MPa (SWT 0.8368), sigma_ar = sqrt(E W_D) = 244.9 MPa
sigma_a = 200.0, sigma_m =    0.0: W_D = 0.5579 MPa (SWT 0.5579), sigma_ar = sqrt(E W_D) = 200.0 MPa
sigma_a = 200.0, sigma_m =  -90.0: W_D = 0.3905 MPa (SWT 0.3068), sigma_ar = sqrt(E W_D) = 167.3 MPa
sigma_a = 150.0, sigma_m = -180.0: W_D = 0.0628 MPa (SWT 0.0000), sigma_ar = sqrt(E W_D) =  67.1 MPa
in-phase, sigma_m = 100: mSWT 0.7879 = SWT 0.7879
```

For $\sigma_m = -90$ MPa, $\sigma_{ar} = \sqrt{(110 + 30)\cdot 200} = 167.3$ MPa (Eq. 12b), against
$\sqrt{110\cdot 200} = 148.3$ MPa for SWT. The compression–compression cycle
($\sigma_{\max} = -30$ MPa) is damaging, as in the paper's Fig. 10.

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models KUJAWSKI2021 --all-variants --results-dir Results/current/guides
```

None of the supplied workbooks has a compressive normal mean on the critical plane, so on all of
them `KUJAWSKI2021` gives exactly the `SWT` damage parameters and lives.

## Comparison capabilities and outputs

- **Life:** master curve on $W_D$ (`Nf_expected`).
- **DP comparison:** supported. $W_D$ is homogeneous of degree 2 in the load; use `Error_load`
  across methods.
- **`_ext`:** available.

## Verification and limitations

[`tests/methods/test_kujawski2021.py`](../../tests/methods/test_kujawski2021.py):

- the stress factor of Eqs. (11a)–(11b), both branches and the compression–compression range;
- Eqs. (12a)–(13b): $\sqrt{E W_D} = \sigma_{ar}$ in uniaxial elastic loading;
- the ratio to SWT under a negative mean;
- identity with `SWT`, case and plane level, without a compressive normal mean.

Not reproduced: Figs. 8–11. They plot data of other authors and print no table.

Limitations:

- uniaxial in the source; the multiaxial plane rule is a repository choice;
- the 1/3 relaxation was observed in strain-controlled tests; it is applied here to the nominal
  mean stress;
- linear-elastic strains.

## References

- D. Kujawski (2021), An interpretation and modification of the SWT function, *Fatigue & Fracture
  of Engineering Materials & Structures* 44(8), 2202–2209. <https://doi.org/10.1111/ffe.13488>
- K. N. Smith, P. Watson, T. H. Topper (1970), A stress–strain function for the fatigue of metals,
  *Journal of Materials* 5(4), 767–778.
