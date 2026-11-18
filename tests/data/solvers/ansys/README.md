# In-House Ansys Biaxial Fatigue Benchmark: Holed Plate

This dataset was generated directly with Ansys MAPDL 2026 R1 running in batch mode on WorkstationUniTN using `hole_plate_fatigue.mac`. It is entirely self-generated and contains no third-party or imported external repository files.

## Model Description
* **Physics**: 2D plane stress ($t = 1.0\,\mathrm{mm}$) with linear elasticity ($E = 200\,\mathrm{GPa}$, $\nu = 0.3$).
* **Element Type**: `PLANE183` (8-node quadratic 2D solid element) with `KEYOPT(3)=3` (plane stress with thickness).
* **Geometry**: Quarter-symmetric plate ($X \ge 0, Y \ge 0$), width $W = 100\,\mathrm{mm}$, height $L = 100\,\mathrm{mm}$, circular notch radius $a = 10\,\mathrm{mm}$ ($d/D = 2a/2W = 0.1$).
* **Boundary Conditions**: Symmetry on $Y = 0$ ($U_y = 0$) and $X = 0$ ($U_x = 0$).

## 2-Step Biaxial Fatigue Cycle
* **Load Step 1 ($t = 1.0\,\mathrm{s}$)**: Uniaxial tension $\sigma_0 = 100\,\mathrm{MPa}$ applied to the top edge ($Y = L$).
  - Notch A $(a, 0)$: Tensile stress concentration $\sigma_y \approx +304.5\,\mathrm{MPa}$, $\sigma_x \approx +5.8\,\mathrm{MPa}$.
  - Notch B $(0, a)$: Compressive notch root $\sigma_x \approx -103.2\,\mathrm{MPa}$, $\sigma_y \approx -9\,\mathrm{MPa}$.
* **Load Step 2 ($t = 2.0\,\mathrm{s}$)**: Uniaxial tension $\sigma_0 = 100\,\mathrm{MPa}$ applied to the right edge ($X = W$), top edge unloaded.
  - Notch A $(a, 0)$: Compressive notch root $\sigma_y \approx -103.2\,\mathrm{MPa}$, $\sigma_x \approx -9\,\mathrm{MPa}$.
  - Notch B $(0, a)$: Tensile stress concentration $\sigma_x \approx +304.5\,\mathrm{MPa}$, $\sigma_y \approx +5.8\,\mathrm{MPa}$.

The principal stress directions rotate by $90^\circ$ between load steps, producing an authentic multiaxial, non-proportional alternating fatigue cycle.

## Analytical Verification
According to Howland (1930) and Peterson's Stress Concentration Factors for a finite-width plate with a central circular hole ($d/D = 0.1$):
$$K_{tn} = 2 + (1 - d/D)^3 = 2.729 \implies K_{tg} = \frac{K_{tn}}{1 - d/D} \approx 3.032 \implies \sigma_{\mathrm{peak}} \approx 303.2\,\mathrm{MPa}$$
The refined FE solution peak is $\sigma_{\mathrm{peak}} = 304.48\,\mathrm{MPa}$, agreeing with Howland's closed-form baseline to within $+0.41\%$.

## Files
* `hole_plate_fatigue.mac`: Complete APDL source script (64 hole divisions, 1.5 mm ligament element sizing).
* `hole_plate_fatigue.out`: MAPDL listing containing global Cartesian `NLIST` coordinates and multi-step `PRNSOL, S, COMP` stress tables.
* `hole_plate_fatigue.rst`: MAPDL binary result file (2,864 corner nodes with stored stress).
* `hole_plate_fatigue.csv`: Base exchange CSV (`multifat-fe-csv 1.0`, 638 kB).
