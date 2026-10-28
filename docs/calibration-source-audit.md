# Initial primary-source calibration audit for Issue #14

Status: source passages reviewed in the three papers supplied by the maintainer on 28 October 2026. Fatemi-Socie is a separate maintainer-prescribed calibration. These findings establish the intended scientific procedures but **do not claim that the existing production code reproduces them**.

## Fatemi-Socie (FS)

**Calibration specification (maintainer-provided, original paper not supplied):** Determine `k_FS` by minimising fatigue-life scatter using compatible fully reversed axial tension-compression and pure torsional fatigue failures (`R = -1`). Preserve the search objective, selected axial/torsional observation IDs, fit parameters, critical-plane definition, bounds, numerical settings and source/batch. Run-outs do not become failures. Mark the original-paper verification status `maintainer_defined`, not `source_verified`.

**Default plane:** evaluate the FS damage parameter on the plane of maximum shear strain amplitude (`gamma_a`, equivalently the shear strain range for a fixed cycle convention); use `FS_ext` only on explicit request to maximise the *entire* FS damage parameter. This is currently how `FatemiSocie.evaluate_on_plane` and `CPModelAdapter` label their two results. The ordinary call must resolve to the first variant without requiring the user to type a suffix.

## Findley (FIN), Findley (1959)

Source: *A Theory for the Effect of Mean Stress on Fatigue of Metals Under Combined Torsion and Axial Load or Bending*, DOI [10.1115/1.4008327](https://doi.org/10.1115/1.4008327).

- Page 1, Eqs. (1)-(2), introduce the alternating shear and maximum normal stress relation with constants `f` and `k` defined **for a given life**. This does not amount to fitting a generic damage-versus-life scatter objective.
- Page 2, Eq. (13), gives the critical-plane direction for a fully reversed cycle through `tan(2 theta) = 1/k`, and Eq. (18) describes combined proportional loading with superimposed static stresses.
- Page 4, Eq. (30), relates the reversed axial (bending) and torsional fatigue strengths `b/t` at the **same reference life** to `k`. The source-derived `author` calibration should use this strength ratio with consistent amplitude definitions and physically admissible `k`; equivalence between bending and axial test conditions must be reviewed before treating them as interchangeable.
- The paper explicitly restricts its derivation to proportional/synchronous combined bending or axial load and torsion with mean stresses in corresponding proportions. This is **not** an independent validation of arbitrary non-proportional 3D paths.

**Default plane:** maximise the entire `tau_a + k_FI sigma_n,max` parameter. The current `Findley` class already sets its selection metric equal to the full damage parameter. Hence the plain `FIN` and explicit `FIN_ext` are numerically equivalent under the existing implementation, although the source and variant selections should be recorded.

The repository's current `k_FI` scatter minimisation is an **alternative repository fit**, not yet the author-defined default. If fatigue limits at matching reference life are absent, state that this source calibration is unavailable rather than silently switching to a generic fit.

## Li et al. (2021), LI

Source: *Multiaxial fatigue life prediction for metals by means of an improved strain energy density-based critical plane criterion*, DOI [10.1016/j.euromechsol.2021.104353](https://doi.org/10.1016/j.euromechsol.2021.104353).

- Page 4, Eq. (3), defines the shear and normal strain-energy contributions and the weight `L`. Critical-plane selection is described as the maximum shear strain energy density plane with the greater normal strain energy density in case of equivalent maxima.
- Page 5, Eqs. (20)-(25), construct the torsion-based fatigue-life relation from the shear strain-life and stress-life coefficients `tau_f`, `gamma_f`, `b0` and `c0`.
- Page 6, Section 2.4, Eqs. (31)-(32), explicitly calibrates `L`: (1) estimate the four torsional coefficients using least squares on fully reversed torsion tests; (2) build the uniaxial fitting quantity `P_F(L)` from the fully reversed axial tests; (3) evaluate `L` from **0 to 10 in steps of 0.01**, choosing the minimum.

The current `L_LI` generic scatter minimisation over combined axial and torsional failures is therefore **not** the prescribed method. The paper uses a source-defined shear-based right-hand life curve, not an arbitrary subsequent damage-parameter master curve. The default requires suitable torsional and axial fatigue data with valid strain/stress components. Record separate torsional curve and axial `L` calibration observation IDs.

The implemented LI selection metric is the shear contribution, but the exact tie-break rule and stress/strain range convention require equation-level reproduction tests before source verification is complete.

## Zhu et al. (2018), MGSE_ZHU and MSWT

Source: *Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials*, DOI [10.1016/j.ijfatigue.2018.03.028](https://doi.org/10.1016/j.ijfatigue.2018.03.028).

### MGSE

- Pages 7-8, Eqs. (7)-(12), introduce `k_MGSE` and derive a relation using reversed axial fatigue properties; the source permits fitting the uniaxial fatigue data or averaging values calculated at designated lives.
- Pages 12-13, Section 4.1, give the practical reference-life selection. For the metals analysed, the authors use **5,000-50,000 cycles** to calculate the average `k_MGSE`. A lower-life interval is used for sintered porous iron. The choice is material-/regime-dependent and should be recorded explicitly.
- The manuscript suggests `k_MGSE = 1` for ductile materials when data or properties are unavailable, but this is an **empirical fallback**, not a data-derived calibrated constant.

A source-defined implementation must state whether it follows the equation-derived, life-interval average or direct-fitting route. The current fixed life interval and numerical procedure are repository settings that require comparison against the source's actual equations, including elasticity versus plasticity assumptions.

### Modified SWT (MSWT)

- Page 8, Eqs. (13)-(15), define the weight `a` and model form; the critical plane maximises the full MSWT damage parameter. Cracking-mode ranges constrain interpretation and the paper does not provide one universally identifiable algebraic formula for `a`.
- Page 14, Section 4.1, recommends identifying `a` by bringing the axial and torsional damage-parameter-versus-life curves as close as possible. The **exact numerical mismatch objective and weighting must be declared as an implementation choice** unless separately confirmed with an original source.
- Using a generic life-scatter objective is one possible `all_eligible` optimisation, but it must not be presented as the manuscript's default procedure without verification.

## Source-verified versus implementation-verified

The primary-source passages above have been checked, but numerical reproduction of their formulas and fitting procedures remains to be done. Do not mark an implementation `source_verified` until equations, data selection, parameter bounds, plane selection and at least one independent reference example agree.

Where two sources prescribe different calibrations for a nominally similar method, record the exact implemented formulation, article, equation and procedure identity instead of assuming a universal rule.

## Immediate follow-up tests

1. Verify that `FS` returns the damage on the maximum `gamma_a` plane and `FS_ext` returns the largest damage across all planes; confirm cases where these differ.
2. Verify that plain `FIN` uses maximisation of `tau_a + k_FI sigma_n,max`, and that requesting `FIN_ext` does not silently alter this result.
3. Reproduce Findley Eq. (30), Li Eqs. (31)-(32), Zhu Eq. (12) and the selected MGSE interval using controlled synthetic or published examples.
4. Verify that `author` mode never silently falls back to the generic scatter objective, and that all fit alternatives remain explicitly labelled with their selected observations and independent validation sets.
5. Verify that every new calibration record is immutable and retained, including default and user-forced calibrations.
