# Feedback on the FatPy nomenclature v01

Notes from mapping the multifat parameters onto `nomenclature_v01.xlsx` (sheets *Material
parameters*, *Experiment parameters*, *Abbreviations*, *checkup*). multifat
(`andreachiocca/multiaxial-fatigue`, FABER WG4.3) is expected to merge into FatPy. The resulting
correspondence is the *FatPy nomenclature* table of [`docs/notation.md`](../notation.md). Most
multifat material quantities have a FatPy key. The points below are what a multiaxial
fatigue-assessment code needs, or found ambiguous.

## 1. Quantities without a key

| Quantity | Why it is needed | Suggested key |
| --- | --- | --- |
| Run-out flag (censored test) | A run-out is a lower bound on life, not a failure. Fits and error statistics must treat it differently. | `runout` (bool) or `test_status` (`failure` / `runout`) |
| Torsional ultimate strength $\tau_b$ | Used by mean-shear corrections (e.g. Gerber in Liu, Xie & Hao 2025) | `ult_shear_strength` |
| Reference life of a fatigue limit | A limit is defined at a number of cycles (e.g. $2\cdot10^6$). Run-outs stopped earlier are not limit states. | `fat_lim_ref_life` |
| Per-component load description | Multiaxial tests need, per component (axial, torsional, radial): amplitude, mean and phase. `stress_amp`, `mean_stress` and `phase_shift` are scalars. | e.g. `stress_amp_ax`, `stress_amp_tor`, `phase_shift_tor` |
| Strain amplitudes | Strain-controlled and strain-based criteria | `strain_amp`, `shear_strain_amp` |
| Method constants | Criteria carry fitted constants (Fatemi–Socie $k$, Findley $k$, …) | a convention such as `<method>_<symbol>`, e.g. `fs_k`, `findley_k` |
| Source of each value | Several sources often report one property. multifat keeps every source row with its DOI. | `reference` per value |

## 2. Ambiguous definitions

- **`fat_lim_ten_0`, `fat_lim_tor_0`.** They don't say whether the repeated (R = 0) limit is
  the amplitude or the maximum stress; the two differ by a factor of 2. multifat stores the
  maximum.
- **`fat_lim_ten_m1` versus `fat_lim_ben_m1`.** The split is useful, and multifat now records it
  per source. A bending limit is usually higher than the push-pull limit because of the stress
  gradient. The note could say which one a criterion's calibration expects.
- **`biaxiality_ratio` (λ).** Is it $\sigma_2/\sigma_1$, or the amplitude ratio
  $\sigma_a/\tau_a$ used in tension–torsion papers? Both are in use.
- **`stress_ratio` (R).** For multiaxial loading it is defined per component; a single R is
  ambiguous.

## 3. Units

- Young's modulus is in **GPa** while the torsion elastic modulus and every stress are in
  **MPa**. A single unit for stresses and moduli would avoid factor-1000 errors; multifat uses
  MPa throughout.
- `hardness` lists the unit as `HV, HB, HRC`. The scale would need its own field.

## 4. Spelling and form of the keys

- **Typo.** "coeffecient" in the *Abbreviations* and *checkup* sheets.
- **`coef` suffix.** The abbreviation list gives `coef`, but keys mix `coef` (`fat_strength_coef`)
  and `coeff` (`cyc_hard_coeff`, `bergmann_coeff`).
- **`cyc_hard_ext`** (cyclic hardening exponent) should probably be `cyc_hard_exp`, by analogy
  with `fat_strength_exp`.
- **Torsion symbols.** The LaTeX symbols of the torsion coefficients are written $\sigma_f't$
  and $\varepsilon_f't$. The usual symbols are $\tau_f'$ and $\gamma_f'$.
- **Invalid Python identifiers.** `fat_lim_tor_A_-1` and `fat_lim_tor_B-1` contain `-`. The
  `m1` suffix of `fat_lim_ten_m1` already solves this.
- **Rows without a final key.** About 30 rows, e.g. `C_FIN`, `W_F-1`, `NF_F-1`, `R_AS`,
  `KD_DEFF_S`, `L_CHAR`, `D_K_TH`, have no description or key yet.

## 5. What multifat does in the meantime

- multifat keeps its own Python names. `docs/notation.md` lists the FatPy key of each one, for
  translation at the merge.
- New multifat names are built from the FatPy abbreviation list when no paper symbol applies.
- Workbooks may record how the axial R = -1 fatigue limit was obtained, through the Summary
  column `Axial fatigue limit loading`. This maps it onto `fat_lim_ten_m1` or `fat_lim_ben_m1`.
