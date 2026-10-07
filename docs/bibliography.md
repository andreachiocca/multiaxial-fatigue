# Bibliography and source mapping

The method names below match the registry. `_ext` variants inherit the base source but change plane selection; see [the audit](method-audit.md). A citation associates a formulation with its source; it does not imply that full text or all validation data were available.

## FS

Fatemi & Socie (1988), A critical plane approach to multiaxial fatigue damage including out-of-phase loading.

Original article: https://doi.org/10.1111/j.1460-2695.1988.tb01169.x

## FIN

Findley (1959), A theory for the effect of mean stress on fatigue of metals under combined torsion and axial load or bending.

Original article: https://doi.org/10.1115/1.4008327

## SWT

Smith, Watson & Topper (1970), A stress-strain function for the fatigue of metals, Journal of Materials 5(4), 767–778.

No DOI is asserted; the full bibliographic citation above identifies the original publication.

## MATAKE

Matake (1977), An explanation on fatigue limit under combined stress.

Original article: https://doi.org/10.1299/jsme1958.20.257

## MCD

McDiarmid (1991), A general criterion for high cycle multiaxial fatigue failure.

Original article: https://doi.org/10.1111/j.1460-2695.1991.tb00673.x

## MGSE_YU

Yu, Zhu, Liu & Liu (2017), Multiaxial fatigue damage parameter and life prediction without any additional material constants, Eq. (11).

Original article: https://doi.org/10.3390/ma10080923

## MGSE_ZHU

Zhu et al. (2018), Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials.

Original article: https://doi.org/10.1016/j.ijfatigue.2018.03.028

## MSWT

Modified SWT form reported by Zhu et al. (2018), Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials.

Original article: https://doi.org/10.1016/j.ijfatigue.2018.03.028

## MKBM

Li et al. (2011), Multiaxial fatigue life prediction for various metallic materials based on the critical plane approach.

Original article: https://doi.org/10.1016/j.ijfatigue.2010.07.003

## LIU2021

Liu, Ran, Wei & Zhang (2021), A critical plane-based multiaxial fatigue life prediction method considering the material sensitivity and the shear stress.

Original article: https://doi.org/10.1016/j.ijpvp.2021.104532

## LI

Li et al. (2021), Multiaxial fatigue life prediction for metals by means of an improved strain energy density-based critical plane criterion.

Original article: https://doi.org/10.1016/j.euromechsol.2021.104353

## GSE

Ince & Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life prediction under proportional and non-proportional loadings.

Original article: https://doi.org/10.1016/j.ijfatigue.2013.10.007

## GSA

Ince & Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life prediction under proportional and non-proportional loadings.

Original article: https://doi.org/10.1016/j.ijfatigue.2013.10.007

## BP

Böhme & Papuga (2024), Advancements in stress-based multiaxial fatigue prediction: A data-driven approach and a new criterion.

Original article: https://doi.org/10.1111/ffe.14281

## CAIM

Böhme, Papuga & Lange (2026), Decoupling amplitude and mean stress effects in multiaxial fatigue: A critical–integral approach.

Original article: https://doi.org/10.1111/ffe.70244

## CARSPA

Carpinteri & Spagnoli (2001), Multiaxial high-cycle fatigue criterion for hard metals.

Original article: https://doi.org/10.1016/S0142-1123(00)00075-X

## DANGVAN

Dang Van (1973), Sur la résistance à la fatigue des métaux, Sciences et Technique de l'Armement 47, 647–722.

No DOI is asserted; the full bibliographic citation above identifies the original publication.

## OTT

Tveit et al. (2024), A continuum approach to multiaxial high-cycle fatigue modeling for ductile metallic materials.

Original article: https://doi.org/10.1016/j.rineng.2024.102171

## SWTD

Kujawski (2014), A deviatoric version of the SWT parameter.

Original article: https://doi.org/10.1016/j.ijfatigue.2013.12.002

## ZHU_EDP

Zhu et al. (2019), A novel energy-based equivalent damage parameter for multiaxial fatigue life prediction.

Original article: https://doi.org/10.1016/j.ijfatigue.2018.11.025

## Additional formulation sources

- GSE/GSA equations were cross-checked in Yu et al. (2017), Eqs. (3)–(4), [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5578289/). MGSE_YU is Eq. (11), not that paper's separate final proposed parameter.
- MSWT lineage: Jiang (2000), *A fatigue criterion for general multiaxial loading*, https://doi.org/10.1046/j.1460-2695.2000.00247.x. The current code implements the scalar form recorded in Zhu (2018), not a full incremental implementation of Jiang (2000).
- BP implementation equations and R=0 estimates: Böhme, Papuga & Lange (2026), https://doi.org/10.1111/ffe.70244. The user-supplied MATLAB file `damage_criteria_CAIM.m` applies specifically to CAIM.
- FS maximum-damage extension: Chiocca (2024), *Closed-form solution for the Fatemi-Socie extended critical plane parameter in case of linear elasticity and proportional loading*, https://doi.org/10.1111/ffe.14153. The repository uses a numerical plane scan, not this closed-form solution.
