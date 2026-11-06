# Notation

Every concept carries one name in the package, and that name is the symbol of the source paper
where a function implements one equation. This page is the single table behind that rule; code
reviews and the `tests/test_notation.py` guard both read it. Arguments that are interface
machinery rather than concepts (a callback, a plotting option, a search accumulator) are recorded
in `tests/data/notation_allow.txt` with the reason they are not vocabulary; that file holds
exactly the names outside this table.

A load cycle is a single-harmonic history

    S(t) = S_mean + S_sin sin(w t) + S_cos cos(w t)
    E(t) = E_mean + E_sin sin(w t) + E_cos cos(w t)

where `S` is the stress tensor and `E` the strain tensor. `E` on its own is Young's modulus; every
strain tensor carries the `E_` prefix, which matches the keys of the importer `H` dictionary.
Critical-plane methods instead receive the two load-step tensors `S0, S1, E0, E1` (the public
interface), from which `S_mean` and `S_sin` follow when no harmonic dictionary is supplied.

`R_list` is a list of plane frames, not a list of stress ratios. Each frame is a rotation matrix
whose third column is the plane normal; `R` is one such frame and `n` / `normals` are the plane
normal(s) themselves.

| Concept | Paper symbol | Python name | Workbook key | Unit |
|---|---|---|---|---|
| Stress tensor, harmonic mean / sine / cosine | M, A, B | `S_mean`, `S_sin`, `S_cos` | imported | MPa |
| Strain tensor, harmonic mean / sine / cosine | - | `E_mean`, `E_sin`, `E_cos` | imported | - |
| Stress tensor | S | `S` | - | MPa |
| Load-step tensors | - | `S0`, `S1`, `E0`, `E1` | - | MPa, - |
| Load-step tensors rotated onto one candidate plane | - | `S0r`, `S1r`, `E0r`, `E1r` | - | MPa, - |
| Harmonic dictionary of the importer (the `H` keys) | - | `harmonics` | imported | - |
| Plane frames, all planes / one plane | - | `R_list`, `R` | - | - |
| Plane set of the search (frames plus the rotated tensors) | - | `planes` | - | - |
| Plane request: maximum shear stress of the history (the tau_max rule) | tau_max | `tau_max` | - | MPa |
| Plane request: plane strain ranges | - | `strain_ranges` | - | - |
| Plane normal, one / many | n | `n`, `normals` | - | - |
| Axial reversed limit (amplitude) | f_-1, s_-1 | `f_m1` | Sigm1 | MPa |
| Torsional reversed limit (amplitude) | t_-1 | `t_m1` | Taum1 | MPa |
| Axial repeated limit (maximum stress) | f_0, s_0 | `f_0` | Sig0 | MPa |
| Torsional repeated limit (maximum stress) | t_0 | `t_0` | Tau0 | MPa |
| Ultimate tensile strength | sigma_f, S_u | `su` | Su | MPa |
| Ultimate torsional strength | tau_b, T_u | `tu` | Tu | MPa |
| Normal stress amplitude / mean on a plane | sigma_a, sigma_m | `sigma_a`, `sigma_m` | - | MPa |
| Shear stress amplitude / mean on a plane | tau_a, tau_m | `tau_a`, `tau_m` | - | MPa |
| Hydrostatic amplitude / mean / maximum | sigma_H,a, sigma_H,m | `sigma_H_a`, `sigma_H_m`, `sigma_H_max` | - | MPa |
| Largest normal stress on a plane | sigma_n,max | `sigma_n_max` | - | MPa |
| Largest normal stress over the two load steps (SWT, MSWT, KUJAWSKI2021) | sigma_max | `sigma_max` | - | MPa |
| Normal stress on a plane at load step 0 / 1 | - | `sigma_z0`, `sigma_z1` | - | MPa |
| Stress components on the free surface | sigma_x, sigma_y, sigma_z, tau_xy | `sigma_x`, `sigma_y`, `sigma_z`, `tau_xy` | - | MPa |
| Mean-stress sensitivity (Dubbel) | M | `M` | - | - |
| MMMH mean and alternating part, Eq. (4) | M_M, R_M | `M_M`, `R_M` | - | MPa |
| Normal / shear stress range on a plane | Delta sigma_n, Delta tau | `d_sigma_n`, `d_tau` | - | MPa |
| Normal / shear strain range on a plane | Delta eps, Delta gamma | `d_eps`, `d_gam` | - | - |
| Shear range components on a plane | - | `d_tau_xz`, `d_tau_yz`, `d_gam_xz`, `d_gam_yz` | - | MPa, - |
| Amplitude of sqrt(J2) | sqrt(J2,a) | `sqrt_J2_a` | - | MPa |
| Young's modulus | E | `E` | E_MPa | MPa |
| Poisson's ratio, elastic / plastic | nu, nu_e, nu_p | `nu`, `nu_e`, `nu_p` | nu | - |
| Weibull exponent / scale stress | m, s_0 | `m`, `s0` | - | -, MPa |
| Damage parameter, selection metric | DP | `DP`, `metric` | - | model |
| Per-plane damage and metric arrays | - | `damage`, `metric` | - | model |
| Number of cycles to failure | N_f | `Nf` | Nf_cycles | cycles |
| Survival probability of a fitted master curve | P_s | `P_s` | - | - |
| Master-curve damage parameter at the observed life | DP_e | `DP_fit` | - | model |
| Weights of tau_a^2, sigma_a^2, tau_a tau_m, sigma_a sigma_m in BP / CAIM | a_BP..d_BP, a_CAIM..d_CAIM | `a_BP`, `b_BP`, `c_BP`, `d_BP`, `a_CAIM`, `b_CAIM`, `c_CAIM`, `d_CAIM` | - | - |
| Coefficients of the stress scale factor surface, Anes et al. (2014) Eq. (31) | a..d, f..i | `a_SSF`, `b_SSF`, `c_SSF`, `d_SSF`, `f_SSF`, `g_SSF`, `h_SSF`, `i_SSF` | a_SSF..i_SSF | MPa^-k |
| Tension/torsion fatigue-limit ratio and its square | kappa, kappa^2 | `kappa`, `kappa_sq` | - | - |
| History of the stress / strain tensor of FE locations over ordered load cases, `(M, N, 3, 3)` | - | `S_hist`, `E_hist` | - | MPa, - |
| Stress / strain range of a history under a range strategy | Delta sigma, Delta eps | `S_range`, `E_range` | - | MPa, - |
| Stable source identity of an FE location (node, element or integration point) | - | `location_id` | - | - |
| Pair of load-case indices a range strategy selects | a, b | `pair` | - | - |
| Axis of a time history and what it measures (`time` or `cycle`) | - | `time`, `axis_kind` | - | s or cycles |
| Elastic / plastic part of the strain history of FE locations, `(M, N, 3, 3)` | eps^e, eps^p | `E_elastic_hist`, `E_plastic_hist` | - | - |
| Temperature history of FE locations, `(M, N)` | T | `T_hist` | - | solver unit |
| Instance (part) and material label of an FE location, with the node label a location is the triple (instance, ID, material) | - | `instance`, `material` | - | - |
| Load function of a unit load case, `S(t) = sum_k f_k(t) S_k` | f_k | `load_functions` | - | - |
| Solver step number of a load case (need not be contiguous) | - | `step` | - | - |
| Kind of FE analysis (`linear_elastic`, `elastic_plastic`) and role of its steps (`time_sample`, `unit_case`) | - | `analysis_type`, `step_role` | - | - |
| Per-case inputs of a native life law | - | `aux` | - | model |
| Fitted constants of a native life law | - | `fit` | - | model |

## Rules

1. A function that implements one published equation uses that equation's symbols, including the
   subscript that the paper prints (`f_m1`, `t_m1`, `sigma_H_a`, `sqrt_J2_a`).
2. A paper symbol that clashes with another concept in the package takes a suffix, never a
   synonym: `nu_e` / `nu_p` for the elastic and plastic Poisson ratios, `S_sin_dev` for the
   deviatoric part of a harmonic coefficient.
3. No single-letter name that the paper does not use. Loop indices (`i`, `k`) and matrix axes
   (`p`, `q`) are the exception.
4. Layout is not part of the vocabulary: black owns spacing, line length and trailing commas.
5. Python source, its docstrings and the tables generated from it spell the name `Boehme`, and
   the guard rejects `Böhme` there. Hand-written documentation keeps the author's spelling
   `Böhme`: this includes the method catalogue (`docs/catalogue/methods.json`, the source of the
   generated `METHODS.md` and the per-family pages), `docs/bibliography.md` and the README.

## Retired synonyms

These names do not appear in the package any more. The guard in `tests/test_notation.py` fails if
one of them is reintroduced as a local name or parameter.

| Retired | Use instead |
|---|---|
| `mean_part`, `avg_part`, `mean_tensor`, `mean_stress` | `S_mean` |
| `sin_part`, `sine_part`, `sin_tensor`, `sine_tensor`, `sin_amp` | `S_sin` |
| `cos_part`, `cosine_part`, `cos_tensor`, `cosine_tensor`, `cos_amp` | `S_cos` |
| `h_amp`, `h_avg`, `h_mean`, `hydro_amp`, `hydro_mean`, `mean_hydrostatic` | `sigma_H_a`, `sigma_H_m` |
| `sig_m1`, `s_m1`, `axial_limit`, `sigma_minus1`, `tension_limit`, `f_lim`, `sig_w` | `f_m1` |
| `torsion_limit`, `tau_minus1`, `t_lim` | `t_m1` |
| `sig0`, `sig1`, `sig_max` | `sigma_z0`, `sigma_z1`, `sigma_max` |
| `sig_n_max` | `sigma_n_max` |
| `dtau_a`, `dtau_b`, `delta_tau`, `DeltaTau` | `d_tau_xz`, `d_tau_yz`, `d_tau` |
| `dgam_a`, `dgam_b`, `delta_gamma` | `d_gam_xz`, `d_gam_yz`, `d_gam` |
| `young`, `youngs`, `young_modulus` | `E` |
| `poisson`, `poisson_e` | `nu`, `nu_e` |
| `dirs`, `directions`, `plane_normals` | `normals` |
| `dmg`, `sel_dmg` | `damage` |
| `met`, `mtr` | `metric` |
| `sel_damage`, `sel_metric`, `sel_normals` | `selection_damage`, `selection_metric`, `selection_normals` |
| `fs_value`, `swt_value`, `cs_max` | `DP` |
| `exponent` (Weibull), `scale_stress` | `m`, `s0` |
| `shape_beta` | `beta` |
| `tensile_strength`, `rot_bending_limit` | `s_ten`, `s_rb` |
| `amp_vals`, `life` (function parameter) | `amplitude`, `Nf` |
| `evaluate_fn` | `evaluate_normals` |
| `mat_name` | `material` |
| `spec`, `src` | `case` |
| `goal` | `target` |
| `sens` | `M` |
| `ratio`, `k2`, `k_ratio` | `kappa`, `kappa_sq` |
| `sig_x`, `sig_y`, `sig_z` | `sigma_x`, `sigma_y`, `sigma_z` |
| `elastic_strain`, `plastic_strain` | `eps_e`, `eps_p` |
| `algorithm_key` | `selection` |
| `tolerance` | `tol` |
| `record_id` | `case_id` |
| `book`, `mat`, `dataset`, `consts` (CLI layer) | `workbook`, `material`, `data`, `params` |

The same rule covers private function aliases: the package has no `_old = _new` lines, so every
helper has exactly one name and the tests import that name.

## FatPy nomenclature

`multifat` is expected to merge into FatPy, the Python
package of COST Action CA23109 (FABER). FatPy's draft nomenclature (`nomenclature_v01`, 2025)
names material and experiment quantities in snake case, built from a shared abbreviation list
(`fat`, `lim`, `ten`, `tor`, `ben`, `coef`, `exp`, `amp`, `eq`, …). The table records the FatPy
key of each multifat name that has one, so that results and workbooks can be translated at the
merge. **The Python names above stay the reference inside multifat**: this is a correspondence,
not a rename.

Units differ for one quantity. multifat keeps every stress and modulus in MPa, while FatPy v01
gives Young's modulus in GPa (but the shear modulus in MPa): $E_\text{FatPy} = E/1000$.

| multifat name | FatPy key (v01) | FatPy unit | Note |
|---|---|---|---|
| `E` | `elastic_modulus` | GPa | multifat: MPa |
| `nu` | `poissons_ratio` | - | |
| `nu_p` | `poissons_ratio_pl` | - | |
| `Sy` | `yield_strength` | MPa | |
| `Su` | `ult_tensile_strength` | MPa | |
| `Tu` | - | - | no FatPy key yet (torsional ultimate strength) |
| `Sigm1` | `fat_lim_ten_m1` or `fat_lim_ben_m1` | MPa | push-pull or bending, per the Summary column `Axial fatigue limit loading` (see `docs/workbook-format.md`); unrecorded sources stay unresolved |
| `Taum1` | `fat_lim_tor_m1` | MPa | |
| `Sig0` | `fat_lim_ten_0` | MPa | multifat stores the maximum stress of the R = 0 cycle; FatPy v01 does not state amplitude or maximum |
| `Tau0` | `fat_lim_tor_0` | MPa | as `Sig0` |
| `sigma_f`, `b_fat` | `fat_strength_coef`, `fat_strength_exp` | MPa, - | |
| `epsilon_f`, `c_fat` | `fat_ductility_coef`, `fat_ductility_exp` | - | |
| `tau_f`, `b0` | `fat_strength_coef_tor`, `fat_strength_exp_tor` | MPa, - | |
| `gamma_f`, `c0` | `fat_ductility_coef_tor`, `fat_ductility_exp_tor` | - | |
| `K_prime`, `n_prime` | `cyc_hard_coeff`, `cyc_hard_ext` | MPa, - | FatPy key spellings as in v01 |
| `Nf` | `tot_fat_life` | cycles | |
| `sigma_a`, `sigma_m` | `stress_amp`, `mean_stress` | MPa | load-level amplitude and mean (multifat also uses them on a plane) |
| `phase_shift` | `phase_shift` | degrees | |
| `k_FS`, `k_FI`, `L_LI`, `k_MGSE`, `a_MSWT`, … | - | - | method constants: no FatPy key yet |

FatPy v01 has no key for a run-out (censored observation) flag, for the torsional ultimate
strength, or for method constants. New multifat names should be built from the FatPy
abbreviation list when no paper symbol applies.
