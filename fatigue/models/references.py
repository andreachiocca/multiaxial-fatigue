"""Bibliography and audit decisions. See docs/method-audit.md for scope."""
IMPLEMENTATION_REVISION = "2026-10-07-audit1"
REFERENCES = {
    "FS": ("Fatemi & Socie (1988), A critical plane approach to multiaxial fatigue damage including out-of-phase loading", "10.1111/j.1460-2695.1988.tb01169.x"),
    "FIN": ("Findley (1959), A theory for the effect of mean stress on fatigue of metals under combined torsion and axial load or bending", "10.1115/1.4008327"),
    "SWT": ("Smith, Watson & Topper (1970), A stress-strain function for the fatigue of metals, Journal of Materials 5(4), 767–778", ""),
    "MATAKE": ("Matake (1977), An explanation on fatigue limit under combined stress", "10.1299/jsme1958.20.257"),
    "MCD": ("McDiarmid (1991), A general criterion for high cycle multiaxial fatigue failure", "10.1111/j.1460-2695.1991.tb00673.x"),
    "MGSE_YU": ("Yu, Zhu, Liu & Liu (2017), Multiaxial fatigue damage parameter and life prediction without any additional material constants, Eq. (11)", "10.3390/ma10080923"),
    "MGSE_ZHU": ("Zhu et al. (2018), Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials", "10.1016/j.ijfatigue.2018.03.028"),
    "MSWT": ("Modified SWT form reported by Zhu et al. (2018), Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials", "10.1016/j.ijfatigue.2018.03.028"),
    "MKBM": ("Li et al. (2011), Multiaxial fatigue life prediction for various metallic materials based on the critical plane approach", "10.1016/j.ijfatigue.2010.07.003"),
    "LIU2021": ("Liu, Ran, Wei & Zhang (2021), A critical plane-based multiaxial fatigue life prediction method considering the material sensitivity and the shear stress", "10.1016/j.ijpvp.2021.104532"),
    "LI": ("Li et al. (2021), Multiaxial fatigue life prediction for metals by means of an improved strain energy density-based critical plane criterion", "10.1016/j.euromechsol.2021.104353"),
    "GSE": ("Ince & Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life prediction under proportional and non-proportional loadings", "10.1016/j.ijfatigue.2013.10.007"),
    "GSA": ("Ince & Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life prediction under proportional and non-proportional loadings", "10.1016/j.ijfatigue.2013.10.007"),
    "BP": ("Böhme & Papuga (2024), Advancements in stress-based multiaxial fatigue prediction: A data-driven approach and a new criterion", "10.1111/ffe.14281"),
    "CAIM": ("Böhme, Papuga & Lange (2026), Decoupling amplitude and mean stress effects in multiaxial fatigue: A critical–integral approach", "10.1111/ffe.70244"),
    "CARSPA": ("Carpinteri & Spagnoli (2001), Multiaxial high-cycle fatigue criterion for hard metals", "10.1016/S0142-1123(00)00075-X"),
    "DANGVAN": ("Dang Van (1973), Sur la résistance à la fatigue des métaux, Sciences et Technique de l'Armement 47, 647–722", ""),
    "OTT": ("Tveit et al. (2024), A continuum approach to multiaxial high-cycle fatigue modeling for ductile metallic materials", "10.1016/j.rineng.2024.102171"),
    "SWTD": ("Kujawski (2014), A deviatoric version of the SWT parameter", "10.1016/j.ijfatigue.2013.12.002"),
    "ZHU_EDP": ("Zhu et al. (2019), A novel energy-based equivalent damage parameter for multiaxial fatigue life prediction", "10.1016/j.ijfatigue.2018.11.025"),
}
DISABLED_METHODS = {
    "CARSPA": "Missing the published critical-plane orientation construction; maximum-shear/max-damage scans are not that criterion.",
    "DANGVAN": "Combines shear amplitude with independently reconstructed hydrostatic endpoints, losing their time relationship and corrupting the trace for out-of-phase histories.",
    "OTT": "Only a proportional cycle reduction; non-proportional extension is unverified and fatigue-ratio calibration silently falls back to m=2 outside its supported range.",
    "SWTD": "Separately sorted stress and strain eigenvalues do not preserve paired physical directions; also ignores harmonic histories.",
    "ZHU_EDP": "An acknowledged proxy: omits elastoplastic work W0 and substitutes a phase proxy for the published non-proportional factor.",
}


def method_reference(name):
    citation, doi = REFERENCES.get(name.removesuffix("_ext"), ("Unregistered method", ""))
    return {"citation": citation, "url": "https://doi.org/"+doi if doi else "",
            "variant": "maximum damage" if name.endswith("_ext") else "standard"}


def ensure_enabled(name):
    base = name.removesuffix("_ext")
    if base in DISABLED_METHODS:
        raise NotImplementedError(f"{name} disabled by numerical audit: {DISABLED_METHODS[base]} See docs/method-audit.md")
