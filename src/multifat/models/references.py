"""Bibliography and audit decisions. See docs/method-audit.md for scope."""

IMPLEMENTATION_REVISION = "2026-10-07-audit1"
REFERENCES = {
    "FS": (
        "Fatemi & Socie (1988), A critical plane approach to multiaxial fatigue damage including out-of-phase loading",
        "10.1111/j.1460-2695.1988.tb01169.x",
    ),
    "FIN": (
        "Findley (1959), A theory for the effect of mean stress on fatigue of metals under combined torsion and axial load or bending",
        "10.1115/1.4008327",
    ),
    "SWT": (
        "Smith, Watson & Topper (1970), A stress-strain function for the fatigue of metals, Journal of Materials 5(4), 767–778",
        "",
    ),
    "MATAKE": (
        "Matake (1977), An explanation on fatigue limit under combined stress",
        "10.1299/jsme1958.20.257",
    ),
    "MCD": (
        "McDiarmid (1991), A general criterion for high cycle multiaxial fatigue failure",
        "10.1111/j.1460-2695.1991.tb00673.x",
    ),
    "MGSE_YU": (
        "Yu, Zhu, Liu & Liu (2017), Multiaxial fatigue damage parameter and life prediction without any additional material constants, Eq. (11)",
        "10.3390/ma10080923",
    ),
    "MGSE_ZHU": (
        "Zhu et al. (2018), Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials",
        "10.1016/j.ijfatigue.2018.03.028",
    ),
    "MSWT": (
        "Modified SWT form reported by Zhu et al. (2018), Evaluation and comparison of critical plane criteria for multiaxial fatigue analysis of ductile and brittle materials",
        "10.1016/j.ijfatigue.2018.03.028",
    ),
    "MKBM": (
        "Li et al. (2011), Multiaxial fatigue life prediction for various metallic materials based on the critical plane approach",
        "10.1016/j.ijfatigue.2010.07.003",
    ),
    "LIU2021": (
        "Liu, Ran, Wei & Zhang (2021), A critical plane-based multiaxial fatigue life prediction method considering the material sensitivity and the shear stress",
        "10.1016/j.ijpvp.2021.104532",
    ),
    "LI": (
        "Li et al. (2021), Multiaxial fatigue life prediction for metals by means of an improved strain energy density-based critical plane criterion",
        "10.1016/j.euromechsol.2021.104353",
    ),
    "GSE": (
        "Ince & Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life prediction under proportional and non-proportional loadings",
        "10.1016/j.ijfatigue.2013.10.007",
    ),
    "GSA": (
        "Ince & Glinka (2014), A generalized fatigue damage parameter for multiaxial fatigue life prediction under proportional and non-proportional loadings",
        "10.1016/j.ijfatigue.2013.10.007",
    ),
    "BP": (
        "Boehme & Papuga (2024), Advancements in stress-based multiaxial fatigue prediction: A data-driven approach and a new criterion",
        "10.1111/ffe.14281",
    ),
    "CAIM": (
        "Boehme, Papuga & Lange (2026), Decoupling amplitude and mean stress effects in multiaxial fatigue: A critical–integral approach",
        "10.1111/ffe.70244",
    ),
    "CARSPA": (
        "Carpinteri & Spagnoli (2001), Multiaxial high-cycle fatigue criterion for hard metals",
        "10.1016/S0142-1123(00)00075-X",
    ),
    "DANGVAN": (
        "Dang Van (1973), Sur la résistance à la fatigue des métaux, Sciences et Technique de l'Armement 47, 647–722",
        "",
    ),
    "OTT": (
        "Tveit et al. (2024), A continuum approach to multiaxial high-cycle fatigue modeling for ductile metallic materials",
        "10.1016/j.rineng.2024.102171",
    ),
    "TVEIT2025": (
        "Tveit, Ottosen, Erduran & Reyes (2025), Continuum mechanics fatigue model for arbitrary ratios of shear and normal stress fatigue strengths, Eqs. (1)-(12), (24)",
        "10.1016/j.mechmat.2025.105452",
    ),
    "SWTD": (
        "Kujawski (2014), A deviatoric version of the SWT parameter",
        "10.1016/j.ijfatigue.2013.12.002",
    ),
    "ZHU_EDP": (
        "Zhu et al. (2019), A novel energy-based equivalent damage parameter for multiaxial fatigue life prediction",
        "10.1016/j.ijfatigue.2018.11.025",
    ),
    "CROSSLAND": (
        "Crossland (1956), Effect of large hydrostatic pressures on the torsional fatigue strength of an alloy steel; as restated by Papadopoulos et al. (1997), A comparative study of multiaxial high-cycle fatigue criteria for metals, Eqs. (44), (65)-(66)",
        "10.1016/S0142-1123(96)00064-3",
    ),
    "SINES": (
        "Sines (1959), Behavior of metals under complex static and alternating stresses; as restated by Papadopoulos et al. (1997), A comparative study of multiaxial high-cycle fatigue criteria for metals, Eqs. (44), (62)-(64)",
        "10.1016/S0142-1123(96)00064-3",
    ),
    "PAPADOPOULOS": (
        "Papadopoulos (2001), Long life fatigue under multiaxial loading, Eqs. (10), (25), (41)",
        "10.1016/S0142-1123(01)00059-7",
    ),
    "PCR": (
        "Papuga & Růžička (2008), Two new multiaxial criteria for high cycle fatigue computation, Eq. (17), Appendix B",
        "10.1016/j.ijfatigue.2007.02.015",
    ),
    "PI": (
        "Papuga & Růžička (2008), Two new multiaxial criteria for high cycle fatigue computation, Eq. (16), Appendix A",
        "10.1016/j.ijfatigue.2007.02.015",
    ),
    "MDC": (
        "Margetin, Ďurka & Chmelko (2016), Multiaxial fatigue criterion based on parameters from torsion and axial S-N curve, Eqs. (4)-(8)",
        "10.3221/IGF-ESIS.37.20",
    ),
    "MWCM": (
        "Susmel & Lazzarin (2002), A bi-parametric Wöhler curve for high cycle multiaxial fatigue assessment, Eqs. (5), (8)-(13)",
        "10.1046/j.1460-2695.2002.00462.x",
    ),
    "VHN": (
        "Vu, Halm & Nadot (2010), Multiaxial fatigue criterion for complex loading based on stress invariants, Eqs. (6)-(21)",
        "10.1016/j.ijfatigue.2009.11.006",
    ),
    "PBP": (
        "Cristofori, Susmel & Tovo (2008), A stress invariant based criterion to estimate fatigue damage under multiaxial loading, Eqs. (4)-(13)",
        "10.1016/j.ijfatigue.2007.11.006",
    ),
    "MMMH": (
        "Mertens & Liebich (2017), A time dependent equivalent stress function for proportional and non-proportional transient loaded and notched metallic components, Eqs. (1), (4), (20)-(25), smooth surface; endurance-limit Haigh line with M of Dubbel (2004) E1",
        "10.1016/j.ijfatigue.2017.06.003",
    ),
    "COIN": (
        "Schmidt & Kraft (2023), A new equivalent stress approach based on complex invariants: The COIN LiWI method, Eqs. (21)-(32)",
        "10.1016/j.ijfatigue.2022.107474",
    ),
    "MOREL": (
        "Morel (1998), A fatigue life prediction method based on a mesoscopic approach in constant amplitude multiaxial loading, Eqs. (9)-(11), (33)-(38)",
        "10.1046/j.1460-2695.1998.00452.x",
    ),
    "EK": (
        "Ellyin & Kujawski (1993), A multiaxial fatigue criterion including mean-stress effect, ASTM STP 1191, Eqs. (5b), (11)-(18), high-cycle elastic form",
        "10.1520/STP24795S",
    ),
    "MMCCM": (
        "Wang & Susmel (2016), The Modified Manson-Coffin Curve Method to estimate fatigue lifetime under complex constant and variable amplitude multiaxial fatigue loading, Eqs. (1)-(8)",
        "10.1016/j.ijfatigue.2015.10.005",
    ),
    "VSE": (
        "Liu (1993), A method based on virtual strain-energy parameters for multiaxial fatigue life prediction, ASTM STP 1191, Eqs. (7)-(11)",
        "10.1520/STP24796S",
    ),
    "RFB": (
        "Robert, Fogue & Bahuaud (1994), Fatigue life prediction under periodical or random multiaxial stress states, ASTM STP 1231, Eqs. (8), (15)-(25), Appendix",
        "10.1520/STP13959S",
    ),
}
DISABLED_METHODS = {
    "CARSPA": "Missing the published critical-plane orientation construction; maximum-shear/max-damage scans are not that criterion.",
    "ZHU_EDP": "Eq. (23) needs the measured non-proportional factor F_NP (Eq. 22, stabilised stress amplitudes) and the uniaxial elastoplastic work dW_0 (cyclic-curve constants H_c, h_c); the workbook has elastic tensors only. The Eq. (19)-(23) kernel is implemented with explicit inputs.",
}


def method_reference(name):
    citation, doi = REFERENCES.get(name.removesuffix("_ext"), ("Unregistered method", ""))
    return {
        "citation": citation,
        "url": "https://doi.org/" + doi if doi else "",
        "variant": "maximum damage" if name.endswith("_ext") else "standard",
    }


def ensure_enabled(name):
    base = name.removesuffix("_ext")
    if base in DISABLED_METHODS:
        raise NotImplementedError(
            f"{name} disabled by numerical audit: {DISABLED_METHODS[base]} See docs/method-audit.md"
        )
