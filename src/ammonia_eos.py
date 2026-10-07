"""Pure-ammonia Peng–Robinson EOS and saturation-pressure routines.

The project requires ammonia to be treated with a real-gas EOS rather than
a black-box ammonia property library.
"""

import math
from dataclasses import dataclass

R = 8.314462618  # J/mol/K
MW_NH3 = 17.031e-3  # kg/mol
TC = 405.56       # K
PC = 11.333e6     # Pa
OMEGA = 0.256


def lee_kesler_pitzer_psat(T_K: float) -> float:
    """Generalized Pitzer/Lee-Kesler vapour-pressure correlation, Pa."""
    Tr = T_K / TC
    if Tr <= 0:
        raise ValueError("Temperature must be positive.")
    if Tr >= 1:
        return PC

    ln_Pr0 = (
        5.92714
        - 6.09648 / Tr
        - 1.28862 * math.log(Tr)
        + 0.169347 * Tr**6
    )
    ln_Pr1 = (
        15.2518
        - 15.6875 / Tr
        - 13.4721 * math.log(Tr)
        + 0.43577 * Tr**6
    )
    return PC * math.exp(ln_Pr0 + OMEGA * ln_Pr1)


def psat_nh3(T_C: float) -> float:
    return lee_kesler_pitzer_psat(T_C + 273.15)


def _pr_parameters(T_K):
    kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA**2
    alpha = (1 + kappa * (1 - math.sqrt(T_K / TC)))**2
    a = 0.45724 * R**2 * TC**2 / PC * alpha
    b = 0.07780 * R * TC / PC
    return a, b


def _z_roots(T_K, P_Pa):
    a, b = _pr_parameters(T_K)
    A = a * P_Pa / (R**2 * T_K**2)
    B = b * P_Pa / (R * T_K)

    coeff = [
        1.0,
        -(1.0 - B),
        A - 3.0 * B**2 - 2.0 * B,
        -(A * B - B**2 - B**3),
    ]
    roots = np_roots_real(coeff)
    return sorted(roots), A, B


def np_roots_real(coeff):
    # Local import keeps this module lightweight.
    import numpy as np
    roots = np.roots(coeff)
    return [float(r.real) for r in roots if abs(float(r.imag)) < 1e-8]


def z_factor(T_K, P_Pa, phase="vapor"):
    roots, _, _ = _z_roots(T_K, P_Pa)
    if not roots:
        raise RuntimeError("No real Peng–Robinson Z root.")
    return max(roots) if phase == "vapor" else min(roots)


def ideal_cp_nh3(T_K):
    # Smooth engineering approximation around the project operating range.
    # Used only for the ideal-gas contribution; departure is supplied by PR.
    return 35.06 + 0.0105 * (T_K - 273.15)  # J/mol/K


def ideal_h_nh3(T_K, T_ref=273.15):
    a = 35.06
    b = 0.0105
    return (a * (T_K - T_ref) + 0.5 * b *
            ((T_K - 273.15)**2 - (T_ref - 273.15)**2)) / MW_NH3


def pr_enthalpy_departure(T_K, P_Pa, phase="vapor"):
    """PR enthalpy departure in kJ/kg."""
    a, b = _pr_parameters(T_K)
    Z = z_factor(T_K, P_Pa, phase)
    _, _, B = _z_roots(T_K, P_Pa)

    # Numerical derivative of a(T) is sufficiently stable for this model.
    dT = max(1e-3, 1e-4 * T_K)
    ap, _ = _pr_parameters(T_K + dT)
    am, _ = _pr_parameters(T_K - dT)
    da_dT = (ap - am) / (2 * dT)

    sqrt2 = math.sqrt(2.0)
    term_log = math.log(
        (Z + (1 + sqrt2) * B) /
        (Z + (1 - sqrt2) * B)
    )
    hdep_molar = (
        R * T_K * (Z - 1.0)
        + (T_K * da_dT - a) /
        (2.0 * sqrt2 * b) * term_log
    )
    return hdep_molar / 1000.0 / MW_NH3


def nh3_enthalpy(T_C, P_kPa, phase="vapor"):
    T_K = T_C + 273.15
    P_Pa = P_kPa * 1000.0
    return ideal_h_nh3(T_K) + pr_enthalpy_departure(T_K, P_Pa, phase)


def nh3_entropy(T_C, P_kPa, phase="vapor"):
    """Approximate entropy from ideal-gas contribution + PR pressure correction."""
    T_K = T_C + 273.15
    P_Pa = P_kPa * 1000.0
    cp = ideal_cp_nh3(T_K)
    s_ideal = (
        cp * math.log(T_K / 273.15)
        - R * math.log(P_Pa / 101325.0)
    ) / 1000.0 / MW_NH3

    # A compact EOS correction based on Z.
    Z = z_factor(T_K, P_Pa, phase)
    s_dep = R * math.log(max(Z, 1e-12)) / 1000.0 / MW_NH3
    return s_ideal + s_dep


def saturation_state(T_C):
    P = psat_nh3(T_C) / 1000.0
    hf = nh3_enthalpy(T_C, P, "liquid")
    hg = nh3_enthalpy(T_C, P, "vapor")
    return {"P_kPa": P, "h_f": hf, "h_g": hg, "h_fg": hg - hf}
