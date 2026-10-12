"""Peng-Robinson equation of state utilities for ammonia."""

import math
import numpy as np
from scipy.optimize import brentq

R = 8.314462618
TC = 405.56
PC = 11.333e6
OMEGA = 0.256
MW = 0.017031


def pr_parameters(T):
    kappa = (
        0.37464 + 1.54226 * OMEGA
        - 0.26992 * OMEGA**2
    )
    alpha = (
        1 + kappa * (1 - math.sqrt(T / TC))
    )**2

    a = 0.45724 * R**2 * TC**2 / PC * alpha
    b = 0.07780 * R * TC / PC
    return a, b


def pr_roots(T, P):
    a, b = pr_parameters(T)
    A = a * P / (R**2 * T**2)
    B = b * P / (R * T)

    roots = np.roots([
        1.0,
        -(1 - B),
        A - 3 * B**2 - 2 * B,
        -(A * B - B**2 - B**3),
    ])

    real = sorted(
        float(z.real) for z in roots
        if abs(z.imag) < 1e-8 and z.real > B
    )

    if not real:
        raise ValueError("No physical EOS root.")

    return real


def ln_fugacity_coefficient(Z, A, B):
    s2 = math.sqrt(2)

    return (
        Z - 1 - math.log(Z - B)
        - A / (2 * s2 * B)
        * math.log(
            (Z + (1 + s2) * B)
            / (Z + (1 - s2) * B)
        )
    )


def saturation_pressure(T):
    """Approximate ammonia saturation pressure in Pa."""
    if not 195.0 < T < TC:
        raise ValueError("Temperature outside supported range.")

    def residual(logP):
        P = math.exp(logP)
        a, b = pr_parameters(T)
        A = a * P / (R**2 * T**2)
        B = b * P / (R * T)

        roots = pr_roots(T, P)

        if len(roots) < 2:
            return float("nan")

        return (
            ln_fugacity_coefficient(roots[-1], A, B)
            - ln_fugacity_coefficient(roots[0], A, B)
        )

    grid = np.linspace(math.log(1e3), math.log(PC * 0.999), 400)
    previous_x = grid[0]
    previous_y = residual(previous_x)

    for x in grid[1:]:
        y = residual(x)

        if math.isfinite(previous_y) and math.isfinite(y):
            if previous_y * y < 0:
                return math.exp(brentq(residual, previous_x, x))

        previous_x, previous_y = x, y

    # Corresponding-states fallback; less accurate than EOS fugacity equality.
    Tr = T / TC
    return PC * math.exp(5.373 * (1 + OMEGA) * (1 - 1 / Tr))


def compressibility(T, P, phase="vapor"):
    roots = pr_roots(T, P)
    return roots[-1] if phase == "vapor" else roots[0]


def ideal_cp_ammonia(T):
    return 35.0 + 0.025 * (T - 273.15)


def ideal_enthalpy(T, Tref=273.15):
    dT = T - Tref
    return 35.0 * dT + 0.0125 * dT**2


def ideal_entropy(T, Tref=273.15):
    return (
        35.0 * math.log(T / Tref)
        + 0.025 * (T - Tref - Tref * math.log(T / Tref))
    )


def departure(T, P, phase="vapor"):
    a, b = pr_parameters(T)
    Z = compressibility(T, P, phase)
    B = b * P / (R * T)

    dT = max(1e-3, T * 1e-5)
    da_dT = (
        pr_parameters(T + dT)[0]
        - pr_parameters(T - dT)[0]
    ) / (2 * dT)

    s2 = math.sqrt(2)
    logterm = math.log(
        (Z + (1 + s2) * B)
        / (Z + (1 - s2) * B)
    )

    hres = (
        R * T * (Z - 1)
        + (T * da_dT - a) / (2 * s2 * b) * logterm
    )

    sres = (
        R * math.log(Z - B)
        + da_dT / (2 * s2 * b) * logterm
    )

    return hres, sres


def enthalpy(T, P, phase="vapor"):
    return (
        ideal_enthalpy(T) + departure(T, P, phase)[0]
    ) / MW / 1000.0


def entropy(T, P, phase="vapor"):
    return (
        ideal_entropy(T)
        - R * math.log(P / 101325.0)
        + departure(T, P, phase)[1]
    ) / MW / 1000.0
