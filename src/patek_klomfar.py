"""Ammonia-water bubble/dew temperature correlation utilities."""

import math

P0 = 2.0e6
T0 = 100.0

BUBBLE = [
    (0, 0, 3.223020), (0, 1, -0.3842060),
    (0, 2, 0.04609650), (0, 3, -0.003789450),
    (0, 4, 0.0001356100), (1, 0, 0.4877550),
    (1, 1, -0.1201080), (1, 2, 0.01061540),
    (2, 3, -0.0005335890), (4, 0, 7.850410),
    (5, 0, -11.59410), (5, 1, -0.05231500),
    (6, 0, 0.04895960), (13, 1, 0.04210590),
]

DEW = [
    (0, 0, 3.240040), (0, 1, -0.3959200),
    (0, 2, -0.04356240), (0, 3, -0.002189430),
    (1, 0, -1.435260), (1, 1, 1.052560),
    (1, 2, -0.07192810), (2, 0, 12.23620),
    (2, 1, -2.243680), (3, 0, -21.07800),
    (3, 1, 1.108340), (4, 0, 14.53990),
    (4, 2, 0.6443120), (5, 0, -2.212640),
    (5, 2, -0.7562660), (6, 0, -1.355290),
    (7, 2, 0.1835410),
]


def _temperature(P_Pa, x, coefficients, dew=False):
    if P_Pa <= 0 or not 0 <= x <= 1:
        raise ValueError("Invalid pressure or composition.")

    p_term = math.log(P0 / P_Pa)
    total = 0.0

    for m, n, coefficient in coefficients:
        exponent = m / 4.0 if dew else m
        total += (
            coefficient * (1 - x)**exponent * p_term**n
        )

    return T0 * total


def bubble_temperature(P_Pa, x_liquid):
    """Return bubble temperature in K."""
    return _temperature(P_Pa, x_liquid, BUBBLE)


def dew_temperature(P_Pa, y_vapor):
    """Return dew temperature in K."""
    return _temperature(P_Pa, y_vapor, DEW, dew=True)


def solution_enthalpy(T_K, x_nh3, weak=False):
    """Approximate solution enthalpy in kJ/kg; requires validation."""
    cp = 4.18 * (1 - x_nh3) + 4.70 * x_nh3
    return cp * (T_K - 273.15)


def vapor_enthalpy(T_K, y_nh3):
    """Approximate vapour enthalpy in kJ/kg; requires validation."""
    cp = 1.90 * y_nh3 + 1.86 * (1 - y_nh3)
    return 1250.0 * y_nh3 + cp * (T_K - 273.15)
