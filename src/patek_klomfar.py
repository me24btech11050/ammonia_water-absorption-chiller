"""Pátek–Klomfar ammonia–water correlations.

Mass-fraction inputs are converted to ammonia mole fraction internally.
The correlations are used in the pressure/temperature/composition range
relevant to the project (up to about 20 bar for the enthalpy fits).
"""

import math
from scipy.optimize import brentq

T0 = 100.0       # K
P0_MPa = 2.0     # MPa
H0_L = 100.0     # kJ/kg, as used with the saturated-liquid correlation
H0_G = 1000.0    # kJ/kg
T0_HG = 324.0    # K


# Bubble/dew coefficients supplied in the project statement.
BUBBLE = [
    (0,0,3.223020e0),(0,1,-3.842060e-1),(0,2,4.609650e-2),
    (0,3,-3.789450e-3),(0,4,1.356100e-4),(1,0,4.877550e-1),
    (1,1,-1.201080e-1),(1,2,1.061540e-2),(2,3,-5.335890e-4),
    (4,0,7.850410e0),(5,0,-1.159410e1),(5,1,-5.231500e-2),
    (6,0,4.895960e0),(13,1,4.210590e-2),
]

DEW = [
    (0,0,3.240040e0),(0,1,-3.959200e-1),(0,2,4.356240e-2),
    (0,3,-2.189430e-3),(1,0,-1.435260e0),(1,1,1.052560e0),
    (1,2,-7.192810e-2),(2,0,1.223620e1),(2,1,-2.243680e0),
    (3,0,-2.01780e1),(3,1,1.108340e0),(4,0,1.453990e1),
    (4,2,6.443120e-1),(5,0,-2.212640e0),(5,2,-7.562660e-1),
    (6,0,-1.355290e0),(7,2,1.835410e-1),
]

# Saturated-liquid enthalpy coefficients:
# hL = h0 sum a_i (T/T0 - 1)^m_i x^n_i
HL = [
    (0,1,-7.61080),(0,4,25.6905),(0,8,-247.092),(0,9,325.952),
    (0,12,-158.854),(0,14,61.9084),(1,0,11.4314),(1,1,1.18157),
    (2,1,2.84179),(3,3,7.41609),(5,3,891.844),(5,4,-1613.09),
    (5,5,622.106),(6,2,-207.588),(6,4,-6.87393),(8,0,3.50716),
]

# Saturated-vapor enthalpy coefficients:
# hg = h0 sum a_i (1-T/T0)^m_i (1-y)^(n_i/4)
HG = [
    (0,0,-68.7393),(0,1,35.0716),(0,2,16.1309),(0,3,-89.1844),
    (0,4,7.41609),(0,5,2.84179),(0,6,1.18157),(0,7,11.4314),
    (1,0,61.9084),(2,1,-158.854),(3,2,-247.092),(3,3,325.952),
    (4,0,25.6905),(4,1,-76.1080),(4,2,30.8482),(4,3,98.8009),
    (5,0,-96.1248),
]


def mass_to_mole_fraction(w):
    """Ammonia mass fraction -> ammonia mole fraction."""
    return (w / 17.031) / ((w / 17.031) + ((1.0 - w) / 18.015))


def mole_to_mass_fraction(x):
    return x * 17.031 / (x * 17.031 + (1.0 - x) * 18.015)


def bubble_temperature(P_kPa, x):
    p = P_kPa / 1000.0
    if not (0 < p <= 2.0):
        raise ValueError("Pátek–Klomfar correlation intended for 0 < P <= 2 MPa.")
    L = math.log(P0_MPa / p)
    return T0 * sum(a * (1.0-x)**m * L**n for m,n,a in BUBBLE)


def dew_temperature(P_kPa, y):
    p = P_kPa / 1000.0
    if not (0 < p <= 2.0):
        raise ValueError("Pátek–Klomfar correlation intended for 0 < P <= 2 MPa.")
    L = math.log(P0_MPa / p)
    return T0 * sum(a * (1.0-y)**(m/4.0) * L**n for m,n,a in DEW)


def vapor_mass_fraction_from_liquid(P_kPa, x_mass):
    """Equilibrium vapor mass fraction from a liquid mass fraction."""
    x = mass_to_mole_fraction(x_mass)
    p = P_kPa / 1000.0
    L = math.log(P0_MPa / p)

    # Pátek–Klomfar vapor composition relation.
    s = sum(a * x**n * L**m for m,n,a in [
        (0,0,19.8022017),(0,1,-11.8092669),(0,6,27.7479980),
        (0,7,-28.8634277),(1,0,-59.1616608),(2,1,5.78091305),
        (2,2,-6.21736743),(3,2,-3421.98402),(4,3,11940.3127),
        (5,4,-24541.3777),(6,5,29159.1865),(7,6,-18.478229),
        (7,7,23.4819434),(8,7,4803.10617)
    ])
    y = 1.0 - math.exp(math.log(max(1.0-x, 1e-14)) * s)
    return mole_to_mass_fraction(max(0.0, min(1.0, y)))


def saturated_liquid_enthalpy(T_K, x_mass):
    x = mass_to_mole_fraction(x_mass)
    tau = T_K / T0 - 1.0
    return H0_L * sum(a * tau**m * x**n for m,n,a in HL)


def saturated_vapor_enthalpy(T_K, y_mass):
    y = mass_to_mole_fraction(y_mass)
    theta = 1.0 - T_K / T0_HG
    return H0_G * sum(a * theta**m * (1.0-y)**(n/4.0)
                      for m,n,a in HG)


def solution_state(P_kPa, T_C, x_mass):
    """Saturated liquid solution state."""
    T_K = T_C + 273.15
    h = saturated_liquid_enthalpy(T_K, x_mass)
    return {"P_kPa": P_kPa, "T_C": T_C, "x_mass": x_mass, "h_kJkg": h}


def solve_liquid_concentration(P_kPa, T_C):
    """Solve bubble-point liquid concentration from P and T."""
    f = lambda x: bubble_temperature(P_kPa, mass_to_mole_fraction(x)) - (T_C + 273.15)
    return brentq(f, 1e-8, 1.0 - 1e-8)
