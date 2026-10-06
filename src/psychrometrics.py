"""
Psychrometric calculations for ME3160 Project 2.

Reference:
T_ref = 0 °C
P_ref = 101.325 kPa

Moist air is treated as an ideal-gas mixture of dry air and
water vapour, following the project statement.
"""

import math


# Atmospheric pressure
P_ATM = 101.325  # kPa

# Specific heats used by the project formulation
CP_DA = 1.006     # kJ/(kg_da K)
CP_V = 1.86       # kJ/(kg_v K)

# Latent heat constant used in the project enthalpy relation
H_FG_0 = 2501.0   # kJ/kg


def saturation_pressure_water(T_C):
    """
    Saturation pressure of water vapour in kPa.

    This function will be completed using the permitted
    water-property calculation method.
    """
    raise NotImplementedError(
        "Water saturation-pressure correlation will be implemented next."
    )


def humidity_ratio_from_rh(T_C, RH, P_kPa=P_ATM):
    """
    Calculate humidity ratio from dry-bulb temperature and RH.

    Parameters
    ----------
    T_C : float
        Dry-bulb temperature [°C].
    RH : float
        Relative humidity as a fraction (0 to 1).
    P_kPa : float
        Atmospheric pressure [kPa].

    Returns
    -------
    float
        Humidity ratio [kg water/kg dry air].
    """

    P_sat = saturation_pressure_water(T_C)
    P_v = RH * P_sat

    return 0.622 * P_v / (P_kPa - P_v)


def relative_humidity_from_humidity_ratio(
    T_C, w, P_kPa=P_ATM
):
    """
    Calculate relative humidity from temperature and humidity ratio.
    """

    P_sat = saturation_pressure_water(T_C)

    P_v = P_kPa * w / (0.622 + w)

    return P_v / P_sat


def moist_air_enthalpy(T_C, w):
    """
    Moist-air specific enthalpy.

    h = Cp_da*T + w*(2501 + Cp_v*T)

    Parameters
    ----------
    T_C : float
        Dry-bulb temperature [°C].
    w : float
        Humidity ratio [kg/kg dry air].

    Returns
    -------
    float
        Enthalpy [kJ/kg dry air].
    """

    return CP_DA * T_C + w * (H_FG_0 + CP_V * T_C)


def humidity_ratio_from_vapor_pressure(P_v, P_kPa=P_ATM):
    """
    Calculate humidity ratio from water-vapour partial pressure.
    """

    return 0.622 * P_v / (P_kPa - P_v)


def vapor_pressure_from_humidity_ratio(w, P_kPa=P_ATM):
    """
    Calculate water-vapour partial pressure from humidity ratio.
    """

    return P_kPa * w / (0.622 + w)


def adiabatic_mixing(
    T1_C,
    RH1,
    T2_C,
    RH2,
    fraction_1,
    fraction_2,
    P_kPa=P_ATM,
):
    """
    Adiabatically mix two moist-air streams.

    Fractions are based on dry-air mass flow.

    Returns
    -------
    dict
        Mixed-air temperature, humidity ratio, enthalpy,
        and vapour pressure.
    """

    w1 = humidity_ratio_from_rh(T1_C, RH1, P_kPa)
    w2 = humidity_ratio_from_rh(T2_C, RH2, P_kPa)

    h1 = moist_air_enthalpy(T1_C, w1)
    h2 = moist_air_enthalpy(T2_C, w2)

    w_mix = (
        fraction_1 * w1
        + fraction_2 * w2
    )

    h_mix = (
        fraction_1 * h1
        + fraction_2 * h2
    )

    # Solve the enthalpy equation for T_mix:
    #
    # h = Cp_da*T + w*(2501 + Cp_v*T)
    #
    # h = (Cp_da + w*Cp_v)*T + w*2501

    T_mix_C = (
        h_mix - w_mix * H_FG_0
    ) / (
        CP_DA + w_mix * CP_V
    )

    P_v_mix = vapor_pressure_from_humidity_ratio(
        w_mix,
        P_kPa
    )

    return {
        "T_C": T_mix_C,
        "w": w_mix,
        "h": h_mix,
        "P_v_kPa": P_v_mix,
    }


if __name__ == "__main__":
    print("Psychrometric module loaded successfully.")
