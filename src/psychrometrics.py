"""Ideal-gas moist-air routines for ME3160 Project 2."""

from CoolProp.CoolProp import PropsSI
import numpy as np

P_ATM = 101.325  # kPa
CP_DA = 1.006    # kJ/(kg_da K)
CP_V = 1.86      # kJ/(kg_v K)
H_FG_0 = 2501.0  # kJ/kg


def saturation_pressure_water(T_C: float) -> float:
    """Water saturation pressure, kPa."""
    return PropsSI("P", "T", float(T_C) + 273.15, "Q", 0, "Water") / 1000.0


def humidity_ratio_from_rh(T_C: float, RH: float, P_kPa: float = P_ATM) -> float:
    """Humidity ratio kg_v/kg_da."""
    if not (0.0 <= RH <= 1.0):
        raise ValueError("RH must be between 0 and 1.")
    ps = saturation_pressure_water(T_C)
    pv = RH * ps
    if pv >= P_kPa:
        raise ValueError("Water-vapour partial pressure must be below total pressure.")
    return 0.622 * pv / (P_kPa - pv)


def vapor_pressure_from_humidity_ratio(w: float, P_kPa: float = P_ATM) -> float:
    return P_kPa * w / (0.622 + w)


def relative_humidity_from_humidity_ratio(T_C: float, w: float,
                                          P_kPa: float = P_ATM) -> float:
    pv = vapor_pressure_from_humidity_ratio(w, P_kPa)
    return pv / saturation_pressure_water(T_C)


def moist_air_enthalpy(T_C: float, w: float) -> float:
    """kJ/kg dry air, using the project reference."""
    return CP_DA * T_C + w * (H_FG_0 + CP_V * T_C)


def temperature_from_h_w(h: float, w: float) -> float:
    return (h - H_FG_0 * w) / (CP_DA + CP_V * w)


def adiabatic_mixing(T1_C, RH1, T2_C, RH2, fraction_1, fraction_2,
                     P_kPa=P_ATM):
    """Adiabatic mixing on a dry-air mass basis."""
    if abs(fraction_1 + fraction_2 - 1.0) > 1e-12:
        raise ValueError("Dry-air fractions must sum to one.")

    w1 = humidity_ratio_from_rh(T1_C, RH1, P_kPa)
    w2 = humidity_ratio_from_rh(T2_C, RH2, P_kPa)
    h1 = moist_air_enthalpy(T1_C, w1)
    h2 = moist_air_enthalpy(T2_C, w2)

    w = fraction_1 * w1 + fraction_2 * w2
    h = fraction_1 * h1 + fraction_2 * h2
    T = temperature_from_h_w(h, w)
    pv = vapor_pressure_from_humidity_ratio(w, P_kPa)

    return {"T_C": T, "w": w, "h": h, "P_v_kPa": pv}


def wet_bulb_temperature(T_db_C, RH, P_kPa=P_ATM):
    """Wet-bulb temperature from an adiabatic-saturation energy balance."""
    from scipy.optimize import brentq

    w1 = humidity_ratio_from_rh(T_db_C, RH, P_kPa)
    h1 = moist_air_enthalpy(T_db_C, w1)

    def residual(Twb):
        ws = humidity_ratio_from_rh(Twb, 1.0, P_kPa)
        # Liquid-water enthalpy referenced consistently at 0 C.
        hw = 4.186 * Twb
        h_out = moist_air_enthalpy(Twb, ws)
        # Approximate adiabatic saturation balance per kg dry air.
        return h_out + (ws - w1) * hw - h1

    lo = -50.0
    hi = T_db_C
    return brentq(residual, lo, hi)


def coil_from_adp(T_mix_C, w_mix, T_supply_C, BF=0.15, P_kPa=P_ATM):
    """Solve coil ADP and outlet humidity from BF relations."""
    if not (0 < BF < 1):
        raise ValueError("BF must be between zero and one.")

    adp = (T_supply_C - BF * T_mix_C) / (1.0 - BF)
    w_adp = humidity_ratio_from_rh(adp, 1.0, P_kPa)
    w_supply = w_adp + BF * (w_mix - w_adp)
    h_mix = moist_air_enthalpy(T_mix_C, w_mix)
    h_supply = moist_air_enthalpy(T_supply_C, w_supply)
    q_coil = h_mix - h_supply

    return {
        "T_ADP_C": adp,
        "w_ADP": w_adp,
        "T_supply_C": T_supply_C,
        "w_supply": w_supply,
        "h_supply": h_supply,
        "q_coil_kJ_per_kg_da": q_coil,
    }


def room_moisture_iteration(T_room_C, internal_latent_W, m_da_kg_s,
                            w_supply, P_kPa=P_ATM, tol=1e-9):
    """
    Steady room moisture balance:
        m_da*(w_room-w_supply) = m_water_generation
    with room RH determined from the converged humidity ratio.
    """
    m_water_kg_s = internal_latent_W / 2_500_000.0
    w_room = w_supply + m_water_kg_s / m_da_kg_s
    rh = relative_humidity_from_humidity_ratio(T_room_C, w_room, P_kPa)
    return {"w_room": w_room, "RH_room": rh, "m_condensate_kg_s": 0.0}
