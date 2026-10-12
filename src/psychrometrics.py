"""Psychrometric calculations for moist air."""

from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

P_ATM = 101325.0


def psat_water(T_C):
    return PropsSI("P", "T", T_C + 273.15, "Q", 0, "Water")


def humidity_ratio(T_C, RH):
    if not 0.0 <= RH <= 1.0:
        raise ValueError("RH must be between 0 and 1.")

    pv = RH * psat_water(T_C)

    if pv >= P_ATM:
        raise ValueError("Vapour pressure exceeds atmospheric pressure.")

    return 0.621945 * pv / (P_ATM - pv)


def moist_air_enthalpy(T_C, w):
    return 1.006 * T_C + w * (2501.0 + 1.86 * T_C)


def wet_bulb_temperature(T_db, RH):
    w = humidity_ratio(T_db, RH)
    h = moist_air_enthalpy(T_db, w)

    def balance(T_wb):
        ws = humidity_ratio(T_wb, 1.0)
        hs = moist_air_enthalpy(T_wb, ws)
        return hs - (ws - w) * 4.186 * T_wb - h

    return brentq(balance, -50.0, T_db)


def mix_air(T1, RH1, T2, RH2, fraction1=0.70):
    w1 = humidity_ratio(T1, RH1)
    w2 = humidity_ratio(T2, RH2)

    h1 = moist_air_enthalpy(T1, w1)
    h2 = moist_air_enthalpy(T2, w2)

    w = fraction1 * w1 + (1 - fraction1) * w2
    h = fraction1 * h1 + (1 - fraction1) * h2
    T = (h - 2501.0 * w) / (1.006 + 1.86 * w)

    return {"T_C": T, "w": w, "h_kJ_kg_da": h}


def coil_from_bypass(T_mix, w_mix, T_supply, bypass_factor=0.15):
    T_adp = (T_supply - bypass_factor * T_mix) / (
        1 - bypass_factor
    )
    w_adp = humidity_ratio(T_adp, 1.0)
    w_supply = w_adp + bypass_factor * (w_mix - w_adp)

    return {"ADP_C": T_adp, "w_supply": w_supply}
