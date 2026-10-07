"""Single-effect NH3-H2O absorption-cycle model for Project 2."""

from dataclasses import dataclass
from .ammonia_eos import psat_nh3, nh3_enthalpy
from .patek_klomfar import (
    solve_liquid_concentration,
    saturated_liquid_enthalpy,
    mole_to_mass_fraction,
    mass_to_mole_fraction,
)
import math


@dataclass
class CycleResult:
    COP: float
    Q_evap_kW: float
    Q_gen_kW: float
    Q_cond_kW: float
    Q_abs_kW: float
    W_pump_kW: float
    m_ref_kg_s: float
    m_strong_kg_s: float
    m_weak_kg_s: float
    x_strong: float
    x_weak: float
    circulation_ratio: float
    P_low_kPa: float
    P_high_kPa: float
    T_gen_C: float
    T_cond_C: float
    T_abs_C: float
    T_evap_C: float
    refrigerant_vapor_m3_s: float
    pump_vol_m3_h: float
    pump_head_kPa: float
    feasible: bool
    notes: str


def _solution_concentration(P_kPa, T_C):
    try:
        return solve_liquid_concentration(P_kPa, T_C)
    except Exception:
        # If the requested temperature lies outside the fitted VLE range,
        # clamp to a numerically stable concentration range.
        lo, hi = 1e-5, 0.99999
        from scipy.optimize import minimize_scalar
        f = lambda x: abs(
            (solve_bubble_temperature(P_kPa, x) - (T_C + 273.15))
        )
        r = minimize_scalar(f, bounds=(lo, hi), method="bounded")
        return r.x


def solve_bubble_temperature(P_kPa, x_mass):
    from .patek_klomfar import bubble_temperature, mass_to_mole_fraction
    return bubble_temperature(P_kPa, mass_to_mole_fraction(x_mass))


def cycle_performance(Q_evap_kW, T_evap_C=2.0, T_cond_C=45.0,
                      T_gen_C=130.0, y_ref=0.995,
                      eta_pump=0.70):
    """
    Engineering single-effect cycle.

    Strong solution: saturated liquid at absorber pressure and absorber T.
    Weak solution: saturated liquid at generator pressure and generator T.
    Refrigerant: pure NH3 between evaporator and condenser.
    """
    P_low = psat_nh3(T_evap_C) / 1000.0
    P_high = psat_nh3(T_cond_C) / 1000.0

    x_s = _solution_concentration(P_low, T_cond_C)
    x_w = _solution_concentration(P_high, T_gen_C)

    if x_s <= x_w + 1e-4:
        x_w = max(0.01, x_s - 0.05)

    denom = max(y_ref - x_w, 1e-8)
    m_weak_per_ref = (x_s - y_ref) / (x_s - x_w)
    m_strong_per_ref = (y_ref - x_w) / (x_s - x_w)

    # The formulas above are signed depending on stream labeling.
    # Use the positive circulation ratio directly.
    f = (y_ref - x_w) / max(x_s - x_w, 1e-8)
    m_strong = f
    m_weak = max(f - 1.0, 1e-8)

    h_ref_evap = nh3_enthalpy(T_evap_C, P_low, "vapor")
    h_ref_cond = nh3_enthalpy(T_cond_C, P_high, "liquid")
    q_e_per_kg = max(h_ref_evap - h_ref_cond, 1.0)

    m_ref = Q_evap_kW / q_e_per_kg

    h_w = saturated_liquid_enthalpy(T_gen_C + 273.15, x_w)
    h_s = saturated_liquid_enthalpy(T_cond_C + 273.15, x_s)

    # Solution pump: incompressible approximation with density ~ 850 kg/m3.
    rho_sol = 850.0
    deltaP = max(P_high - P_low, 1.0)
    vdot = (m_strong * m_ref) / rho_sol
    Wpump = vdot * deltaP / eta_pump  # kW because kPa*m3/s = kW

    # Component heat balances.
    m_s = m_strong * m_ref
    m_w = m_weak * m_ref
    Q_gen = max(m_ref * h_ref_cond + m_w * h_w
                - m_s * h_s + Wpump, 0.0)
    Q_cond = m_ref * (nh3_enthalpy(T_cond_C, P_high, "vapor")
                      - h_ref_cond)
    Q_abs = max(Q_gen + Q_evap_kW + Wpump - Q_cond, 0.0)

    COP = Q_evap_kW / max(Q_gen, 1e-9)

    # Ideal-gas estimate of ammonia vapour volumetric flow.
    Rspec = 8.314462618 / MW_NH3
    v_ref = Rspec * (T_evap_C + 273.15) / (P_low * 1000.0)
    vdot_ref = m_ref * v_ref

    feasible = (
        P_high <= 2000.0
        and T_cond_C >= 25.0
        and x_s > x_w
        and 0.45 <= COP <= 0.8
    )

    notes = ""
    if P_high > 2000:
        notes += "High-pressure vessel limit exceeded. "
    if T_cond_C < 25:
        notes += "Condenser temperature below required 25 C. "

    return CycleResult(
        COP=COP,
        Q_evap_kW=Q_evap_kW,
        Q_gen_kW=Q_gen,
        Q_cond_kW=Q_cond,
        Q_abs_kW=Q_abs,
        W_pump_kW=Wpump,
        m_ref_kg_s=m_ref,
        m_strong_kg_s=m_s,
        m_weak_kg_s=m_w,
        x_strong=x_s,
        x_weak=x_w,
        circulation_ratio=m_strong,
        P_low_kPa=P_low,
        P_high_kPa=P_high,
        T_gen_C=T_gen_C,
        T_cond_C=T_cond_C,
        T_abs_C=T_cond_C,
        T_evap_C=T_evap_C,
        refrigerant_vapor_m3_s=vdot_ref,
        pump_vol_m3_h=vdot * 3600.0,
        pump_head_kPa=deltaP / eta_pump,
        feasible=feasible,
        notes=notes,
    )
