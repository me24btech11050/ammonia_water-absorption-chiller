"""Preliminary single-effect ammonia-water absorption cycle."""

from dataclasses import dataclass
from src.ammonia_eos import saturation_pressure


@dataclass
class CycleResult:
    Q_evap_kW: float
    Q_gen_kW: float
    Q_cond_kW: float
    Q_abs_kW: float
    W_pump_kW: float
    COP: float
    P_high_bar: float
    P_low_bar: float
    T_gen_C: float
    T_cond_C: float
    x_strong: float
    x_weak: float
    m_refrigerant_kg_s: float
    m_strong_kg_s: float
    m_weak_kg_s: float


def solve_cycle(
    Q_evap_kW,
    T_evap_C=2.0,
    T_cond_C=40.0,
    T_gen_C=130.0,
    refrigerant_purity=0.995,
):
    """Return a preliminary cycle estimate."""
    if Q_evap_kW <= 0:
        raise ValueError("Cooling duty must be positive.")

    P_low = saturation_pressure(T_evap_C + 273.15)
    P_high = saturation_pressure(T_cond_C + 273.15)

    if P_high / 1e5 > 20:
        raise ValueError("High-side pressure exceeds 20 bar.")

    # Preliminary estimates; replace with validated property calculations.
    refrigerant_effect_kJ_kg = max(
        1050.0, 1250.0 - 2.0 * (T_cond_C - 40.0)
    )
    m_ref = Q_evap_kW / refrigerant_effect_kJ_kg

    x_strong = min(0.92, max(0.55, 0.78 - 0.0015 * (T_cond_C - 35)))
    x_weak = min(0.65, max(0.25, 0.48 - 0.0010 * (T_gen_C - 120)))

    if x_strong <= x_weak:
        raise ValueError("Strong-solution concentration must exceed weak.")

    # Approximate mass balance; concentrations must ultimately be
    # defined consistently on a mass or mole basis.
    m_weak = m_ref * (x_strong - refrigerant_purity) / (
        x_strong - x_weak
    )
    m_strong = m_weak + m_ref

    pump_work = 0.02 * m_strong
    cop = max(
        0.35,
        min(0.65, 0.58 - 0.004 * (T_cond_C - 40) + 0.001 * (T_gen_C - 130)),
    )

    Q_gen = Q_evap_kW / cop
    Q_cond = Q_gen + Q_evap_kW + pump_work - 0.08 * Q_gen
    Q_abs = max(
        0.0, Q_gen + Q_evap_kW + pump_work - Q_cond
    )

    return CycleResult(
        Q_evap_kW=Q_evap_kW,
        Q_gen_kW=Q_gen,
        Q_cond_kW=Q_cond,
        Q_abs_kW=Q_abs,
        W_pump_kW=pump_work,
        COP=cop,
        P_high_bar=P_high / 1e5,
        P_low_bar=P_low / 1e5,
        T_gen_C=T_gen_C,
        T_cond_C=T_cond_C,
        x_strong=x_strong,
        x_weak=x_weak,
        m_refrigerant_kg_s=m_ref,
        m_strong_kg_s=m_strong,
        m_weak_kg_s=m_weak,
    )


def required_generator_temperature(P_high_bar):
    """Preliminary estimate only; not an equilibrium calculation."""
    return 75.0 + 3.0 * P_high_bar
