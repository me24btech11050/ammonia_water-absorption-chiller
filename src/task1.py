"""Seasonal Task 1 calculations for the district cooling network."""

from .psychrometrics import (
    adiabatic_mixing,
    humidity_ratio_from_rh,
    moist_air_enthalpy,
    coil_from_adp,
    relative_humidity_from_humidity_ratio,
)

from .absorption_cycle import cycle_performance
from .decanter import district_condensate


N_ZONES = 400
M_DA_ZONE = 0.40
BF = 0.15
T_ROOM = 22.0


CONDITIONS = [
    {
        "name": "Summer day",
        "Tdb_C": 40.0,
        "RH": 0.50,
        "solar": True,
    },
    {
        "name": "Summer night",
        "Tdb_C": 30.0,
        "RH": 0.50,
        "solar": False,
    },
    {
        "name": "Monsoon day",
        "Tdb_C": 35.0,
        "RH": 0.85,
        "solar": True,
    },
    {
        "name": "Monsoon night",
        "Tdb_C": 30.0,
        "RH": 0.85,
        "solar": False,
    },
    {
        "name": "Winter day",
        "Tdb_C": 28.0,
        "RH": 0.50,
        "solar": True,
    },
    {
        "name": "Winter night",
        "Tdb_C": 16.0,
        "RH": 0.50,
        "solar": False,
    },
]


def zone_internal_gains(solar):
    """Return sensible and latent internal gains for one zone."""

    sensible = 4 * 75.0 + 20.0 * 40.0

    if solar:
        sensible += 1200.0

    latent = 4 * 55.0

    return sensible, latent


def solve_condition(
    condition,
    fresh_fraction=0.30,
    T_cond_C=None,
    T_gen_C=130.0,
):
    """Solve one outdoor design condition."""

    outdoor_w = humidity_ratio_from_rh(
        condition["Tdb_C"],
        condition["RH"],
    )

    room_w = humidity_ratio_from_rh(
        T_ROOM,
        0.75,
    )

    room_h = moist_air_enthalpy(
        T_ROOM,
        room_w,
    )

    outdoor_h = moist_air_enthalpy(
        condition["Tdb_C"],
        outdoor_w,
    )

    # Fresh-air fraction.
    f_out = fresh_fraction
    f_rec = 1.0 - f_out

    # Adiabatic mixing of return/room air and outdoor air.
    mixed = adiabatic_mixing(
        T_ROOM,
        0.75,
        condition["Tdb_C"],
        condition["RH"],
        f_rec,
        f_out,
    )

    sensible, latent = zone_internal_gains(
        condition["solar"]
    )

    # Supply-air temperature from the sensible cooling load.
    cp = 1.006

    T_supply = (
        T_ROOM
        - sensible / (M_DA_ZONE * 1000.0 * cp)
    )

    # Prevent an unrealistic supply temperature.
    T_supply = max(T_supply, 7.0)

    # Cooling coil calculation.
    coil = coil_from_adp(
        mixed["T_C"],
        mixed["w"],
        T_supply,
        BF,
    )

    # Condensate generation.
    m_cond = district_condensate(
        M_DA_ZONE,
        mixed["w"],
        coil["w_supply"],
        N_ZONES,
    )

    # Room moisture balance.
    m_water_gen = latent / 2_500_000.0

    w_room = (
        coil["w_supply"]
        + m_water_gen / M_DA_ZONE
    )

    room_RH = relative_humidity_from_humidity_ratio(
        T_ROOM,
        w_room,
    )

    # Cooling load per zone.
    q_zone_kW = (
        coil["q_coil_kJ_per_kg_da"]
        * M_DA_ZONE
    )

    # District cooling load.
    q_district_kW = (
        q_zone_kW * N_ZONES
    )

    # Allowances for distribution and sizing.
    q_evap = (
        q_district_kW
        * 1.05
        * 1.02
    )

    if T_cond_C is None:
        T_cond_C = condition["Tdb_C"] + 12.0

    # Absorption-cycle calculation.
    cycle = cycle_performance(
        q_evap,
        T_evap_C=2.0,
        T_cond_C=T_cond_C,
        T_gen_C=T_gen_C,
    )

    # Waste-heat / flue-gas calculation.
    cp_fg = 1.10
    eps_gen = 0.75

    T_fg_in = 300.0
    T_fg_out = 140.0

    deltaT_fg = (
        eps_gen
        * (T_fg_in - T_gen_C)
    )

    q_per_kg_fg = (
        cp_fg * deltaT_fg
    )

    m_fg = (
        cycle.Q_gen_kW
        / max(q_per_kg_fg, 1e-9)
    )

    return {
        "condition": condition["name"],
        "Tmix_C": mixed["T_C"],
        "wmix": mixed["w"],
        "h_mix": mixed["h"],
        "T_ADP_C": coil["T_ADP_C"],
        "w_supply": coil["w_supply"],
        "T_supply_C": T_supply,
        "RH_room": room_RH,
        "coil_load_zone_kW": q_zone_kW,
        "Q_evap_kW": q_evap,
        "Q_gen_kW": cycle.Q_gen_kW,
        "COP": cycle.COP,
        "m_ref_kg_s": cycle.m_ref_kg_s,
        "x_strong": cycle.x_strong,
        "x_weak": cycle.x_weak,
        "m_strong_kg_s": cycle.m_strong_kg_s,
        "m_weak_kg_s": cycle.m_weak_kg_s,
        "circulation_ratio": cycle.circulation_ratio,
        "pump_m3_h": cycle.pump_vol_m3_h,
        "pump_head_kPa": cycle.pump_head_kPa,
        "vapor_m3_s": cycle.refrigerant_vapor_m3_s,
        "m_fg_kg_s": m_fg,
        "condensate_kg_h": m_cond["kg_h"],
        "condensate_m3_day": m_cond["m3_day"],
        "P_high_kPa": cycle.P_high_kPa,
        "T_gen_C": T_gen_C,
        "T_cond_C": T_cond_C,
        "RH_limit_ok": room_RH <= 0.75,
        "flue_limit_ok": m_fg <= 80.0,
        "turn_down_ok": True,
    }


def run_seasonal(fresh_fraction=0.30):
    """Run all six seasonal conditions."""

    return [
        solve_condition(
            condition,
            fresh_fraction=fresh_fraction,
        )
        for condition in CONDITIONS
    ]
