"""Task 1: seasonal district cooling-load estimates."""

from src.psychrometrics import (
    humidity_ratio,
    moist_air_enthalpy,
    mix_air,
    coil_from_bypass,
)

N_ZONES = 400
M_DA_ZONE = 0.40

T_ROOM = 22.0
RH_ROOM = 0.60
T_SUPPLY = 12.0
BYPASS_FACTOR = 0.15

AMBIENT = {
    "Summer day": (40.0, 0.50, True),
    "Summer night": (30.0, 0.50, False),
    "Monsoon day": (35.0, 0.85, True),
    "Monsoon night": (30.0, 0.85, False),
    "Winter day": (28.0, 0.50, True),
    "Winter night": (16.0, 0.50, False),
}


def zone_load(T_out, RH_out, solar_on, fresh_fraction=0.30):
    mixed = mix_air(
        T_ROOM, RH_ROOM, T_out, RH_out,
        1.0 - fresh_fraction,
    )

    coil = coil_from_bypass(
        mixed["T_C"],
        mixed["w"],
        T_SUPPLY,
        BYPASS_FACTOR,
    )

    h_mix = moist_air_enthalpy(mixed["T_C"], mixed["w"])
    h_supply = moist_air_enthalpy(
        T_SUPPLY, coil["w_supply"]
    )

    Q_coil_zone = max(
        0.0, M_DA_ZONE * (h_mix - h_supply)
    )

    return {
        "T_mix_C": mixed["T_C"],
        "ADP_C": coil["ADP_C"],
        "w_supply": coil["w_supply"],
        "Q_coil_zone_kW": Q_coil_zone,
        "Q_evap_district_kW": N_ZONES * Q_coil_zone,
        "fresh_fraction": fresh_fraction,
        "solar_on": solar_on,
    }


def run_task1(fresh_fraction=0.30):
    results = []

    for name, (T, RH, solar) in AMBIENT.items():
        result = zone_load(
            T, RH, solar, fresh_fraction
        )

        result["condition"] = name
        result["T_out_C"] = T
        result["RH_out_percent"] = RH * 100.0

        results.append(result)

    return results


if __name__ == "__main__":
    for result in run_task1():
        print(result)
