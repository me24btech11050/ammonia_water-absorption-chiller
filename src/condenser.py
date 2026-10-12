"""Task 2: seasonal condenser comparison and preliminary sizing."""

import csv
import math
from pathlib import Path

from src.psychrometrics import wet_bulb_temperature
from src.absorption_cycle import (
    solve_cycle,
    required_generator_temperature,
)

U_AIR = 30.0
U_EVAP = 600.0

AMBIENT = {
    "Summer day": (40.0, 0.50),
    "Summer night": (30.0, 0.50),
    "Monsoon day": (35.0, 0.85),
    "Monsoon night": (30.0, 0.85),
    "Winter day": (28.0, 0.50),
    "Winter night": (16.0, 0.50),
}


def calculate_lmtd(dT1, dT2):
    if dT1 <= 0 or dT2 <= 0:
        raise ValueError("Temperature differences must be positive.")

    if math.isclose(dT1, dT2, rel_tol=1e-9):
        return dT1

    return (dT1 - dT2) / math.log(dT1 / dT2)


def size_condenser(
    Q_kW,
    sink_in_C,
    approach_K,
    U_W_m2K,
    cp_kJ_kgK,
    sink_rise_K,
):
    if Q_kW <= 0:
        raise ValueError("Heat duty must be positive.")

    if approach_K <= sink_rise_K:
        raise ValueError("Approach must exceed sink temperature rise.")

    T_cond = sink_in_C + approach_K
    T_out = sink_in_C + sink_rise_K

    LMTD = calculate_lmtd(
        T_cond - sink_in_C,
        T_cond - T_out,
    )

    UA = Q_kW * 1000.0 / LMTD
    area = UA / U_W_m2K
    mass_flow = Q_kW / (cp_kJ_kgK * sink_rise_K)

    return {
        "T_cond_C": T_cond,
        "T_sink_out_C": T_out,
        "LMTD_K": LMTD,
        "UA_W_K": UA,
        "area_m2": area,
        "cooling_flow_kg_s": mass_flow,
    }


def compare_condensers(Q_evap_kW=5000.0):
    """Compare both condenser types at the six ambient conditions."""
    rows = []

    for name, (T_db, RH) in AMBIENT.items():
        T_wb = wet_bulb_temperature(T_db, RH)

        options = {
            "air": size_condenser(
                Q_evap_kW, T_db, 12.0, U_AIR, 1.006, 10.0
            ),
            "evaporative": size_condenser(
                Q_evap_kW, T_wb, 8.0, U_EVAP, 4.186, 5.0
            ),
        }

        for condenser_type, design in options.items():
            try:
                cycle = solve_cycle(
                    Q_evap_kW,
                    T_cond_C=design["T_cond_C"],
                )

                design.update({
                    "Q_gen_kW": cycle.Q_gen_kW,
                    "Q_cond_kW": cycle.Q_cond_kW,
                    "COP": cycle.COP,
                    "P_high_bar": cycle.P_high_bar,
                    "T_gen_required_C":
                        required_generator_temperature(cycle.P_high_bar),
                    "pressure_ok": cycle.P_high_bar <= 20.0,
                    "minimum_Tcond_ok":
                        design["T_cond_C"] >= 25.0,
                })

            except (ValueError, OverflowError):
                design.update({
                    "Q_gen_kW": None,
                    "Q_cond_kW": None,
                    "COP": None,
                    "P_high_bar": None,
                    "T_gen_required_C": None,
                    "pressure_ok": False,
                    "minimum_Tcond_ok":
                        design["T_cond_C"] >= 25.0,
                })

            rows.append({
                "condition": name,
                "condenser_type": condenser_type,
                "T_db_C": T_db,
                "RH_percent": RH * 100.0,
                "T_wb_C": T_wb,
                **design,
            })

    return rows


def save_results(rows, filename="results/tables/condenser_task2.csv"):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved results to {path}")


if __name__ == "__main__":
    rows = compare_condensers()
    save_results(rows)

    for row in rows:
        print(
            row["condition"],
            row["condenser_type"],
            "area =", round(row["area_m2"], 2), "m2",
            "LMTD =", round(row["LMTD_K"], 2), "K",
        )
