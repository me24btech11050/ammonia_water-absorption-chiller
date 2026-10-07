"""Condenser sizing and off-design routines."""

import math
from dataclasses import dataclass
from .psychrometrics import saturation_pressure_water, wet_bulb_temperature
from .absorption_cycle import cycle_performance


@dataclass
class CondenserDesign:
    kind: str
    governing_condition: str
    approach_K: float
    U_W_m2K: float
    area_m2: float
    UA_W_K: float
    heat_duty_kW: float


def lmtd(delta1, delta2):
    if delta1 <= 0 or delta2 <= 0:
        raise ValueError("Temperature differences must be positive.")
    if abs(delta1 - delta2) < 1e-9:
        return delta1
    return (delta1 - delta2) / math.log(delta1 / delta2)


def compare_condenser_types(conditions, Q_by_condition):
    rows = []
    for c in conditions:
        twb = wet_bulb_temperature(c["Tdb_C"], c["RH"])
        q = Q_by_condition[c["name"]]
        # Representative sizing approaches requested by the project.
        Tcond_air = c["Tdb_C"] + 12.0
        Tcond_evap = twb + 8.0

        dt_air = max(Tcond_air - c["Tdb_C"], 1.0)
        dt_evap = max(Tcond_evap - twb, 1.0)

        rows.append({
            "condition": c["name"],
            "Tdb_C": c["Tdb_C"],
            "Twb_C": twb,
            "Q_kW": q,
            "A_air_m2": q * 1000.0 / (30.0 * dt_air),
            "A_evap_m2": q * 1000.0 / (600.0 * dt_evap),
        })
    return rows


def size_condenser(kind, governing, Q_kW, Tdb_C, Twb_C,
                   approach_K):
    if kind == "air":
        U = 30.0
        T_sink = Tdb_C
    elif kind == "evaporative":
        U = 600.0
        T_sink = Twb_C
    else:
        raise ValueError("kind must be 'air' or 'evaporative'")

    Tcond = T_sink + approach_K
    # Simple condensing-side LMTD with an assumed 5 K cooling-medium rise.
    dt1 = Tcond - T_sink
    dt2 = Tcond - (T_sink + 5.0)
    if dt2 <= 0:
        dt2 = dt1
    L = lmtd(dt1, dt2)
    UA = Q_kW * 1000.0 / L
    A = UA / U

    return CondenserDesign(
        kind=kind,
        governing_condition=governing,
        approach_K=approach_K,
        U_W_m2K=U,
        area_m2=A,
        UA_W_K=UA,
        heat_duty_kW=Q_kW,
    )
