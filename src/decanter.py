"""Preliminary condensate collection calculations."""

from src.psychrometrics import humidity_ratio, mix_air

N_ZONES = 400
ZONES_PER_RISER = 20
M_DA_ZONE = 0.40
BYPASS_FACTOR = 0.15
COIL_ADP_C = 7.0


def condensate_rate_kg_s(m_da, w_in, w_out):
    """Condensate rate from the dry-air moisture balance."""
    return max(0.0, m_da * (w_in - w_out))


def estimate_condensate(
    T_out,
    RH_out,
    T_room=22.0,
    RH_room=0.60,
    fresh_fraction=0.30,
):
    """Estimate condensate using a specified coil ADP."""
    mixed = mix_air(
        T_room,
        RH_room,
        T_out,
        RH_out,
        1.0 - fresh_fraction,
    )

    w_adp = humidity_ratio(COIL_ADP_C, 1.0)
    w_supply = (
        w_adp + BYPASS_FACTOR * (mixed["w"] - w_adp)
    )

    per_zone = condensate_rate_kg_s(
        M_DA_ZONE, mixed["w"], w_supply
    )

    district = per_zone * N_ZONES

    return {
        "T_out_C": T_out,
        "RH_out": RH_out,
        "condensate_zone_kg_h": per_zone * 3600,
        "condensate_district_kg_h": district * 3600,
        "district_m3_day": district * 86400 / 1000,
        "number_of_risers": N_ZONES // ZONES_PER_RISER,
    }


if __name__ == "__main__":
    cases = [
        ("Summer day", 40.0, 0.50),
        ("Summer night", 30.0, 0.50),
        ("Monsoon day", 35.0, 0.85),
        ("Monsoon night", 30.0, 0.85),
        ("Winter day", 28.0, 0.50),
        ("Winter night", 16.0, 0.50),
    ]

    for name, T, RH in cases:
        result = estimate_condensate(T, RH)
        print(name, result)
