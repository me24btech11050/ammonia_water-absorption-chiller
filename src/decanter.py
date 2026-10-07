"""Condensate/decanter calculations."""

def condensate_rate_per_zone(m_da_kg_s, w_mix, w_supply):
    """kg/s of condensate at the cooling coil."""
    return max(m_da_kg_s * (w_mix - w_supply), 0.0)


def district_condensate(m_da_kg_s, w_mix, w_supply, n_zones=400):
    mdot = condensate_rate_per_zone(m_da_kg_s, w_mix, w_supply)
    kg_h = mdot * 3600.0 * n_zones
    m3_day = kg_h * 24.0 / 1000.0
    return {"kg_h": kg_h, "m3_day": m3_day}


def riser_count(n_zones=400, zones_per_riser=20):
    return (n_zones + zones_per_riser - 1) // zones_per_riser
