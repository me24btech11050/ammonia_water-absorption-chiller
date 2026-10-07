"""Automated validation tests required by the project statement."""

from src.psychrometrics import (
    humidity_ratio_from_rh,
    moist_air_enthalpy,
    adiabatic_mixing,
    wet_bulb_temperature,
)


def test_basic_enthalpy():
    T = 22.0
    w = 0.01242
    h = moist_air_enthalpy(T, w)
    assert abs(h - 53.71) < 0.05


def test_psychrometric_mixing():
    room = adiabatic_mixing(
        22.0, 0.75,
        40.0, 0.50,
        0.70, 0.30
    )

    assert abs(room["w"] - 0.01575) < 0.00008
    assert abs(room["h"] - 67.84) < 0.10
    assert abs(room["T_C"] - 27.48) < 0.10


def test_summer_day_wet_bulb():
    twb = wet_bulb_temperature(40.0, 0.50)
    assert 29.0 < twb < 32.0


def test_monsoon_day_wet_bulb():
    twb = wet_bulb_temperature(35.0, 0.85)
    assert 31.0 < twb < 34.0
