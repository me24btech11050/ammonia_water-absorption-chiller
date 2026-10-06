from src.psychrometrics import moist_air_enthalpy


def test_basic_enthalpy():
    T = 22.0
    w = 0.01242

    h = moist_air_enthalpy(T, w)

    expected = 53.71

    assert abs(h - expected) < 0.05
