"""Basic psychrometric tests."""

import unittest

from src.psychrometrics import (
    humidity_ratio,
    moist_air_enthalpy,
    mix_air,
    wet_bulb_temperature,
)


class TestPsychrometrics(unittest.TestCase):

    def test_room_air(self):
        w = humidity_ratio(22.0, 0.75)
        h = moist_air_enthalpy(22.0, w)

        self.assertAlmostEqual(w, 0.01242, delta=0.0001)
        self.assertAlmostEqual(h, 53.71, delta=0.15)

    def test_outdoor_air(self):
        w = humidity_ratio(40.0, 0.50)
        h = moist_air_enthalpy(40.0, w)

        self.assertAlmostEqual(w, 0.02352, delta=0.0001)
        self.assertAlmostEqual(h, 100.82, delta=0.15)

    def test_air_mixing(self):
        result = mix_air(
            22.0, 0.75,
            40.0, 0.50,
            0.70,
        )

        self.assertAlmostEqual(
            result["w"], 0.01575, delta=0.00008
        )

    def test_wet_bulb(self):
        twb = wet_bulb_temperature(40.0, 0.50)

        self.assertLess(twb, 40.0)
        self.assertGreater(twb, 0.0)


if __name__ == "__main__":
    unittest.main()
