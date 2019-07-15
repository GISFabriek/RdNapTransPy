# ***********************************************************************
# Author           : Willem A. Ligtendag, De GISFabriek
# Created          : 07-14-2019
#
# Last Modified By : Willem A. Ligtendag, De GISFabriek
# Last Modified On : 07-15-2019
# ***********************************************************************
# Python PORT from C version of RDNAPTRANS
# ***********************************************************************
import unittest
import transformer

DeltaRdPlaces = 3
DeltaAnglePlaces = 8
DeltaHPlaces = 3


class TestRdNapTrans(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._items = []
        cls._items.append(("Texel", (53.160753042, 4.824761912, 42.8614), (117380.1200, 575040.3400, 1.0000)))
        cls._items.append(("Noord-Groningen", (53.419482050, 6.776726674, 42.3586), (247380.5600, 604580.7800, 2.0000)))
        cls._items.append(("Amersfoort", (52.155172897, 5.387203657, 43.2551), (155000.0000, 463000.0000, 0.0000)))
        cls._items.append(("Amersfoort 100m", (52.155172910, 5.387203658, 143.2551),
                         (155000.0000, 463000.0000, 100.0000)))
        cls._items.append(("Zeeuws-Vlaanderen", (51.368607152, 3.397588595, 47.4024), (16460.9100, 377380.2300, 3.0000)))
        cls._items.append(("Zuid-Limburg", (50.792584916, 5.773795548, 245.9478), (182260.4500, 311480.6700, 200.0000)))
        cls._items.append(("Maasvlakte", (51.947393898, 4.072887101, 47.5968), (64640.8900, 440700.0101, 4.0000)))
        cls._items.append(("outside", (48.843030210, 8.723260235, 52.0289), (400000.2300, 100000.4500, 5.0000)))
        cls._items.append(("no_rd&geoid", (50.687420392, 4.608971813, 51.6108), (100000.6700, 300000.8900, 6.0000)))
        cls._items.append(("no_geoid", (51.136825197, 4.601375361, 50.9672), (100000.6700, 350000.8900, 6.0000)))
        cls._items.append(("no_rd", (52.482440839, 4.268403889, 49.9436), (79000.0100, 500000.2300, 7.0000)))
        cls._items.append(("edge_rd", (51.003976532, 3.891247830, 52.7427), (50000.4500, 335999.6700, 8.0000)))

    def test_etrs2rdnap(self):
        for item in self._items:
            name, geograpic, cartesian = item
            phi, lmbd, h = geograpic
            expected_x, expected_y, expected_nap = cartesian
            x_rd, y_rd, nap = transformer.etrs2rdnap(phi, lmbd, h)
            self.assertAlmostEqual(x_rd, expected_x, DeltaRdPlaces)
            self.assertAlmostEqual(y_rd, expected_y, DeltaRdPlaces)
            self.assertAlmostEqual(nap, expected_nap, DeltaHPlaces)

    def test_rdnap2etrs(self):
        for item in self._items:
            name, geograpic, cartesian = item
            x, y, z = cartesian
            expected_phi, expected_lmbd, expected_h = geograpic
            phi, lmbd, h = transformer.rdnap2etrs(x, y, z)
            self.assertAlmostEqual(phi, expected_phi, DeltaAnglePlaces)
            self.assertAlmostEqual(lmbd, expected_lmbd, DeltaRdPlaces)
            self.assertAlmostEqual(h, expected_h, DeltaHPlaces)


if __name__ == '__main__':
    unittest.main()
