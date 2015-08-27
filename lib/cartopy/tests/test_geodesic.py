# (C) British Crown Copyright 2015, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

import unittest

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal
from nose.tools import assert_equal

from cartopy import geodesic


class TestGeodesic(unittest.TestCase):

    def setUp(self):
        """
        Data sampled from GeographicLib Test Data for Geodesics at:
        http://geographiclib.sourceforge.net/html/geodesic.html#testgeod

        """
        self.geod = geodesic.Geodesic()

        # Fill a 10 by 7 numpy array with starting lons, lats, azimuths; ending
        # lons, lats and azimuths and distances to travel.

        self.data = np.array([[0.0000000000, 36.5300423550, 176.1258751622,
                          5.7623446947, -48.1642707791, 175.3343083163,
                          9398502.0434687007],
                         [0.0000000000, 20.8766024619, 6.9012827094,
                          163.9792202999, 64.2764863397, 165.0440144913,
                          10462971.2273696996],
                         [0.0000000000, 59.7405712203, 80.9569174535,
                          80.1969954660, 30.9857449391, 144.4488137288,
                          6549489.1863671001],
                         [0.0000000000, 38.6508883588, 18.3455177945,
                          23.5931524958, 66.3457305181, 37.7145989984,
                          3425212.4767990001],
                         [0.0000000000, 23.2214345509, 165.5720618611,
                          148.3625110902, -68.8453788967, 39.2692310682,
                          14506511.2971898001],
                         [0.0000000000, 31.2989275984, 155.7723493796,
                          93.8764112107, -69.2776346668, 98.5250397385,
                          13370814.5013951007],
                         [0.0000000000, 49.6823298563, 1.0175398481,
                          5.3554086646, 83.8681965431, 6.1667605618,
                          3815028.2543704999],
                         [0.0000000000, 32.7651878215, 98.6494285944,
                          70.3527194957, 2.4777491770, 123.5999412794,
                          8030520.7178932996],
                         [0.0000000000, 46.3648067071, 94.9148631993,
                          56.5676529172, 25.2581951337, 130.4405565458,
                          5485075.9286326999],
                         [0.0000000000, 33.7321188396, 147.9041907517,
                          33.1346935645, -26.3211288531, 150.4502346224,
                          7512675.5414637001]])

    def test_dir(self):
        geod_dir = self.geod.direct(self.data[:, :2], self.data[:, 2],
                                                      self.data[:, 6])
        assert_array_almost_equal(geod_dir, self.data[:, 3:6], decimal=5)

    def test_inverse(self):
        geod_inv = self.geod.inverse(self.data[:, :2], self.data[:, 3:5])
        assert_array_almost_equal(geod_inv, self.data[:, [6, 2, 5]], decimal=5)

    def test_circle(self):
        geod_circle = self.geod.circle(40, 50, 500000, n_samples=3)
        assert_almost_equal(geod_circle,
                            np.array([[40., 54.49349757],
                                      [34.23766162, 47.60355349],
                                      [45.76233838, 47.60355349]]), decimal=5)

    def test_repr(self):
	expected = '<Geodesic: radius=6378137.000, flattening=1/298.257>'
	assert_equal(expected, repr(self.geod))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
