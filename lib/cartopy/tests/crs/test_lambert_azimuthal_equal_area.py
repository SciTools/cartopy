# (C) British Crown Copyright 2015 - 2016, Met Office
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
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

import unittest

import numpy as np
from numpy.testing import assert_almost_equal
from nose.tools import assert_equal

import cartopy.crs as ccrs


class TestLambertAzimuthalEqualArea(unittest.TestCase):
    def test_default(self):
        crs = ccrs.LambertAzimuthalEqualArea()
        expected = ('+ellps=WGS84 +proj=laea +lon_0=0.0 '
                    '+lat_0=0.0 +x_0=0.0 +y_0=0.0 +no_defs')
        assert_equal(crs.proj4_init, expected)

        assert_almost_equal(np.array(crs.x_limits),
                            [-12755636.1863, 12755636.1863],
                            decimal=4)
        assert_almost_equal(np.array(crs.y_limits),
                            [-12727770.598700099, 12727770.598700099],
                            decimal=4)

    def test_eccentric_globe(self):
        globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                           ellipse=None)
        crs = ccrs.LambertAzimuthalEqualArea(globe=globe)
        expected = ('+a=1000 +b=500 +proj=laea +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +no_defs')
        assert_equal(crs.proj4_init, expected)

        assert_almost_equal(np.array(crs.x_limits),
                            [-1999.9, 1999.9], decimal=1)
        assert_almost_equal(np.array(crs.y_limits),
                            [-1380.17298647, 1380.17298647], decimal=4)

    def test_offset(self):
        crs = ccrs.LambertAzimuthalEqualArea()
        crs_offset = ccrs.LambertAzimuthalEqualArea(false_easting=1234,
                                                    false_northing=-4321)
        expected = ('+ellps=WGS84 +proj=laea +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=1234 +y_0=-4321 +no_defs')
        assert_equal(crs_offset.proj4_init, expected)
        assert_equal(tuple(np.array(crs.x_limits) + 1234),
                     crs_offset.x_limits)
        assert_equal(tuple(np.array(crs.y_limits) - 4321),
                     crs_offset.y_limits)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
