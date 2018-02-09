# (C) British Crown Copyright 2018, Met Office
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
"""
Tests for the Equidistant Conic coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal

import cartopy.crs as ccrs


class TestEquidistantConic(object):
    def test_default(self):
        eqdc = ccrs.EquidistantConic()
        expected = ('+ellps=WGS84 +proj=eqdc +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +lat_1=20.0 +lat_2=50.0 +no_defs')
        assert eqdc.proj4_init == expected

        assert_almost_equal(np.array(eqdc.x_limits),
                            (-22784919.35600352, 22784919.35600352),
                            decimal=7)
        assert_almost_equal(np.array(eqdc.y_limits),
                            (-10001965.729313632, 17558791.85156368),
                            decimal=7)

    def test_eccentric_globe(self):
        globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                           ellipse=None)
        eqdc = ccrs.EquidistantConic(globe=globe)
        expected = ('+a=1000 +b=500 +proj=eqdc +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +lat_1=20.0 +lat_2=50.0 +no_defs')
        assert eqdc.proj4_init == expected

        assert_almost_equal(np.array(eqdc.x_limits),
                            (-3016.869847713461, 3016.869847713461),
                            decimal=7)
        assert_almost_equal(np.array(eqdc.y_limits),
                            (-1216.6029342241113, 2511.0574375797723),
                            decimal=7)

    def test_eastings(self):
        eqdc_offset = ccrs.EquidistantConic(false_easting=1234,
                                            false_northing=-4321)

        expected = ('+ellps=WGS84 +proj=eqdc +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=1234 +y_0=-4321 +lat_1=20.0 +lat_2=50.0 +no_defs')
        assert eqdc_offset.proj4_init == expected

    def test_standard_parallels(self):
        eqdc = ccrs.EquidistantConic(standard_parallels=(13, 37))
        expected = ('+ellps=WGS84 +proj=eqdc +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +lat_1=13 +lat_2=37 +no_defs')
        assert eqdc.proj4_init == expected

        eqdc = ccrs.EquidistantConic(standard_parallels=(13, ))
        expected = ('+ellps=WGS84 +proj=eqdc +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +lat_1=13 +no_defs')
        assert eqdc.proj4_init == expected

        eqdc = ccrs.EquidistantConic(standard_parallels=13)
        expected = ('+ellps=WGS84 +proj=eqdc +lon_0=0.0 +lat_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +lat_1=13 +no_defs')
        assert eqdc.proj4_init == expected

    def test_sphere_transform(self):
        # USGS Professional Paper 1395, pg 298
        globe = ccrs.Globe(semimajor_axis=1.0, semiminor_axis=1.0,
                           ellipse=None)
        lat_1 = 29.5
        lat_2 = 45.5
        eqdc = ccrs.EquidistantConic(central_longitude=-96.0,
                                     central_latitude=23,
                                     standard_parallels=(lat_1, lat_2),
                                     globe=globe)
        geodetic = eqdc.as_geodetic()

        expected = ('+a=1.0 +b=1.0 +proj=eqdc +lon_0=-96.0 +lat_0=23 '
                    '+x_0=0.0 +y_0=0.0 +lat_1=29.5 +lat_2=45.5 +no_defs')
        assert eqdc.proj4_init == expected

        assert_almost_equal(np.array(eqdc.x_limits),
                            (-3.520038619089038, 3.520038619089038),
                            decimal=7)
        assert_almost_equal(np.array(eqdc.y_limits),
                            (-1.9722220547535922, 2.7066811021065535),
                            decimal=7)

        result = eqdc.transform_point(-75.0, 35.0, geodetic)

        assert_almost_equal(result, (0.2952057, 0.2424021), decimal=7)

    def test_ellipsoid_transform(self):
        # USGS Professional Paper 1395, p9 299--300
        globe = ccrs.Globe(semimajor_axis=6378206.4,
                           flattening=1 - np.sqrt(1 - 0.00676866),
                           ellipse=None)
        lat_1 = 29.5
        lat_2 = 45.5
        eqdc = ccrs.EquidistantConic(central_latitude=23.0,
                                     central_longitude=-96.0,
                                     standard_parallels=(lat_1, lat_2),
                                     globe=globe)
        geodetic = eqdc.as_geodetic()

        expected = ('+a=6378206.4 +f=0.003390076308689371 +proj=eqdc '
                    '+lon_0=-96.0 +lat_0=23.0 +x_0=0.0 +y_0=0.0 '
                    '+lat_1=29.5 +lat_2=45.5 +no_defs')
        assert eqdc.proj4_init == expected

        assert_almost_equal(np.array(eqdc.x_limits),
                            (-22421870.719894886, 22421870.719894886),
                            decimal=7)
        assert_almost_equal(np.array(eqdc.y_limits),
                            (-12546277.778958388, 17260638.403203618),
                            decimal=7)

        result = eqdc.transform_point(-75.0, 35.0, geodetic)

        assert_almost_equal(result, (1885051.9, 1540507.6), decimal=1)
