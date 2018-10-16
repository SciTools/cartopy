# (C) British Crown Copyright 2016 - 2018, Met Office
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

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=sinu', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


class TestSinusoidal(object):
    def test_default(self):
        crs = ccrs.Sinusoidal()
        other_args = {'ellps=WGS84', 'lon_0=0.0', 'x_0=0.0', 'y_0=0.0'}
        check_proj4_params(crs, other_args)

        assert_almost_equal(np.array(crs.x_limits),
                            [-20037508.3428, 20037508.3428],
                            decimal=4)
        assert_almost_equal(np.array(crs.y_limits),
                            [-10001965.7293, 10001965.7293],
                            decimal=4)

    def test_eccentric_globe(self):
        globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                           ellipse=None)
        crs = ccrs.Sinusoidal(globe=globe)
        other_args = {'a=1000', 'b=500', 'lon_0=0.0', 'x_0=0.0', 'y_0=0.0'}
        check_proj4_params(crs, other_args)

        assert_almost_equal(np.array(crs.x_limits),
                            [-3141.59, 3141.59], decimal=2)
        assert_almost_equal(np.array(crs.y_limits),
                            [-1216.60, 1216.60], decimal=2)

    def test_offset(self):
        crs = ccrs.Sinusoidal()
        crs_offset = ccrs.Sinusoidal(false_easting=1234,
                                     false_northing=-4321)
        other_args = {'ellps=WGS84', 'lon_0=0.0', 'x_0=1234', 'y_0=-4321'}
        check_proj4_params(crs_offset, other_args)
        assert tuple(np.array(crs.x_limits) + 1234) == crs_offset.x_limits
        assert tuple(np.array(crs.y_limits) - 4321) == crs_offset.y_limits

    @pytest.mark.parametrize('lon', [-10.0, 10.0])
    def test_central_longitude(self, lon):
        crs = ccrs.Sinusoidal(central_longitude=lon)
        other_args = {'ellps=WGS84', 'lon_0={}'.format(lon),
                      'x_0=0.0', 'y_0=0.0'}
        check_proj4_params(crs, other_args)

        assert_almost_equal(np.array(crs.x_limits),
                            [-20037508.3428, 20037508.3428],
                            decimal=4)
        assert_almost_equal(np.array(crs.y_limits),
                            [-10001965.7293, 10001965.7293],
                            decimal=4)

    def test_MODIS(self):
        # Testpoints verified with MODLAND Tile Calculator
        # http://landweb.nascom.nasa.gov/cgi-bin/developer/tilemap.cgi
        # Settings: Sinusoidal, Global map coordinates, Forward mapping
        crs = ccrs.Sinusoidal.MODIS
        lons = np.array([-180, -50, 40, 180])
        lats = np.array([-89.999, 30, 20, 89.999])
        expected_x = np.array([-349.33, -4814886.99,
                               4179566.79, 349.33])
        expected_y = np.array([-10007443.48, 3335851.56,
                               2223901.04, 10007443.48])
        assert_almost_equal(crs.transform_points(crs.as_geodetic(),
                                                 lons, lats),
                            np.c_[expected_x, expected_y, [0, 0, 0, 0]],
                            decimal=2)
