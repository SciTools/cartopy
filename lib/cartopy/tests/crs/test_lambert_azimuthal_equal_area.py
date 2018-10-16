# (C) British Crown Copyright 2015 - 2018, Met Office
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
    expected = other_args | {'proj=laea', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


class TestLambertAzimuthalEqualArea(object):
    def test_default(self):
        crs = ccrs.LambertAzimuthalEqualArea()
        other_args = {'ellps=WGS84', 'lon_0=0.0', 'lat_0=0.0', 'x_0=0.0',
                      'y_0=0.0'}
        check_proj4_params(crs, other_args)

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
        other_args = {'a=1000', 'b=500', 'lon_0=0.0', 'lat_0=0.0', 'x_0=0.0',
                      'y_0=0.0'}
        check_proj4_params(crs, other_args)

        assert_almost_equal(np.array(crs.x_limits),
                            [-1999.9, 1999.9], decimal=1)
        assert_almost_equal(np.array(crs.y_limits),
                            [-1380.17298647, 1380.17298647], decimal=4)

    def test_offset(self):
        crs = ccrs.LambertAzimuthalEqualArea()
        crs_offset = ccrs.LambertAzimuthalEqualArea(false_easting=1234,
                                                    false_northing=-4321)
        other_args = {'ellps=WGS84', 'lon_0=0.0', 'lat_0=0.0', 'x_0=1234',
                      'y_0=-4321'}
        check_proj4_params(crs_offset, other_args)
        assert tuple(np.array(crs.x_limits) + 1234) == crs_offset.x_limits
        assert tuple(np.array(crs.y_limits) - 4321) == crs_offset.y_limits

    @pytest.mark.parametrize("latitude", [-90, 90])
    def test_extrema(self, latitude):
        crs = ccrs.LambertAzimuthalEqualArea(central_latitude=latitude)
        other_args = {'ellps=WGS84', 'lon_0=0.0', 'lat_0={}'.format(latitude),
                      'x_0=0.0', 'y_0=0.0'}
        check_proj4_params(crs, other_args)
