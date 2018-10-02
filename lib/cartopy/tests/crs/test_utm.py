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
Tests for the UTM coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=utm', 'no_defs', 'units=m'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


@pytest.mark.parametrize('south', [False, True])
def test_default(south):
    zone = 1  # Limits are fixed, so don't bother checking other zones.
    utm = ccrs.UTM(zone, southern_hemisphere=south)
    other_args = {'ellps=WGS84', 'zone={}'.format(zone)}
    if south:
        other_args |= {'south'}
    check_proj4_params(utm, other_args)

    assert_almost_equal(np.array(utm.x_limits),
                        [-250000, 1250000])
    assert_almost_equal(np.array(utm.y_limits),
                        [-10000000,  25000000])


def test_ellipsoid_transform():
    # USGS Professional Paper 1395, pp 269 - 271
    globe = ccrs.Globe(ellipse='clrk66')
    utm = ccrs.UTM(zone=18, globe=globe)
    geodetic = utm.as_geodetic()

    other_args = {'ellps=clrk66', 'zone=18'}
    check_proj4_params(utm, other_args)

    assert_almost_equal(np.array(utm.x_limits),
                        [-250000, 1250000])
    assert_almost_equal(np.array(utm.y_limits),
                        [-10000000,  25000000])

    result = utm.transform_point(-73.5, 40.5, geodetic)
    assert_almost_equal(result, np.array([127106.5 + 500000, 4484124.4]),
                        decimal=1)

    inverse_result = geodetic.transform_point(result[0], result[1], utm)
    assert_almost_equal(inverse_result, [-73.5, 40.5])
