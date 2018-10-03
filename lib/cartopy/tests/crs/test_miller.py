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
Tests for the Miller coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=mill', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


def test_default():
    mill = ccrs.Miller()
    other_args = {'a=57.29577951308232', 'lon_0=0.0'}
    check_proj4_params(mill, other_args)

    assert_almost_equal(np.array(mill.x_limits),
                        [-180, 180])
    assert_almost_equal(np.array(mill.y_limits),
                        [-131.9758172, 131.9758172])


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    mill = ccrs.Miller(central_longitude=lon)
    other_args = {'a=57.29577951308232', 'lon_0={}'.format(lon)}
    check_proj4_params(mill, other_args)

    assert_almost_equal(np.array(mill.x_limits),
                        [-180, 180])
    assert_almost_equal(np.array(mill.y_limits),
                        [-131.9758172, 131.9758172])


def test_grid():
    # USGS Professional Paper 1395, p 89, Table 14
    globe = ccrs.Globe(semimajor_axis=1.0, ellipse=None)
    mill = ccrs.Miller(central_longitude=0.0, globe=globe)
    geodetic = mill.as_geodetic()

    other_args = {'a=1.0', 'lon_0=0.0'}
    check_proj4_params(mill, other_args)

    assert_almost_equal(np.array(mill.x_limits),
                        [-3.14159265, 3.14159265])
    assert_almost_equal(np.array(mill.y_limits),
                        [-2.3034125, 2.3034125])

    lats, lons = np.mgrid[0:91:5, 0:91:10].reshape((2, -1))
    expected_x = np.deg2rad(lons)
    expected_y = np.array([
        2.30341, 2.04742, 1.83239, 1.64620, 1.48131, 1.33270, 1.19683, 1.07113,
        0.95364, 0.84284, 0.73754, 0.63674, 0.53962, 0.44547, 0.35369, 0.26373,
        0.17510, 0.08734, 0.00000,
    ])[::-1].repeat(10)

    result = mill.transform_points(geodetic, lons, lats)
    assert_almost_equal(result[:, 0], expected_x, decimal=5)
    assert_almost_equal(result[:, 1], expected_y, decimal=5)


def test_sphere_transform():
    # USGS Professional Paper 1395, pp 287 - 288
    globe = ccrs.Globe(semimajor_axis=1.0, ellipse=None)
    mill = ccrs.Miller(central_longitude=0.0, globe=globe)
    geodetic = mill.as_geodetic()

    other_args = {'a=1.0', 'lon_0=0.0'}
    check_proj4_params(mill, other_args)

    assert_almost_equal(np.array(mill.x_limits),
                        [-3.14159265, 3.14159265])
    assert_almost_equal(np.array(mill.y_limits),
                        [-2.3034125, 2.3034125])

    result = mill.transform_point(-75.0, 50.0, geodetic)
    assert_almost_equal(result, [-1.3089969, 0.9536371])

    inverse_result = geodetic.transform_point(result[0], result[1], mill)
    assert_almost_equal(inverse_result, [-75.0, 50.0])
