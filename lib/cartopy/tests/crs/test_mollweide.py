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
Tests for the Mollweide coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=moll', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


def test_default():
    moll = ccrs.Mollweide()
    other_args = {'a=6378137.0', 'lon_0=0'}
    check_proj4_params(moll, other_args)

    assert_almost_equal(np.array(moll.x_limits),
                        [-18040095.6961473, 18040095.6961473])
    assert_almost_equal(np.array(moll.y_limits),
                        [-9020047.8480736, 9020047.8480736])


def test_offset():
    crs = ccrs.Mollweide()
    crs_offset = ccrs.Mollweide(false_easting=1234, false_northing=-4321)
    other_args = {'a=6378137.0', 'lon_0=0', 'x_0=1234', 'y_0=-4321'}
    check_proj4_params(crs_offset, other_args)
    assert tuple(np.array(crs.x_limits) + 1234) == crs_offset.x_limits
    assert tuple(np.array(crs.y_limits) - 4321) == crs_offset.y_limits


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    moll = ccrs.Mollweide(central_longitude=lon)
    other_args = {'a=6378137.0', 'lon_0={}'.format(lon)}
    check_proj4_params(moll, other_args)

    assert_almost_equal(np.array(moll.x_limits),
                        [-18040095.6961473, 18040095.6961473],
                        decimal=5)
    assert_almost_equal(np.array(moll.y_limits),
                        [-9020047.8480736, 9020047.8480736])


def test_grid():
    # USGS Professional Paper 1395, pg 252, Table 42
    globe = ccrs.Globe(ellipse=None,
                       semimajor_axis=0.5**0.5, semiminor_axis=0.5**0.5)
    moll = ccrs.Mollweide(globe=globe)
    geodetic = moll.as_geodetic()

    other_args = {'a=0.7071067811865476', 'b=0.7071067811865476', 'lon_0=0'}
    check_proj4_params(moll, other_args)

    assert_almost_equal(np.array(moll.x_limits),
                        [-2, 2])
    assert_almost_equal(np.array(moll.y_limits),
                        [-1, 1])

    lats = np.arange(0, 91, 5)[::-1]
    lons = np.full_like(lats, 90)
    result = moll.transform_points(geodetic, lons, lats)

    expected_x = np.array([
        0.00000, 0.20684, 0.32593, 0.42316, 0.50706, 0.58111, 0.64712, 0.70617,
        0.75894, 0.80591, 0.84739, 0.88362, 0.91477, 0.94096, 0.96229, 0.97882,
        0.99060, 0.99765, 1.00000,
    ])
    assert_almost_equal(result[:, 0], expected_x, decimal=5)

    expected_y = np.array([
        1.00000, 0.97837, 0.94539, 0.90606, 0.86191, 0.81382, 0.76239, 0.70804,
        0.65116, 0.59204, 0.53097, 0.46820, 0.40397, 0.33850, 0.27201, 0.20472,
        0.13681, 0.06851, 0.00000,
    ])
    assert_almost_equal(result[:, 1], expected_y, decimal=5)


def test_sphere_transform():
    # USGS Professional Paper 1395, pg 367
    globe = ccrs.Globe(semimajor_axis=1.0, semiminor_axis=1.0,
                       ellipse=None)
    moll = ccrs.Mollweide(central_longitude=-90.0,
                          globe=globe)
    geodetic = moll.as_geodetic()

    other_args = {'a=1.0', 'b=1.0', 'lon_0=-90.0'}
    check_proj4_params(moll, other_args)

    assert_almost_equal(np.array(moll.x_limits),
                        [-2.8284271247461903, 2.8284271247461903],
                        decimal=2)
    assert_almost_equal(np.array(moll.y_limits),
                        [-1.4142135623730951, 1.4142135623730951])

    result = moll.transform_point(-75.0, -50.0, geodetic)
    assert_almost_equal(result, [0.1788845, -0.9208758])

    inverse_result = geodetic.transform_point(result[0], result[1], moll)
    assert_almost_equal(inverse_result, [-75.0, -50.0])
