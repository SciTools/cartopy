# (C) British Crown Copyright 2013 - 2018, Met Office
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

from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=merc', 'units=m', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


def test_default():
    crs = ccrs.Mercator()

    other_args = {'ellps=WGS84', 'lon_0=0.0', 'x_0=0.0', 'y_0=0.0'}
    check_proj4_params(crs, other_args)
    assert_almost_equal(crs.boundary.bounds,
                        [-20037508, -15496571, 20037508, 18764656], decimal=0)


def test_eccentric_globe():
    globe = ccrs.Globe(semimajor_axis=10000, semiminor_axis=5000,
                       ellipse=None)
    crs = ccrs.Mercator(globe=globe, min_latitude=-40, max_latitude=40)
    other_args = {'a=10000', 'b=5000', 'lon_0=0.0', 'x_0=0.0', 'y_0=0.0'}
    check_proj4_params(crs, other_args)

    assert_almost_equal(crs.boundary.bounds,
                        [-31415.93, -2190.5, 31415.93, 2190.5], decimal=2)

    assert_almost_equal(crs.x_limits, [-31415.93, 31415.93], decimal=2)
    assert_almost_equal(crs.y_limits, [-2190.5, 2190.5], decimal=2)


def test_equality():
    default = ccrs.Mercator()
    crs = ccrs.Mercator(min_latitude=0)
    crs2 = ccrs.Mercator(min_latitude=0)

    # Check the == and != operators.
    assert crs == crs2
    assert crs != default
    assert hash(crs) != hash(default)
    assert hash(crs) == hash(crs2)


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    crs = ccrs.Mercator(central_longitude=lon)
    other_args = {'ellps=WGS84', 'lon_0={}'.format(lon), 'x_0=0.0', 'y_0=0.0'}
    check_proj4_params(crs, other_args)

    assert_almost_equal(crs.boundary.bounds,
                        [-20037508, -15496570, 20037508, 18764656], decimal=0)


def test_latitude_true_scale():
    lat_ts = 20.0
    crs = ccrs.Mercator(latitude_true_scale=lat_ts)
    other_args = {'ellps=WGS84', 'lon_0=0.0', 'x_0=0.0', 'y_0=0.0',
                  'lat_ts={}'.format(lat_ts)}
    check_proj4_params(crs, other_args)

    assert_almost_equal(crs.boundary.bounds,
                        [-18836475, -14567718, 18836475, 17639917], decimal=0)


def test_easting_northing():
    false_easting = 1000000
    false_northing = -2000000
    crs = ccrs.Mercator(false_easting=false_easting,
                        false_northing=false_northing)
    other_args = {'ellps=WGS84', 'lon_0=0.0', 'x_0={}'.format(false_easting),
                  'y_0={}'.format(false_northing)}
    check_proj4_params(crs, other_args)

    assert_almost_equal(crs.boundary.bounds,
                        [-19037508, -17496571, 21037508, 16764656], decimal=0)


def test_scale_factor():
    # Should be same as lat_ts=20 for a sphere
    scale_factor = 0.939692620786
    crs = ccrs.Mercator(scale_factor=scale_factor,
                        globe=ccrs.Globe(ellipse='sphere'))
    other_args = {'ellps=sphere', 'lon_0=0.0', 'x_0=0.0', 'y_0=0.0',
                  'k_0={:.12f}'.format(scale_factor)}
    check_proj4_params(crs, other_args)

    assert_almost_equal(crs.boundary.bounds,
                        [-18808021, -14585266, 18808021, 17653216], decimal=0)
