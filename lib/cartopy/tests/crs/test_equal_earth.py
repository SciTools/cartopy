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
Tests for the Equal Earth coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


pytestmark = pytest.mark.skipif(ccrs.PROJ4_VERSION < (5, 2, 0),
                                reason='Proj is too old.')


def check_proj_params(crs, other_args):
    expected = other_args | {'proj=eqearth', 'no_defs'}
    proj_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == proj_params


def test_default():
    eqearth = ccrs.EqualEarth()
    other_args = {'ellps=WGS84', 'lon_0=0'}
    check_proj_params(eqearth, other_args)

    assert_almost_equal(eqearth.x_limits,
                        [-17243959.0622169, 17243959.0622169])
    assert_almost_equal(eqearth.y_limits,
                        [-8392927.59846646, 8392927.59846646])
    # Expected aspect ratio from the paper.
    assert_almost_equal(np.diff(eqearth.x_limits) / np.diff(eqearth.y_limits),
                        2.05458, decimal=5)


def test_offset():
    crs = ccrs.EqualEarth()
    crs_offset = ccrs.EqualEarth(false_easting=1234, false_northing=-4321)
    other_args = {'ellps=WGS84', 'lon_0=0', 'x_0=1234', 'y_0=-4321'}
    check_proj_params(crs_offset, other_args)
    assert tuple(np.array(crs.x_limits) + 1234) == crs_offset.x_limits
    assert tuple(np.array(crs.y_limits) - 4321) == crs_offset.y_limits


def test_eccentric_globe():
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                       ellipse=None)
    eqearth = ccrs.EqualEarth(globe=globe)
    other_args = {'a=1000', 'b=500', 'lon_0=0'}
    check_proj_params(eqearth, other_args)

    assert_almost_equal(eqearth.x_limits,
                        [-2248.43664092550, 2248.43664092550])
    assert_almost_equal(eqearth.y_limits,
                        [-1094.35228122148, 1094.35228122148])
    # Expected aspect ratio from the paper.
    assert_almost_equal(np.diff(eqearth.x_limits) / np.diff(eqearth.y_limits),
                        2.05458, decimal=5)


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    eqearth = ccrs.EqualEarth(central_longitude=lon)
    other_args = {'ellps=WGS84', 'lon_0={}'.format(lon)}
    check_proj_params(eqearth, other_args)

    assert_almost_equal(eqearth.x_limits,
                        [-17243959.0622169, 17243959.0622169], decimal=5)
    assert_almost_equal(eqearth.y_limits,
                        [-8392927.59846646, 8392927.59846646])
    # Expected aspect ratio from the paper.
    assert_almost_equal(np.diff(eqearth.x_limits) / np.diff(eqearth.y_limits),
                        2.05458, decimal=5)
