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
Tests for the InterruptedGoodeHomolosine coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=igh', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


def test_default():
    igh = ccrs.InterruptedGoodeHomolosine()
    other_args = {'ellps=WGS84', 'lon_0=0'}
    check_proj4_params(igh, other_args)

    assert_almost_equal(np.array(igh.x_limits),
                        [-20037508.3427892, 20037508.3427892])
    assert_almost_equal(np.array(igh.y_limits),
                        [-8683259.7164347, 8683259.7164347])


def test_eccentric_globe():
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                       ellipse=None)
    igh = ccrs.InterruptedGoodeHomolosine(globe=globe)
    other_args = {'a=1000', 'b=500', 'lon_0=0'}
    check_proj4_params(igh, other_args)

    assert_almost_equal(np.array(igh.x_limits),
                        [-3141.5926536, 3141.5926536])
    assert_almost_equal(np.array(igh.y_limits),
                        [-1361.410035, 1361.410035])


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    igh = ccrs.InterruptedGoodeHomolosine(central_longitude=lon)
    other_args = {'ellps=WGS84', 'lon_0={}'.format(lon)}
    check_proj4_params(igh, other_args)

    assert_almost_equal(np.array(igh.x_limits),
                        [-20037508.3427892, 20037508.3427892],
                        decimal=5)
    assert_almost_equal(np.array(igh.y_limits),
                        [-8683259.7164347, 8683259.7164347])
