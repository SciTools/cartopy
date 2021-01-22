# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
"""
Tests for the InterruptedGoodeHomolosine coordinate system.

"""

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_default():
    igh = ccrs.InterruptedGoodeHomolosine()
    other_args = {'ellps=WGS84', 'lon_0=0'}
    check_proj_params('igh', igh, other_args)

    assert_almost_equal(np.array(igh.x_limits),
                        [-20037508.3427892, 20037508.3427892])
    assert_almost_equal(np.array(igh.y_limits),
                        [-8683259.7164347, 8683259.7164347])


def test_eccentric_globe():
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                       ellipse=None)
    igh = ccrs.InterruptedGoodeHomolosine(globe=globe)
    other_args = {'a=1000', 'b=500', 'lon_0=0'}
    check_proj_params('igh', igh, other_args)

    assert_almost_equal(np.array(igh.x_limits),
                        [-3141.5926536, 3141.5926536])
    assert_almost_equal(np.array(igh.y_limits),
                        [-1361.410035, 1361.410035])


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    igh = ccrs.InterruptedGoodeHomolosine(central_longitude=lon)
    other_args = {'ellps=WGS84', 'lon_0={}'.format(lon)}
    check_proj_params('igh', igh, other_args)

    assert_almost_equal(np.array(igh.x_limits),
                        [-20037508.3427892, 20037508.3427892],
                        decimal=5)
    assert_almost_equal(np.array(igh.y_limits),
                        [-8683259.7164347, 8683259.7164347])
