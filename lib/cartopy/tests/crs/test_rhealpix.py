# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Tests for the RHEALPix projection.

"""
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_defaults():
    crs = ccrs.RHEALPix()
    expected = {'ellps=WGS84', 'lon_0=0', 'north_square=0', 'south_square=0'}
    check_proj_params('rhealpix', crs, expected)


def test_square_positions():
    crs = ccrs.RHEALPix(north_square=1, south_square=2)
    expected = {'ellps=WGS84', 'lon_0=0', 'north_square=1', 'south_square=2'}
    check_proj_params('rhealpix', crs, expected)


def test_invalid_square_positions():
    with pytest.raises(ValueError, match='north_square must be'):
        ccrs.RHEALPix(north_square=4)
    with pytest.raises(ValueError, match='south_square must be'):
        ccrs.RHEALPix(south_square=-1)


def test_central_longitude():
    crs = ccrs.RHEALPix(north_square=2, central_longitude=-124.8)
    expected = {'ellps=WGS84', 'lon_0=-124.8', 'north_square=2',
                'south_square=0'}
    check_proj_params('rhealpix', crs, expected)


def test_eccentric_globe():
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500, ellipse=None)
    crs = ccrs.RHEALPix(globe=globe)
    expected = {'a=1000', 'b=500', 'lon_0=0', 'north_square=0',
                'south_square=0'}
    check_proj_params('rhealpix', crs, expected)

    assert_almost_equal(np.array(crs.x_limits), [-3141.5927, 3141.5927],
                        decimal=4)
    assert_almost_equal(np.array(crs.y_limits), [-2356.1945, 2356.1945],
                        decimal=4)
