# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Tests for the HEALPix projection.

"""
import numpy as np
from numpy.testing import assert_almost_equal

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_defaults():
    crs = ccrs.HEALPix()
    expected = {'ellps=WGS84', 'lon_0=0'}
    check_proj_params('healpix', crs, expected)


def test_central_longitude():
    crs = ccrs.HEALPix(central_longitude=124.8)
    expected = {'ellps=WGS84', 'lon_0=124.8'}
    check_proj_params('healpix', crs, expected)


def test_eccentric_globe():
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500, ellipse=None)
    crs = ccrs.HEALPix(globe=globe)
    expected = {'a=1000', 'b=500', 'lon_0=0'}
    check_proj_params('healpix', crs, expected)

    assert_almost_equal(np.array(crs.x_limits), [-3141.5927, 3141.5927],
                        decimal=4)
    assert_almost_equal(np.array(crs.y_limits), [-1570.7963, 1570.7963],
                        decimal=4)
