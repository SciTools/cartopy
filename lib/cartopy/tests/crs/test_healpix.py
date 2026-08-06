# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Tests for the HEALPix projection.

"""
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
