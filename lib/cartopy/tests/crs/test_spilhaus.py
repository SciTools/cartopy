# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Tests for the HEALPix projection.
"""
import pytest

import cartopy.crs as ccrs
from .helpers import check_proj_params


common_arg = {
    'azi=40.17823482',
    'datum=WGS84',
    'ellps=WGS84',
    'k_0=1.414213562373095',
    'lat_0=-49.56371678',
    'lon_0=66.94970198',
    'no_defs',
    'units=m',
}
def test_defaults():
    crs = ccrs.Spilhaus()
    expected = {'rot=45','x_0=0.0','y_0=0.0'} | common_arg
    check_proj_params('spilhaus', crs, expected)

@pytest.mark.parametrize("orientation",[0,1,2,3])
def test_greenland_at(orientation):
    crs = ccrs.Spilhaus(greenland_at = orientation)
    expected = {f'rot={orientation*90+45}','x_0=0.0','y_0=0.0'} | common_arg
    check_proj_params('spilhaus', crs, expected)

@pytest.mark.parametrize("orientation",['some random string',4,2.5])
def test_disallowed_orientation(orientation):
    with pytest.raises(ValueError):
        ccrs.Spilhaus(greenland_at = orientation)
