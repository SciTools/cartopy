# (C) British Crown Copyright 2014 - 2018, Met Office
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
Tests for the Rotated Geodetic coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=ob_tran', 'o_proj=latlon',
                             'to_meter=0.0174532925199433', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


def test_default():
    geos = ccrs.RotatedGeodetic(30, 15, 27)
    other_args = {'datum=WGS84', 'ellps=WGS84', 'lon_0=210', 'o_lat_p=15',
                  'o_lon_p=27'}
    check_proj4_params(geos, other_args)
