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
Tests for the Transverse Mercator projection, including OSGB and OSNI.

"""

from __future__ import (absolute_import, division, print_function)

import cartopy.crs as ccrs
from .helpers import check_proj_params


common_other_args = {'o_proj=latlon', 'to_meter=0.0174532925199433'}


def test_default():
    geos = ccrs.RotatedPole(60, 50, 80)
    other_args = common_other_args | {'ellps=WGS84', 'lon_0=240', 'o_lat_p=50',
                                      'o_lon_p=80'}
    check_proj_params('ob_tran', geos, other_args)
