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
Tests for RectangularHealpix projection.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
import pytest

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_defaults():
    crs = ccrs.RectangularHealpix()
    expected = {'ellps=WGS84', 'lon_0=0', 'north_square=0', 'south_square=0'}
    check_proj_params('rhealpix', crs, expected)


def test_square_positions():
    crs = ccrs.RectangularHealpix(north_square=1, south_square=2)
    expected = {'ellps=WGS84', 'lon_0=0', 'north_square=1', 'south_square=2'}
    check_proj_params('rhealpix', crs, expected)


def test_central_longitude():
    crs = ccrs.RectangularHealpix(north_square=2, central_longitude=-124.8)
    expected = {'ellps=WGS84', 'lon_0=-124.8', 'north_square=2', 'south_square=0'}
    check_proj_params('rhealpix', crs, expected)
