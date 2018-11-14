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

from __future__ import (absolute_import, division, print_function)

from datetime import datetime

import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
from cartopy.tests.mpl import ImageTesting


@ImageTesting(['nightshade_platecarree'])
def test_nightshade_image():
    # Test the actual creation of the image
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    dt = datetime(2018, 11, 10, 0, 0)
    ax.set_global()
    ax.add_feature(Nightshade(dt, alpha=0.75))
