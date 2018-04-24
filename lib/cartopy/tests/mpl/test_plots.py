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

from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs


def test_empty_plot():
    """Test making a plot with empty arrays."""
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Mercator())
    ax.plot([], [], transform=ccrs.PlateCarree())
    fig.savefig(BytesIO())


def test_triplot_bbox_tight():
    """Test triplot with a tight bbox (#1060)."""
    x = np.degrees([-0.101, -0.090, -0.069])
    y = np.degrees([0.872, 0.883, 0.888])
    triangles = np.asarray([[0, 1, 2]])

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.OSGB())
    ax.triplot(x, y, triangles, transform=ccrs.Geodetic())
    fig.savefig(BytesIO(), bbox_inches='tight')
