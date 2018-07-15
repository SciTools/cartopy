# (C) British Crown Copyright 2016 - 2018, Met Office
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

import matplotlib.pyplot as plt
from matplotlib.testing.decorators import cleanup
import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy.interpolate import NearestNDInterpolator
from scipy.signal import convolve2d

import cartopy.crs as ccrs


@cleanup
def test_contour_plot_bounds():
    x = np.linspace(-2763217.0, 2681906.0, 200)
    y = np.linspace(-263790.62, 3230840.5, 130)
    data = np.hypot(*np.meshgrid(x, y)) / 2e5

    proj_lcc = ccrs.LambertConformal(central_longitude=-95,
                                     central_latitude=25,
                                     standard_parallels=[25])
    ax = plt.axes(projection=proj_lcc)
    ax.contourf(x, y, data, levels=np.arange(0, 40, 1))
    assert_array_almost_equal(ax.get_extent(),
                              np.array([x[0], x[-1], y[0], y[-1]]))


@cleanup
def test_contour_linear_ring():
    """Test contourf with a section that only has 3 points."""
    ax = plt.axes([0.01, 0.05, 0.898, 0.85], projection=ccrs.Mercator(),
                  aspect='equal')
    ax.set_extent([-99.6, -89.0, 39.8, 45.5])

    xbnds = ax.get_xlim()
    ybnds = ax.get_ylim()
    ll = ccrs.Geodetic().transform_point(xbnds[0], ybnds[0], ax.projection)
    ul = ccrs.Geodetic().transform_point(xbnds[0], ybnds[1], ax.projection)
    ur = ccrs.Geodetic().transform_point(xbnds[1], ybnds[1], ax.projection)
    lr = ccrs.Geodetic().transform_point(xbnds[1], ybnds[0], ax.projection)
    xi = np.linspace(min(ll[0], ul[0]), max(lr[0], ur[0]), 100)
    yi = np.linspace(min(ll[1], ul[1]), max(ul[1], ur[1]), 100)
    xi, yi = np.meshgrid(xi, yi)
    nn = NearestNDInterpolator((np.arange(-94, -85), np.arange(36, 45)),
                               np.arange(9))
    vals = nn(xi, yi)
    lons = xi
    lats = yi
    window = np.ones((6, 6))
    vals = convolve2d(vals, window / window.sum(), mode='same',
                      boundary='symm')
    ax.contourf(lons, lats, vals, np.arange(9), transform=ccrs.PlateCarree())

    plt.draw()
