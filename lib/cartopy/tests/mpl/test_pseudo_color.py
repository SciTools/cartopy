# (C) British Crown Copyright 2013 - 2017, Met Office
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

import io
import matplotlib.pyplot as plt
try:
    from unittest import mock
except ImportError:
    import mock
from nose.tools import assert_equal
import numpy as np

import cartopy.crs as ccrs


def test_pcolormesh_fully_masked():
    data = np.ma.masked_all((30, 40))

    # Check that a fully masked data array doesn't trigger a pcolor call.
    with mock.patch('cartopy.mpl.geoaxes.GeoAxes.pcolor') as pcolor:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.pcolormesh(np.linspace(-90, 90, 40), np.linspace(0, 360, 30), data)
        assert_equal(pcolor.call_count, 0, ("pcolor shouldn't have been "
                                            "called, but was."))
        plt.close()


def test_pcolormesh_partially_masked():
    data = np.ma.masked_all((30, 40))
    data[0:100] = 10

    # Check that a partially masked data array does trigger a pcolor call.
    with mock.patch('cartopy.mpl.geoaxes.GeoAxes.pcolor') as pcolor:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.pcolormesh(np.linspace(-90, 90, 40), np.linspace(0, 360, 30), data)
        assert_equal(pcolor.call_count, 1, ("pcolor should have been "
                                            "called exactly once."))
        plt.close()


def test_pcolormesh_invisible():
    data = np.zeros((3, 3))

    # Check that a fully invisible mesh doesn't fail.
    with mock.patch('cartopy.mpl.geoaxes.GeoAxes.pcolor') as pcolor:
        ax = plt.axes(projection=ccrs.Orthographic())
        ax.pcolormesh(np.linspace(-75, 75, 3), np.linspace(105, 255, 3), data,
                      transform=ccrs.PlateCarree())
        assert_equal(pcolor.call_count, 0, ("pcolor shouldn't have been "
                                            "called, but was."))
        plt.close()


def test_savefig_tight():
    nx, ny = 36, 18
    xbnds = np.linspace(0, 360, nx, endpoint=True)
    ybnds = np.linspace(-90, 90, ny, endpoint=True)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    ax = plt.subplot(211, projection=ccrs.Robinson())
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
