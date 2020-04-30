# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import io
from unittest import mock

import matplotlib.pyplot as plt
import numpy as np
import pytest

import cartopy.crs as ccrs
from cartopy.tests.mpl import MPL_VERSION


def test_pcolormesh_fully_masked():
    data = np.ma.masked_all((30, 40))

    # Check that a fully masked data array doesn't trigger a pcolor call.
    with mock.patch('cartopy.mpl.geoaxes.GeoAxes.pcolor') as pcolor:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.pcolormesh(np.linspace(-90, 90, 40), np.linspace(0, 360, 30), data)
        assert pcolor.call_count == 0, ("pcolor shouldn't have been called, "
                                        "but was.")
        plt.close()


def test_pcolormesh_partially_masked():
    data = np.ma.masked_all((30, 40))
    data[0:100] = 10

    # Check that a partially masked data array does trigger a pcolor call.
    with mock.patch('cartopy.mpl.geoaxes.GeoAxes.pcolor') as pcolor:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.pcolormesh(np.linspace(-90, 90, 40), np.linspace(0, 360, 30), data)
        assert pcolor.call_count == 1, ("pcolor should have been called "
                                        "exactly once.")
        plt.close()


def test_pcolormesh_invisible():
    data = np.zeros((3, 3))

    # Check that a fully invisible mesh doesn't fail.
    with mock.patch('cartopy.mpl.geoaxes.GeoAxes.pcolor') as pcolor:
        ax = plt.axes(projection=ccrs.Orthographic())
        ax.pcolormesh(np.linspace(-75, 75, 3), np.linspace(105, 255, 3), data,
                      transform=ccrs.PlateCarree())
        assert pcolor.call_count == 0, ("pcolor shouldn't have been called, "
                                        "but was.")
        plt.close()


@pytest.mark.skipif(MPL_VERSION < '2.1.0', reason='Matplotlib is broken.')
def test_savefig_tight():
    nx, ny = 36, 18
    xbnds = np.linspace(0, 360, nx, endpoint=True)
    ybnds = np.linspace(-90, 90, ny, endpoint=True)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    plt.subplot(211, projection=ccrs.Robinson())
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()
