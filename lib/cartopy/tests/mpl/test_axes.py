# (C) British Crown Copyright 2011 - 2018, Met Office
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

from matplotlib.testing.decorators import cleanup
import matplotlib.path as mpath
import matplotlib.pyplot as plt
try:
    from unittest import mock
except ImportError:
    import mock
import numpy as np
import pytest

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import InterProjectionTransform, GeoAxes
from cartopy.tests.mpl.test_caching import CallCounter


class TestNoSpherical(object):
    def setup_method(self):
        self.ax = plt.axes(projection=ccrs.PlateCarree())
        self.data = np.arange(12).reshape((3, 4))

    def teardown_method(self):
        plt.clf()
        plt.close()

    def test_contour(self):
        with pytest.raises(ValueError):
            self.ax.contour(self.data, transform=ccrs.Geodetic())

    def test_contourf(self):
        with pytest.raises(ValueError):
            self.ax.contourf(self.data, transform=ccrs.Geodetic())

    def test_pcolor(self):
        with pytest.raises(ValueError):
            self.ax.pcolor(self.data, transform=ccrs.Geodetic())

    def test_pcolormesh(self):
        with pytest.raises(ValueError):
            self.ax.pcolormesh(self.data, transform=ccrs.Geodetic())


def test_transform_PlateCarree_shortcut():
    src = ccrs.PlateCarree(central_longitude=0)
    target = ccrs.PlateCarree(central_longitude=180)

    # of the 3 paths, 2 of them cannot be short-cutted.
    pth1 = mpath.Path([[0.5, 0], [10, 10]])
    pth2 = mpath.Path([[0.5, 91], [10, 10]])
    pth3 = mpath.Path([[-0.5, 0], [10, 10]])

    trans = InterProjectionTransform(src, target)

    counter = CallCounter(target, 'project_geometry')

    with counter:
        trans.transform_path(pth1)
        # pth1 should allow a short-cut.
        assert counter.count == 0

    with counter:
        trans.transform_path(pth2)
        assert counter.count == 1

    with counter:
        trans.transform_path(pth3)
        assert counter.count == 2


class Test_InterProjectionTransform():
    def pc_2_pc(self):
        return InterProjectionTransform(
            ccrs.PlateCarree(), ccrs.PlateCarree())

    def pc_2_rob(self):
        return InterProjectionTransform(ccrs.PlateCarree(), ccrs.Robinson())

    def rob_2_rob_shifted(self):
        return InterProjectionTransform(
            ccrs.Robinson(), ccrs.Robinson(central_longitude=0))

    def test_eq(self):
        assert self.pc_2_pc() == self.pc_2_pc()
        assert self.pc_2_rob() == self.pc_2_rob()
        assert self.rob_2_rob_shifted() == self.rob_2_rob_shifted()

        assert not self.pc_2_rob() == self.rob_2_rob_shifted()
        assert not self.pc_2_pc() == 'not a transform obj'

    def test_ne(self):
        assert not self.pc_2_pc() != self.pc_2_pc()
        print(self.pc_2_pc() != self.pc_2_rob())
        assert self.pc_2_pc() != self.pc_2_rob()


class Test_Axes_add_geometries():
    def teardown_method(self):
        plt.close()

    @mock.patch('cartopy.mpl.geoaxes.GeoAxes.add_feature')
    @mock.patch('cartopy.feature.ShapelyFeature')
    def test_styler_kwarg(self, ShapelyFeature, add_feature_method):
        ax = GeoAxes(plt.figure(), [0, 0, 1, 1],
                     map_projection=ccrs.Robinson())
        ax.add_geometries(mock.sentinel.geometries, mock.sentinel.crs,
                          styler=mock.sentinel.styler, wibble='wobble')

        ShapelyFeature.assert_called_once_with(
            mock.sentinel.geometries, mock.sentinel.crs, wibble='wobble')

        add_feature_method.assert_called_once_with(
            ShapelyFeature(), styler=mock.sentinel.styler)


@cleanup
def test_geoaxes_subplot():
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    assert str(ax.__class__) == "<class 'cartopy.mpl.geoaxes.GeoAxesSubplot'>"
