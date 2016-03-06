# (C) British Crown Copyright 2011 - 2016, Met Office
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

import unittest


from matplotlib.testing.decorators import cleanup
import matplotlib.path as mpath
import matplotlib.pyplot as plt
from nose.tools import assert_equal
import numpy as np


import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import InterProjectionTransform
from .test_caching import CallCounter


class TestNoSpherical(unittest.TestCase):
    def setUp(self):
        self.ax = plt.axes(projection=ccrs.PlateCarree())
        self.data = np.arange(12).reshape((3, 4))

    def tearDown(self):
        plt.clf()
        plt.close()

    def test_contour(self):
        with self.assertRaises(ValueError):
            self.ax.contour(self.data, transform=ccrs.Geodetic())

    def test_contourf(self):
        with self.assertRaises(ValueError):
            self.ax.contourf(self.data, transform=ccrs.Geodetic())

    def test_pcolor(self):
        with self.assertRaises(ValueError):
            self.ax.pcolor(self.data, transform=ccrs.Geodetic())

    def test_pcolormesh(self):
        with self.assertRaises(ValueError):
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
        assert_equal(counter.count, 0)

    with counter:
        trans.transform_path(pth2)
        assert_equal(counter.count, 1)

    with counter:
        trans.transform_path(pth3)
        assert_equal(counter.count, 2)


@cleanup
def test_geoaxes_subplot():
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    assert_equal(str(ax.__class__),
                 "<class 'cartopy.mpl.geoaxes.GeoAxesSubplot'>")


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
