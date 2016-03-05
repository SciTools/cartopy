# (C) British Crown Copyright 2015 - 2016, Met Office
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
try:
    from unittest import mock
except ImportError:
    import mock

from nose.tools import raises, assert_equal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import cleanup

import cartopy.crs as ccrs


# Note, other tests for quiver exist in test_mpl_integration.

class TestQuiverShapes(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(-60, 42.5, 10)
        self.y = np.linspace(30, 72.5, 7)
        self.x2d, self.y2d = np.meshgrid(self.x, self.y)
        self.u = np.cos(np.deg2rad(self.y2d))
        self.v = np.cos(2. * np.deg2rad(self.x2d))
        self.rp = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
        self.pc = ccrs.PlateCarree()
        self.ax = plt.axes(projection=self.pc)

    @cleanup
    def test_quiver_transform_xyuv_1d(self):
        with mock.patch('matplotlib.axes.Axes.quiver') as patch:
            self.ax.quiver(self.x2d.ravel(), self.y2d.ravel(),
                           self.u.ravel(), self.v.ravel(), transform=self.rp)
        args, kwargs = patch.call_args
        assert_equal(len(args), 5)
        assert_equal(sorted(kwargs.keys()), [u'transform'])
        shapes = [arg.shape for arg in args[1:]]
        # Assert that all the shapes have been broadcast.
        assert_equal(shapes, [(70, )] * 4)

    @cleanup
    def test_quiver_transform_xy_1d_uv_2d(self):
        with mock.patch('matplotlib.axes.Axes.quiver') as patch:
            self.ax.quiver(self.x, self.y, self.u, self.v, transform=self.rp)
        args, kwargs = patch.call_args
        assert_equal(len(args), 5)
        assert_equal(sorted(kwargs.keys()), [u'transform'])
        shapes = [arg.shape for arg in args[1:]]
        # Assert that all the shapes have been broadcast.
        assert_equal(shapes, [(7, 10)] * 4)

    @raises(ValueError)
    def test_quiver_transform_xy_2d_uv_1d(self):
        self.ax.quiver(self.x2d, self.y2d,
                       self.u.ravel(), self.v.ravel(), transform=self.rp)

    @raises(ValueError)
    def test_quiver_transform_inconsistent_shape(self):
        self.ax.quiver(self.x, self.y,
                       self.u.ravel(), self.v.ravel(), transform=self.rp)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
