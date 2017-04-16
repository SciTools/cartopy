# (C) British Crown Copyright 2015, Met Office
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
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

from nose.tools import raises

import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs


class TestQuiverShapes(object):

    @classmethod
    def setup_class(cls):
        cls.x = np.arange(-60, 42.5, 2.5)
        cls.y = np.arange(30, 72.5, 2.5)
        cls.x2d, cls.y2d = np.meshgrid(cls.x, cls.y)
        cls.u = np.cos(np.deg2rad(cls.y2d))
        cls.v = np.cos(2. * np.deg2rad(cls.x2d))
        cls.rp = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
        cls.pc = ccrs.PlateCarree()

    def test_quiver_transform_xyuv_1d(self):
        fig, ax = plt.subplots(subplot_kw=dict(projection=self.pc))
        ax.quiver(self.x2d.ravel(), self.y2d.ravel(),
                  self.u.ravel(), self.v.ravel(), transform=self.rp)

    def test_quiver_transform_xyuv_2d(self):
        fig, ax = plt.subplots(subplot_kw=dict(projection=self.pc))
        ax.quiver(self.x2d, self.y2d, self.u, self.v, transform=self.rp)

    def test_quiver_transform_xy_1d_uv_2d(self):
        fig, ax = plt.subplots(subplot_kw=dict(projection=self.pc))
        ax.quiver(self.x, self.y, self.u, self.v, transform=self.rp)

    @raises(ValueError)
    def test_quiver_transform_xy_2d_uv_1d(self):
        fig, ax = plt.subplots(subplot_kw=dict(projection=self.pc))
        ax.quiver(self.x2d, self.y2d,
                  self.u.ravel(), self.v.ravel(), transform=self.rp)

    @raises(ValueError)
    def test_quiver_transform_inconsistent_shape(self):
        fig, ax = plt.subplots(subplot_kw=dict(projection=self.pc))
        ax.quiver(self.x, self.y,
                  self.u.ravel(), self.v.ravel(), transform=self.rp)
