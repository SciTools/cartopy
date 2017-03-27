# (C) British Crown Copyright 2017, Met Office
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

import cartopy.crs as ccrs
import numpy as np
import unittest
from nose.tools import assert_true


class TestTransformPoints:
    def __init__(self):
        x = np.array([-68, 142.5])
        y = np.array([-72.5, -77.1])

        self.x2d, self.y2d = np.meshgrid(x, y)

        self.src_proj = ccrs.PlateCarree()

    def test_transform_to_orthographic(self):
        target_proj = ccrs.Orthographic()
        proj_xyz = target_proj.transform_points(self.src_proj,
                                                self.x2d, self.y2d)
        # target_lon_0 = (np.max(self.x2d) - np.min(self.x2d))/2
        # target_lat_0 = (np.max(self.y2d) - np.min(self.y2d))/2
        # target_proj.proj4_params.update(lon_0=target_lon_0, lat_0=target_lat_0)
        assert_true(np.inf not in proj_xyz)

    def test_transform_to_transverse_mercator(self):
        target_proj = ccrs.TransverseMercator()
        proj_xyz = target_proj.transform_points(self.src_proj,
                                                self.x2d, self.y2d)
        assert_true(np.inf not in proj_xyz)

    def test_transform_to_gnomonic(self):
        target_proj = ccrs.Gnomonic()
        proj_xyz = target_proj.transform_points(self.src_proj,
                                                self.x2d, self.y2d)
        assert_true(np.inf not in proj_xyz)

    def test_transform_to_geostationary(self):
        target_proj = ccrs.Geostationary()
        proj_xyz = target_proj.transform_points(self.src_proj,
                                                self.x2d, self.y2d)
        assert_true(np.inf not in proj_xyz)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
