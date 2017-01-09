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




# TODO: Decide how I want to test this:
# Graphics test? (Probably not).
# Test for nans or infs in transformed arrays? (Maybe too contrived)
# Some kind of test of the geoaxes (like in test_contour.py) - Check mpl for options:
# contains_point()?
# get_visible()?
# has_data()? (unlikely, but worth a try...)

class TestTransformPoints(unittest.TestCase):
    def __init__(self):
        x = np.linspace(0, 360, 60).astype(np.float32)
        y = np.linspace(-90, 90, 90).astype(np.float32)

        self.x2d, self.y2d = np.meshgrid(x, y)

        self.src_proj = ccrs.PlateCarree()
    def test_transform_to_orthographic(self):
        target_proj = ccrs.Orthographic()
        proj_xyz = target_proj.transform_points(self.src_proj,
                                                self.x2d, self.y2d)
        assert_true(np.inf not in proj_xyz)

    # def test_transform_to_transverse_mercator():
    #
    #
    #
    #
    # def test_transform_to_gnomonic():
    #
    #
    #
    # def test_transform_to_geostationary():




if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
