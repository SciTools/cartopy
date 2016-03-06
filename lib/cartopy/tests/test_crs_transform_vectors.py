# (C) British Crown Copyright 2013 - 2016, Met Office
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
import warnings

import numpy as np
from numpy.testing import assert_array_almost_equal

import cartopy.crs as ccrs


class TestTransformVectors(unittest.TestCase):

    def test_transform(self):
        # Test some simple vectors to make sure they are transformed
        # correctly.
        rlons = np.array([-90., 0, 90., 180.])
        rlats = np.array([0., 0., 0., 0.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        # transform grid eastward vectors
        ut, vt = target_proj.transform_vectors(src_proj,
                                               rlons,
                                               rlats,
                                               np.ones([4]),
                                               np.zeros([4]))
        assert_array_almost_equal(ut, np.array([0, 1, 0, -1]), decimal=2)
        assert_array_almost_equal(vt, np.array([-1, 0, 1, 0]), decimal=2)
        # transform grid northward vectors
        ut, vt = target_proj.transform_vectors(src_proj,
                                               rlons,
                                               rlats,
                                               np.zeros([4]),
                                               np.ones([4]))
        assert_array_almost_equal(ut, np.array([1, 0, -1, 0]), decimal=2)
        assert_array_almost_equal(vt, np.array([0, 1, 0, -1]), decimal=2)
        # transform grid north-eastward vectors
        ut, vt = target_proj.transform_vectors(src_proj,
                                               rlons,
                                               rlats,
                                               np.ones([4]),
                                               np.ones([4]))
        assert_array_almost_equal(ut, np.array([1, 1, -1, -1]), decimal=2)
        assert_array_almost_equal(vt, np.array([-1, 1, 1, -1]), decimal=2)

    def test_transform_and_inverse(self):
        # Check a full circle transform back to the native projection.
        x = np.arange(-60, 42.5, 2.5)
        y = np.arange(30, 72.5, 2.5)
        x2d, y2d = np.meshgrid(x, y)
        u = np.cos(np.deg2rad(y2d))
        v = np.cos(2. * np.deg2rad(x2d))
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        proj_xyz = target_proj.transform_points(src_proj, x2d, y2d)
        xt, yt = proj_xyz[..., 0], proj_xyz[..., 1]
        ut, vt = target_proj.transform_vectors(src_proj, x2d, y2d, u, v)
        utt, vtt = src_proj.transform_vectors(target_proj, xt, yt, ut, vt)
        assert_array_almost_equal(u, utt, decimal=4)
        assert_array_almost_equal(v, vtt, decimal=4)

    def test_invalid_input_domain(self):
        # If an input coordinate is outside the input projection domain
        # we should be able to handle it correctly.
        rlon = np.array([270.])
        rlat = np.array([0.])
        u = np.array([1.])
        v = np.array([0.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([0]), decimal=2)
        assert_array_almost_equal(vt, np.array([-1]), decimal=2)

    def test_invalid_x_domain(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection x-domain it should be handled correctly as long as
        # it is not a corner point.
        rlon = np.array([180.])
        rlat = np.array([0.])
        u = np.array([1.])
        v = np.array([0.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([-1]), decimal=2)
        assert_array_almost_equal(vt, np.array([0.]), decimal=2)

    def test_invalid_y_domain(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection y-domain it should be handled correctly as long as
        # it is not a corner point.
        rlon = np.array([0.])
        rlat = np.array([90.])
        u = np.array([0.])
        v = np.array([1.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([0.]), decimal=2)
        assert_array_almost_equal(vt, np.array([1.]), decimal=2)

    def test_invalid_xy_domain_corner(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection x and y-domain it should be handled correctly.
        rlon = np.array([180.])
        rlat = np.array([90.])
        u = np.array([1.])
        v = np.array([1.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([0.]), decimal=2)
        assert_array_almost_equal(vt, np.array([-2**.5]), decimal=2)

    def test_invalid_x_domain_corner(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection x-domain and is a corner point, it may be handled
        # incorrectly and a warning should be raised.
        rlon = np.array([180.])
        rlat = np.array([90.])
        u = np.array([1.])
        v = np.array([-1.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with self.assertRaises(UserWarning):
                ut, vt = target_proj.transform_vectors(
                    src_proj, rlon, rlat, u, v)

    def test_invalid_y_domain_corner(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection y-domain and is a corner point, it may be handled
        # incorrectly and a warning should be raised.
        rlon = np.array([180.])
        rlat = np.array([90.])
        u = np.array([-1.])
        v = np.array([1.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with self.assertRaises(UserWarning):
                ut, vt = target_proj.transform_vectors(
                    src_proj, rlon, rlat, u, v)
