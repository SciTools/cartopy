# (C) British Crown Copyright 2013 - 2018, Met Office
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
"""
Tests for the Transverse Mercator projection, including OSGB and OSNI.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

import cartopy.crs as ccrs


class TestTransverseMercator(object):
    def setup_class(self):
        self.point_a = (-3.474083, 50.727301)
        self.point_b = (0.5, 50.5)
        self.src_crs = ccrs.PlateCarree()

    def test_default(self):
        proj = ccrs.TransverseMercator()
        res = proj.transform_point(*self.point_a, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(res, (-245269.53180633,
                                                   5627508.74354959))
        res = proj.transform_point(*self.point_b, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(res, (35474.63566645,
                                                   5596583.41949901))

    def test_osgb_vals(self):
        proj = ccrs.TransverseMercator(central_longitude=-2,
                                       central_latitude=49,
                                       scale_factor=0.9996012717,
                                       false_easting=400000,
                                       false_northing=-100000,
                                       globe=ccrs.Globe(datum='OSGB36',
                                                        ellipse='airy'))
        res = proj.transform_point(*self.point_a, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(res, (295971.28667707,
                                                   93064.27666368))
        res = proj.transform_point(*self.point_b, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(res, (577274.98380140,
                                                   69740.49227181))

    def test_nan(self):
        proj = ccrs.TransverseMercator()
        res = proj.transform_point(0.0, float('nan'), src_crs=self.src_crs)
        assert np.all(np.isnan(res))
        res = proj.transform_point(float('nan'), 0.0, src_crs=self.src_crs)
        assert np.all(np.isnan(res))


class TestOSGB(object):
    def setup_class(self):
        self.point_a = (-3.474083, 50.727301)
        self.point_b = (0.5, 50.5)
        self.src_crs = ccrs.PlateCarree()
        self.nan = float('nan')

    def test_default(self):
        proj = ccrs.OSGB()
        res = proj.transform_point(*self.point_a, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(res, (295971.28667707,
                                                   93064.27666368))
        res = proj.transform_point(*self.point_b, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(res, (577274.98380140,
                                                   69740.49227181))

    def test_nan(self):
        proj = ccrs.OSGB()
        res = proj.transform_point(0.0, float('nan'), src_crs=self.src_crs)
        assert np.all(np.isnan(res))
        res = proj.transform_point(float('nan'), 0.0, src_crs=self.src_crs)
        assert np.all(np.isnan(res))


class TestOSNI(object):
    def setup_class(self):
        self.point_a = (-6.826286, 54.725116)
        self.src_crs = ccrs.PlateCarree()
        self.nan = float('nan')

    def test_default(self):
        proj = ccrs.OSNI()
        res = proj.transform_point(*self.point_a, src_crs=self.src_crs)
        np.testing.assert_array_almost_equal(
            res, (275614.26762651594, 386984.206429612),
            decimal=0 if ccrs.PROJ4_VERSION < (5, 0, 0) else 6)

    def test_nan(self):
        proj = ccrs.OSNI()
        res = proj.transform_point(0.0, float('nan'), src_crs=self.src_crs)
        assert np.all(np.isnan(res))
        res = proj.transform_point(float('nan'), 0.0, src_crs=self.src_crs)
        assert np.all(np.isnan(res))
