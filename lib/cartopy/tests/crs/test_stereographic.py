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

import numpy as np
from numpy.testing import assert_almost_equal
from nose.tools import assert_equal

import cartopy.crs as ccrs


class TestStereographic(unittest.TestCase):
    def test_default(self):
        stereo = ccrs.Stereographic()
        expected = ('+ellps=WGS84 +proj=stere +lat_0=0.0 '
                    '+lon_0=0.0 +x_0=0.0 +y_0=0.0 +no_defs')
        assert_equal(stereo.proj4_init, expected)

        assert_almost_equal(np.array(stereo.x_limits),
                            [-5e7, 5e7], decimal=4)
        assert_almost_equal(np.array(stereo.y_limits),
                            [-5e7, 5e7], decimal=4)

    def test_eccentric_globe(self):
        globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                           ellipse=None)
        stereo = ccrs.Stereographic(globe=globe)
        expected = ('+a=1000 +b=500 +proj=stere +lat_0=0.0 +lon_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +no_defs')
        assert_equal(stereo.proj4_init, expected)

        # The limits in this test are sensible values, but are by no means
        # a "correct" answer - they mean that plotting the crs results in a
        # reasonable map.
        assert_almost_equal(np.array(stereo.x_limits),
                            [-7839.27971444, 7839.27971444], decimal=4)
        assert_almost_equal(np.array(stereo.y_limits),
                            [-3932.82587779, 3932.82587779], decimal=4)

    def test_true_scale(self):
        # The "true_scale_latitude" parameter to Stereographic appears
        # meaningless. This test just ensures that the correct proj4
        # string is being created. (#339)
        stereo = ccrs.Stereographic(true_scale_latitude=10)
        expected = ('+ellps=WGS84 +proj=stere +lat_0=0.0 +lon_0=0.0 '
                    '+x_0=0.0 +y_0=0.0 +lat_ts=10 +no_defs')
        assert_equal(stereo.proj4_init, expected)

    def test_eastings(self):
        stereo = ccrs.Stereographic()
        stereo_offset = ccrs.Stereographic(false_easting=1234,
                                           false_northing=-4321)

        expected = ('+ellps=WGS84 +proj=stere +lat_0=0.0 +lon_0=0.0 '
                    '+x_0=1234 +y_0=-4321 +no_defs')
        assert_equal(stereo_offset.proj4_init, expected)
        assert_equal(tuple(np.array(stereo.x_limits) + 1234),
                     stereo_offset.x_limits)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
