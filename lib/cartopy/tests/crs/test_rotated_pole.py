# (C) British Crown Copyright 2014 - 2016, Met Office
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
Tests for the Rotated Geodetic coordinate system.

"""

from __future__ import (absolute_import, division, print_function)

import unittest

from numpy.testing import assert_almost_equal
from nose.tools import assert_equal

import cartopy.crs as ccrs


class TestRotatedGeodetic(unittest.TestCase):
    def check_proj4_params(self, crs, expected):
        pro4_params = sorted(crs.proj4_init.split(' +'))
        assert_equal(expected, pro4_params)

    def test_default(self):
        geos = ccrs.RotatedGeodetic(30, 15, 27)
        expected = ['+datum=WGS84', 'ellps=WGS84', 'lon_0=210', 'no_defs',
                    'o_lat_p=15', 'o_lon_p=27', 'o_proj=latlon',
                    'proj=ob_tran', 'to_meter=0.0174532925199433']
        self.check_proj4_params(geos, expected)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
