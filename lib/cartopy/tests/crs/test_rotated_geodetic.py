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
Tests for the Transverse Mercator projection, including OSGB and OSNI.

"""

from __future__ import (absolute_import, division, print_function)

import unittest

from numpy.testing import assert_almost_equal
from nose.tools import assert_equal

import cartopy.crs as ccrs


class TestRotatedPole(unittest.TestCase):
    def check_proj4_params(self, crs, expected):
        pro4_params = sorted(crs.proj4_init.split(' +'))
        assert_equal(expected, pro4_params)

    def test_default(self):
        geos = ccrs.RotatedPole(60, 50, 80)
        expected = ['+ellps=WGS84', 'lon_0=240', 'no_defs', 'o_lat_p=50',
                    'o_lon_p=80', 'o_proj=latlon', 'proj=ob_tran',
                    'to_meter=0.0174532925199433']
        self.check_proj4_params(geos, expected)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
