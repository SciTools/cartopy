# (C) British Crown Copyright 2013, Met Office
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
"""
Tests for the Transverse Mercator projection, including OSGB and OSNI.

"""

import unittest

from numpy.testing import assert_almost_equal
from nose.tools import assert_equal

import cartopy.crs as ccrs


class TestGeostationary(unittest.TestCase):
    def check_proj4_params(self, crs, expected):
        pro4_params = sorted(crs.proj4_init.split(' +'))
        assert_equal(expected, pro4_params)

    def test_default(self):
        geos = ccrs.Geostationary()
        expected = ['+ellps=WGS84', 'h=35785831', 'lat_0=0', 'lon_0=0.0',
                    'no_defs', 'proj=geos', 'units=m', 'x_0=0', 'y_0=0']
        self.check_proj4_params(geos, expected)

        assert_almost_equal(geos.boundary.bounds,
                            (-5364970.19679699, -5370680.80015303,
                             5372584.78443894, 5370680.80015303),
                            decimal=4)

    def test_eccentric_globe(self):
        globe = ccrs.Globe(semimajor_axis=10000, semiminor_axis=5000,
                           ellipse=None)
        geos = ccrs.Geostationary(satellite_height=50000,
                                  globe=globe)
        expected = ['+a=10000', 'b=5000', 'h=50000', 'lat_0=0', 'lon_0=0.0',
                    'no_defs', 'proj=geos', 'units=m', 'x_0=0', 'y_0=0']
        self.check_proj4_params(geos, expected)

        assert_almost_equal(geos.boundary.bounds,
                            (-8245.7306, -4531.3879, 8257.4339, 4531.3879),
                            decimal=4)

    def test_eastings(self):
        geos = ccrs.Geostationary(false_easting=5000000,
                                  false_northing=-125000,)
        expected = ['+ellps=WGS84', 'h=35785831', 'lat_0=0', 'lon_0=0.0',
                    'no_defs', 'proj=geos', 'units=m', 'x_0=5000000',
                    'y_0=-125000']
        self.check_proj4_params(geos, expected)

        assert_almost_equal(geos.boundary.bounds,
                            (-364970.19679699, -5495680.80015303,
                             10372584.78443894, 5245680.80015303),
                            decimal=4)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
