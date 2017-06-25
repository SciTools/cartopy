# (C) British Crown Copyright 2013 - 2017, Met Office
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
Tests for the Geostationary projection.

"""

from __future__ import (absolute_import, division, print_function)

from numpy.testing import assert_almost_equal

import cartopy.crs as ccrs


# Note: code here is now shared with the NearsidePerspective test.
def check_proj4_params(crs, expected):
    pro4_params = sorted(crs.proj4_init.split(' +'))
    assert expected == pro4_params


class GeostationaryTestsMixin(object):
    test_class = ccrs.Geostationary
    expected_proj_name = 'geos'

    def adjust_expected_params(self, expected):
        # Only for Geostationary do we expect the sweep parameter
        if self.expected_proj_name == 'geos':
            expected.insert(-3, 'sweep=y')

    def test_default(self):
        geos = self.test_class()
        expected = ['+ellps=WGS84', 'h=35785831', 'lat_0=0.0', 'lon_0=0.0',
                    'no_defs',
                    'proj={}'.format(self.expected_proj_name),
                    'units=m', 'x_0=0', 'y_0=0']
        self.adjust_expected_params(expected)

        check_proj4_params(geos, expected)

        assert_almost_equal(geos.boundary.bounds,
                            (-5434177.81588539, -5434177.81588539,
                             5434177.81588539, 5434177.81588539),
                            decimal=4)

    def test_eccentric_globe(self):
        globe = ccrs.Globe(semimajor_axis=10000, semiminor_axis=5000,
                           ellipse=None)
        geos = self.test_class(satellite_height=50000,
                               globe=globe)
        expected = ['+a=10000', 'b=5000', 'h=50000', 'lat_0=0.0', 'lon_0=0.0',
                    'no_defs',
                    'proj={}'.format(self.expected_proj_name),
                    'units=m', 'x_0=0', 'y_0=0']
        self.adjust_expected_params(expected)

        check_proj4_params(geos, expected)

        assert_almost_equal(geos.boundary.bounds,
                            (-8372.4040, -4171.5043, 8372.4040, 4171.5043),
                            decimal=4)

    def test_eastings(self):
        geos = self.test_class(false_easting=5000000,
                               false_northing=-125000,)
        expected = ['+ellps=WGS84', 'h=35785831', 'lat_0=0.0', 'lon_0=0.0',
                    'no_defs',
                    'proj={}'.format(self.expected_proj_name),
                    'units=m', 'x_0=5000000',
                    'y_0=-125000']
        self.adjust_expected_params(expected)

        check_proj4_params(geos, expected)

        assert_almost_equal(geos.boundary.bounds,
                            (-434177.81588539, -5559177.81588539,
                             10434177.81588539, 5309177.81588539),
                            decimal=4)


class TestGeostationary(GeostationaryTestsMixin, object):
    def test_sweep(self):
        geos = ccrs.Geostationary(sweep_axis='x')
        expected = ['+ellps=WGS84', 'h=35785831', 'lat_0=0.0', 'lon_0=0.0',
                    'no_defs', 'proj=geos', 'sweep=x',
                    'units=m', 'x_0=0', 'y_0=0']

        check_proj4_params(geos, expected)

        pt = geos.transform_point(-60, 25, ccrs.PlateCarree())

        assert_almost_equal(pt,
                            (-4529521.6442, 2437479.4195),
                            decimal=4)
