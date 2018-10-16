# (C) British Crown Copyright 2016 - 2018, Met Office
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
Tests for the NearsidePerspective projection.

"""

from __future__ import (absolute_import, division, print_function)

from numpy.testing import assert_almost_equal

from cartopy.tests.crs.test_geostationary import (GeostationaryTestsMixin,
                                                  check_proj4_params)

from cartopy.crs import NearsidePerspective


class TestEquatorialDefault(GeostationaryTestsMixin, object):
    # Check that it behaves just like Geostationary, in the absence of a
    # central_latitude parameter.
    test_class = NearsidePerspective
    expected_proj_name = 'nsper'


class TestOwnSpecifics(object):
    def test_central_latitude(self):
        # Check the effect of the added 'central_latitude' key.
        geos = NearsidePerspective(central_latitude=53.7)
        other_args = {'ellps=WGS84', 'h=35785831', 'lat_0=53.7', 'lon_0=0.0',
                      'x_0=0', 'y_0=0'}
        check_proj4_params('nsper', geos, other_args)

        assert_almost_equal(geos.boundary.bounds,
                            (-5434177.81588539, -5434177.81588539,
                             5434177.81588539, 5434177.81588539),
                            decimal=4)
