# (C) British Crown Copyright 2011 - 2013, Met Office
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

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal, assert_not_equal, assert_true

import cartopy.crs as ccrs


def test_defaults():
    crs = ccrs.LambertConformal()
    assert_equal(crs.proj4_init, ('+ellps=WGS84 +proj=lcc +lon_0=-96.0 '
                                  '+lat_0=39.0 +x_0=0.0 +y_0=0.0 +lat_1=33 '
                                  '+lat_2=45 +no_defs'))


def test_default_with_cutoff():
    crs = ccrs.LambertConformal(cutoff=-80)
    crs2 = ccrs.LambertConformal(cutoff=-80)
    default = ccrs.LambertConformal()

    assert_equal(crs.proj4_init, ('+ellps=WGS84 +proj=lcc +lon_0=-96.0 '
                                  '+lat_0=39.0 +x_0=0.0 +y_0=0.0 +lat_1=33 '
                                  '+lat_2=45 +no_defs'))

    # Check the behaviour of !=, == and (not ==) for the different cutoffs.
    assert_equal(crs, crs2)
    assert_true(crs != default)
    assert_not_equal(crs, default)

    assert_not_equal(hash(crs), hash(default))
    assert_equal(hash(crs), hash(crs2))

    assert_array_almost_equal(crs.y_limits,
                              (-49788019.81822971, 30793476.084826108))


def test_specific_lambert():
    # This projection comes from EPSG Projection 3034 - ETRS89 / ETRS-LCC.
    crs = ccrs.LambertConformal(central_longitude=10,
                                secant_latitudes=(35, 65),
                                central_latitude=52,
                                false_easting=4000000,
                                false_northing=2800000,
                                globe=ccrs.Globe(ellipse='GRS80'))
    assert_equal(crs.proj4_init, ('+ellps=GRS80 +proj=lcc +lon_0=10 '
                                  '+lat_0=52 +x_0=4000000 +y_0=2800000 '
                                  '+lat_1=35 +lat_2=65 +no_defs'))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
