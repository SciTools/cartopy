# (C) British Crown Copyright 2011 - 2018, Met Office
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

from numpy.testing import assert_array_almost_equal
import pytest

import cartopy.crs as ccrs


def check_proj4_params(crs, other_args):
    expected = other_args | {'proj=lcc', 'no_defs'}
    pro4_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == pro4_params


def test_defaults():
    crs = ccrs.LambertConformal()
    other_args = {'ellps=WGS84', 'lon_0=-96.0', 'lat_0=39.0', 'x_0=0.0',
                  'y_0=0.0', 'lat_1=33', 'lat_2=45'}
    check_proj4_params(crs, other_args)


def test_default_with_cutoff():
    crs = ccrs.LambertConformal(cutoff=-80)
    crs2 = ccrs.LambertConformal(cutoff=-80)
    default = ccrs.LambertConformal()

    other_args = {'ellps=WGS84', 'lon_0=-96.0', 'lat_0=39.0', 'x_0=0.0',
                  'y_0=0.0', 'lat_1=33', 'lat_2=45'}
    check_proj4_params(crs, other_args)

    # Check the behaviour of !=, == and (not ==) for the different cutoffs.
    assert crs == crs2
    assert crs != default

    assert hash(crs) != hash(default)
    assert hash(crs) == hash(crs2)

    assert_array_almost_equal(crs.y_limits,
                              (-49788019.81831982, 30793476.08487709))


def test_specific_lambert():
    # This projection comes from EPSG Projection 3034 - ETRS89 / ETRS-LCC.
    crs = ccrs.LambertConformal(central_longitude=10,
                                standard_parallels=(35, 65),
                                central_latitude=52,
                                false_easting=4000000,
                                false_northing=2800000,
                                globe=ccrs.Globe(ellipse='GRS80'))
    other_args = {'ellps=GRS80', 'lon_0=10', 'lat_0=52',
                  'x_0=4000000', 'y_0=2800000', 'lat_1=35', 'lat_2=65'}
    check_proj4_params(crs, other_args)


class Test_LambertConformal_standard_parallels(object):
    def test_single_value(self):
        crs = ccrs.LambertConformal(standard_parallels=[1.])
        other_args = {'ellps=WGS84', 'lon_0=-96.0', 'lat_0=39.0',
                      'x_0=0.0', 'y_0=0.0', 'lat_1=1.0'}
        check_proj4_params(crs, other_args)

    def test_no_parallel(self):
        with pytest.raises(ValueError, message='1 or 2 standard parallels'):
            ccrs.LambertConformal(standard_parallels=[])

    def test_too_many_parallel(self):
        with pytest.raises(ValueError, message='1 or 2 standard parallels'):
            ccrs.LambertConformal(standard_parallels=[1, 2, 3])

    def test_single_spole(self):
        s_pole_crs = ccrs.LambertConformal(standard_parallels=[-1.])
        assert_array_almost_equal(s_pole_crs.x_limits,
                                  (-19840440, 19840440.),
                                  decimal=0)
        assert_array_almost_equal(s_pole_crs.y_limits,
                                  (-370239953, -8191953),
                                  decimal=0)

    def test_single_npole(self):
        n_pole_crs = ccrs.LambertConformal(standard_parallels=[1.])
        assert_array_almost_equal(n_pole_crs.x_limits,
                                  (-20222156, 20222156),
                                  decimal=0)
        assert_array_almost_equal(n_pole_crs.y_limits,
                                  (-8164817, 360848720),
                                  decimal=0)
