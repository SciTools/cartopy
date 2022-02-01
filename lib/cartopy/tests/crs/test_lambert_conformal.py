# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from numpy.testing import assert_array_almost_equal
import pytest

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_defaults():
    crs = ccrs.LambertConformal()
    other_args = {'ellps=WGS84', 'lon_0=-96.0', 'lat_0=39.0', 'x_0=0.0',
                  'y_0=0.0', 'lat_1=33', 'lat_2=45'}
    check_proj_params('lcc', crs, other_args)


def test_default_with_cutoff():
    crs = ccrs.LambertConformal(cutoff=-80)
    crs2 = ccrs.LambertConformal(cutoff=-80)
    default = ccrs.LambertConformal()

    other_args = {'ellps=WGS84', 'lon_0=-96.0', 'lat_0=39.0', 'x_0=0.0',
                  'y_0=0.0', 'lat_1=33', 'lat_2=45'}
    check_proj_params('lcc', crs, other_args)

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
    check_proj_params('lcc', crs, other_args)


class Test_LambertConformal_standard_parallels:
    def test_single_value(self):
        crs = ccrs.LambertConformal(standard_parallels=[1.])
        other_args = {'ellps=WGS84', 'lon_0=-96.0', 'lat_0=39.0',
                      'x_0=0.0', 'y_0=0.0', 'lat_1=1.0'}
        check_proj_params('lcc', crs, other_args)

    def test_no_parallel(self):
        with pytest.raises(ValueError, match='1 or 2 standard parallels'):
            ccrs.LambertConformal(standard_parallels=[])

    def test_too_many_parallel(self):
        with pytest.raises(ValueError, match='1 or 2 standard parallels'):
            ccrs.LambertConformal(standard_parallels=[1, 2, 3])

    def test_single_spole(self):
        s_pole_crs = ccrs.LambertConformal(standard_parallels=[-1.])
        assert_array_almost_equal(s_pole_crs.x_limits,
                                  (-19939660, 19939660),
                                  decimal=0)
        assert_array_almost_equal(s_pole_crs.y_limits,
                                  (-735590302, -8183795),
                                  decimal=0)

    def test_single_npole(self):
        n_pole_crs = ccrs.LambertConformal(standard_parallels=[1.])
        assert_array_almost_equal(n_pole_crs.x_limits,
                                  (-20130569, 20130569),
                                  decimal=0)
        assert_array_almost_equal(n_pole_crs.y_limits,
                                  (-8170229, 726200683),
                                  decimal=0)
