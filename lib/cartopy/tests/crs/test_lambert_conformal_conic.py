from __future__ import (absolute_import, division, print_function)

from numpy.testing import assert_array_almost_equal
import pytest

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_defaults():
    conus = ccrs.LambertConformalConic()
    other_args = {'ellps=intl', 'lon_0=-97.5', 'lat_0=38.5',
                  'lat_1=38.5', 'lat_2=38.5', 'x_0=0.0', 'y_0=0.0', 'units=m'}
    check_proj_params('lcc', conus, other_args)


def test_specific_region():
    aus_kwargs = dict(
        central_longitude=134.489,
        central_latitude=-23.993,
        x_limits=(-2.25e6, 2.25e6), y_limits=(-1.9e6, 1.9e6)
    )
    aus = ccrs.LambertConformalConic(**aus_kwargs)
    aus2 = ccrs.LambertConformalConic(**aus_kwargs)
    default = ccrs.LambertConformalConic()
    other_args = {'ellps=intl', 'lon_0=134.489', 'lat_0=-23.993',
                  'lat_1=-23.993', 'lat_2=-23.993' 'units=m',
                  'x_0=0.0', 'y_0=0.0'}
    check_proj_params('lcc', aus, other_args)
    assert aus == aus2
    assert aus != default
    assert hash(aus) != hash(default)
    assert hash(aus) == hash(aus2)
    assert_array_almost_equal(aus.x_limits, (-2250000.0000001, 2249999.999999))


def test_too_many_standard_parallels():
    with pytest.raises(ValueError, match='1 or 2 standard parallels'):
        ccrs.LambertConformalConic(standard_parallels=(23.5, 37, 42.5))
