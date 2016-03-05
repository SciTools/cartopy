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

from __future__ import (absolute_import, division, print_function)

try:
    from unittest.mock import Mock
except ImportError:
    from mock import Mock
from nose.tools import assert_equal
try:
    from nose.tools import assert_raises_regex
except ImportError:
    from nose.tools import assert_raises_regexp as assert_raises_regex
from matplotlib.axes import Axes

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter


def test_LatitudeFormatter_bad_axes():
    formatter = LatitudeFormatter()
    formatter.axis = Mock(axes=Mock(Axes, projection=ccrs.PlateCarree()))
    message = 'This formatter can only be used with cartopy axes.'
    with assert_raises_regex(TypeError, message):
        formatter(0)


def test_LatitudeFormatter_bad_projection():
    formatter = LatitudeFormatter()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=ccrs.Orthographic()))
    message = 'This formatter cannot be used with non-rectangular projections.'
    with assert_raises_regex(TypeError, message):
        formatter(0)


def test_LongitudeFormatter_bad_axes():
    formatter = LongitudeFormatter()
    formatter.axis = Mock(axes=Mock(Axes, projection=ccrs.PlateCarree()))
    message = 'This formatter can only be used with cartopy axes.'
    with assert_raises_regex(TypeError, message):
        formatter(0)


def test_LongitudeFormatter_bad_projection():
    formatter = LongitudeFormatter()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=ccrs.Orthographic()))
    message = 'This formatter cannot be used with non-rectangular projections.'
    with assert_raises_regex(TypeError, message):
        formatter(0)


def test_LatitudeFormatter():
    formatter = LatitudeFormatter()
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'90\u00B0S', u'60\u00B0S', u'30\u00B0S', u'0\u00B0',
                u'30\u00B0N', u'60\u00B0N', u'90\u00B0N']
    assert_equal(result, expected)


def test_LatitudeFormatter_degree_symbol():
    formatter = LatitudeFormatter(degree_symbol='')
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'90S', u'60S', u'30S', u'0',
                u'30N', u'60N', u'90N']
    assert_equal(result, expected)


def test_LatitudeFormatter_number_format():
    formatter = LatitudeFormatter(number_format='.2f')
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'90.00\u00B0S', u'60.00\u00B0S', u'30.00\u00B0S',
                u'0.00\u00B0', u'30.00\u00B0N', u'60.00\u00B0N',
                u'90.00\u00B0N']
    assert_equal(result, expected)


def test_LatitudeFormatter_mercator():
    formatter = LatitudeFormatter()
    p = ccrs.Mercator()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-15496570.739707904, -8362698.548496634,
                  -3482189.085407435, 0.0, 3482189.085407435,
                  8362698.548496634, 15496570.739707898]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'80\u00B0S', u'60\u00B0S', u'30\u00B0S', u'0\u00B0',
                u'30\u00B0N', u'60\u00B0N', u'80\u00B0N']
    assert_equal(result, expected)


def test_LatitudeFormatter_small_numbers():
    formatter = LatitudeFormatter(number_format='.7f')
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [40.1275150, 40.1275152, 40.1275154]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'40.1275150\u00B0N', u'40.1275152\u00B0N',
                u'40.1275154\u00B0N']
    assert_equal(result, expected)


def test_LongitudeFormatter_central_longitude_0():
    formatter = LongitudeFormatter(dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180\u00B0W', u'120\u00B0W', u'60\u00B0W', u'0\u00B0',
                u'60\u00B0E', u'120\u00B0E', u'180\u00B0E']
    assert_equal(result, expected)


def test_LongitudeFormatter_central_longitude_180():
    formatter = LongitudeFormatter(zero_direction_label=True)
    p = ccrs.PlateCarree(central_longitude=180)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'0\u00B0E', u'60\u00B0E', u'120\u00B0E', u'180\u00B0',
                u'120\u00B0W', u'60\u00B0W', u'0\u00B0W']
    assert_equal(result, expected)


def test_LongitudeFormatter_central_longitude_120():
    formatter = LongitudeFormatter()
    p = ccrs.PlateCarree(central_longitude=120)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'60\u00B0W', u'0\u00B0', u'60\u00B0E', u'120\u00B0E',
                u'180\u00B0', u'120\u00B0W', u'60\u00B0W']
    assert_equal(result, expected)


def test_LongitudeFormatter_degree_symbol():
    formatter = LongitudeFormatter(degree_symbol='',
                                   dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180W', u'120W', u'60W', u'0', u'60E', u'120E', u'180E']
    assert_equal(result, expected)


def test_LongitudeFormatter_number_format():
    formatter = LongitudeFormatter(number_format='.2f',
                                   dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180.00\u00B0W', u'120.00\u00B0W', u'60.00\u00B0W',
                u'0.00\u00B0', u'60.00\u00B0E', u'120.00\u00B0E',
                u'180.00\u00B0E']
    assert_equal(result, expected)


def test_LongitudeFormatter_mercator():
    formatter = LongitudeFormatter(dateline_direction_label=True)
    p = ccrs.Mercator()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-20037508.342783064, -13358338.895188706,
                  -6679169.447594353, 0.0, 6679169.447594353,
                  13358338.895188706, 20037508.342783064]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180\u00B0W', u'120\u00B0W', u'60\u00B0W', u'0\u00B0',
                u'60\u00B0E', u'120\u00B0E', u'180\u00B0E']
    assert_equal(result, expected)


def test_LongitudeFormatter_small_numbers_0():
    formatter = LongitudeFormatter(number_format='.7f')
    p = ccrs.PlateCarree(central_longitude=0)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-17.1142343, -17.1142340, -17.1142337]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'17.1142343\u00B0W', u'17.1142340\u00B0W',
                u'17.1142337\u00B0W']
    assert_equal(result, expected)


def test_LongitudeFormatter_small_numbers_180():
    formatter = LongitudeFormatter(zero_direction_label=True,
                                   number_format='.7f')
    p = ccrs.PlateCarree(central_longitude=180)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-17.1142343, -17.1142340, -17.1142337]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'162.8857657\u00B0E', u'162.8857660\u00B0E',
                u'162.8857663\u00B0E']
    assert_equal(result, expected)
