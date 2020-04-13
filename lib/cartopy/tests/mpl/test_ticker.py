# (C) British Crown Copyright 2014 - 2020, Met Office
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
import matplotlib.pyplot as plt
import pytest
import numpy as np

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import (LatitudeFormatter, LongitudeFormatter,
                                LatitudeLocator, LongitudeLocator)

ONE_MIN = 1 / 60.
ONE_SEC = 1 / 3600.


def test_LatitudeFormatter_bad_projection():
    formatter = LatitudeFormatter()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=ccrs.Orthographic()))
    match = r'This formatter cannot be used with non-rectangular projections\.'
    with pytest.raises(TypeError, match=match):
        formatter(0)


def test_LongitudeFormatter_bad_projection():
    formatter = LongitudeFormatter()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=ccrs.Orthographic()))
    match = r'This formatter cannot be used with non-rectangular projections\.'
    with pytest.raises(TypeError, match=match):
        formatter(0)


def test_LatitudeFormatter():
    formatter = LatitudeFormatter()
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'90\u00B0S', u'60\u00B0S', u'30\u00B0S', u'0\u00B0',
                u'30\u00B0N', u'60\u00B0N', u'90\u00B0N']
    assert result == expected


def test_LatitudeFormatter_degree_symbol():
    formatter = LatitudeFormatter(degree_symbol='')
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'90S', u'60S', u'30S', u'0',
                u'30N', u'60N', u'90N']
    assert result == expected


def test_LatitudeFormatter_number_format():
    formatter = LatitudeFormatter(number_format='.2f', dms=False)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'90.00\u00B0S', u'60.00\u00B0S', u'30.00\u00B0S',
                u'0.00\u00B0', u'30.00\u00B0N', u'60.00\u00B0N',
                u'90.00\u00B0N']
    assert result == expected


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
    assert result == expected


def test_LatitudeFormatter_small_numbers():
    formatter = LatitudeFormatter(number_format='.7f', dms=False)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [40.1275150, 40.1275152, 40.1275154]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'40.1275150\u00B0N', u'40.1275152\u00B0N',
                u'40.1275154\u00B0N']
    assert result == expected


def test_LongitudeFormatter_central_longitude_0():
    formatter = LongitudeFormatter(dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180\u00B0W', u'120\u00B0W', u'60\u00B0W', u'0\u00B0',
                u'60\u00B0E', u'120\u00B0E', u'180\u00B0E']
    assert result == expected


def test_LongitudeFormatter_central_longitude_180():
    formatter = LongitudeFormatter(zero_direction_label=True)
    p = ccrs.PlateCarree(central_longitude=180)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'0\u00B0E', u'60\u00B0E', u'120\u00B0E', u'180\u00B0',
                u'120\u00B0W', u'60\u00B0W', u'0\u00B0W']
    assert result == expected


def test_LongitudeFormatter_central_longitude_120():
    formatter = LongitudeFormatter()
    p = ccrs.PlateCarree(central_longitude=120)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'60\u00B0W', u'0\u00B0', u'60\u00B0E', u'120\u00B0E',
                u'180\u00B0', u'120\u00B0W', u'60\u00B0W']
    assert result == expected


def test_LongitudeFormatter_degree_symbol():
    formatter = LongitudeFormatter(degree_symbol='',
                                   dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180W', u'120W', u'60W', u'0', u'60E', u'120E', u'180E']
    assert result == expected


def test_LongitudeFormatter_number_format():
    formatter = LongitudeFormatter(number_format='.2f', dms=False,
                                   dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'180.00\u00B0W', u'120.00\u00B0W', u'60.00\u00B0W',
                u'0.00\u00B0', u'60.00\u00B0E', u'120.00\u00B0E',
                u'180.00\u00B0E']
    assert result == expected


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
    assert result == expected


def test_LongitudeFormatter_small_numbers_0():
    formatter = LongitudeFormatter(number_format='.7f', dms=False)
    p = ccrs.PlateCarree(central_longitude=0)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-17.1142343, -17.1142340, -17.1142337]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'17.1142343\u00B0W', u'17.1142340\u00B0W',
                u'17.1142337\u00B0W']
    assert result == expected


def test_LongitudeFormatter_small_numbers_180():
    formatter = LongitudeFormatter(zero_direction_label=True, dms=False,
                                   number_format='.7f')
    p = ccrs.PlateCarree(central_longitude=180)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-17.1142343, -17.1142340, -17.1142337]
    result = [formatter(tick) for tick in test_ticks]
    expected = [u'162.8857657\u00B0E', u'162.8857660\u00B0E',
                u'162.8857663\u00B0E']
    assert result == expected


@pytest.mark.parametrize("test_ticks,expected",
                         [pytest.param([-3.75, -3.5],
                                       [u"3\u00B0W45'", u"3\u00B0W30'"],
                                       id='minutes_no_hide'),
                          pytest.param([-3.5, -3.],
                                       [u"30'", u"3\u00B0W"],
                                       id='minutes_hide'),
                          pytest.param([-3. - 2 * ONE_MIN - 30 * ONE_SEC],
                                       [u"3\u00B0W2'30''"],
                                       id='seconds'),
                          ])
def test_LongitudeFormatter_minutes_seconds(test_ticks, expected):
    formatter = LongitudeFormatter(dms=True, auto_hide=True)
    formatter.set_locs(test_ticks)
    result = [formatter(tick) for tick in test_ticks]
    assert result == expected


@pytest.mark.parametrize("test_ticks,expected",
                         [pytest.param([-3.75, -3.5],
                                       [u"3\u00B0S45'", u"3\u00B0S30'"],
                                       id='minutes_no_hide'),
                          ])
def test_LatitudeFormatter_minutes_seconds(test_ticks, expected):
    formatter = LatitudeFormatter(dms=True, auto_hide=True)
    formatter.set_locs(test_ticks)
    result = [formatter(tick) for tick in test_ticks]
    assert result == expected


@pytest.mark.parametrize("cls,letter",
                         [(LongitudeFormatter, 'E'),
                          (LatitudeFormatter, 'N')])
def test_lonlatformatter_non_geoaxes(cls, letter):
    ticks = [2, 2.5]
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot([0, 10], [0, 1])
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(cls(degree_symbol='', dms=False))
    fig.canvas.draw()
    ticklabels = [t.get_text() for t in ax.get_xticklabels()]
    assert ticklabels == ['{:g}{}'.format(v, letter) for v in ticks]
    plt.close()


@pytest.mark.parametrize("cls,vmin,vmax,expected",
                         [pytest.param(LongitudeLocator, -180, 180,
                                       [-180., -120., -60., 0.,
                                        60., 120., 180.],
                                       id='lon_large'),
                          pytest.param(LatitudeLocator, -180, 180,
                                       [-90.0, -60.0, -30.0, 0.0,
                                        30.0, 60.0, 90.0],
                                       id='lat_large'),
                          pytest.param(LongitudeLocator, -10, 0,
                                       [-10.5, -9., -7.5, -6., -4.5,
                                        -3., -1.5, 0.],
                                       id='lon_medium'),
                          pytest.param(LongitudeLocator, -1, 0,
                                       np.array([-60., -50., -40., -30.,
                                                 -20., -10.,   0.]) / 60,
                                       id='lon_small'),
                          pytest.param(LongitudeLocator, 0, 2 * ONE_MIN,
                                       np.array([0., 18., 36., 54., 72., 90.,
                                                 108., 126.]) / 3600,
                                       id='lon_tiny'),
                          ])
def test_LongitudeLocator(cls, vmin, vmax, expected):
    locator = cls(dms=True)
    result = locator.tick_values(vmin, vmax)
    np.testing.assert_allclose(result, expected)
