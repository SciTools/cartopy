# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from unittest.mock import Mock

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
    expected = ['90\u00B0S', '60\u00B0S', '30\u00B0S', '0\u00B0',
                '30\u00B0N', '60\u00B0N', '90\u00B0N']
    assert result == expected


def test_LatitudeFormatter_direction_label():
    formatter = LatitudeFormatter(direction_label=False)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['-90\u00B0', '-60\u00B0', '-30\u00B0', '0\u00B0',
                '30\u00B0', '60\u00B0', '90\u00B0']
    assert result == expected


def test_LatitudeFormatter_degree_symbol():
    formatter = LatitudeFormatter(degree_symbol='')
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['90S', '60S', '30S', '0',
                '30N', '60N', '90N']
    assert result == expected


def test_LatitudeFormatter_number_format():
    formatter = LatitudeFormatter(number_format='.2f', dms=False)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-90, -60, -30, 0, 30, 60, 90]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['90.00\u00B0S', '60.00\u00B0S', '30.00\u00B0S',
                '0.00\u00B0', '30.00\u00B0N', '60.00\u00B0N',
                '90.00\u00B0N']
    assert result == expected


def test_LatitudeFormatter_mercator():
    formatter = LatitudeFormatter()
    p = ccrs.Mercator()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-15496570.739707904, -8362698.548496634,
                  -3482189.085407435, 0.0, 3482189.085407435,
                  8362698.548496634, 15496570.739707898]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['80\u00B0S', '60\u00B0S', '30\u00B0S', '0\u00B0',
                '30\u00B0N', '60\u00B0N', '80\u00B0N']
    assert result == expected


def test_LatitudeFormatter_small_numbers():
    formatter = LatitudeFormatter(number_format='.7f', dms=False)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [40.1275150, 40.1275152, 40.1275154]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['40.1275150\u00B0N', '40.1275152\u00B0N',
                '40.1275154\u00B0N']
    assert result == expected


def test_LongitudeFormatter_direction_label():
    formatter = LongitudeFormatter(direction_label=False,
                                   dateline_direction_label=True,
                                   zero_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['-180\u00B0', '-120\u00B0', '-60\u00B0', '0\u00B0',
                '60\u00B0', '120\u00B0', '180\u00B0']
    assert result == expected


def test_LongitudeFormatter_central_longitude_0():
    formatter = LongitudeFormatter(dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['180\u00B0W', '120\u00B0W', '60\u00B0W', '0\u00B0',
                '60\u00B0E', '120\u00B0E', '180\u00B0E']
    assert result == expected


def test_LongitudeFormatter_central_longitude_180():
    formatter = LongitudeFormatter(zero_direction_label=True)
    p = ccrs.PlateCarree(central_longitude=180)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['0\u00B0E', '60\u00B0E', '120\u00B0E', '180\u00B0',
                '120\u00B0W', '60\u00B0W', '0\u00B0W']
    assert result == expected


def test_LongitudeFormatter_central_longitude_120():
    formatter = LongitudeFormatter()
    p = ccrs.PlateCarree(central_longitude=120)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['60\u00B0W', '0\u00B0', '60\u00B0E', '120\u00B0E',
                '180\u00B0', '120\u00B0W', '60\u00B0W']
    assert result == expected


def test_LongitudeFormatter_degree_symbol():
    formatter = LongitudeFormatter(degree_symbol='',
                                   dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['180W', '120W', '60W', '0', '60E', '120E', '180E']
    assert result == expected


def test_LongitudeFormatter_number_format():
    formatter = LongitudeFormatter(number_format='.2f', dms=False,
                                   dateline_direction_label=True)
    p = ccrs.PlateCarree()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-180, -120, -60, 0, 60, 120, 180]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['180.00\u00B0W', '120.00\u00B0W', '60.00\u00B0W',
                '0.00\u00B0', '60.00\u00B0E', '120.00\u00B0E',
                '180.00\u00B0E']
    assert result == expected


def test_LongitudeFormatter_mercator():
    formatter = LongitudeFormatter(dateline_direction_label=True)
    p = ccrs.Mercator()
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-20037508.342783064, -13358338.895188706,
                  -6679169.447594353, 0.0, 6679169.447594353,
                  13358338.895188706, 20037508.342783064]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['180\u00B0W', '120\u00B0W', '60\u00B0W', '0\u00B0',
                '60\u00B0E', '120\u00B0E', '180\u00B0E']
    assert result == expected


def test_LongitudeFormatter_small_numbers_0():
    formatter = LongitudeFormatter(number_format='.7f', dms=False)
    p = ccrs.PlateCarree(central_longitude=0)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-17.1142343, -17.1142340, -17.1142337]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['17.1142343\u00B0W', '17.1142340\u00B0W',
                '17.1142337\u00B0W']
    assert result == expected


def test_LongitudeFormatter_small_numbers_180():
    formatter = LongitudeFormatter(zero_direction_label=True, dms=False,
                                   number_format='.7f')
    p = ccrs.PlateCarree(central_longitude=180)
    formatter.axis = Mock(axes=Mock(GeoAxes, projection=p))
    test_ticks = [-17.1142343, -17.1142340, -17.1142337]
    result = [formatter(tick) for tick in test_ticks]
    expected = ['162.8857657\u00B0E', '162.8857660\u00B0E',
                '162.8857663\u00B0E']
    assert result == expected


@pytest.mark.parametrize("test_ticks,expected",
                         [pytest.param([-3.75, -3.5],
                                       ["3\u00B045'W", "3\u00B030'W"],
                                       id='minutes_no_hide'),
                          pytest.param([-3.5, -3.],
                                       ["30'", "3\u00B0W"],
                                       id='minutes_hide'),
                          pytest.param([-3. - 2 * ONE_MIN - 30 * ONE_SEC],
                                       ["3\u00B02'30''W"],
                                       id='seconds'),
                          ])
def test_LongitudeFormatter_minutes_seconds(test_ticks, expected):
    formatter = LongitudeFormatter(dms=True, auto_hide=True)
    formatter.set_locs(test_ticks)
    result = [formatter(tick) for tick in test_ticks]
    assert result == expected


@pytest.mark.parametrize("test_ticks,expected",
                         [pytest.param([-3.75, -3.5],
                                       ["-3\u00B045'", "-3\u00B030'"],
                                       id='minutes_no_hide'),
                          pytest.param([-3.5, -3.],
                                       ["30'", "-3\u00B0"],
                                       id='minutes_hide'),
                          pytest.param([-3. - 2 * ONE_MIN - 30 * ONE_SEC],
                                       ["-3\u00B02'30''"],
                                       id='seconds'),
                          ])
def test_LongitudeFormatter_minutes_seconds_direction_label(test_ticks,
                                                            expected):
    formatter = LongitudeFormatter(
        dms=True, auto_hide=True, direction_label=False)
    formatter.set_locs(test_ticks)
    result = [formatter(tick) for tick in test_ticks]
    assert result == expected


@pytest.mark.parametrize("test_ticks,expected",
                         [pytest.param([-3.75, -3.5],
                                       ["3\u00B045'S", "3\u00B030'S"],
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
    assert ticklabels == [f'{v:g}{letter}' for v in ticks]
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
                                                 -20., -10., 0.]) / 60,
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


def test_lonlatformatter_decimal_point():
    xticker = LongitudeFormatter(decimal_point=',', number_format='0.2f')
    yticker = LatitudeFormatter(decimal_point=',', number_format='0.2f')
    assert xticker(-10) == "10,00°W"
    assert yticker(-10) == "10,00°S"


def test_lonlatformatter_cardinal_labels():
    xticker = LongitudeFormatter(cardinal_labels={'west': 'O'})
    yticker = LatitudeFormatter(cardinal_labels={'south': 'South'})
    assert xticker(-10) == "10°O"
    assert yticker(-10) == "10°South"
