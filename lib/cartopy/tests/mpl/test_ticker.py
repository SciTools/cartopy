# (C) British Crown Copyright 2014, Met Office
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
from nose.tools import raises
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from cartopy.tests.mpl import ImageTesting


def _run_test(projection, xticks, yticks, xformatter, yformatter):
    ax = plt.axes(projection=projection)
    ax.set_global()
    ax.coastlines()
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(xformatter)
    ax.yaxis.set_major_formatter(yformatter)


@ImageTesting(['ticks_central_longitude_0'])
def test_central_longitude_0():
    _run_test(ccrs.PlateCarree(central_longitude=0),
              [-180, -120, -60, 0, 60, 120, 180],
              [-90, -60, -30, 0, 30, 60, 90],
              LongitudeFormatter(dateline_direction_label=True),
              LatitudeFormatter())


@ImageTesting(['ticks_central_longitude_180'])
def test_central_longitude_180():
    _run_test(ccrs.PlateCarree(central_longitude=180),
              [0, 60, 120, 180, 240, 300, 360],
              [-90, -60, -30, 0, 30, 60, 90],
              LongitudeFormatter(zero_direction_label=True),
              LatitudeFormatter())


@ImageTesting(['ticks_central_longitude_120'])
def test_central_longitude_120():
    _run_test(ccrs.PlateCarree(central_longitude=120),
              [-60, 0, 60, 120, 180, 240, 300],
              [-90, -60, -30, 0, 30, 60, 90],
              LongitudeFormatter(degree_symbol='', number_format='.2f'),
              LatitudeFormatter(degree_symbol='', number_format='.2f'))


@ImageTesting(['ticks_mercator'])
def test_mercator():
    _run_test(ccrs.Mercator(),
              [-180, -120, -60, 0, 60, 120, 180],
              [-80, -60, -30, 0, 30, 60, 80],
              LongitudeFormatter(dateline_direction_label=True),
              LatitudeFormatter())


@raises(TypeError)
def test__GeoFormatter_invalid_axes():
    ax = plt.axes()
    try:
        ax.xaxis.set_major_formatter(LongitudeFormatter)
    finally:
        plt.close()


@raises(TypeError)
def test__GeoFormatter_invalid_projection():
    ax = plt.axes(projection=ccrs.Stereographic())
    try:
        ax.xaxis.set_major_formatter(LongitudeFormatter)
    finally:
        plt.close()
