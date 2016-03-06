# (C) British Crown Copyright 2011 - 2016, Met Office
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

import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker
import nose.tools

import cartopy.crs as ccrs
from cartopy.tests.mpl import ImageTesting


def _format_lat(val, i):
    if val > 0:
        return '%.0fN' % val
    elif val < 0:
        return '%.0fS' % abs(val)
    else:
        return '0'


def _format_lon(val, i):
    # Apply periodic boundary conditions, with an almost equal test on 180 lon.
    while val > 180:
        val -= 360
    while val < -180:
        val += 360
    if abs(abs(val) - 180.) <= 1e-06 or val == 0:
        return '%.0f' % abs(val)
    elif val > 0:
        return '%.0fE' % val
    elif val < 0:
        return '%.0fW' % abs(val)


test_fn_suffix = '' if mpl.__version__ >= '1.5' else '_pre_mpl_1.5'


@ImageTesting(['xticks_no_transform' + test_fn_suffix], tolerance=0.5)
def test_set_xticks_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lon))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lat))
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xticks([-135, -45, 45, 135], minor=True)


@ImageTesting(['xticks_cylindrical' + test_fn_suffix], tolerance=0.5)
def test_set_xticks_cylindrical():
    ax = plt.axes(projection=ccrs.Mercator(
                  min_latitude=-85.,
                  max_latitude=85.,
                  globe=ccrs.Globe(semimajor_axis=math.degrees(1))))
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lon))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lat))
    ax.set_xticks([-180, -90, 0, 90, 180], crs=ccrs.PlateCarree())
    ax.set_xticks([-135, -45, 45, 135], minor=True, crs=ccrs.PlateCarree())


def test_set_xticks_non_cylindrical():
    ax = plt.axes(projection=ccrs.Orthographic())
    with nose.tools.assert_raises(RuntimeError):
        ax.set_xticks([-180, -90, 0, 90, 180], crs=ccrs.Geodetic())
    with nose.tools.assert_raises(RuntimeError):
        ax.set_xticks([-135, -45, 45, 135], minor=True, crs=ccrs.Geodetic())
    plt.close()


@ImageTesting(['yticks_no_transform' + test_fn_suffix], tolerance=0.5)
def test_set_yticks_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lon))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lat))
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_yticks([-75, -45, 15, 45, 75], minor=True)


@ImageTesting(['yticks_cylindrical' + test_fn_suffix], tolerance=0.5)
def test_set_yticks_cylindrical():
    ax = plt.axes(projection=ccrs.Mercator(
                  min_latitude=-85.,
                  max_latitude=85.,
                  globe=ccrs.Globe(semimajor_axis=math.degrees(1))))
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lon))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(_format_lat))
    ax.set_yticks([-60, -30, 0, 30, 60], crs=ccrs.PlateCarree())
    ax.set_yticks([-75, -45, 15, 45, 75], minor=True, crs=ccrs.PlateCarree())


def test_set_yticks_non_cylindrical():
    ax = plt.axes(projection=ccrs.Orthographic())
    with nose.tools.assert_raises(RuntimeError):
        ax.set_yticks([-60, -30, 0, 30, 60], crs=ccrs.Geodetic())
    with nose.tools.assert_raises(RuntimeError):
        ax.set_yticks([-75, -45, 15, 45, 75], minor=True, crs=ccrs.Geodetic())
    plt.close()


@ImageTesting(['xyticks' + test_fn_suffix], tolerance=0.5)
def test_set_xyticks():
    fig = plt.figure(figsize=(10, 10))
    projections = (ccrs.PlateCarree(),
                   ccrs.Mercator(globe=ccrs.Globe(
                       semimajor_axis=math.degrees(1))),
                   ccrs.TransverseMercator())
    x = -3.275024
    y = 50.753998
    for i, prj in enumerate(projections, 1):
        ax = fig.add_subplot(3, 1, i, projection=prj)
        ax.set_extent([-12.5, 4, 49, 60], ccrs.Geodetic())
        ax.coastlines('110m')
        p, q = prj.transform_point(x, y, ccrs.Geodetic())
        ax.set_xticks([p])
        ax.set_yticks([q])


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
