# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.tests.mpl import ImageTesting


@pytest.mark.natural_earth
@ImageTesting(['xticks_no_transform'])
def test_set_xticks_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xticks([-135, -45, 45, 135], minor=True)


@pytest.mark.natural_earth
@ImageTesting(['xticks_cylindrical'])
def test_set_xticks_cylindrical():
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=-85, max_latitude=85))
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    ax.set_xticks([-180, -90, 0, 90, 180], crs=ccrs.PlateCarree())
    ax.set_xticks([-135, -45, 45, 135], minor=True, crs=ccrs.PlateCarree())


def test_set_xticks_non_cylindrical():
    ax = plt.axes(projection=ccrs.Orthographic())
    with pytest.raises(RuntimeError):
        ax.set_xticks([-180, -90, 0, 90, 180], crs=ccrs.Geodetic())
    with pytest.raises(RuntimeError):
        ax.set_xticks([-135, -45, 45, 135], minor=True, crs=ccrs.Geodetic())
    plt.close()


@pytest.mark.natural_earth
@ImageTesting(['yticks_no_transform'])
def test_set_yticks_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('110m')
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_yticks([-75, -45, -15, 15, 45, 75], minor=True)


@pytest.mark.natural_earth
@ImageTesting(['yticks_cylindrical'])
def test_set_yticks_cylindrical():
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=-85, max_latitude=85))
    ax.coastlines('110m')
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    ax.set_yticks([-60, -30, 0, 30, 60], crs=ccrs.PlateCarree())
    ax.set_yticks([-75, -45, -15, 15, 45, 75], minor=True,
                  crs=ccrs.PlateCarree())


def test_set_yticks_non_cylindrical():
    ax = plt.axes(projection=ccrs.Orthographic())
    with pytest.raises(RuntimeError):
        ax.set_yticks([-60, -30, 0, 30, 60], crs=ccrs.Geodetic())
    with pytest.raises(RuntimeError):
        ax.set_yticks([-75, -45, -15, 15, 45, 75], minor=True,
                      crs=ccrs.Geodetic())
    plt.close()


@pytest.mark.natural_earth
@ImageTesting(['xyticks'])
def test_set_xyticks():
    fig = plt.figure(figsize=(10, 10))
    projections = (ccrs.PlateCarree(),
                   ccrs.Mercator(),
                   ccrs.TransverseMercator(approx=False))
    x = -3.275024
    y = 50.753998
    for i, prj in enumerate(projections, 1):
        ax = fig.add_subplot(3, 1, i, projection=prj)
        ax.set_extent([-12.5, 4, 49, 60], ccrs.Geodetic())
        ax.coastlines('110m')
        p, q = prj.transform_point(x, y, ccrs.Geodetic())
        ax.set_xticks([p])
        ax.set_yticks([q])
