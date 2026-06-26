# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
from cartopy.mpl import _MPL_311
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='xticks_no_transform.png', style='mpl20',
                               tolerance=5.20 if not _MPL_311 else 0.5)
def test_set_xticks_no_transform():
    plt.rc('figure.subplot', left=0.125, right=0.9, bottom=0.1, top=0.9)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xticks([-135, -45, 45, 135], minor=True)
    return fig


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='xticks_cylindrical.png', style='mpl20',
                               tolerance=5.04 if not _MPL_311 else 0.5)
def test_set_xticks_cylindrical():
    plt.rc('figure.subplot', left=0.125, right=0.9, bottom=0.1, top=0.9)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=-85, max_latitude=85))
    ax.coastlines('110m')
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    ax.set_xticks([-180, -90, 0, 90, 180], crs=ccrs.PlateCarree())
    ax.set_xticks([-135, -45, 45, 135], minor=True, crs=ccrs.PlateCarree())
    return fig


def test_set_xticks_non_cylindrical():
    ax = plt.axes(projection=ccrs.Orthographic())
    with pytest.raises(RuntimeError):
        ax.set_xticks([-180, -90, 0, 90, 180], crs=ccrs.Geodetic())
    with pytest.raises(RuntimeError):
        ax.set_xticks([-135, -45, 45, 135], minor=True, crs=ccrs.Geodetic())


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='yticks_no_transform.png', style='mpl20',
                               tolerance=4.90 if not _MPL_311 else 0.5)
def test_set_yticks_no_transform():
    plt.rc('figure.subplot', left=0.125, right=0.9, bottom=0.1, top=0.9)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('110m')
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_yticks([-75, -45, -15, 15, 45, 75], minor=True)
    return fig


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='yticks_cylindrical.png', style='mpl20',
                               tolerance=4.89 if not _MPL_311 else 0.5)
def test_set_yticks_cylindrical():
    plt.rc('figure.subplot', left=0.125, right=0.9, bottom=0.1, top=0.9)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=-85, max_latitude=85))
    ax.coastlines('110m')
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
    ax.set_yticks([-60, -30, 0, 30, 60], crs=ccrs.PlateCarree())
    ax.set_yticks([-75, -45, -15, 15, 45, 75], minor=True,
                  crs=ccrs.PlateCarree())
    return fig


def test_set_yticks_non_cylindrical():
    ax = plt.axes(projection=ccrs.Orthographic())
    with pytest.raises(RuntimeError):
        ax.set_yticks([-60, -30, 0, 30, 60], crs=ccrs.Geodetic())
    with pytest.raises(RuntimeError):
        ax.set_yticks([-75, -45, -15, 15, 45, 75], minor=True,
                      crs=ccrs.Geodetic())


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='xyticks.png', style='mpl20',
                               tolerance=5.42 if not _MPL_311 else 0.5)
def test_set_xyticks():
    plt.rc('figure.subplot', left=0.125, right=0.9, bottom=0.1, top=0.9)
    plt.rc('axes.formatter', limits=(-7, 7))
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
    return fig
