# (C) British Crown Copyright 2011 - 2019, Met Office
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

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
try:
    from unittest import mock
except ImportError:
    import mock
import pytest

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

from cartopy.tests.mpl import MPL_VERSION, ImageTesting


TEST_PROJS = [
    ccrs.PlateCarree,
    ccrs.AlbersEqualArea,
    ccrs.AzimuthalEquidistant,
    ccrs.LambertConformal,
    ccrs.LambertCylindrical,
    ccrs.Mercator,
    ccrs.Miller,
    ccrs.Mollweide,
    ccrs.Orthographic,
    ccrs.Robinson,
    ccrs.Sinusoidal,
    ccrs.Stereographic,
    ccrs.TransverseMercator,
    #     ccrs.UTM,
    ccrs.InterruptedGoodeHomolosine,
    ccrs.RotatedPole,
    ccrs.OSGB,
    ccrs.EuroPP,
    #     ccrs.Geostationary,  # same as orthographic
    #     ccrs.NearsidePerspective,  # same as orthographic
    ccrs.Gnomonic,
    ccrs.LambertAzimuthalEqualArea,
    ccrs.NorthPolarStereo,
    #     ccrs.OSNI,  # fails in 0.17.0
    ccrs.SouthPolarStereo,
]

RP = ccrs.RotatedPole(pole_longitude=180.0,
                      pole_latitude=36.0,
                      central_rotated_longitude=-106.0,
                      globe=ccrs.Globe(semimajor_axis=6370000,
                                       semiminor_axis=6370000))


@pytest.mark.natural_earth
@ImageTesting(['gridliner1'])
def test_gridliner():
    ny, nx = 2, 4

    plt.figure(figsize=(10, 10))

    ax = plt.subplot(nx, ny, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(linestyle=':')

    ax = plt.subplot(nx, ny, 2, projection=ccrs.OSGB())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(linestyle=':')

    ax = plt.subplot(nx, ny, 3, projection=ccrs.OSGB())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.PlateCarree(), color='blue', linestyle='-')
    ax.gridlines(ccrs.OSGB(), linestyle=':')

    ax = plt.subplot(nx, ny, 4, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.NorthPolarStereo(), alpha=0.5,
                 linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 5, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    osgb = ccrs.OSGB()
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb, linestyle=':')

    ax = plt.subplot(nx, ny, 6, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(alpha=0.5, linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 7, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines()
    osgb = ccrs.OSGB()
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb, linestyle=':')

    ax = plt.subplot(nx, ny, 8,
                     projection=ccrs.Robinson(central_longitude=135))
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.PlateCarree(), alpha=0.5, linewidth=1.5, linestyle='-')

    delta = 1.5e-2
    plt.subplots_adjust(left=0 + delta, right=1 - delta,
                        top=1 - delta, bottom=0 + delta)


def test_gridliner_specified_lines():
<<<<<<< HEAD
    xs = [0, 60, 120, 180, 240, 360]
    ys = [-90, -60, -30, 0, 30, 60, 90]
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl = GeoAxes.gridlines(ax, xlocs=xs, ylocs=ys)
=======
    meridians = [0, 60, 120, 180, 240, 360]
    parallels = [-90, -60, -30, 0, 30, 60, 90]

    def mpl_connext(*args):
        pass

    canvas = mock.Mock(mpl_connext=mpl_connext)
    fig = mock.Mock(spec=matplotlib.figure.Figure, canvas=canvas)
    ax = mock.Mock(_gridliners=[], spec=GeoAxes, figure=fig)
    gl = GeoAxes.gridlines(ax, xlocs=meridians, ylocs=parallels)
>>>>>>> 45d0479abed27f73c5bb8bac35036772a3faa773
    assert isinstance(gl.xlocator, mticker.FixedLocator)
    assert isinstance(gl.ylocator, mticker.FixedLocator)
    assert gl.xlocator.tick_values(None, None).tolist() == meridians
    assert gl.ylocator.tick_values(None, None).tolist() == parallels


# The tolerance on this test is particularly high because of the high number
# of text objects. A new testing strategy is needed for this kind of test.
if MPL_VERSION >= '2.0':
    grid_label_image = 'gridliner_labels'
    grid_label_inline_image = 'gridliner_labels_inline'
    grid_label_inline_usa_image = 'gridliner_labels_inline_usa'
else:
    grid_label_image = 'gridliner_labels_1.5'
    grid_label_inline_image = 'gridliner_labels_inline_1.5'
    grid_label_inline_usa_image = 'gridliner_labels_inline_usa_1.5'


@pytest.mark.natural_earth
@ImageTesting([grid_label_image])
def test_grid_labels():
    fig = plt.figure(figsize=(10, 10))

    crs_pc = ccrs.PlateCarree()
    crs_merc = ccrs.Mercator()

    ax = fig.add_subplot(3, 2, 1, projection=crs_pc)
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    # Check that adding labels to Mercator gridlines gives an error.
    # (Currently can only label PlateCarree gridlines.)
    ax = fig.add_subplot(3, 2, 2,
                         projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines()

    ax.set_title('Known bug')
    gl = ax.gridlines(crs=crs_pc, draw_labels=True)
    gl.top_labels = False
    gl.left_labels = False
    gl.xlines = False

    ax = fig.add_subplot(3, 2, 3, projection=crs_merc)
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.xlabel_style = gl.ylabel_style = {'size': 9}

<<<<<<< HEAD
    # Check that labelling the gridlines on an OSGB plot gives an error.
    # (Currently can only draw these on PlateCarree or Mercator plots.)
    ax = fig.add_subplot(3, 2, 4, projection=crs_osgb)
    ax.coastlines()
    ax.remove()

    ax = fig.add_subplot(3, 2, 4, projection=crs_pc)
=======
    ax = plt.subplot(3, 2, 4, projection=crs_pc)
>>>>>>> 45d0479abed27f73c5bb8bac35036772a3faa773
    ax.coastlines()
    gl = ax.gridlines(
        crs=crs_pc, linewidth=2, color='gray', alpha=0.5, linestyle=':')
    gl.bottom_labels = True
    gl.right_labels = True
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-180, -45, 45, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red'}
    gl.xpadding = 10
    gl.ypadding = 15

    # trigger a draw at this point and check the appropriate artists are
    # populated on the gridliner instance
    fig.canvas.draw()

    assert len(gl.bottom_label_artists) == 4
    assert len(gl.top_label_artists) == 0
    assert len(gl.left_label_artists) == 0
    assert len(gl.right_label_artists) != 0
    assert len(gl.xline_artists) == 0

    ax = fig.add_subplot(3, 2, 5, projection=crs_pc)
    ax.set_extent([-20, 10.0, 45.0, 70.0])
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    ax = fig.add_subplot(3, 2, 6, projection=crs_merc)
    ax.set_extent([-20, 10.0, 45.0, 70.0], crs=crs_pc)
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.rotate_labels = False
    gl.xlabel_style = gl.ylabel_style = {'size': 9}

    # Increase margins between plots to stop them bumping into one another.
<<<<<<< HEAD
    plt.subplots_adjust(wspace=0.25, hspace=0.25)


@pytest.mark.natural_earth
@ImageTesting([grid_label_inline_image])
def test_grid_labels_inline():
    plt.figure(figsize=(35, 30))
    for i, proj in enumerate(TEST_PROJS, 1):
        if isinstance(proj(), ccrs.RotatedPole):
            ax = plt.subplot(7, 4, i, projection=RP)
        else:
            ax = plt.subplot(7, 4, i, projection=proj())
        ax.gridlines(draw_labels=True, auto_inline=True)
        ax.coastlines()
        ax.set_title(proj, y=1.08)


@pytest.mark.natural_earth
@ImageTesting([grid_label_inline_usa_image])
def test_grid_labels_inline_usa():
    top = 49.3457868  # north lat
    left = -124.7844079  # west long
    right = -66.9513812  # east long
    bottom = 24.7433195  # south lat
    plt.figure(figsize=(35, 30))
    for i, proj in enumerate(TEST_PROJS, 1):
        if isinstance(proj(), ccrs.RotatedPole):
            ax = plt.subplot(7, 4, i, projection=RP)
        else:
            ax = plt.subplot(7, 4, i, projection=proj())
        try:
            ax.set_extent([left, right, bottom, top],
                          crs=ccrs.PlateCarree())
        except Exception:
            pass
        ax.set_title(proj, y=1.08)
        ax.gridlines(draw_labels=True, auto_inline=True, clip_on=True)
        ax.coastlines()
=======
    plt.subplots_adjust(wspace=0.25, hspace=0.25, top=.98, left=.07,
                        bottom=0.02, right=0.93)
>>>>>>> 45d0479abed27f73c5bb8bac35036772a3faa773
