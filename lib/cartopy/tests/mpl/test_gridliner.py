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

import matplotlib as mpl
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
try:
    from unittest import mock
except ImportError:
    import mock
from nose.tools import assert_raises
import numpy as np

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

from cartopy.tests.mpl import ImageTesting


@ImageTesting(['gridliner1'])
def test_gridliner():
    ny, nx = 2, 4

    plt.figure(figsize=(10, 10))

    ax = plt.subplot(nx, ny, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines()

    ax = plt.subplot(nx, ny, 2, projection=ccrs.OSGB())
    ax.set_global()
    ax.coastlines()
    ax.gridlines()

    ax = plt.subplot(nx, ny, 3, projection=ccrs.OSGB())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.PlateCarree(), color='blue', linestyle='-')
    ax.gridlines(ccrs.OSGB())

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
    ax.gridlines(osgb)

    ax = plt.subplot(nx, ny, 6, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines()
    ax.gridlines(alpha=0.5, linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 7, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines()
    osgb = ccrs.OSGB()
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb)

    ax = plt.subplot(nx, ny, 8,
                     projection=ccrs.Robinson(central_longitude=135))
    ax.set_global()
    ax.coastlines()
    ax.gridlines(ccrs.PlateCarree(), alpha=0.5, linewidth=1.5, linestyle='-')

    delta = 1.5e-2
    plt.subplots_adjust(left=0 + delta, right=1 - delta,
                        top=1 - delta, bottom=0 + delta)


def test_gridliner_specified_lines():
    xs = [0, 60, 120, 180, 240, 360]
    ys = [-90, -60, -30, 0, 30, 60, 90]
    ax = mock.Mock(_gridliners=[], spec=GeoAxes)
    gl = GeoAxes.gridlines(ax, xlocs=xs, ylocs=ys)
    assert isinstance(gl.xlocator, mticker.FixedLocator)
    assert isinstance(gl.ylocator, mticker.FixedLocator)
    assert gl.xlocator.tick_values(None, None).tolist() == xs
    assert gl.ylocator.tick_values(None, None).tolist() == ys


# The tolerance on this test is particularly high because of the high number
# of text objects. A new testing strategy is needed for this kind of test.
@ImageTesting(['gridliner_labels'
               if mpl.__version__ >= '1.5' else
               'gridliner_labels_pre_mpl_1.5'])
def test_grid_labels():
    plt.figure(figsize=(8, 10))

    crs_pc = ccrs.PlateCarree()
    crs_merc = ccrs.Mercator()
    crs_osgb = ccrs.OSGB()

    ax = plt.subplot(3, 2, 1, projection=crs_pc)
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    # Check that adding labels to Mercator gridlines gives an error.
    # (Currently can only label PlateCarree gridlines.)
    ax = plt.subplot(3, 2, 2,
                     projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines()
    with assert_raises(TypeError):
        ax.gridlines(crs=crs_merc, draw_labels=True)

    ax.set_title('Known bug')
    gl = ax.gridlines(crs=crs_pc, draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlines = False

    ax = plt.subplot(3, 2, 3, projection=crs_merc)
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    # Check that labelling the gridlines on an OSGB plot gives an error.
    # (Currently can only draw these on PlateCarree or Mercator plots.)
    ax = plt.subplot(3, 2, 4, projection=crs_osgb)
    ax.coastlines()
    with assert_raises(TypeError):
        ax.gridlines(draw_labels=True)

    ax = plt.subplot(3, 2, 4, projection=crs_pc)
    ax.coastlines()
    gl = ax.gridlines(
        crs=crs_pc, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_bottom = True
    gl.ylabels_right = True
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-180, -45, 45, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red'}

    # trigger a draw at this point and check the appropriate artists are
    # populated on the gridliner instance
    FigureCanvasAgg(plt.gcf()).draw()

    assert len(gl.xlabel_artists) == 4
    assert len(gl.ylabel_artists) == 5
    assert len(gl.ylabel_artists) == 5
    assert len(gl.xline_artists) == 0

    ax = plt.subplot(3, 2, 5, projection=crs_pc)
    ax.set_extent([-20, 10.0, 45.0, 70.0])
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    ax = plt.subplot(3, 2, 6, projection=crs_merc)
    ax.set_extent([-20, 10.0, 45.0, 70.0], crs=crs_pc)
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    # Increase margins between plots to stop them bumping into one another.
    plt.subplots_adjust(wspace=0.25, hspace=0.25)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
