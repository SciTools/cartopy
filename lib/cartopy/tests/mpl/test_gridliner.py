# (C) British Crown Copyright 2011 - 2014, Met Office
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

import unittest
import warnings

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from nose.tools import assert_raises

import cartopy.crs as ccrs
from cartopy.tests.mpl import ImageTesting
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER


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


# The tolerance on this test is particularly high because of the high number
# of text objects. A new testing stratergy is needed for this kind of test.
@ImageTesting(['gridliner_labels'], tolerance=2.6)
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
    plt.draw()
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
