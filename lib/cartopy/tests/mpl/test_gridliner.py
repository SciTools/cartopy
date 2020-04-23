# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from __future__ import (absolute_import, division, print_function)

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pytest

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeLocator, LongitudeFormatter
from cartopy.mpl.gridliner import (
    LATITUDE_FORMATTER, LONGITUDE_FORMATTER,
    classic_locator, classic_formatter)


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
    ccrs.InterruptedGoodeHomolosine,
    (ccrs.RotatedPole,
     dict(pole_longitude=180.0,
          pole_latitude=36.0,
          central_rotated_longitude=-106.0,
          globe=ccrs.Globe(semimajor_axis=6370000,
                           semiminor_axis=6370000))),
    (ccrs.OSGB, dict(approx=False)),
    ccrs.EuroPP,
    ccrs.Geostationary,
    ccrs.NearsidePerspective,
    ccrs.Gnomonic,
    ccrs.LambertAzimuthalEqualArea,
    ccrs.NorthPolarStereo,
    (ccrs.OSNI, dict(approx=False)),
    ccrs.SouthPolarStereo,
]


@pytest.mark.natural_earth
@ImageTesting(['gridliner1'],
              # Robinson projection is slightly better in Proj 6+.
              tolerance=0.7 if ccrs.PROJ4_VERSION >= (6, 0, 0) else 0.5)
def test_gridliner():
    ny, nx = 2, 4

    plt.figure(figsize=(10, 10))

    ax = plt.subplot(nx, ny, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(linestyle=':')

    ax = plt.subplot(nx, ny, 2, projection=ccrs.OSGB(approx=False))
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(linestyle=':')

    ax = plt.subplot(nx, ny, 3, projection=ccrs.OSGB(approx=False))
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(ccrs.PlateCarree(), color='blue', linestyle='-')
    ax.gridlines(ccrs.OSGB(approx=False), linestyle=':')

    ax = plt.subplot(nx, ny, 4, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(ccrs.NorthPolarStereo(), alpha=0.5,
                 linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 5, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(resolution="110m")
    osgb = ccrs.OSGB(approx=False)
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb, linestyle=':')

    ax = plt.subplot(nx, ny, 6, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(alpha=0.5, linewidth=1.5, linestyle='-')

    ax = plt.subplot(nx, ny, 7, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines(resolution="110m")
    osgb = ccrs.OSGB(approx=False)
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb, linestyle=':')

    ax = plt.subplot(nx, ny, 8,
                     projection=ccrs.Robinson(central_longitude=135))
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(ccrs.PlateCarree(), alpha=0.5, linewidth=1.5, linestyle='-')

    delta = 1.5e-2
    plt.subplots_adjust(left=0 + delta, right=1 - delta,
                        top=1 - delta, bottom=0 + delta)


def test_gridliner_specified_lines():
    meridians = [0, 60, 120, 180, 240, 360]
    parallels = [-90, -60, -30, 0, 30, 60, 90]

    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl = GeoAxes.gridlines(ax, xlocs=meridians, ylocs=parallels)
    assert isinstance(gl.xlocator, mticker.FixedLocator)
    assert isinstance(gl.ylocator, mticker.FixedLocator)
    assert gl.xlocator.tick_values(None, None).tolist() == meridians
    assert gl.ylocator.tick_values(None, None).tolist() == parallels


# The tolerance on these tests are particularly high because of the high number
# of text objects. A new testing strategy is needed for this kind of test.
grid_label_tol = grid_label_inline_tol = grid_label_inline_usa_tol = 0.5
if MPL_VERSION >= '2.0':
    grid_label_image = 'gridliner_labels'
    if ccrs.PROJ4_VERSION < (4, 9, 3):
        # A 0-longitude label is missing on older Proj versions.
        grid_label_tol = 1.8
    grid_label_inline_image = 'gridliner_labels_inline'
    grid_label_inline_usa_image = 'gridliner_labels_inline_usa'
    if ccrs.PROJ4_VERSION == (4, 9, 1):
        # AzimuthalEquidistant was previously broken.
        grid_label_inline_tol = 7.9
        grid_label_inline_usa_tol = 7.7
    elif ccrs.PROJ4_VERSION < (5, 0, 0):
        # Stereographic was previously broken.
        grid_label_inline_tol = 6.4
        grid_label_inline_usa_tol = 4.0
else:
    grid_label_image = 'gridliner_labels_1.5'
    grid_label_tol = 1.8
    grid_label_inline_image = 'gridliner_labels_inline_1.5'
    grid_label_inline_usa_image = 'gridliner_labels_inline_usa_1.5'
    if ccrs.PROJ4_VERSION >= (5, 0, 0):
        # Stereographic was fixed, but test image was not updated.
        grid_label_inline_tol = 7.9
        grid_label_inline_usa_tol = 7.9
    elif ccrs.PROJ4_VERSION >= (4, 9, 2):
        # AzimuthalEquidistant was fixed, but test image was not updated.
        grid_label_inline_tol = 5.4
        grid_label_inline_usa_tol = 7.2
if (5, 0, 0) <= ccrs.PROJ4_VERSION < (5, 1, 0):
    # Several projections are broken in these versions, so not plotted.
    grid_label_inline_tol += 5.1
    grid_label_inline_usa_tol += 5.5
elif (6, 0, 0) <= ccrs.PROJ4_VERSION:
    # Better Robinson projection causes some text movement.
    grid_label_inline_tol += 1.2


@pytest.mark.natural_earth
@ImageTesting([grid_label_image], tolerance=grid_label_tol)
def test_grid_labels():
    fig = plt.figure(figsize=(10, 10))

    crs_pc = ccrs.PlateCarree()
    crs_merc = ccrs.Mercator()

    ax = fig.add_subplot(3, 2, 1, projection=crs_pc)
    ax.coastlines(resolution="110m")
    ax.gridlines(draw_labels=True)

    # Check that adding labels to Mercator gridlines gives an error.
    # (Currently can only label PlateCarree gridlines.)
    ax = fig.add_subplot(3, 2, 2,
                         projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines(resolution="110m")

    ax.set_title('Known bug')
    gl = ax.gridlines(crs=crs_pc, draw_labels=True)
    gl.top_labels = False
    gl.left_labels = False
    gl.xlines = False

    ax = fig.add_subplot(3, 2, 3, projection=crs_merc)
    ax.coastlines(resolution="110m")
    gl = ax.gridlines(draw_labels=True)
    gl.xlabel_style = gl.ylabel_style = {'size': 9}

    ax = plt.subplot(3, 2, 4, projection=crs_pc)
    ax.coastlines(resolution="110m")
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
    ax.coastlines(resolution="110m")
    ax.gridlines(draw_labels=True)

    ax = fig.add_subplot(3, 2, 6, projection=crs_merc)
    ax.set_extent([-20, 10.0, 45.0, 70.0], crs=crs_pc)
    ax.coastlines(resolution="110m")
    gl = ax.gridlines(draw_labels=True)
    gl.rotate_labels = False
    gl.xlabel_style = gl.ylabel_style = {'size': 9}

    # Increase margins between plots to stop them bumping into one another.
    plt.subplots_adjust(wspace=0.25, hspace=0.25)


@pytest.mark.natural_earth
@ImageTesting([grid_label_inline_image], tolerance=grid_label_inline_tol)
def test_grid_labels_inline():
    plt.figure(figsize=(35, 35))
    for i, proj in enumerate(TEST_PROJS, 1):
        if isinstance(proj, tuple):
            proj, kwargs = proj
        else:
            kwargs = {}
        ax = plt.subplot(7, 4, i, projection=proj(**kwargs))
        if (ccrs.PROJ4_VERSION[:2] == (5, 0) and
                proj in (ccrs.Orthographic, ccrs.AlbersEqualArea,
                         ccrs.Geostationary, ccrs.NearsidePerspective)):
            # Above projections are broken, so skip labels.
            # Add gridlines anyway to minimize image differences.
            ax.gridlines()
        else:
            ax.gridlines(draw_labels=True, auto_inline=True)
        ax.coastlines(resolution="110m")
        ax.set_title(proj, y=1.075)
    plt.subplots_adjust(wspace=0.35, hspace=0.35)


@pytest.mark.natural_earth
@ImageTesting([grid_label_inline_usa_image],
              tolerance=grid_label_inline_usa_tol)
def test_grid_labels_inline_usa():
    top = 49.3457868  # north lat
    left = -124.7844079  # west long
    right = -66.9513812  # east long
    bottom = 24.7433195  # south lat
    plt.figure(figsize=(35, 35))
    for i, proj in enumerate(TEST_PROJS, 1):
        if isinstance(proj, tuple):
            proj, kwargs = proj
        else:
            kwargs = {}
        ax = plt.subplot(7, 4, i, projection=proj(**kwargs))
        try:
            ax.set_extent([left, right, bottom, top],
                          crs=ccrs.PlateCarree())
        except Exception:
            pass
        ax.set_title(proj, y=1.075)
        if (ccrs.PROJ4_VERSION[:2] == (5, 0) and
                proj in (ccrs.Orthographic, ccrs.AlbersEqualArea,
                         ccrs.Geostationary, ccrs.NearsidePerspective)):
            # Above projections are broken, so skip labels.
            # Add gridlines anyway to minimize image differences.
            ax.gridlines()
        else:
            ax.gridlines(draw_labels=True, auto_inline=True, clip_on=True)
        ax.coastlines(resolution="110m")
    plt.subplots_adjust(wspace=0.35, hspace=0.35)


@pytest.mark.parametrize(
    "proj,gcrs,xloc,xfmt,xloc_expected,xfmt_expected",
    [
       (ccrs.PlateCarree(), ccrs.PlateCarree(),
        [10, 20], None, mticker.FixedLocator, LongitudeFormatter),
       (ccrs.PlateCarree(), ccrs.Mercator(),
        [10, 20], None, mticker.FixedLocator, classic_formatter),
       (ccrs.PlateCarree(), ccrs.PlateCarree(),
        mticker.MaxNLocator(nbins=9), None,
        mticker.MaxNLocator, LongitudeFormatter),
       (ccrs.PlateCarree(), ccrs.Mercator(),
        mticker.MaxNLocator(nbins=9), None,
        mticker.MaxNLocator, classic_formatter),
       (ccrs.PlateCarree(), ccrs.PlateCarree(),
        None, None, LongitudeLocator, LongitudeFormatter),
       (ccrs.PlateCarree(), ccrs.Mercator(),
        None, None, classic_locator.__class__, classic_formatter),
       (ccrs.PlateCarree(), ccrs.PlateCarree(),
        None, mticker.StrMethodFormatter('{x}'),
        LongitudeLocator, mticker.StrMethodFormatter),
     ])
def test_gridliner_default_fmtloc(
        proj, gcrs, xloc, xfmt, xloc_expected, xfmt_expected):
    plt.figure()
    ax = plt.subplot(111, projection=proj)
    gl = ax.gridlines(crs=gcrs, draw_labels=False, xlocs=xloc, xformatter=xfmt)
    plt.close()
    assert isinstance(gl.xlocator, xloc_expected)
    assert isinstance(gl.xformatter, xfmt_expected)
