# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pytest
from shapely.geos import geos_version

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import (LATITUDE_FORMATTER, LONGITUDE_FORMATTER,
                                   classic_formatter, classic_locator)
from cartopy.mpl.ticker import LongitudeFormatter, LongitudeLocator


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
    pytest.param((ccrs.InterruptedGoodeHomolosine, dict(emphasis='land')),
                 id='InterruptedGoodeHomolosine'),
    pytest.param(
        (ccrs.RotatedPole,
         dict(pole_longitude=180.0, pole_latitude=36.0,
              central_rotated_longitude=-106.0)),
        id='RotatedPole'),
    pytest.param((ccrs.OSGB, dict(approx=False)), id='OSGB'),
    ccrs.EuroPP,
    ccrs.Geostationary,
    ccrs.NearsidePerspective,
    ccrs.Gnomonic,
    ccrs.LambertAzimuthalEqualArea,
    ccrs.NorthPolarStereo,
    pytest.param((ccrs.OSNI, dict(approx=False)), id='OSNI'),
    ccrs.SouthPolarStereo,
]


@pytest.mark.natural_earth
# Robinson projection is slightly better in Proj 6+.
@pytest.mark.mpl_image_compare(filename='gridliner1.png', tolerance=0.7)
def test_gridliner():
    ny, nx = 2, 4

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot(nx, ny, 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(linestyle=':')

    ax = fig.add_subplot(nx, ny, 2, projection=ccrs.OSGB(approx=False))
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(linestyle=':')

    ax = fig.add_subplot(nx, ny, 3, projection=ccrs.OSGB(approx=False))
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(ccrs.PlateCarree(), color='blue', linestyle='-')
    ax.gridlines(ccrs.OSGB(approx=False), linestyle=':')

    ax = fig.add_subplot(nx, ny, 4, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(ccrs.NorthPolarStereo(), alpha=0.5,
                 linewidth=1.5, linestyle='-')

    ax = fig.add_subplot(nx, ny, 5, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines(resolution="110m")
    osgb = ccrs.OSGB(approx=False)
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb, linestyle=':')

    ax = fig.add_subplot(nx, ny, 6, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(alpha=0.5, linewidth=1.5, linestyle='-')

    ax = fig.add_subplot(nx, ny, 7, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    ax.coastlines(resolution="110m")
    osgb = ccrs.OSGB(approx=False)
    ax.set_extent(tuple(osgb.x_limits) + tuple(osgb.y_limits), crs=osgb)
    ax.gridlines(osgb, linestyle=':')

    ax = fig.add_subplot(nx, ny, 8,
                         projection=ccrs.Robinson(central_longitude=135))
    ax.set_global()
    ax.coastlines(resolution="110m")
    ax.gridlines(ccrs.PlateCarree(), alpha=0.5, linewidth=1.5, linestyle='-')

    delta = 1.5e-2
    fig.subplots_adjust(left=0 + delta, right=1 - delta,
                        top=1 - delta, bottom=0 + delta)
    return fig


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
grid_label_tol = 3.9


@pytest.mark.skipif(geos_version == (3, 9, 0), reason="GEOS intersection bug")
@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='gridliner_labels.png',
                               tolerance=grid_label_tol)
def test_grid_labels():
    fig = plt.figure(figsize=(10, 10))

    crs_pc = ccrs.PlateCarree()
    crs_merc = ccrs.Mercator()

    ax = fig.add_subplot(3, 2, 1, projection=crs_pc)
    ax.coastlines(resolution="110m")
    ax.gridlines(draw_labels=True)

    # Check that adding labels to Mercator gridlines gives an error.
    # (Currently can only label PlateCarree gridlines.)
    ax = fig.add_subplot(
        3, 2, 2, projection=ccrs.PlateCarree(central_longitude=180))
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

    ax = fig.add_subplot(3, 2, 4, projection=crs_pc)
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
    assert len([
        lb for lb in gl.bottom_label_artists if lb.get_visible()]) == 4
    assert len([
        lb for lb in gl.top_label_artists if lb.get_visible()]) == 0
    assert len([
        lb for lb in gl.left_label_artists if lb.get_visible()]) == 0
    assert len([
        lb for lb in gl.right_label_artists if lb.get_visible()]) != 0
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
    fig.subplots_adjust(wspace=0.25, hspace=0.25)

    return fig


@pytest.mark.skipif(geos_version == (3, 9, 0), reason="GEOS intersection bug")
@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='gridliner_labels_tight.png',
                               tolerance=2.92)
def test_grid_labels_tight():
    # Ensure tight layout accounts for gridlines
    fig = plt.figure(figsize=(7, 5))

    crs_pc = ccrs.PlateCarree()
    crs_merc = ccrs.Mercator()

    ax = fig.add_subplot(2, 2, 1, projection=crs_pc)
    ax.coastlines(resolution="110m")
    ax.gridlines(draw_labels=True)

    ax = fig.add_subplot(2, 2, 2, projection=crs_merc)
    ax.coastlines(resolution="110m")
    ax.gridlines(draw_labels=True)

    # Matplotlib tight layout is also incorrect if cartopy fails
    # to adjust aspect ratios first. Relevant when aspect ratio has
    # changed due to set_extent.
    ax = fig.add_subplot(2, 2, 3, projection=crs_pc)
    ax.set_extent([-20, 10.0, 45.0, 70.0])
    ax.coastlines(resolution="110m")
    ax.gridlines(draw_labels=True)

    ax = fig.add_subplot(2, 2, 4, projection=crs_merc)
    ax.set_extent([-20, 10.0, 45.0, 70.0], crs=crs_pc)
    ax.coastlines(resolution="110m")
    gl = ax.gridlines(draw_labels=True)
    gl.rotate_labels = False

    # Apply tight layout
    fig.tight_layout()

    # Ensure gridliners were drawn
    for ax in fig.axes:
        for gl in ax._gridliners:
            assert hasattr(gl, '_drawn') and gl._drawn

    return fig


@pytest.mark.skipif(geos_version == (3, 9, 0), reason="GEOS intersection bug")
@pytest.mark.natural_earth
@pytest.mark.parametrize('proj', TEST_PROJS)
@pytest.mark.mpl_image_compare(style='mpl20')
def test_grid_labels_inline(proj):
    fig = plt.figure()
    if isinstance(proj, tuple):
        proj, kwargs = proj
    else:
        kwargs = {}
    ax = fig.add_subplot(projection=proj(**kwargs))
    ax.gridlines(draw_labels=True, auto_inline=True)
    ax.coastlines(resolution="110m")
    return fig


@pytest.mark.skipif(geos_version == (3, 9, 0), reason="GEOS intersection bug")
@pytest.mark.natural_earth
@pytest.mark.parametrize('proj', TEST_PROJS)
@pytest.mark.mpl_image_compare(style='mpl20', tolerance=0.79)
def test_grid_labels_inline_usa(proj):
    top = 49.3457868  # north lat
    left = -124.7844079  # west long
    right = -66.9513812  # east long
    bottom = 24.7433195  # south lat
    fig = plt.figure()
    if isinstance(proj, tuple):
        proj, kwargs = proj
    else:
        kwargs = {}
    ax = fig.add_subplot(projection=proj(**kwargs))
    try:
        ax.set_extent([left, right, bottom, top],
                      crs=ccrs.PlateCarree())
    except Exception:
        pytest.skip('Projection does not support changing extent')
    ax.gridlines(draw_labels=True, auto_inline=True, clip_on=True)
    ax.coastlines(resolution="110m")
    return fig


@pytest.mark.skipif(geos_version == (3, 9, 0), reason="GEOS intersection bug")
@pytest.mark.mpl_image_compare(filename='gridliner_labels_bbox_style.png',
                               tolerance=grid_label_tol)
def test_gridliner_labels_bbox_style():
    top = 49.3457868  # north lat
    left = -124.7844079  # west long
    right = -66.9513812  # east long
    bottom = 24.7433195  # south lat

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.coastlines(resolution="110m")
    ax.set_extent([left, right, bottom, top],
                  crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True)

    gl.labels_bbox_style = {
        "pad": 0,
        "visible": True,
        "facecolor": "white",
        "edgecolor": "black",
        "boxstyle": "round, pad=0.2",
    }

    return fig


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
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    gl = ax.gridlines(crs=gcrs, draw_labels=False, xlocs=xloc, xformatter=xfmt)
    assert isinstance(gl.xlocator, xloc_expected)
    assert isinstance(gl.xformatter, xfmt_expected)


def test_gridliner_line_limits():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
    ax.set_global()
    # Test a single value passed in which represents (-lim, lim)
    xlim, ylim = 125, 75
    gl = ax.gridlines(xlim=xlim, ylim=ylim)
    fig.canvas.draw_idle()

    paths = gl.xline_artists[0].get_paths() + gl.yline_artists[0].get_paths()
    for path in paths:
        assert (np.min(path.vertices, axis=0) >= (-xlim, -ylim)).all()
        assert (np.max(path.vertices, axis=0) <= (xlim, ylim)).all()

    # Test a pair of values passed in which represents the min, max
    xlim = (-125, 150)
    ylim = (50, 70)
    gl = ax.gridlines(xlim=xlim, ylim=ylim)
    fig.canvas.draw_idle()

    paths = gl.xline_artists[0].get_paths() + gl.yline_artists[0].get_paths()
    for path in paths:
        assert (np.min(path.vertices, axis=0) >= (xlim[0], ylim[0])).all()
        assert (np.max(path.vertices, axis=0) <= (xlim[1], ylim[1])).all()


@pytest.mark.parametrize(
    "draw_labels, result",
    [
        (True,
         {'left': ['40°N'],
          'right': ['40°N', '50°N'],
          'top': ['70°E', '100°E', '130°E'],
          'bottom': ['100°E']}),
        (False,
         {'left': [],
          'right': [],
          'top': [],
          'bottom': []}),
        (['top', 'left'],
         {'left': ['40°N'],
          'right': [],
          'top': ['70°E', '100°E', '130°E'],
          'bottom': []}),
        ({'top': 'x', 'right': 'y'},
         {'left': [],
          'right': ['40°N', '50°N'],
          'top': ['70°E', '100°E', '130°E'],
          'bottom': []}),
        ({'left': 'x'},
         {'left': ['70°E'],
          'right': [],
          'top': [],
          'bottom': []}),
        ({'top': 'y'},
         {'left': [],
          'right': [],
          'top': ['50°N'],
          'bottom': []}),
    ])
def test_gridliner_draw_labels_param(draw_labels, result):
    fig = plt.figure()
    lambert_crs = ccrs.LambertConformal(central_longitude=105)
    ax = fig.add_subplot(projection=lambert_crs)
    ax.set_extent([75, 130, 18, 54], crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=draw_labels, rotate_labels=False, dms=True,
                      x_inline=False, y_inline=False)
    gl.xlocator = mticker.FixedLocator([70, 100, 130])
    gl.ylocator = mticker.FixedLocator([40, 50])
    fig.canvas.draw()
    res = {}
    for loc in 'left', 'right', 'top', 'bottom':
        artists = getattr(gl, f'{loc}_label_artists')
        res[loc] = [a.get_text() for a in artists if a.get_visible()]
    assert res == result


def test_gridliner_formatter_kwargs():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([-80, -40.0, 10.0, -30.0])
    gl = ax.gridlines(draw_labels=True, dms=False,
                      formatter_kwargs=dict(cardinal_labels={'west': 'O'}))
    fig.canvas.draw()
    labels = [a.get_text() for a in gl.bottom_label_artists if a.get_visible()]
    assert labels == ['75°O', '70°O', '65°O', '60°O', '55°O', '50°O', '45°O']
