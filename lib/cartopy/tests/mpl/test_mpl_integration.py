# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import re
import warnings

import numpy as np
import matplotlib.pyplot as plt
from packaging.version import parse as parse_version
import pytest

import cartopy.crs as ccrs

from cartopy.tests.mpl import MPL_VERSION, ImageTesting


# This is due to a change in MPL 3.5 contour line paths changing
# ever so slightly.
contour_tol = 2.25
@pytest.mark.natural_earth
@ImageTesting(['global_contour_wrap'], style='mpl20',
              tolerance=contour_tol)
def test_global_contour_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contour(x, y, data, transform=ccrs.PlateCarree())


@pytest.mark.natural_earth
@ImageTesting(['global_contour_wrap'], style='mpl20',
              tolerance=contour_tol)
def test_global_contour_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contour(x, y, data)


@pytest.mark.natural_earth
@ImageTesting(['global_contourf_wrap'])
def test_global_contourf_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contourf(x, y, data, transform=ccrs.PlateCarree())


@pytest.mark.natural_earth
@ImageTesting(['global_contourf_wrap'])
def test_global_contourf_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contourf(x, y, data)


@pytest.mark.natural_earth
@ImageTesting(['global_pcolor_wrap'])
def test_global_pcolor_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))[:-1, :-1]
    plt.pcolor(x, y, data, transform=ccrs.PlateCarree())


@pytest.mark.natural_earth
@ImageTesting(['global_pcolor_wrap'])
def test_global_pcolor_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))[:-1, :-1]
    plt.pcolor(x, y, data)


@pytest.mark.natural_earth
@ImageTesting(['global_scatter_wrap'])
def test_global_scatter_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    # By default the coastline feature will be drawn after patches.
    # By setting zorder we can ensure our scatter points are drawn
    # after the coastlines.
    ax.coastlines(zorder=0)
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.scatter(x, y, c=data, transform=ccrs.PlateCarree())


@pytest.mark.natural_earth
@ImageTesting(['global_scatter_wrap'])
def test_global_scatter_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(zorder=0)
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.scatter(x, y, c=data)


@pytest.mark.natural_earth
@ImageTesting(['global_hexbin_wrap'],
              tolerance=2 if MPL_VERSION < parse_version('3.2') else 0.5)
def test_global_hexbin_wrap():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(zorder=2)
    x, y = np.meshgrid(np.arange(-179, 181), np.arange(-90, 91))
    data = np.sin(np.sqrt(x**2 + y**2))
    plt.hexbin(
        x.flatten(),
        y.flatten(),
        C=data.flatten(),
        gridsize=20,
        zorder=1,
    )


@pytest.mark.natural_earth
@ImageTesting(['global_hexbin_wrap'],
              tolerance=2 if MPL_VERSION < parse_version('3.2') else 0.5)
def test_global_hexbin_wrap_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(zorder=2)
    x, y = np.meshgrid(np.arange(0, 360), np.arange(-90, 91))
    # wrap values so to match x values from test_global_hexbin_wrap
    x_wrap = np.where(x >= 180, x-360, x)
    data = np.sin(np.sqrt(x_wrap**2 + y**2))
    plt.hexbin(
        x.flatten(),
        y.flatten(),
        C=data.flatten(),
        gridsize=20,
        zorder=1,
    )


@ImageTesting(['global_map'],
              tolerance=0.55)
def test_global_map():
    plt.axes(projection=ccrs.Robinson())
#    ax.coastlines()
#    ax.gridlines(5)

    plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())

    plt.plot([-0.08, 132], [51.53, 43.17], color='red',
             transform=ccrs.PlateCarree())

    plt.plot([-0.08, 132], [51.53, 43.17], color='blue',
             transform=ccrs.Geodetic())


@pytest.mark.filterwarnings("ignore:Unable to determine extent")
@pytest.mark.natural_earth
@ImageTesting(['simple_global'])
def test_simple_global():
    plt.axes(projection=ccrs.PlateCarree())
    plt.gca().coastlines()
    # produces a global map, despite not having needed to set the limits


@pytest.mark.filterwarnings("ignore:Unable to determine extent")
@pytest.mark.natural_earth
@ImageTesting(['multiple_projections5'],
              tolerance=2.05)
def test_multiple_projections():

    projections = [ccrs.PlateCarree(),
                   ccrs.Robinson(),
                   ccrs.RotatedPole(pole_latitude=45, pole_longitude=180),
                   ccrs.OSGB(approx=True),
                   ccrs.TransverseMercator(approx=True),
                   ccrs.Mercator(min_latitude=-85., max_latitude=85.),
                   ccrs.LambertCylindrical(),
                   ccrs.Miller(),
                   ccrs.Gnomonic(),
                   ccrs.Stereographic(),
                   ccrs.NorthPolarStereo(),
                   ccrs.SouthPolarStereo(),
                   ccrs.Orthographic(),
                   ccrs.Mollweide(),
                   ccrs.InterruptedGoodeHomolosine(emphasis='land'),
                   ccrs.EckertI(),
                   ccrs.EckertII(),
                   ccrs.EckertIII(),
                   ccrs.EckertIV(),
                   ccrs.EckertV(),
                   ccrs.EckertVI(),
                   ]

    rows = np.ceil(len(projections) / 5).astype(int)

    fig = plt.figure(figsize=(10, 2 * rows))
    for i, prj in enumerate(projections, 1):
        ax = fig.add_subplot(rows, 5, i, projection=prj)

        ax.set_global()

        ax.coastlines(resolution="110m")

        plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
        plt.plot([-0.08, 132], [51.53, 43.17], color='red',
                 transform=ccrs.PlateCarree())
        plt.plot([-0.08, 132], [51.53, 43.17], color='blue',
                 transform=prj.as_geodetic())


@pytest.mark.natural_earth
@ImageTesting(['multiple_projections520'], tolerance=0.65)
def test_multiple_projections_520():
    # Test projections added in Proj 5.2.0.

    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.EqualEarth())

    ax.set_global()

    ax.coastlines()

    ax.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())

    ax.plot([-0.08, 132], [51.53, 43.17], color='red',
            transform=ccrs.PlateCarree())

    ax.plot([-0.08, 132], [51.53, 43.17], color='blue',
            transform=ccrs.Geodetic())


def test_cursor_values():
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    x, y = np.array([-969100., -4457000.])
    r = ax.format_coord(x, y)
    assert (r.encode('ascii', 'ignore') ==
            b'-9.691e+05, -4.457e+06 (50.716617N, 12.267069W)')

    ax = plt.axes(projection=ccrs.PlateCarree())
    x, y = np.array([-181.5, 50.])
    r = ax.format_coord(x, y)
    assert (r.encode('ascii', 'ignore') ==
            b'-181.5, 50 (50.000000N, 178.500000E)')

    ax = plt.axes(projection=ccrs.Robinson())
    x, y = np.array([16060595.2, 2363093.4])
    r = ax.format_coord(x, y)
    assert re.search(b'1.606e\\+07, 2.363e\\+06 '
                     b'\\(22.09[0-9]{4}N, 173.70[0-9]{4}E\\)',
                     r.encode('ascii', 'ignore'))

    plt.close()


@pytest.mark.natural_earth
@ImageTesting(['natural_earth_interface'], tolerance=1.21)
def test_axes_natural_earth_interface():
    rob = ccrs.Robinson()

    ax = plt.axes(projection=rob)

    with warnings.catch_warnings(record=True) as all_warnings:
        warnings.simplefilter('always')

        ax.natural_earth_shp('rivers_lake_centerlines', edgecolor='black',
                             facecolor='none')
        ax.natural_earth_shp('lakes', facecolor='blue')

    assert len(all_warnings) == 2
    for warning in all_warnings:
        msg = str(warning.message)
        assert 'deprecated' in msg
        assert 'add_feature' in msg


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_global_wrap1'], tolerance=1.27)
def test_pcolormesh_global_with_wrap1():
    # make up some realistic data with bounds (such as data from the UM)
    nx, ny = 36, 18
    xbnds = np.linspace(0, 360, nx, endpoint=True)
    ybnds = np.linspace(-90, 90, ny, endpoint=True)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(212, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


def test_pcolormesh_get_array_with_mask():
    # make up some realistic data with bounds (such as data from the UM)
    nx, ny = 36, 18
    xbnds = np.linspace(0, 360, nx, endpoint=True)
    ybnds = np.linspace(-90, 90, ny, endpoint=True)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    c = plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    assert c._wrapped_collection_fix is not None, \
        'No pcolormesh wrapping was done when it should have been.'

    assert np.array_equal(data.ravel(), c.get_array()), \
        'Data supplied does not match data retrieved in wrapped case'

    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    # Case without wrapping
    nx, ny = 36, 18
    xbnds = np.linspace(-60, 60, nx, endpoint=True)
    ybnds = np.linspace(-80, 80, ny, endpoint=True)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data2 = data[:-1, :-1]

    ax = plt.subplot(212, projection=ccrs.PlateCarree())
    c = plt.pcolormesh(xbnds, ybnds, data2, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    assert getattr(c, "_wrapped_collection_fix", None) is None, \
        'pcolormesh wrapping was done when it should not have been.'

    assert np.array_equal(data2.ravel(), c.get_array()), \
        'Data supplied does not match data retrieved in unwrapped case'


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_global_wrap2'], tolerance=1.87)
def test_pcolormesh_global_with_wrap2():
    # make up some realistic data with bounds (such as data from the UM)
    nx, ny = 36, 18
    xbnds, xstep = np.linspace(0, 360, nx - 1, retstep=True, endpoint=True)
    ybnds, ystep = np.linspace(-90, 90, ny - 1, retstep=True, endpoint=True)
    xbnds -= xstep / 2
    ybnds -= ystep / 2
    xbnds = np.append(xbnds, xbnds[-1] + xstep)
    ybnds = np.append(ybnds, ybnds[-1] + ystep)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(212, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_global_wrap3'], tolerance=1.42)
def test_pcolormesh_global_with_wrap3():
    nx, ny = 33, 17
    xbnds = np.linspace(-1.875, 358.125, nx, endpoint=True)
    ybnds = np.linspace(91.25, -91.25, ny, endpoint=True)
    xbnds, ybnds = np.meshgrid(xbnds, ybnds)

    data = np.exp(np.sin(np.deg2rad(xbnds)) + np.cos(np.deg2rad(ybnds)))

    # this step is not necessary, but makes the plot even harder to do (i.e.
    # it really puts cartopy through its paces)
    ybnds = np.append(ybnds, ybnds[:, 1:2], axis=1)
    xbnds = np.append(xbnds, xbnds[:, 1:2] + 360, axis=1)
    data = np.ma.concatenate([data, data[:, 0:1]], axis=1)

    data = data[:-1, :-1]
    data = np.ma.masked_greater(data, 2.6)

    ax = plt.subplot(311, projection=ccrs.PlateCarree(-45))
    c = plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                       snap=False)
    assert c._wrapped_collection_fix is not None, \
        'No pcolormesh wrapping was done when it should have been.'

    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(312, projection=ccrs.PlateCarree(-1.87499952))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(313, projection=ccrs.Robinson(-2))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_global_wrap3'], tolerance=1.42)
def test_pcolormesh_set_array_with_mask():
    """Testing that set_array works with masked arrays properly."""
    nx, ny = 33, 17
    xbnds = np.linspace(-1.875, 358.125, nx, endpoint=True)
    ybnds = np.linspace(91.25, -91.25, ny, endpoint=True)
    xbnds, ybnds = np.meshgrid(xbnds, ybnds)

    data = np.exp(np.sin(np.deg2rad(xbnds)) + np.cos(np.deg2rad(ybnds)))

    # this step is not necessary, but makes the plot even harder to do (i.e.
    # it really puts cartopy through its paces)
    ybnds = np.append(ybnds, ybnds[:, 1:2], axis=1)
    xbnds = np.append(xbnds, xbnds[:, 1:2] + 360, axis=1)
    data = np.ma.concatenate([data, data[:, 0:1]], axis=1)

    data = data[:-1, :-1]
    data = np.ma.masked_greater(data, 2.6)
    norm = plt.Normalize(np.min(data), np.max(data))
    bad_data = np.ones(data.shape)
    # Start with the opposite mask and then swap back in the set_array call
    bad_data_mask = np.ma.array(bad_data, mask=~data.mask)

    ax = plt.subplot(311, projection=ccrs.PlateCarree(-45))
    c = plt.pcolormesh(xbnds, ybnds, bad_data,
                       norm=norm, transform=ccrs.PlateCarree(),
                       snap=False)
    c.set_array(data.ravel())
    assert c._wrapped_collection_fix is not None, \
        'No pcolormesh wrapping was done when it should have been.'

    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(312, projection=ccrs.PlateCarree(-1.87499952))
    c2 = plt.pcolormesh(xbnds, ybnds, bad_data_mask,
                        norm=norm, transform=ccrs.PlateCarree(),
                        snap=False)
    c2.set_array(data.ravel())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(313, projection=ccrs.Robinson(-2))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_global_wrap3'], tolerance=1.42)
def test_pcolormesh_set_clim_with_mask():
    """Testing that set_clim works with masked arrays properly."""
    nx, ny = 33, 17
    xbnds = np.linspace(-1.875, 358.125, nx, endpoint=True)
    ybnds = np.linspace(91.25, -91.25, ny, endpoint=True)
    xbnds, ybnds = np.meshgrid(xbnds, ybnds)

    data = np.exp(np.sin(np.deg2rad(xbnds)) + np.cos(np.deg2rad(ybnds)))

    # this step is not necessary, but makes the plot even harder to do (i.e.
    # it really puts cartopy through its paces)
    ybnds = np.append(ybnds, ybnds[:, 1:2], axis=1)
    xbnds = np.append(xbnds, xbnds[:, 1:2] + 360, axis=1)
    data = np.ma.concatenate([data, data[:, 0:1]], axis=1)

    data = data[:-1, :-1]
    data = np.ma.masked_greater(data, 2.6)

    bad_initial_norm = plt.Normalize(-100, 100)

    ax = plt.subplot(311, projection=ccrs.PlateCarree(-45))
    c = plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                       norm=bad_initial_norm, snap=False)
    assert c._wrapped_collection_fix is not None, \
        'No pcolormesh wrapping was done when it should have been.'

    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(312, projection=ccrs.PlateCarree(-1.87499952))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(313, projection=ccrs.Robinson(-2))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree(),
                   snap=False)
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    # Fix clims on c so that test passes
    c.set_clim(data.min(), data.max())


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_limited_area_wrap'], tolerance=1.82)
def test_pcolormesh_limited_area_wrap():
    # make up some realistic data with bounds (such as data from the UM's North
    # Atlantic Europe model)
    nx, ny = 22, 36
    xbnds = np.linspace(311.91998291, 391.11999512, nx, endpoint=True)
    ybnds = np.linspace(-23.59000015, 24.81000137, ny, endpoint=True)
    x, y = np.meshgrid(xbnds, ybnds)
    data = ((np.sin(np.deg2rad(x))) / 10. + np.exp(np.cos(np.deg2rad(y))))
    data = data[:-1, :-1]

    rp = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)

    plt.figure(figsize=(10, 6))

    ax = plt.subplot(221, projection=ccrs.PlateCarree())
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Spectral',
                   snap=False)
    ax.coastlines()

    ax = plt.subplot(222, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Spectral',
                   snap=False)
    ax.coastlines()
    ax.set_global()

    # draw the same plot, only more zoomed in, and using the 2d versions
    # of the coordinates (just to test that 1d and 2d are both suitably
    # being fixed)
    ax = plt.subplot(223, projection=ccrs.PlateCarree())
    plt.pcolormesh(x, y, data, transform=rp, cmap='Spectral',
                   snap=False)
    ax.coastlines()
    ax.set_extent([-70, 0, 0, 80])

    ax = plt.subplot(224, projection=rp)
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Spectral',
                   snap=False)
    ax.coastlines()


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_single_column_wrap'], tolerance=0.7)
def test_pcolormesh_single_column_wrap():
    # Check a wrapped mesh like test_pcolormesh_limited_area_wrap, but only use
    # a single column, which could break depending on how wrapping is
    # determined.
    ny = 36
    xbnds = np.array([360.9485619, 364.71999105])
    ybnds = np.linspace(-23.59000015, 24.81000137, ny, endpoint=True)
    x, y = np.meshgrid(xbnds, ybnds)
    data = ((np.sin(np.deg2rad(x))) / 10. + np.exp(np.cos(np.deg2rad(y))))
    data = data[:-1, :-1]

    rp = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)

    plt.figure(figsize=(10, 6))

    ax = plt.subplot(111, projection=ccrs.PlateCarree(180))
    # TODO: Remove snap when updating this image
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Spectral',
                   snap=False)
    ax.coastlines()
    ax.set_global()


def test_pcolormesh_diagonal_wrap():
    # Check that a cell with the top edge on one side of the domain
    # and the bottom edge on the other gets wrapped properly
    xs = [[160, 170], [190, 200]]
    ys = [[-10, -10], [10, 10]]
    zs = [[0]]

    ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(xs, ys, zs)

    # And that the wrapped_collection is added
    assert hasattr(mesh, "_wrapped_collection_fix")


def test_pcolormesh_nan_wrap():
    # Check that data with nan's as input still creates
    # the proper number of pcolor cells and those aren't
    # masked in the process.
    xs, ys = np.meshgrid([120, 160, 200], [-30, 0, 30])
    data = np.ones((2, 2)) * np.nan

    ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(xs, ys, data)
    pcolor = getattr(mesh, "_wrapped_collection_fix")
    assert len(pcolor.get_paths()) == 2


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_goode_wrap'])
def test_pcolormesh_goode_wrap():
    # global data on an Interrupted Goode Homolosine projection
    # shouldn't spill outside projection boundary
    x = np.linspace(0, 360, 73)
    y = np.linspace(-87.5, 87.5, 36)
    X, Y = np.meshgrid(*[np.deg2rad(c) for c in (x, y)])
    Z = np.cos(Y) + 0.375 * np.sin(2. * X)
    Z = Z[:-1, :-1]
    ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine(emphasis='land'))
    ax.coastlines()
    # TODO: Remove snap when updating this image
    ax.pcolormesh(x, y, Z, transform=ccrs.PlateCarree(), snap=False)


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_mercator_wrap'], tolerance=0.93)
def test_pcolormesh_mercator_wrap():
    x = np.linspace(0, 360, 73)
    y = np.linspace(-87.5, 87.5, 36)
    X, Y = np.meshgrid(*[np.deg2rad(c) for c in (x, y)])
    Z = np.cos(Y) + 0.375 * np.sin(2. * X)
    Z = Z[:-1, :-1]
    ax = plt.axes(projection=ccrs.Mercator())
    ax.coastlines()
    ax.pcolormesh(x, y, Z, transform=ccrs.PlateCarree(), snap=False)


@pytest.mark.natural_earth
@ImageTesting(['pcolormesh_mercator_wrap'], tolerance=0.93)
def test_pcolormesh_wrap_set_array():
    x = np.linspace(0, 360, 73)
    y = np.linspace(-87.5, 87.5, 36)
    X, Y = np.meshgrid(*[np.deg2rad(c) for c in (x, y)])
    Z = np.cos(Y) + 0.375 * np.sin(2. * X)
    Z = Z[:-1, :-1]
    ax = plt.axes(projection=ccrs.Mercator())
    norm = plt.Normalize(np.min(Z), np.max(Z))
    ax.coastlines()
    # Start off with bad data
    coll = ax.pcolormesh(x, y, np.ones(Z.shape), norm=norm,
                         transform=ccrs.PlateCarree(), snap=False)
    # Now update the plot with the set_array method
    coll.set_array(Z.ravel())


@pytest.mark.parametrize('shading, input_size, expected', [
    pytest.param('auto', 3, 4, id='auto same size'),
    pytest.param('auto', 4, 4, id='auto input larger'),
    pytest.param('nearest', 3, 4, id='nearest same size'),
    pytest.param('nearest', 4, 4, id='nearest input larger'),
    pytest.param('flat', 4, 4, id='flat input larger'),
    pytest.param('gouraud', 3, 3, id='gouraud same size')
])
def test_pcolormesh_shading(shading, input_size, expected):
    # Testing that the coordinates are all broadcast as expected with
    # the various shading options
    # The data shape is (3, 3) and we are changing the input shape
    # based upon that
    ax = plt.axes(projection=ccrs.PlateCarree())

    x = np.arange(input_size)
    y = np.arange(input_size)
    d = np.zeros((3, 3))

    coll = ax.pcolormesh(x, y, d, shading=shading)
    # We can use coll.get_coordinates() once MPL >= 3.5 is required
    # For now, we use the private variable for testing
    assert coll._coordinates.shape == (expected, expected, 2)


@pytest.mark.natural_earth
@ImageTesting(['quiver_plate_carree'])
def test_quiver_plate_carree():
    x = np.arange(-60, 42.5, 2.5)
    y = np.arange(30, 72.5, 2.5)
    x2d, y2d = np.meshgrid(x, y)
    u = np.cos(np.deg2rad(y2d))
    v = np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    plt.figure(figsize=(6, 6))
    # plot on native projection
    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution="110m")
    ax.quiver(x, y, u, v, mag)
    # plot on a different projection
    ax = plt.subplot(212, projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.quiver(x, y, u, v, mag, transform=ccrs.PlateCarree())


@pytest.mark.natural_earth
@ImageTesting(['quiver_rotated_pole'])
def test_quiver_rotated_pole():
    nx, ny = 22, 36
    x = np.linspace(311.91998291, 391.11999512, nx, endpoint=True)
    y = np.linspace(-23.59000015, 24.81000137, ny, endpoint=True)
    x2d, y2d = np.meshgrid(x, y)
    u = np.cos(np.deg2rad(y2d))
    v = -2. * np.cos(2. * np.deg2rad(y2d)) * np.sin(np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    rp = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
    plot_extent = [x[0], x[-1], y[0], y[-1]]
    # plot on native projection
    plt.figure(figsize=(6, 6))
    ax = plt.subplot(211, projection=rp)
    ax.set_extent(plot_extent, crs=rp)
    ax.coastlines()
    ax.quiver(x, y, u, v, mag)
    # plot on different projection
    ax = plt.subplot(212, projection=ccrs.PlateCarree())
    ax.set_extent(plot_extent, crs=rp)
    ax.coastlines()
    ax.quiver(x, y, u, v, mag, transform=rp)


@pytest.mark.natural_earth
@ImageTesting(['quiver_regrid'], tolerance=1.3)
def test_quiver_regrid():
    x = np.arange(-60, 42.5, 2.5)
    y = np.arange(30, 72.5, 2.5)
    x2d, y2d = np.meshgrid(x, y)
    u = np.cos(np.deg2rad(y2d))
    v = np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.quiver(x, y, u, v, mag, transform=ccrs.PlateCarree(),
              regrid_shape=30)


@pytest.mark.natural_earth
@ImageTesting(['quiver_regrid_with_extent'], tolerance=0.54)
def test_quiver_regrid_with_extent():
    x = np.arange(-60, 42.5, 2.5)
    y = np.arange(30, 72.5, 2.5)
    x2d, y2d = np.meshgrid(x, y)
    u = np.cos(np.deg2rad(y2d))
    v = np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    target_extent = [-3e6, 2e6, -6e6, -2.5e6]
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.quiver(x, y, u, v, mag, transform=ccrs.PlateCarree(),
              regrid_shape=10, target_extent=target_extent)


@pytest.mark.natural_earth
@ImageTesting(['barbs_plate_carree'])
def test_barbs():
    x = np.arange(-60, 45, 5)
    y = np.arange(30, 75, 5)
    x2d, y2d = np.meshgrid(x, y)
    u = 40 * np.cos(np.deg2rad(y2d))
    v = 40 * np.cos(2. * np.deg2rad(x2d))
    plot_extent = [-60, 40, 30, 70]
    plt.figure(figsize=(6, 6))
    # plot on native projection
    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution="110m")
    ax.barbs(x, y, u, v, length=4, linewidth=.25)
    # plot on a different projection
    ax = plt.subplot(212, projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution="110m")
    ax.barbs(x, y, u, v, transform=ccrs.PlateCarree(), length=4, linewidth=.25)


@pytest.mark.natural_earth
@ImageTesting(['barbs_regrid'])
def test_barbs_regrid():
    x = np.arange(-60, 42.5, 2.5)
    y = np.arange(30, 72.5, 2.5)
    x2d, y2d = np.meshgrid(x, y)
    u = 40 * np.cos(np.deg2rad(y2d))
    v = 40 * np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.barbs(x, y, u, v, mag, transform=ccrs.PlateCarree(),
             length=4, linewidth=.4, regrid_shape=20)


@pytest.mark.natural_earth
@ImageTesting(['barbs_regrid_with_extent'], tolerance=0.54)
def test_barbs_regrid_with_extent():
    x = np.arange(-60, 42.5, 2.5)
    y = np.arange(30, 72.5, 2.5)
    x2d, y2d = np.meshgrid(x, y)
    u = 40 * np.cos(np.deg2rad(y2d))
    v = 40 * np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    target_extent = [-3e6, 2e6, -6e6, -2.5e6]
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.barbs(x, y, u, v, mag, transform=ccrs.PlateCarree(),
             length=4, linewidth=.25, regrid_shape=10,
             target_extent=target_extent)


@pytest.mark.natural_earth
@ImageTesting(['barbs_1d'])
def test_barbs_1d():
    x = np.array([20., 30., -17., 15.])
    y = np.array([-1., 35., 11., 40.])
    u = np.array([23., -18., 2., -11.])
    v = np.array([5., -4., 19., 11.])
    plot_extent = [-21, 40, -5, 45]
    plt.figure(figsize=(6, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution="110m")
    ax.barbs(x, y, u, v, transform=ccrs.PlateCarree(),
             length=8, linewidth=1, color='#7f7f7f')


@pytest.mark.natural_earth
@ImageTesting(['barbs_1d_transformed'])
def test_barbs_1d_transformed():
    x = np.array([20., 30., -17., 15.])
    y = np.array([-1., 35., 11., 40.])
    u = np.array([23., -18., 2., -11.])
    v = np.array([5., -4., 19., 11.])
    plot_extent = [-20, 31, -5, 45]
    plt.figure(figsize=(6, 5))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.barbs(x, y, u, v, transform=ccrs.PlateCarree(),
             length=8, linewidth=1, color='#7f7f7f')


@pytest.mark.natural_earth
@ImageTesting(
    ['streamplot'], style='mpl20',
    tolerance=(42 if MPL_VERSION.release[:2] < (3, 2) else
               9.77 if MPL_VERSION.release[:2] < (3, 5) else 0.5))
def test_streamplot():
    x = np.arange(-60, 42.5, 2.5)
    y = np.arange(30, 72.5, 2.5)
    x2d, y2d = np.meshgrid(x, y)
    u = np.cos(np.deg2rad(y2d))
    v = np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.streamplot(x, y, u, v, transform=ccrs.PlateCarree(),
                  density=(1.5, 2), color=mag, linewidth=2*mag)
