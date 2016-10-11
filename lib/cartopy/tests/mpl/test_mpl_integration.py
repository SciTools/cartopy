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
import warnings

from nose.tools import assert_equal
try:
    from nose.tools import assert_regex
except ImportError:
    from nose.tools import assert_regexp_matches as assert_regex
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import six

import cartopy.crs as ccrs

from cartopy.tests.mpl import ImageTesting


_ROB_TOL = 0.5 if ccrs.PROJ4_VERSION < (4, 9) else 0.1


@ImageTesting(['global_contour_wrap'])
def test_global_contour_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contour(x, y, data, transform=ccrs.PlateCarree())


@ImageTesting(['global_contour_wrap'])
def test_global_contour_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contour(x, y, data)


@ImageTesting(['global_contourf_wrap'])
def test_global_contourf_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contourf(x, y, data, transform=ccrs.PlateCarree())


@ImageTesting(['global_contourf_wrap'])
def test_global_contourf_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contourf(x, y, data)


global_pcolor_wrap = ('global_pcolor_wrap'
                      if mpl.__version__ >= '1.5' else
                      'global_pcolor_wrap_pre_mpl_1.5')


@ImageTesting([global_pcolor_wrap])
def test_global_pcolor_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.pcolor(x, y, data, transform=ccrs.PlateCarree())


@ImageTesting([global_pcolor_wrap])
def test_global_pcolor_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.pcolor(x, y, data)


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


@ImageTesting(['global_scatter_wrap'])
def test_global_scatter_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(zorder=0)
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.scatter(x, y, c=data)


@ImageTesting(['global_map'],
              tolerance=16 if ccrs.PROJ4_VERSION < (4, 9) else 0.1)
def test_global_map():
    ax = plt.axes(projection=ccrs.Robinson())
#    ax.coastlines()
#    ax.gridlines(5)

    plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())

    plt.plot([-0.08, 132], [51.53, 43.17], color='red',
             transform=ccrs.PlateCarree())

    plt.plot([-0.08, 132], [51.53, 43.17], color='blue',
             transform=ccrs.Geodetic())


@ImageTesting(['simple_global'])
def test_simple_global():
    plt.axes(projection=ccrs.PlateCarree())
    plt.gca().coastlines()
    # produces a global map, despite not having needed to set the limits


@ImageTesting(['multiple_projections1'])
def test_multiple_projections():

    projections = [ccrs.PlateCarree(),
                   ccrs.Robinson(),
                   ccrs.RotatedPole(pole_latitude=45, pole_longitude=180),
                   ccrs.OSGB(),
                   ccrs.TransverseMercator(),
                   ccrs.Mercator(
                       globe=ccrs.Globe(semimajor_axis=math.degrees(1)),
                       min_latitude=-85., max_latitude=85.),
                   ccrs.LambertCylindrical(),
                   ccrs.Miller(),
                   ccrs.Gnomonic(),
                   ccrs.Stereographic(),
                   ccrs.NorthPolarStereo(),
                   ccrs.SouthPolarStereo(),
                   ccrs.Orthographic(),
                   ccrs.Mollweide(),
                   ccrs.InterruptedGoodeHomolosine(),
                   ]

    fig = plt.figure(figsize=(10, 10))
    for i, prj in enumerate(projections, 1):
        ax = fig.add_subplot(5, 5, i, projection=prj)

        ax.set_global()

        ax.coastlines()

        plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())

        plt.plot([-0.08, 132], [51.53, 43.17], color='red',
                 transform=ccrs.PlateCarree())

        plt.plot([-0.08, 132], [51.53, 43.17], color='blue',
                 transform=ccrs.Geodetic())


def test_cursor_values():
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    x, y = np.array([-969100.]), np.array([-4457000.])
    r = ax.format_coord(x, y)
    assert_equal(r.encode('ascii', 'ignore'),
                 six.b('-9.691e+05, -4.457e+06 (50.716617N, 12.267069W)'))

    ax = plt.axes(projection=ccrs.PlateCarree())
    x, y = np.array([-181.5]), np.array([50.])
    r = ax.format_coord(x, y)
    assert_equal(r.encode('ascii', 'ignore'),
                 six.b('-181.5, 50 (50.000000N, 178.500000E)'))

    ax = plt.axes(projection=ccrs.Robinson())
    x, y = np.array([16060595.2]), np.array([2363093.4])
    r = ax.format_coord(x, y)
    assert_regex(r.encode('ascii', 'ignore'),
                 six.b('1.606e\\+07, 2.363e\\+06 '
                       '\\(22.09[0-9]{4}N, 173.70[0-9]{4}E\\)'))

    plt.close()


@ImageTesting(['natural_earth_interface'], tolerance=_ROB_TOL)
def test_axes_natural_earth_interface():
    rob = ccrs.Robinson()

    ax = plt.axes(projection=rob)

    with warnings.catch_warnings(record=True) as all_warnings:
        warnings.simplefilter('always')

        ax.natural_earth_shp('rivers_lake_centerlines', edgecolor='black',
                             facecolor='none')
        ax.natural_earth_shp('lakes', facecolor='blue')

    assert_equal(len(all_warnings), 2)
    for warning in all_warnings:
        msg = str(warning.message)
        assert 'deprecated' in msg
        assert 'add_feature' in msg


@ImageTesting(['pcolormesh_global_wrap1'])
def test_pcolormesh_global_with_wrap1():
    # make up some realistic data with bounds (such as data from the UM)
    nx, ny = 36, 18
    xbnds = np.linspace(0, 360, nx, endpoint=True)
    ybnds = np.linspace(-90, 90, ny, endpoint=True)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(212, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


@ImageTesting(['pcolormesh_global_wrap2'])
def test_pcolormesh_global_with_wrap2():
    # make up some realistic data with bounds (such as data from the UM)
    nx, ny = 36, 18
    xbnds, xstep = np.linspace(0, 360, nx - 1, retstep=True, endpoint=True)
    ybnds, ystep = np.linspace(-90, 90, nx - 1, retstep=True, endpoint=True)
    xbnds -= xstep / 2
    ybnds -= ystep / 2
    xbnds = np.append(xbnds, xbnds[-1] + xstep)
    ybnds = np.append(ybnds, ybnds[-1] + ystep)

    x, y = np.meshgrid(xbnds, ybnds)
    data = np.exp(np.sin(np.deg2rad(x)) + np.cos(np.deg2rad(y)))
    data = data[:-1, :-1]

    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(212, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


@ImageTesting(['pcolormesh_global_wrap3'], tolerance=_ROB_TOL)
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
    c = plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    assert c._wrapped_collection_fix is not None, \
        'No pcolormesh wrapping was done when it should have been.'

    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(312, projection=ccrs.PlateCarree(-1.87499952))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible

    ax = plt.subplot(313, projection=ccrs.Robinson(-2))
    plt.pcolormesh(xbnds, ybnds, data, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()  # make sure everything is visible


@ImageTesting(['pcolormesh_limited_area_wrap'])
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
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Set1')
    ax.coastlines()

    ax = plt.subplot(222, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Set1')
    ax.coastlines()
    ax.set_global()

    # draw the same plot, only more zoomed in, and using the 2d versions
    # of the coordinates (just to test that 1d and 2d are both suitably
    # being fixed)
    ax = plt.subplot(223, projection=ccrs.PlateCarree(180))
    plt.pcolormesh(x, y, data, transform=rp, cmap='Set1')
    ax.coastlines()
    ax.set_extent([-70, 0, 0, 80])

    ax = plt.subplot(224, projection=rp)
    plt.pcolormesh(xbnds, ybnds, data, transform=rp, cmap='Set1')
    ax.coastlines()


@ImageTesting(['pcolormesh_goode_wrap'])
def test_pcolormesh_goode_wrap():
    # global data on an Interrupted Goode Homolosine projection
    # shouldn't spill outside projection boundary
    x = np.linspace(0, 360, 73)
    y = np.linspace(-87.5, 87.5, 36)
    X, Y = np.meshgrid(*[np.deg2rad(c) for c in (x, y)])
    Z = np.cos(Y) + 0.375 * np.sin(2. * X)
    Z = Z[:-1, :-1]
    ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
    ax.coastlines()
    ax.pcolormesh(x, y, Z, transform=ccrs.PlateCarree())


@ImageTesting(['pcolormesh_mercator_wrap'])
def test_pcolormesh_mercator_wrap():
    x = np.linspace(0, 360, 73)
    y = np.linspace(-87.5, 87.5, 36)
    X, Y = np.meshgrid(*[np.deg2rad(c) for c in (x, y)])
    Z = np.cos(Y) + 0.375 * np.sin(2. * X)
    Z = Z[:-1, :-1]
    ax = plt.axes(projection=ccrs.Mercator())
    ax.coastlines()
    ax.pcolormesh(x, y, Z, transform=ccrs.PlateCarree())


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
    ax.coastlines()
    ax.quiver(x, y, u, v, mag)
    # plot on a different projection
    ax = plt.subplot(212, projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.quiver(x, y, u, v, mag, transform=ccrs.PlateCarree())


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


@ImageTesting(['quiver_regrid'])
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


@ImageTesting(['quiver_regrid_with_extent'])
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


@ImageTesting(['barbs_plate_carree'])
def test_barbs():
    x = np.arange(-60, 45, 5)
    y = np.arange(30, 75, 5)
    x2d, y2d = np.meshgrid(x, y)
    u = 40 * np.cos(np.deg2rad(y2d))
    v = 40 * np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    plot_extent = [-60, 40, 30, 70]
    plt.figure(figsize=(6, 6))
    # plot on native projection
    ax = plt.subplot(211, projection=ccrs.PlateCarree())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.barbs(x, y, u, v, length=4, linewidth=.25)
    # plot on a different projection
    ax = plt.subplot(212, projection=ccrs.NorthPolarStereo())
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.barbs(x, y, u, v, transform=ccrs.PlateCarree(), length=4, linewidth=.25)


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


@ImageTesting(['barbs_regrid_with_extent'])
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
    ax.coastlines()
    ax.barbs(x, y, u, v, transform=ccrs.PlateCarree(),
             length=8, linewidth=1, color='#7f7f7f')


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


@ImageTesting(['streamplot'
               if mpl.__version__ >= '1.5' else
               'streamplot_pre_mpl_1.5'])
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


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
