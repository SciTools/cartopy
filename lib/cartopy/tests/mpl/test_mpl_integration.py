# (C) British Crown Copyright 2011 - 2012, Met Office
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

from nose.tools import assert_equal
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs

from cartopy.tests.mpl import ImageTesting


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


@ImageTesting(['global_pcolor_wrap'])
def test_global_pcolor_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.pcolor(x, y, data, transform=ccrs.PlateCarree())


@ImageTesting(['global_pcolor_wrap'])
def test_global_pcolor_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.pcolor(x, y, data)


@ImageTesting(['global_scatter_wrap'])
def test_global_scatter_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.scatter(x, y, c=data, transform=ccrs.PlateCarree())


@ImageTesting(['global_scatter_wrap'])
def test_global_scatter_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.scatter(x, y, c=data)


@ImageTesting(['global_map'])
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
                   ccrs.Mercator(),
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

        ax.set_title(prj.__class__.__name__)

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
    assert_equal(r.encode('ascii', 'ignore'), '-9.691e+05, -4.457e+06 (50.716617N, 12.267069W)')
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    x, y = np.array([-181.5]), np.array([50.])
    r = ax.format_coord(x, y)
    assert_equal(r.encode('ascii', 'ignore'), '-181.5, 50 (50.000000N, 178.500000E)')

    ax = plt.axes(projection=ccrs.Robinson())
    x, y = np.array([16060595.2]), np.array([2363093.4])
    r = ax.format_coord(x, y)
    assert_equal(r.encode('ascii', 'ignore'), '1.606e+07, 2.363e+06 (22.095524N, 173.709136E)')

    plt.close()


@ImageTesting(['natural_earth_interface'])
def test_axes_natural_earth_interface():
    rob = ccrs.Robinson()
    
    ax = plt.axes(projection=rob)
    
    ax.natural_earth_shp('rivers-lake-centerlines', edgecolor='black', facecolor='none')
    ax.natural_earth_shp('lakes', facecolor='blue')


if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
