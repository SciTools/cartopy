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

import numpy as np
from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt

import cartopy.crs as ccrs

from cartopy.tests.mpl import image_comparison


@image_comparison(baseline_images=['global_contour_wrap'])
def test_global_contour_wrap_new_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contourf(x, y, data, transform=ccrs.PlateCarree())


@image_comparison(baseline_images=['global_contour_wrap'])
def test_global_contour_wrap_no_transform():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    x, y = np.meshgrid(np.linspace(0, 360), np.linspace(-90, 90))
    data = np.sin(np.sqrt(x ** 2 + y ** 2))
    plt.contourf(x, y, data)


@image_comparison(baseline_images=['global_map'])
def test_global_map():
    ax = plt.axes(projection=ccrs.Robinson())
#    ax.coastlines()
#    ax.gridlines(5)

    plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())

    plt.plot([-0.08, 132], [51.53, 43.17], color='red',
             transform=ccrs.PlateCarree())

    plt.plot([-0.08, 132], [51.53, 43.17], color='blue',
             transform=ccrs.Geodetic())


@image_comparison(baseline_images=['multiple_projections1'])
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

#
#@image_comparison(baseline_images=['image_transform'])
#def test_image_transforms():
#    plt.subplot(131, projection=ccrs.PlateCarree())
#

if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
