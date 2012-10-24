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

import operator
import unittest

import matplotlib
import matplotlib.pyplot as plt
import numpy

from cartopy.tests.mpl import ImageTesting
import cartopy.crs as ccrs
import cartopy.img_transform


class TestRegrid(unittest.TestCase):
    def test_array_dims(self):
        # Source data
        source_nx = 100
        source_ny = 100
        source_x = numpy.linspace(-180.0, 180.0, source_nx).astype(numpy.float64)
        source_y = numpy.linspace(-90, 90.0, source_ny).astype(numpy.float64)
        source_x, source_y = numpy.meshgrid(source_x, source_y)
        data = numpy.arange(source_nx * source_ny,
                            dtype=numpy.int32).reshape(source_ny, source_nx)
        source_cs = ccrs.Geodetic()

        # Target grid
        target_nx = 23
        target_ny = 45
        target_proj = ccrs.PlateCarree()
        target_x, target_y, extent = cartopy.img_transform.mesh_projection(target_proj,
                                                                           target_nx,
                                                                           target_ny)

        
        # Perform regrid
        new_array = cartopy.img_transform.regrid(data, source_x, source_y, source_cs,
                           target_proj, target_x, target_y)
        
        # Check dimensions of return array
        self.assertEqual(new_array.shape, target_x.shape)
        self.assertEqual(new_array.shape, target_y.shape)
        self.assertEqual(new_array.shape, (target_ny, target_nx))

    def test_different_dims(self):
        # Source data
        source_nx = 100
        source_ny = 100
        source_x = numpy.linspace(-180.0, 180.0, source_nx).astype(numpy.float64)
        source_y = numpy.linspace(-90, 90.0, source_ny).astype(numpy.float64)
        source_x, source_y = numpy.meshgrid(source_x, source_y)
        data = numpy.arange(source_nx * source_ny,
                            dtype=numpy.int32).reshape(source_ny, source_nx)
        source_cs = ccrs.Geodetic()

        # Target grids (different shapes)
        target_x_shape = (23, 45)
        target_y_shape = (23, 44)
        target_x = numpy.arange(reduce(operator.mul, target_x_shape)).reshape(target_x_shape).astype(numpy.float64)
        target_y = numpy.arange(reduce(operator.mul, target_y_shape)).reshape(target_y_shape).astype(numpy.float64)
        target_proj = ccrs.PlateCarree()
        
        # Attempt regrid
        with self.assertRaises(ValueError):
            new_array = cartopy.img_transform.regrid(data, source_x, source_y, source_cs,
                               target_proj, target_x, target_y)


@ImageTesting(['regrid_blue_marble'])
def test_regrid_blue_marble_img():
    # Source data
    filename = '/data/local/dataZoo/cartography/raster/blue_marble_720_360.png'
    nx = 720
    ny = 360
    source_proj = ccrs.PlateCarree()
    source_x, source_y, source_extent = cartopy.img_transform.mesh_projection(source_proj, nx, ny)
    data = plt.imread(filename)
    # Flip vertically to match source_x/source_y orientation
    data = data[::-1]

    # Target grid
    target_nx = 300
    target_ny = 300
    target_proj = ccrs.InterruptedGoodeHomolosine()
    target_x, target_y, target_extent = cartopy.img_transform.mesh_projection(target_proj,
                                                                              target_nx,
                                                                              target_ny)

    # Perform regrid
    new_array = cartopy.img_transform.regrid(data, source_x, source_y, source_proj,
                                             target_proj, target_x, target_y)

    # Plot
    fig = plt.figure(figsize=(10, 10))
    gs = matplotlib.gridspec.GridSpec(nrows=4, ncols=1, hspace=1.5, wspace=0.5)
    # Set up axes and title
    ax = plt.subplot(gs[0], frameon=False, projection=target_proj)
    plt.imshow(new_array, origin='lower', extent=target_extent)
    ax.coastlines()
    # Plot each colour slice (tests masking)
    cmaps = {'red': 'Reds', 'green': 'Greens', 'blue': 'Blues'}
    for i, colour in enumerate(['red', 'green', 'blue']):
        ax = plt.subplot(gs[i + 1], frameon=False, projection=target_proj)
        ax.set_title(colour)
        plt.imshow(new_array[:, :, i], extent=target_extent, origin='lower',
                       cmap=cmaps[colour])
        ax.coastlines()

    # Tighten up layout
    gs.tight_layout(plt.gcf())
        

if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
