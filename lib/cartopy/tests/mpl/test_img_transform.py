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

import operator
import os
import unittest

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from cartopy import config
from cartopy.tests.mpl import ImageTesting
import cartopy.crs as ccrs
import cartopy.img_transform as im_trans
from functools import reduce


class TestRegrid(unittest.TestCase):
    def test_array_dims(self):
        # Source data
        source_nx = 100
        source_ny = 100
        source_x = np.linspace(-180.0,
                               180.0,
                               source_nx).astype(np.float64)
        source_y = np.linspace(-90, 90.0, source_ny).astype(np.float64)
        source_x, source_y = np.meshgrid(source_x, source_y)
        data = np.arange(source_nx * source_ny,
                         dtype=np.int32).reshape(source_ny, source_nx)
        source_cs = ccrs.Geodetic()

        # Target grid
        target_nx = 23
        target_ny = 45
        target_proj = ccrs.PlateCarree()
        target_x, target_y, extent = im_trans.mesh_projection(target_proj,
                                                              target_nx,
                                                              target_ny)

        # Perform regrid
        new_array = im_trans.regrid(data, source_x, source_y, source_cs,
                                    target_proj, target_x, target_y)

        # Check dimensions of return array
        self.assertEqual(new_array.shape, target_x.shape)
        self.assertEqual(new_array.shape, target_y.shape)
        self.assertEqual(new_array.shape, (target_ny, target_nx))

    def test_different_dims(self):
        # Source data
        source_nx = 100
        source_ny = 100
        source_x = np.linspace(-180.0, 180.0,
                               source_nx).astype(np.float64)
        source_y = np.linspace(-90, 90.0,
                               source_ny).astype(np.float64)
        source_x, source_y = np.meshgrid(source_x, source_y)
        data = np.arange(source_nx * source_ny,
                         dtype=np.int32).reshape(source_ny, source_nx)
        source_cs = ccrs.Geodetic()

        # Target grids (different shapes)
        target_x_shape = (23, 45)
        target_y_shape = (23, 44)
        target_x = np.arange(reduce(operator.mul, target_x_shape),
                             dtype=np.float64).reshape(target_x_shape)
        target_y = np.arange(reduce(operator.mul, target_y_shape),
                             dtype=np.float64).reshape(target_y_shape)
        target_proj = ccrs.PlateCarree()

        # Attempt regrid
        with self.assertRaises(ValueError):
            im_trans.regrid(data, source_x, source_y, source_cs,
                            target_proj, target_x, target_y)


@ImageTesting(['regrid_image'], tolerance=0.4)
def test_regrid_image():
    # Source data
    fname = os.path.join(config["repo_data_dir"], 'raster', 'natural_earth',
                         '50-natural-earth-1-downsampled.png')
    nx = 720
    ny = 360
    source_proj = ccrs.PlateCarree()
    source_x, source_y, _ = im_trans.mesh_projection(source_proj, nx, ny)
    data = plt.imread(fname)
    # Flip vertically to match source_x/source_y orientation
    data = data[::-1]

    # Target grid
    target_nx = 300
    target_ny = 300
    target_proj = ccrs.InterruptedGoodeHomolosine()
    target_x, target_y, target_extent = im_trans.mesh_projection(target_proj,
                                                                 target_nx,
                                                                 target_ny)

    # Perform regrid
    new_array = im_trans.regrid(data, source_x, source_y, source_proj,
                                target_proj, target_x, target_y)

    # Plot
    fig = plt.figure(figsize=(10, 10))
    gs = matplotlib.gridspec.GridSpec(nrows=4, ncols=1,
                                      hspace=1.5, wspace=0.5)
    # Set up axes and title
    ax = plt.subplot(gs[0], frameon=False, projection=target_proj)
    plt.imshow(new_array, origin='lower', extent=target_extent)
    ax.coastlines()
    # Plot each color slice (tests masking)
    cmaps = {'red': 'Reds', 'green': 'Greens', 'blue': 'Blues'}
    for i, color in enumerate(['red', 'green', 'blue']):
        ax = plt.subplot(gs[i + 1], frameon=False, projection=target_proj)
        plt.imshow(new_array[:, :, i], extent=target_extent, origin='lower',
                   cmap=cmaps[color])
        ax.coastlines()

    # Tighten up layout
    gs.tight_layout(plt.gcf())


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
