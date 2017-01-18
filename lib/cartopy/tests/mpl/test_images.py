# (C) British Crown Copyright 2011 - 2017, Met Office
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

import os
import types

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import shapely.geometry as sgeom

from cartopy import config
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

from cartopy.tests.mpl import ImageTesting
import cartopy.tests.test_img_nest as ctest_nest
import cartopy.tests.test_img_tiles as ctest_tiles


NATURAL_EARTH_IMG = os.path.join(config["repo_data_dir"],
                                 'raster', 'natural_earth',
                                 '50-natural-earth-1-downsampled.png')
REGIONAL_IMG = os.path.join(config['repo_data_dir'], 'raster', 'sample',
                            'Miriam.A2012270.2050.2km.jpg')


# We have an exceptionally large tolerance for the web_tiles test.
# The basemap changes on a regular basis (for seasons) and we really only
# care that it is putting images onto the map which are roughly correct.
@ImageTesting(['web_tiles'], tolerance=2)
def test_web_tiles():
    extent = [-15, 0.1, 50, 60]
    target_domain = sgeom.Polygon([[extent[0], extent[1]],
                                   [extent[2], extent[1]],
                                   [extent[2], extent[3]],
                                   [extent[0], extent[3]],
                                   [extent[0], extent[1]]])
    map_prj = cimgt.GoogleTiles().crs

    ax = plt.subplot(2, 2, 1, projection=map_prj)
    gt = cimgt.GoogleTiles()
    gt._image_url = types.MethodType(ctest_tiles.GOOGLE_IMAGE_URL_REPLACEMENT,
                                     gt)
    img, extent, origin = gt.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=gt.crs,
              interpolation='bilinear', origin=origin)
    ax.coastlines(color='white')

    ax = plt.subplot(2, 2, 2, projection=map_prj)
    qt = cimgt.QuadtreeTiles()
    img, extent, origin = qt.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=qt.crs,
              interpolation='bilinear', origin=origin)
    ax.coastlines(color='white')

    ax = plt.subplot(2, 2, 3, projection=map_prj)
    osm = cimgt.OSM()
    img, extent, origin = osm.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=osm.crs,
              interpolation='bilinear', origin=origin)
    ax.coastlines()


@ImageTesting(['image_merge'])
def test_image_merge():
    # tests the basic image merging functionality
    tiles = []
    for i in range(1, 4):
        for j in range(0, 3):
            tiles.append((i, j, 2))

    gt = cimgt.GoogleTiles()
    gt._image_url = types.MethodType(ctest_tiles.GOOGLE_IMAGE_URL_REPLACEMENT,
                                     gt)
    images_to_merge = []
    for tile in tiles:
        img, extent, origin = gt.get_image(tile)
        img = np.array(img)
        x = np.linspace(extent[0], extent[1], img.shape[1], endpoint=False)
        y = np.linspace(extent[2], extent[3], img.shape[0], endpoint=False)
        images_to_merge.append([img, x, y, origin])

    img, extent, origin = cimgt._merge_tiles(images_to_merge)
    ax = plt.axes(projection=gt.crs)
    ax.set_global()
    ax.coastlines()
    plt.imshow(img, origin=origin, extent=extent, alpha=0.5)


@ImageTesting(['imshow_natural_earth_ortho'])
def test_imshow():
    source_proj = ccrs.PlateCarree()
    img = plt.imread(NATURAL_EARTH_IMG)
    # Convert the image to a byte array, rather than float, which is the
    # form that JPG images would be loaded with imread.
    img = (img * 255).astype('uint8')
    ax = plt.axes(projection=ccrs.Orthographic())
    ax.imshow(img, origin='upper', transform=source_proj,
              extent=[-180, 180, -90, 90])


@ImageTesting(['imshow_regional_projected'])
def test_imshow_projected():
    source_proj = ccrs.PlateCarree()
    img_extent = (-120.67660000000001, -106.32104523100001,
                  13.2301484511245, 30.766899999999502)
    img = plt.imread(REGIONAL_IMG)
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(img_extent, crs=source_proj)
    ax.coastlines(resolution='50m')
    ax.imshow(img, extent=img_extent, origin='upper', transform=source_proj)


@ImageTesting(['imshow_natural_earth_ortho'], tolerance=0.7)
def test_stock_img():
    ax = plt.axes(projection=ccrs.Orthographic())
    ax.stock_img()


@ImageTesting(['imshow_natural_earth_ortho'])
def test_pil_Image():
    img = Image.open(NATURAL_EARTH_IMG)
    source_proj = ccrs.PlateCarree()
    ax = plt.axes(projection=ccrs.Orthographic())
    ax.imshow(img, origin='upper', transform=source_proj,
              extent=[-180, 180, -90, 90])


@ImageTesting(['imshow_natural_earth_ortho'], tolerance=0.7)
def test_background_img():
    ax = plt.axes(projection=ccrs.Orthographic())
    ax.background_img(name='ne_shaded', resolution='low')


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
