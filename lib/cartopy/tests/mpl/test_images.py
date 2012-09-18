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

from nose.tools import assert_equal, assert_raises
import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

from cartopy.tests.mpl import image_comparison
import cartopy.tests.test_img_nest as ctest_nest


@image_comparison(baseline_images=['web_tiles'])
def test_web_tiles():
    extent = [-15, 00, 50, 60]
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]],
                                              [extent[2], extent[1]],
                                              [extent[2], extent[3]],
                                              [extent[0], extent[3]],
                                              [extent[0], extent[1]]])

    ax = plt.subplot(3, 2, 1, projection=ccrs.Mercator())
    gt = cimgt.GoogleTiles()
    img, extent, origin = gt.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=ccrs.Mercator(),
              interpolation='bilinear', origin=origin)
    ax.coastlines(color='white')


    ax = plt.subplot(3, 2, 2, projection=ccrs.Mercator())
    qt = cimgt.QuadtreeTiles()
    img, extent, origin = qt.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=ccrs.Mercator(),
              interpolation='bilinear', origin=origin)
    ax.coastlines(color='white')


    ax = plt.subplot(3, 2, 3, projection=ccrs.Mercator())
    mq_osm = cimgt.MapQuestOSM()
    img, extent, origin = mq_osm.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=ccrs.Mercator(),
              interpolation='bilinear', origin=origin)
    ax.coastlines()


    ax = plt.subplot(3, 2, 4, projection=ccrs.Mercator())
    mq_oa = cimgt.MapQuestOpenAerial()
    img, extent, origin = mq_oa.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=ccrs.Mercator(),
              interpolation='bilinear', origin=origin)
    ax.coastlines()


    ax = plt.subplot(3, 2, 5, projection=ccrs.Mercator())
    osm = cimgt.OSM()
    img, extent, origin = osm.image_for_domain(target_domain, 1)
    ax.imshow(np.array(img), extent=extent, transform=ccrs.Mercator(),
              interpolation='bilinear', origin=origin)
    ax.coastlines()



@image_comparison(baseline_images=['image_nest'])
def test_image_nest():
    nest_z0_z1 = ctest_nest.gen_nest()

    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_xlim(-45, 45)
    ax.set_ylim(-45, 90)
    ax.coastlines()
    ax.add_image(nest_z0_z1, 'aerial z1 test')


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
