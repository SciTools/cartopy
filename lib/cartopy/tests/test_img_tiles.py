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
from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt
import shapely.geometry

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt


def test_google_wts():
    gt = cimgt.GoogleTiles()

    extent = [-15, 00, 50, 60]
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]],
                                              [extent[2], extent[1]],
                                              [extent[2], extent[3]],
                                              [extent[0], extent[3]],
                                              [extent[0], extent[1]]])

    with assert_raises(AssertionError):
        list(gt.find_images(target_domain, -1))
    assert_equal(tuple(gt.find_images(target_domain, 0)), ((0, 0, 0),))
    assert_equal(tuple(gt.find_images(target_domain, 2)), ((1, 1, 2), (2, 1, 2)))

    assert_equal(list(gt.subtiles((0, 0, 0))), [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)])
    assert_equal(list(gt.subtiles((1, 0, 1))), [(2, 0, 2), (2, 1, 2), (3, 0, 2), (3, 1, 2)])

    with assert_raises(AssertionError):
        gt.tileextent((0, 1, 0))

    assert_equal(gt.tileextent((0, 0, 0)), (-180.0, 180.0, 179.41035067677481, -179.41035067677487))
    assert_equal(gt.tileextent((2, 0, 2)), (0.0, 90.0, 179.41035067677481, 89.705175338387392))
    assert_equal(gt.tileextent((0, 2, 2)), (-180.0, -90.0, -2.8421709430404007e-14, -89.70517533838742))
    assert_equal(gt.tileextent((2, 2, 2)), (0.0, 90.0, -2.8421709430404007e-14, -89.70517533838742))
    assert_equal(gt.tileextent((8, 9, 4)), (0.0, 22.5, -22.426293834596891, -44.852587669193753))


def test_quadtree_wts():
    qt = cimgt.QuadtreeTiles()

    extent = [-15, 00, 50, 60]
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]],
                                              [extent[2], extent[1]],
                                              [extent[2], extent[3]],
                                              [extent[0], extent[3]],
                                              [extent[0], extent[1]]])

    with assert_raises(ValueError):
        list(qt.find_images(target_domain, 0))

    assert_equal(qt.tms_to_quadkey((1, 1, 1)), '1')
    assert_equal(qt.quadkey_to_tms('1'), (1, 1, 1))

    assert_equal(qt.tms_to_quadkey((8, 9, 4)), '1220')
    assert_equal(qt.quadkey_to_tms('1220'), (8, 9, 4))

    assert_equal(tuple(qt.find_images(target_domain, 1)), ('0', '1'))
    assert_equal(tuple(qt.find_images(target_domain, 2)), ('03', '12'))

    assert_equal(list(qt.subtiles('0')), ['00', '01', '02', '03'])
    assert_equal(list(qt.subtiles('11')), ['110', '111', '112', '113'])

    with assert_raises(ValueError):
        qt.tileextent('4')

    assert_equal(qt.tileextent(''), (-180.0, 180.0, 179.4103506767748, -179.41035067677487))
    assert_equal(qt.tileextent(qt.tms_to_quadkey((2, 0, 2), google=True)),
                 (0.0, 90.0, 179.41035067677481, 89.705175338387392))
    assert_equal(qt.tileextent(qt.tms_to_quadkey((0, 2, 2), google=True)),
                 (-180.0, -90.0, -2.8421709430404007e-14, -89.70517533838742))
    assert_equal(qt.tileextent(qt.tms_to_quadkey((2, 2, 2), google=True)),
                 (0.0, 90.0, -2.8421709430404007e-14, -89.70517533838742))
    assert_equal(qt.tileextent(qt.tms_to_quadkey((8, 9, 4), google=True)),
                 (0.0, 22.5, -22.426293834596891, -44.852587669193753))


def image_comparison(baseline_images=None, extensions=('png',), tol=1e-3):
    # changes the mpl default to only use PNGs
    return mpl_image_comparison(baseline_images, extensions, tol)


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


    ax.coastlines()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
