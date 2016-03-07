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

import types

from nose.tools import assert_equal, assert_raises
import numpy as np
from numpy.testing import assert_array_almost_equal as assert_arr_almost
import shapely.geometry as sgeom

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt


#: Maps Google tile coordinates to native mercator coordinates as defined
#: by https://goo.gl/pgJi.
KNOWN_EXTENTS = {(0, 0, 0): (-20037508.342789244, 20037508.342789244,
                             -20037508.342789244, 20037508.342789244),
                 (2, 0, 2): (0., 10018754.17139462,
                             10018754.17139462, 20037508.342789244),
                 (0, 2, 2): (-20037508.342789244, -10018754.171394622,
                             -10018754.171394622, 0),
                 (2, 2, 2): (0, 10018754.17139462,
                             -10018754.171394622, 0),
                 (8, 9, 4): (0, 2504688.542848654,
                             -5009377.085697312, -2504688.542848654),
                 }


def GOOGLE_IMAGE_URL_REPLACEMENT(self, tile):
    url = ('https://chart.googleapis.com/chart?chst=d_text_outline&'
           'chs=256x256&chf=bg,s,00000055&chld=FFFFFF|16|h|000000|b||||'
           'Google:%20%20(' + str(tile[0]) + ',' + str(tile[1]) + ')'
           '|Zoom%20' + str(tile[2]) + '||||||______________________'
           '______')
    return url


def test_google_tile_styles():
    """
    Tests that setting the Google Maps tile style works as expected.

    This is essentially just assures information is properly propagated through
    the class structure.
    """
    reference_url = ("https://mts0.google.com/vt/lyrs={style}@177000000&hl=en"
                     "&src=api&x=1&y=2&z=3&s=G")
    tile = ["1", "2", "3"]

    # Default is street.
    gt = cimgt.GoogleTiles()
    url = gt._image_url(tile)
    assert_equal(reference_url.format(style="m"), url)

    # Street
    gt = cimgt.GoogleTiles(style="street")
    url = gt._image_url(tile)
    assert_equal(reference_url.format(style="m"), url)

    # Satellite
    gt = cimgt.GoogleTiles(style="satellite")
    url = gt._image_url(tile)
    assert_equal(reference_url.format(style="s"), url)

    # Terrain
    gt = cimgt.GoogleTiles(style="terrain")
    url = gt._image_url(tile)
    assert_equal(reference_url.format(style="t"), url)

    # Streets only
    gt = cimgt.GoogleTiles(style="only_streets")
    url = gt._image_url(tile)
    assert_equal(reference_url.format(style="h"), url)

    # Exception is raised if unknown style is passed.
    with assert_raises(ValueError):
        cimgt.GoogleTiles(style="random_style")


def test_google_wts():
    gt = cimgt.GoogleTiles()

    ll_target_domain = sgeom.box(-15, 50, 0, 60)
    multi_poly = gt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

    with assert_raises(AssertionError):
        list(gt.find_images(target_domain, -1))
    assert_equal(tuple(gt.find_images(target_domain, 0)),
                 ((0, 0, 0),))
    assert_equal(tuple(gt.find_images(target_domain, 2)),
                 ((1, 1, 2), (2, 1, 2)))

    assert_equal(list(gt.subtiles((0, 0, 0))),
                 [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)])
    assert_equal(list(gt.subtiles((1, 0, 1))),
                 [(2, 0, 2), (2, 1, 2), (3, 0, 2), (3, 1, 2)])

    with assert_raises(AssertionError):
        gt.tileextent((0, 1, 0))

    assert_arr_almost(gt.tileextent((0, 0, 0)), KNOWN_EXTENTS[(0, 0, 0)])
    assert_arr_almost(gt.tileextent((2, 0, 2)), KNOWN_EXTENTS[(2, 0, 2)])
    assert_arr_almost(gt.tileextent((0, 2, 2)), KNOWN_EXTENTS[(0, 2, 2)])
    assert_arr_almost(gt.tileextent((2, 2, 2)), KNOWN_EXTENTS[(2, 2, 2)])
    assert_arr_almost(gt.tileextent((8, 9, 4)), KNOWN_EXTENTS[(8, 9, 4)])


def test_tile_bbox_y0_at_south_pole():
    tms = cimgt.MapQuestOpenAerial()

    # Check the y0_at_north_pole keywords returns the appropriate bounds.
    assert_arr_almost(tms.tile_bbox(8, 6, 4, y0_at_north_pole=False),
                      np.array(KNOWN_EXTENTS[(8, 9, 4)]).reshape([2, 2]))


def test_tile_find_images():
    gt = cimgt.GoogleTiles()
    # Test the find_images method on a GoogleTiles instance.
    ll_target_domain = sgeom.box(-10, 50, 10, 60)
    multi_poly = gt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

    assert_equal([(7, 4, 4), (7, 5, 4), (8, 4, 4), (8, 5, 4)],
                 list(gt.find_images(target_domain, 4)))


def test_image_for_domain():
    gt = cimgt.GoogleTiles()
    gt._image_url = types.MethodType(GOOGLE_IMAGE_URL_REPLACEMENT, gt)

    ll_target_domain = sgeom.box(-10, 50, 10, 60)
    multi_poly = gt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

    _, extent, _ = gt.image_for_domain(target_domain, 6)

    ll_extent = ccrs.Geodetic().transform_points(gt.crs,
                                                 np.array(extent[:2]),
                                                 np.array(extent[2:]))
    assert_arr_almost(ll_extent[:, :2],
                      [[-11.25, 48.92249926],
                       [11.25, 61.60639637]])


def test_quadtree_wts():
    qt = cimgt.QuadtreeTiles()

    ll_target_domain = sgeom.box(-15, 50, 0, 60)
    multi_poly = qt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

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

    assert_arr_almost(qt.tileextent(''), KNOWN_EXTENTS[(0, 0, 0)])
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((2, 0, 2), google=True)),
                      KNOWN_EXTENTS[(2, 0, 2)])
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((0, 2, 2), google=True)),
                      KNOWN_EXTENTS[(0, 2, 2)])
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((2, 0, 2), google=True)),
                      KNOWN_EXTENTS[(2, 0, 2)])
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((2, 2, 2), google=True)),
                      KNOWN_EXTENTS[(2, 2, 2)])
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((8, 9, 4), google=True)),
                      KNOWN_EXTENTS[(8, 9, 4)])


def test_mapbox_tiles():
    token = 'foo'
    map_id = 'bar'
    tile = [0, 1, 2]
    exp_url = 'https://api.tiles.mapbox.com/v4/bar/2/0/1.png?access_token=foo'

    mapbox_sample = cimgt.MapboxTiles(token, map_id)
    url_str = mapbox_sample._image_url(tile)
    assert_equal(url_str, exp_url)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
