<<<<<<< HEAD
# Copyright Cartopy Contributors
=======
# (C) British Crown Copyright 2011 - 2020, Met Office
>>>>>>> Fix licence headers
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import hashlib
import json
import os
import types

import numpy as np
from numpy.testing import assert_array_almost_equal as assert_arr_almost
from PIL import Image
import pytest
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
if ccrs.PROJ4_VERSION == (5, 0, 0):
    KNOWN_EXTENTS = {
        (0, 0, 0): (-20037508.342789244, 20037508.342789244,
                    -19994827.892149, 19994827.892149),
        (2, 0, 2): (0, 10018754.171395,
                    9997413.946075, 19994827.892149),
        (0, 2, 2): (-20037508.342789244, -10018754.171394622,
                    -9997413.946075, 0),
        (2, 2, 2): (0, 10018754.171395,
                    -9997413.946075, 0),
        (8, 9, 4): (0, 2504688.542849,
                    -4998706.973037, -2499353.486519),
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
    assert reference_url.format(style="m") == url

    # Street
    gt = cimgt.GoogleTiles(style="street")
    url = gt._image_url(tile)
    assert reference_url.format(style="m") == url

    # Satellite
    gt = cimgt.GoogleTiles(style="satellite")
    url = gt._image_url(tile)
    assert reference_url.format(style="s") == url

    # Terrain
    gt = cimgt.GoogleTiles(style="terrain")
    url = gt._image_url(tile)
    assert reference_url.format(style="t") == url

    # Streets only
    gt = cimgt.GoogleTiles(style="only_streets")
    url = gt._image_url(tile)
    assert reference_url.format(style="h") == url

    # Exception is raised if unknown style is passed.
    with pytest.raises(ValueError):
        cimgt.GoogleTiles(style="random_style")


def test_google_wts():
    gt = cimgt.GoogleTiles()

    ll_target_domain = sgeom.box(-15, 50, 0, 60)
    multi_poly = gt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

    with pytest.raises(AssertionError):
        list(gt.find_images(target_domain, -1))
    assert (tuple(gt.find_images(target_domain, 0)) ==
                 ((0, 0, 0),))
    assert (tuple(gt.find_images(target_domain, 2)) ==
                 ((1, 1, 2), (2, 1, 2)))

    assert (list(gt.subtiles((0, 0, 0))) ==
            [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)])
    assert (list(gt.subtiles((1, 0, 1))) ==
            [(2, 0, 2), (2, 1, 2), (3, 0, 2), (3, 1, 2)])

    with pytest.raises(AssertionError):
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

    assert (list(gt.find_images(target_domain, 4)) ==
            [(7, 4, 4), (7, 5, 4), (8, 4, 4), (8, 5, 4)])


@pytest.mark.network
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
    if ccrs.PROJ4_VERSION == (5, 0, 0):
        assert_arr_almost(ll_extent[:, :2],
                          [[-11.25, 49.033955],
                           [11.25, 61.687101]])
    else:
        assert_arr_almost(ll_extent[:, :2],
                          [[-11.25, 48.92249926],
                           [11.25, 61.60639637]])


def test_quadtree_wts():
    qt = cimgt.QuadtreeTiles()

    ll_target_domain = sgeom.box(-15, 50, 0, 60)
    multi_poly = qt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

    with pytest.raises(ValueError):
        list(qt.find_images(target_domain, 0))

    assert qt.tms_to_quadkey((1, 1, 1)) == '1'
    assert qt.quadkey_to_tms('1') == (1, 1, 1)

    assert qt.tms_to_quadkey((8, 9, 4)) == '1220'
    assert qt.quadkey_to_tms('1220') == (8, 9, 4)

    assert tuple(qt.find_images(target_domain, 1)) == ('0', '1')
    assert tuple(qt.find_images(target_domain, 2)) == ('03', '12')

    assert list(qt.subtiles('0')) == ['00', '01', '02', '03']
    assert list(qt.subtiles('11')) == ['110', '111', '112', '113']

    with pytest.raises(ValueError):
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


def test_mapbox_tiles_api_url():
    token = 'foo'
    map_name = 'bar'
    tile = [0, 1, 2]
    exp_url = ('https://api.mapbox.com/v4/mapbox.bar'
               '/2/0/1.png?access_token=foo')

    mapbox_sample = cimgt.MapboxTiles(token, map_name)
    url_str = mapbox_sample._image_url(tile)
    assert url_str == exp_url


def test_mapbox_style_tiles_api_url():
    token = 'foo'
    username = 'baz'
    map_id = 'bar'
    tile = [0, 1, 2]
    exp_url = ('https://api.mapbox.com/styles/v1/'
               'baz/bar/tiles/256/2/0/1'
               '?access_token=foo')

    mapbox_sample = cimgt.MapboxStyleTiles(token, username, map_id)
    url_str = mapbox_sample._image_url(tile)
    assert url_str == exp_url


def test_ordnance_survey_tile_styles():
    """
    Tests that setting the Ordnance Survey tile style works as expected.

    This is essentially just assures information is properly propagated through
    the class structure.
    """
    dummy_apikey = "None"

    ref_url = ('https://api2.ordnancesurvey.co.uk/'
               'mapping_api/v1/service/wmts?'
               'key=None&height=256&width=256&tilematrixSet=EPSG%3A3857&'
               'version=1.0.0&style=true&layer={layer}%203857&'
               'SERVICE=WMTS&REQUEST=GetTile&format=image%2Fpng&'
               'TileMatrix=EPSG%3A3857%3A{z}&TileRow={y}&TileCol={x}')
    tile = ["1", "2", "3"]

    # Default is Road.
    os = cimgt.OrdnanceSurvey(dummy_apikey)
    url = os._image_url(tile)
    assert url == ref_url.format(layer="Road",
                                 z=tile[2], y=tile[1], x=tile[0])

    for layer in ['Outdoor', 'Light', 'Night', 'Leisure']:
        os = cimgt.OrdnanceSurvey(dummy_apikey, layer=layer)
        url = os._image_url(tile)
        assert url == ref_url.format(layer=layer,
                                     z=tile[2], y=tile[1], x=tile[0])

    # Exception is raised if unknown style is passed.
    with pytest.raises(ValueError):
        cimgt.OrdnanceSurvey(dummy_apikey, layer="random_style")


@pytest.mark.network
def test_ordnance_survey_get_image():
    # In order to test fetching map images from OS
    # an API key needs to be provided
    try:
        api_key = os.environ['ORDNANCE_SURVEY_API_KEY']
    except KeyError:
        pytest.skip('ORDNANCE_SURVEY_API_KEY environment variable is unset.')

    os1 = cimgt.OrdnanceSurvey(api_key, layer="Outdoor")
    os2 = cimgt.OrdnanceSurvey(api_key, layer="Night")

    tile = (500, 300, 10)

    img1, extent1, _ = os1.get_image(tile)
    img2, extent2, _ = os2.get_image(tile)

    # Different images for different layers
    assert img1 != img2

    # The extent is the same though
    assert extent1 == extent2


@pytest.mark.network
def test_cache(tmpdir):
    tmpdir_str = tmpdir.strpath

    # Fetch tiles and save them in the cache
    gt = cimgt.GoogleTiles(cache_path=tmpdir_str)
    gt._image_url = types.MethodType(GOOGLE_IMAGE_URL_REPLACEMENT, gt)

    ll_target_domain = sgeom.box(-10, 50, 10, 60)
    multi_poly = gt.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
    target_domain = multi_poly.geoms[0]

    img_init, _, _ = gt.image_for_domain(target_domain, 6)

    # Define expected results
    x_y_z_h = [
        (30, 18, 6, '3f4aad7c536b7d43a16b891283c6ec02'),
        (30, 19, 6, '308d0a772f08e4e0ee6285ca8595913d'),
        (30, 20, 6, 'bb613f149ca6cec329591bd3f50e5a9f'),
        (30, 21, 6, '4f1bb78eb1fe225b77dde6ee1cc33afe'),
        (31, 18, 6, 'e6e2343ff5c6a3f3632a57e673b4a155'),
        (31, 19, 6, '943ec0bc1f55b3b4e7b5c6a17313a6c8'),
        (31, 20, 6, 'e835fd11f30d9976f8aa722205d47af8'),
        (31, 21, 6, 'f9b54ce321ce2a351480f440f70a31eb'),
        (32, 18, 6, 'dd065d5502411e5752dac7bd596cca02'),
        (32, 19, 6, '641b162bc052785d6e1270a49aa42a10'),
        (32, 20, 6, '1de7433ab6f6ce272f269a46c214ff08'),
        (32, 21, 6, 'faa04e22f983e842496fdee07bd67435'),
        (33, 18, 6, 'c097aedf3d78bbe898950244e5e1b82c'),
        (33, 19, 6, '606473c53c86eeda09b06f729f04246d'),
        (33, 20, 6, '541bcd6b12609c2423f0bf141f7f6705'),
        (33, 21, 6, '8ea3e1caff873c7947d2ef4f7f353ae4')
    ]
    base_url = (
        'https://chart.googleapis.com/chart?chst=d_text_outline&chs=256x256&'
        'chf=bg,s,00000055&chld=FFFFFF|16|h|000000|b||||Google:%20%20({},{})'
        '|Zoom%20{}||||||____________________________'
    )

    # Check the results
    files = [i for i in os.listdir(tmpdir_str) if i != "files"]
    hashes = {
        f:
        hashlib.md5(
            open(os.path.join(tmpdir_str, f), mode="rb").read()
        ).hexdigest()
        for f in files
    }
    url_ids = json.load(open(os.path.join(tmpdir_str, "files")))

    assert sorted(hashes.values()) == sorted([
        h for x, y, z, h in x_y_z_h
    ])

    assert sorted(url_ids.keys()) == sorted([
        base_url.format(x, y, z)
        for x, y, z, h in x_y_z_h
    ])

    url_hashes = {url: hashes[uid + ".tiff"] for url, uid in url_ids.items()}
    assert url_hashes == {
        base_url.format(x, y, z): h
        for x, y, z, h in x_y_z_h
    }

    # Update images in cache (all white)
    for f in files:
        filename = os.path.join(tmpdir_str, f)
        img = Image.new('RGB', (255, 255), "white")
        img.save(filename)

    gt_cache = cimgt.GoogleTiles(cache_path=tmpdir_str)
    gt_cache._image_url = types.MethodType(
        GOOGLE_IMAGE_URL_REPLACEMENT, gt_cache)
    img_cache, _, _ = gt_cache.image_for_domain(target_domain, 6)

    # Check that the new image_for_domain() call used cached images
    assert gt_cache.cache.keys() == gt.cache.keys()
    assert (img_cache == 255).all()
