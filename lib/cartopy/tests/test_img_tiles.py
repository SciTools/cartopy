# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import hashlib
import os
import types

import numpy as np
from numpy.testing import assert_array_almost_equal as assert_arr_almost
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
    x_y_z_f_h = [
        (30, 18, 6, '30_18_6.npy', '5f4bcb9e2d21931ad67086a96b0ea679'),
        (30, 19, 6, '30_19_6.npy', '5c26cd7585f869e6d9a6d32530c4ff62'),
        (30, 20, 6, '30_20_6.npy', 'e0861f5af0c3ede62c25213a482b69bf'),
        (30, 21, 6, '30_21_6.npy', '562e0ed47fb5b89b04e26db712ddd224'),
        (31, 18, 6, '31_18_6.npy', '0b8d038ad535d1d0b9671d11e1fed688'),
        (31, 19, 6, '31_19_6.npy', '403bee164411405daea04be26f53af82'),
        (31, 20, 6, '31_20_6.npy', '9a65b7852c4fb7cf38df8c14ded37cea'),
        (31, 21, 6, '31_21_6.npy', '8431309ba8f05e76fd1df47807fd1134'),
        (32, 18, 6, '32_18_6.npy', 'f11e04a071a60266cd197fc46b301f32'),
        (32, 19, 6, '32_19_6.npy', 'ba956e8daf68968d97ce04994d8a5be2'),
        (32, 20, 6, '32_20_6.npy', 'e269f50f4cc79858e4f31bd39bdbad07'),
        (32, 21, 6, '32_21_6.npy', '76f23795f0724af4522a78e9c9d0fa9c'),
        (33, 18, 6, '33_18_6.npy', '132ec5f64985b783c79facdd85de1927'),
        (33, 19, 6, '33_19_6.npy', 'e4e855d2376cf2faa20a8d7da4114a61'),
        (33, 20, 6, '33_20_6.npy', '02009bce9adab5f05a64f4bed3fa5218'),
        (33, 21, 6, '33_21_6.npy', '309beaef09160ce1d849ba77a2d79246')
    ]

    # Check the results
    cache_dir = os.path.join(tmpdir_str, "GoogleTiles")
    files = [i for i in os.listdir(cache_dir)]
    hashes = {
        f:
        hashlib.md5(
            open(os.path.join(cache_dir, f), mode="rb").read()
        ).hexdigest()
        for f in files
    }

    assert sorted(files) == [f for x, y, z, f, h in x_y_z_f_h]
    assert set(files) == gt.cache

    assert sorted(hashes.values()) == sorted([
        h for x, y, z, f, h in x_y_z_f_h
    ])

    # Update images in cache (all white)
    for f in files:
        filename = os.path.join(cache_dir, f)
        img = np.load(filename, allow_pickle=True)
        img.fill(255)
        np.save(filename, img, allow_pickle=True)

    gt_cache = cimgt.GoogleTiles(cache_path=tmpdir_str)
    gt_cache._image_url = types.MethodType(
        GOOGLE_IMAGE_URL_REPLACEMENT, gt_cache)
    img_cache, _, _ = gt_cache.image_for_domain(target_domain, 6)

    # Check that the new image_for_domain() call used cached images
    assert gt_cache.cache == gt.cache
    assert (img_cache == 255).all()
