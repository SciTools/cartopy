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

import io
import os
import shutil
import warnings

from nose.tools import assert_equal, assert_in, assert_true
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from PIL import Image
import shapely.geometry as sgeom
from six.moves import cPickle as pickle

from cartopy import config
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.io.img_nest as cimg_nest
import cartopy.tests as tests


#: An integer version which should be increased if the test data needs
#: to change in some way.
_TEST_DATA_VERSION = 1
_TEST_DATA_DIR = os.path.join(config["data_dir"],
                              'wmts', 'aerial')
#: A global to determine whether the test data has already been made available
#: in this session.
_TEST_DATA_AVAILABLE = False


def test_world_files():
    func = cimg_nest.Img.world_files
    fname = 'one'
    expected = ['one.w', 'one.W', 'ONE.w', 'ONE.W']
    assert_equal(func(fname), expected)

    fname = 'one.png'
    expected = ['one.pngw', 'one.pgw', 'one.PNGW', 'one.PGW',
                'ONE.pngw', 'ONE.pgw', 'ONE.PNGW', 'ONE.PGW']
    assert_equal(func(fname), expected)

    fname = '/one.png'
    expected = ['/one.pngw', '/one.pgw', '/one.PNGW', '/one.PGW',
                '/ONE.pngw', '/ONE.pgw', '/ONE.PNGW', '/ONE.PGW']
    assert_equal(func(fname), expected)

    fname = '/one/two.png'
    expected = ['/one/two.pngw', '/one/two.pgw',
                '/one/two.PNGW', '/one/two.PGW',
                '/one/TWO.pngw', '/one/TWO.pgw',
                '/one/TWO.PNGW', '/one/TWO.PGW']
    assert_equal(func(fname), expected)

    fname = '/one/two/THREE.png'
    expected = ['/one/two/THREE.pngw', '/one/two/THREE.pgw',
                '/one/two/THREE.PNGW', '/one/two/THREE.PGW',
                '/one/two/three.pngw', '/one/two/three.pgw',
                '/one/two/three.PNGW', '/one/two/three.PGW']
    assert_equal(func(fname), expected)


def _save_world(fname, args):
    _world = ('{x_pix_size}\n'
              '{y_rotation}\n'
              '{x_rotation}\n'
              '{y_pix_size}\n'
              '{x_center}\n'
              '{y_center}\n')
    with open(fname, 'w') as fh:
        fh.write(_world.format(**args))


def test_intersect():
    with tests.temp_dir() as base_dir:
        # Zoom level zero.
        # File 1: Parent space of all images.
        z_0_dir = os.path.join(base_dir, 'z_0')
        os.mkdir(z_0_dir)
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=1, y_center=1)
        im = Image.new('RGB', (50, 50))
        fname = os.path.join(z_0_dir, 'p0.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_0_dir, 'p0.tif')
        im.save(fname)

        # Zoom level one.
        # File 1: complete containment within p0.
        z_1_dir = os.path.join(base_dir, 'z_1')
        os.mkdir(z_1_dir)
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=21, y_center=21)
        im = Image.new('RGB', (30, 30))
        fname = os.path.join(z_1_dir, 'p1.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_1_dir, 'p1.tif')
        im.save(fname)

        # Zoom level two.
        # File 1: intersect right edge with p1 left edge.
        z_2_dir = os.path.join(base_dir, 'z_2')
        os.mkdir(z_2_dir)
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=6, y_center=21)
        im = Image.new('RGB', (5, 5))
        fname = os.path.join(z_2_dir, 'p2-1.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_2_dir, 'p2-1.tif')
        im.save(fname)
        # File 2: intersect upper right corner with p1
        #         lower left corner.
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=6, y_center=6)
        im = Image.new('RGB', (5, 5))
        fname = os.path.join(z_2_dir, 'p2-2.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_2_dir, 'p2-2.tif')
        im.save(fname)
        # File 3: complete containment within p1.
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=41, y_center=41)
        im = Image.new('RGB', (5, 5))
        fname = os.path.join(z_2_dir, 'p2-3.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_2_dir, 'p2-3.tif')
        im.save(fname)
        # File 4: overlap with p1 right edge.
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=76, y_center=61)
        im = Image.new('RGB', (5, 5))
        fname = os.path.join(z_2_dir, 'p2-4.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_2_dir, 'p2-4.tif')
        im.save(fname)
        # File 5: overlap with p1 bottom right corner.
        world = dict(x_pix_size=2, y_rotation=0, x_rotation=0,
                     y_pix_size=2, x_center=76, y_center=76)
        im = Image.new('RGB', (5, 5))
        fname = os.path.join(z_2_dir, 'p2-5.tfw')
        _save_world(fname, world)
        fname = os.path.join(z_2_dir, 'p2-5.tif')
        im.save(fname)

        # Provided in reverse order in order to test the area sorting.
        items = [('dummy-z-2', z_2_dir),
                 ('dummy-z-1', z_1_dir),
                 ('dummy-z-0', z_0_dir)]
        nic = cimg_nest.NestedImageCollection.from_configuration('dummy',
                                                                 None,
                                                                 items)

        names = [collection.name for collection in nic._collections]
        zoom_levels = ['dummy-z-0', 'dummy-z-1', 'dummy-z-2']
        assert_true(names, zoom_levels)

        # Check all images are loaded.
        for zoom, expected_image_count in zip(zoom_levels, [1, 1, 5]):
            images = nic._collections_by_name[zoom].images
            assert_equal(len(images), expected_image_count)

        # Check the image ancestry.
        zoom_levels = ['dummy-z-0', 'dummy-z-1']
        assert_equal(sorted(k[0] for k in nic._ancestry.keys()),
                     zoom_levels)

        expected = [('dummy-z-0', ['p1.tif']),
                    ('dummy-z-1', ['p2-3.tif', 'p2-4.tif', 'p2-5.tif'])]
        for zoom, image_names in expected:
            key = [k for k in nic._ancestry.keys() if k[0] == zoom][0]
            ancestry = nic._ancestry[key]
            fnames = sorted([os.path.basename(item[1].filename)
                             for item in ancestry])
            assert_equal(image_names, fnames)

        # Check image retrieval for specific domain.
        items = [(sgeom.box(20, 20, 80, 80), 3),
                 (sgeom.box(20, 20, 75, 75), 1),
                 (sgeom.box(40, 40, 85, 85), 3)]
        for domain, expected in items:
            result = [image for image in nic.find_images(domain,
                                                         'dummy-z-2')]
            assert_equal(len(result), expected)


def _tile_from_img(img):
    """
    Turns an img into the appropriate x, y, z tile based on its filename.

    Imgs have a filename attribute which is something
    like "lib/cartopy/data/wmts/aerial/z_0/x_0_y0.png"

    """
    _, z = os.path.basename(os.path.dirname(img.filename)).split('_')
    xy, _ = os.path.splitext(os.path.basename(img.filename))
    _, x, _, y = xy.split('_')
    return int(x), int(y), int(z)


class RoundedImg(cimg_nest.Img):
    @staticmethod
    def world_file_extent(*args, **kwargs):
        """
        Takes account for the fact that the image tiles are stored with
        imprecise tfw files.

        """
        extent, pix_size = cimg_nest.Img.world_file_extent(*args, **kwargs)
        # round the extent
        extent = tuple(round(v, 4) for v in extent)
        pix_size = tuple(round(v, 4) for v in pix_size)
        return extent, pix_size


def test_nest():
    crs = cimgt.GoogleTiles().crs
    z0 = cimg_nest.ImageCollection('aerial z0 test', crs)
    z0.scan_dir_for_imgs(os.path.join(_TEST_DATA_DIR, 'z_0'),
                         glob_pattern='*.png', img_class=RoundedImg)

    z1 = cimg_nest.ImageCollection('aerial z1 test', crs)
    z1.scan_dir_for_imgs(os.path.join(_TEST_DATA_DIR, 'z_1'),
                         glob_pattern='*.png', img_class=RoundedImg)

    z2 = cimg_nest.ImageCollection('aerial z2 test', crs)
    z2.scan_dir_for_imgs(os.path.join(_TEST_DATA_DIR, 'z_2'),
                         glob_pattern='*.png', img_class=RoundedImg)

    # make sure all the images from z1 are contained by the z0 image. The
    # only reason this might occur is if the tfw files are handling
    # floating point values badly
    for img in z1.images:
        if not z0.images[0].bbox().contains(img.bbox()):
            raise IOError('The test images aren\'t all "contained" by the '
                          'z0 images, the nest cannot possibly work.\n '
                          'img {!s} not contained by {!s}\nExtents: {!s}; '
                          '{!s}'.format(img, z0.images[0], img.extent,
                                        z0.images[0].extent))
    nest_z0_z1 = cimg_nest.NestedImageCollection('aerial test',
                                                 crs,
                                                 [z0, z1])

    nest = cimg_nest.NestedImageCollection('aerial test', crs, [z0, z1, z2])

    z0_key = ('aerial z0 test', z0.images[0])

    assert_true(z0_key in nest_z0_z1._ancestry.keys())
    assert_equal(len(nest_z0_z1._ancestry), 1)

    # check that it has figured out that all the z1 images are children of
    # the only z0 image
    for img in z1.images:
        key = ('aerial z0 test', z0.images[0])
        assert_in(('aerial z1 test', img), nest_z0_z1._ancestry[key])

    x1_y0_z1, = [img for img in z1.images
                 if img.filename.endswith('z_1/x_1_y_0.png')]

    assert_equal((1, 0, 1), _tile_from_img(x1_y0_z1))

    assert_equal([(2, 0, 2), (2, 1, 2), (3, 0, 2), (3, 1, 2)],
                 sorted([_tile_from_img(img) for z, img in
                         nest.subtiles(('aerial z1 test', x1_y0_z1))]))

    nest_from_config = gen_nest()
    # check that the the images in the nest from configuration are the
    # same as those created by hand.
    for name in nest_z0_z1._collections_by_name.keys():
        for img in nest_z0_z1._collections_by_name[name].images:
            collection = nest_from_config._collections_by_name[name]
            assert_in(img, collection.images)

    assert_equal(nest_z0_z1._ancestry, nest_from_config._ancestry)

    # check that a nest can be pickled and unpickled easily.
    s = io.BytesIO()
    pickle.dump(nest_z0_z1, s)
    s.seek(0)
    nest_z0_z1_from_pickle = pickle.load(s)

    assert_equal(nest_z0_z1._ancestry,
                 nest_z0_z1_from_pickle._ancestry)


def test_img_pickle_round_trip():
    """Check that __getstate__ for Img instances is working correctly."""

    img = cimg_nest.Img('imaginary file', (0, 1, 2, 3), 'lower', (1, 2))
    img_from_pickle = pickle.loads(pickle.dumps(img))
    assert_equal(img, img_from_pickle)
    assert_equal(hasattr(img_from_pickle, '_bbox'), True)


def requires_wmts_data(function):
    """
    A decorator which ensures that the WMTS data is available for
    use in testing.

    """
    aerial = cimgt.MapQuestOpenAerial()

    # get web tiles upto 3 zoom levels deep
    tiles = [(0, 0, 0)]
    for tile in aerial.subtiles((0, 0, 0)):
        tiles.append(tile)
    for tile in tiles[1:]:
        for sub_tile in aerial.subtiles(tile):
            tiles.append(sub_tile)

    fname_template = os.path.join(_TEST_DATA_DIR, 'z_{}', 'x_{}_y_{}.png')

    if not os.path.isdir(_TEST_DATA_DIR):
        os.makedirs(_TEST_DATA_DIR)

    data_version_fname = os.path.join(_TEST_DATA_DIR, 'version.txt')

    test_data_version = None
    try:
        with open(data_version_fname, 'r') as fh:
            test_data_version = int(fh.read().strip())
    except IOError:
        pass
    finally:
        if test_data_version != _TEST_DATA_VERSION:
            warnings.warn('WMTS test data is out of date, regenerating at '
                          '{}.'.format(_TEST_DATA_DIR))
            shutil.rmtree(_TEST_DATA_DIR)
            os.makedirs(_TEST_DATA_DIR)
            with open(data_version_fname, 'w') as fh:
                fh.write(str(_TEST_DATA_VERSION))

    # Download the tiles.
    for tile in tiles:
        x, y, z = tile
        fname = fname_template.format(z, x, y)
        if not os.path.exists(fname):
            if not os.path.isdir(os.path.dirname(fname)):
                os.makedirs(os.path.dirname(fname))

            img, extent, _ = aerial.get_image(tile)
            nx, ny = 256, 256
            x_rng = extent[1] - extent[0]
            y_rng = extent[3] - extent[2]

            pix_size_x = x_rng / nx
            pix_size_y = y_rng / ny

            upper_left_center = (extent[0] + pix_size_x / 2,
                                 extent[2] + pix_size_y / 2)

            pgw_fname = fname[:-4] + '.pgw'
            pgw_keys = {'x_pix_size': np.float64(pix_size_x),
                        'y_rotation': 0,
                        'x_rotation': 0,
                        'y_pix_size': np.float64(pix_size_y),
                        'x_center': np.float64(upper_left_center[0]),
                        'y_center': np.float64(upper_left_center[1]),
                        }
            _save_world(pgw_fname, pgw_keys)
            img.save(fname)

    global _TEST_DATA_AVAILABLE
    _TEST_DATA_AVAILABLE = True

    return function


@requires_wmts_data
def test_find_images():
    z2_dir = os.path.join(_TEST_DATA_DIR, 'z_2')
    img_fname = os.path.join(z2_dir, 'x_2_y_0.png')
    world_file_fname = os.path.join(z2_dir, 'x_2_y_0.pgw')
    img = RoundedImg.from_world_file(img_fname, world_file_fname)

    assert_equal(img.filename, img_fname)
    assert_array_almost_equal(img.extent,
                              (0., 10018754.17139462,
                               10018754.17139462, 20037508.342789244),
                              decimal=4)
    assert_equal(img.origin, 'lower')
    assert_array_equal(img, np.array(Image.open(img.filename)))
    assert_equal(img.pixel_size, (39135.7585, 39135.7585))


@requires_wmts_data
def gen_nest():
    from_config = cimg_nest.NestedImageCollection.from_configuration

    files = [['aerial z0 test', os.path.join(_TEST_DATA_DIR, 'z_0')],
             ['aerial z1 test', os.path.join(_TEST_DATA_DIR, 'z_1')],
             ]

    crs = cimgt.GoogleTiles().crs

    nest_z0_z1 = from_config('aerial test',
                             crs, files, glob_pattern='*.png',
                             img_class=RoundedImg)
    return nest_z0_z1


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
