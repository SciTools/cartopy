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

from __future__ import division

import io
import cPickle as pickle
import os

from nose.tools import assert_equal, assert_raises
import numpy as np
from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt
import shapely.geometry

import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.io.img_nest as cimg_nest


_dname = os.path.dirname
# XXX be more clever in the data directory so that users can define a setting.
_TEST_DATA_DIR = os.path.join(_dname(_dname(__file__)), 'data', 'wmts', 'aerial')


def _tile_from_img(img):
        """
        Turns an img into the appropriate x, y, z tile based on it filename alone.

        Imgs have a filename attribute which is something
        like "lib/cartopy/data/wmts/aerial/z_0/x_0_y0.png"

        """
        _, z = os.path.dirname(img.filename).split('_')
        xy, _ = os.path.splitext(os.path.basename(img.filename))
        _, x, _, y = xy.split('_')
        return int(x), int(y), int(z)


class RoundedImageCollection(cimg_nest.ImageCollection):
    """Takes account for the fact that the image tiles are stored with imprecise tfw files."""
    def _extent_finalize(self, extent, fname):
        return tuple(round(num, 1) for num in extent)


def gen_nest():
    gen_test_data()
    nest_z0_z1_from_config = cimg_nest.NestedImageCollection.from_configuration('aerial test',
                                                                     ccrs.Mercator(),
                                                                     [['aerial z0 test', os.path.join(_TEST_DATA_DIR, 'z_0')],
                                                                      ['aerial z1 test', os.path.join(_TEST_DATA_DIR, 'z_1')],
                                                                      ],
                                                                     glob_pattern='*.png',
                                                                     img_collection_cls=RoundedImageCollection,
                                                                     )
    return nest_z0_z1_from_config


def test_nest():
    nest_z0_z1_from_config = gen_nest()

    z0 = RoundedImageCollection('aerial z0 test', ccrs.Mercator())
    z0.scan_dir_for_imgs(os.path.join(_TEST_DATA_DIR, 'z_0'), glob_pattern='*.png')

    z1 = RoundedImageCollection('aerial z1 test', ccrs.Mercator())
    z1.scan_dir_for_imgs(os.path.join(_TEST_DATA_DIR, 'z_1'), glob_pattern='*.png')

    z2 = RoundedImageCollection('aerial z2 test', ccrs.Mercator())
    z2.scan_dir_for_imgs(os.path.join(_TEST_DATA_DIR, 'z_2'), glob_pattern='*.png')

    # make sure all the images from z1 are contained by the z0 image. The only reason this
    # might occur is if the tfw files are handling floating point values badly
    for img in z1.images:
        if not z0.images[0].bbox().contains(img.bbox()):
            raise IOError("The test images aren't all \"contained\" by the z0 images, the nest cannot "
                          "possibly work.\n "
                          'img %s not contained by %s\nExtents: %s; %s' % (img, z0.images[0],
                                                                           img.extent, z0.images[0].extent))
    nest_z0_z1 = cimg_nest.NestedImageCollection('aerial test',
                                                  ccrs.Mercator(),
                                                  [z0, z1])

    nest = cimg_nest.NestedImageCollection('aerial test',
                                            ccrs.Mercator(),
                                            [z0, z1, z2])

    z0_key = ('aerial z0 test', z0.images[0])

    assert z0_key in nest_z0_z1._ancestry.keys()
    assert_equal(len(nest_z0_z1._ancestry), 1)

    # check that it has figured out that all the z1 images are children of the only z0 image
    for img in z1.images:
        assert ('aerial z1 test', img) in nest_z0_z1._ancestry[('aerial z0 test', z0.images[0])]

    x1_y0_z1, = [img for img in z1.images if img.filename.endswith('z_1/x_1_y_0.png')]

    assert_equal((1, 0, 1), _tile_from_img(x1_y0_z1))

    assert_equal([(2, 0, 2), (2, 1, 2), (3, 0, 2), (3, 1, 2)],
                 sorted([_tile_from_img(img) for z, img in nest.subtiles(('aerial z1 test', x1_y0_z1))]))

    # check that the the images in the nest from configuration are the same as those created by hand.
    for name in nest_z0_z1._collections_by_name.keys():
        for img in nest_z0_z1._collections_by_name[name].images:
            assert img in nest_z0_z1_from_config._collections_by_name[name].images

    assert_equal(nest_z0_z1._ancestry, nest_z0_z1_from_config._ancestry)

    # check that a nest can be pickled and unpickled easily.
    s = io.BytesIO()
    pickle.dump(nest_z0_z1, s)
    s.seek(0)
    nest_z0_z1_from_pickle = pickle.load(s)

    assert_equal(nest_z0_z1._ancestry, nest_z0_z1_from_pickle._ancestry)


def gen_test_data():
    aerial = cimgt.MapQuestOpenAerial()

    # get web tiles upto 3 zoom levels deep
    tiles = [(0, 0, 0)]
    for tile in aerial.subtiles((0, 0, 0)):
        tiles.append(tile)
    for tile in tiles[1:]:
        for sub_tile in aerial.subtiles(tile):
            tiles.append(sub_tile)

    data_dir = os.path.join(_TEST_DATA_DIR, 'z_{}', 'x_{}_y_{}.png')

    # download the tiles
    for tile in tiles:
        x, y, z = tile
        fname = data_dir.format(z, x, y)
        if not os.path.exists(fname):
            if not os.path.isdir(os.path.dirname(fname)):
                os.makedirs(os.path.dirname(fname))

            img, extent, _ = aerial.get_image(tile)
            nx, ny = 256, 256
            x_rng = extent[1] - extent[0]
            y_rng = extent[3] - extent[2]

            pix_size_x = x_rng / nx
            pix_size_y = y_rng / ny
            upper_left_center = extent[0] + pix_size_x/2, (extent[2] + pix_size_y/2)

            tfw_fname = fname[:-4] + '.tfw'
            open(tfw_fname, 'w').write('{}\n{}\n{}\n{}\n{}\n{}'.
                                       format(np.float32(pix_size_x),
                                              0, 0,
                                              np.float32(pix_size_y),
                                              np.float32(upper_left_center[0]),
                                              np.float32(upper_left_center[1]),
                                              )
                                       )
            img.save(fname)


if __name__ == '__main__':
#    import nose
#    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
    test_nest()