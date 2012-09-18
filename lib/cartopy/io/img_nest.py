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


import collections
import os
import itertools
import glob

import numpy
from PIL import Image
from shapely.geometry import box


class Img(collections.namedtuple('Img', ['filename', 'extent', 'origin', 'pixel_size'])):
    """
    Represents a simple geolocated image.

    Note: API is likely to change in the future to include a CRS.

    """
    def __new__(cls, *args, **kwargs):
        # ensure any lists given as args or kwargs are turned into tuples.
        new_args = []
        for item in args:
            if isinstance(item, list):
                item = tuple(item)
            new_args.append(item)
        new_kwargs = {}
        for k, item in kwargs.items():
            if isinstance(item, list):
                item = tuple(item)
            new_kwargs[k] = item
        return super(Img, cls).__new__(cls, *new_args, **new_kwargs)

    def __init__(self, *args, **kwargs):
        r = super(Img, self).__init__(self, *args, **kwargs)
        self._bbox = None

    def bbox(self):
        if self._bbox is None:
            x0, x1, y0, y1 = self.extent
            self._bbox = box(x0, y0, x1, y1)
        return self._bbox


class ImageCollection(object):
    """
    Represents a collection of images at the same logical level (typically zoom level).
    """
    def __init__(self, name, crs, images=None):
        self.name = name
        self.crs = crs
        self.images = images or []

    def find_images(self, target_domain):
        """
        Find the images which exist in this collection which intersect with the target domain.

        Target domain is a shapely polygon in native coordinates.

        """
        for fname, extent, domain, origin in self.images:
            if target_domain.intersects(domain):
                yield fname, extent, domain, crs

    def scan_dir_for_imgs(self, directory, glob_pattern='*.tif'):
        """
        Scan the given directory for images with associated tfw files.
        Append any found results to self.images

        .. note:: Does not recursively look through directories.

        """
        imgs = glob.glob(os.path.join(directory, glob_pattern))

        # maps tilename to [bbox_poly, [children]]
        tiles = {}
        for img in imgs:
            dirname, fname = os.path.split(img)
            orig_tfw = fname[:-3] + 'tfw'
            for tfw in [orig_tfw, orig_tfw.upper()]:
                tfw = os.path.join(dirname, tfw)
                if os.path.exists(tfw):
                    break
            else:
                raise ValueError('Image %s has no tfw' % img)

            lines = open(tfw).readlines()
            pix_size = [float(lines[0]), float(lines[3])]
            pix_rotation = [float(lines[1]), float(lines[2])]
            assert pix_rotation == [0., 0.], 'Rotated pixels not currently supported. Image: %s' % img
            ul_corner = [float(lines[4]), float(lines[5])]

            im = Image.open(img, 'r')
            shape = im.size
            min_x, max_x = sorted([ul_corner[0] - pix_size[0] / 2.,
                                   ul_corner[0] + pix_size[0] * shape[0] - pix_size[0] / 2.])
            min_y, max_y = sorted([ul_corner[1] - pix_size[1] / 2.,
                                   ul_corner[1] + pix_size[1] * shape[1] - pix_size[1] / 2.])
            extent = (min_x, max_x, min_y, max_y)
            self.images.append(Img(img, self._extent_finalize(extent, img), 'lower', tuple(pix_size)))

    def _extent_finalize(self, extent, filename):
        """The final extent of an image is passed through this method. This
        is a good place to implement rounding or some other processing."""
        return extent


class NestedImageCollection(object):
    """
    Represents a complex nest of ImageCollections.

    On construction, the image collections are scanned for ancestry, leading
    to fast image finding capabilities.

    A complex (and time consuming to create) NestedImageCollection instance can
    be saved as a pickle file and subsequently be (quickly) restored.

    There is a simplified creation interface for NestedImageCollection
    ``from_configuration`` for more detail.

    """
    def __init__(self, name, crs, collections, _ancestry=None):
        # NOTE: all collections must have the same crs.
        _collection_names = [collection.name for collection in collections]
        assert len(_collection_names) == len(collections), \
               'The collections must have unique names.'

        self.name = name
        self.crs = crs
        self._collections_by_name = {collection.name: collection for collection in collections}
        self._collections = collections
        self._ancestry = {}
        """maps (collection name, image) to a list of children (collection name, image)."""

        tiles = {}

        if _ancestry is not None:
            self._ancestry = _ancestry
        else:
            for parent_collection, collection in itertools.izip(collections, collections[1:]):
                for parent_image in parent_collection.images:
                    for image in collection.images:
                        if self._is_parent(parent_image, image):
                            # append the child image to the parent's ancestry
                            self._ancestry.setdefault((parent_collection.name, parent_image),
                                                      []).append((collection.name, image))

            # TODO check that the ancestry is in a good state (i.e. that each
            # collection has child images)

    @staticmethod
    def _is_parent(potential_parent_image, image):
        """Returns whether the given Image is the parent of image. Used by __init__."""
        return potential_parent_image.bbox().contains(image.bbox())

    def image_for_domain(self, target_domain, target_z):
        # XXX Copied from cartopy.io.img_tiles
        if target_z not in self._collections_by_name:
            # TODO: Handle integer depths also?
            raise ValueError('%s is not one of the possible collections.' % target_z)

        tiles = []
        for tile in self.find_images(target_domain, target_z):
            try:
                img, extent, origin = self.get_image(tile)
            except IOError:
                continue

            img = numpy.array(img)

            x = numpy.linspace(extent[0], extent[1], img.shape[1], endpoint=False)
            y = numpy.linspace(extent[2], extent[3], img.shape[0], endpoint=False)
            tiles.append([numpy.array(img), x, y, origin])

        from cartopy.io.img_tiles import _merge_tiles
        img, extent, origin = _merge_tiles(tiles)
        return img, extent, origin

    def find_images(self, target_domain, target_z, start_tiles=None):
        # XXX Copied from cartopy.io.img_tiles
        if target_z not in self._collections_by_name:
            # TODO: Handle integer depths also?
            raise ValueError('%s is not one of the possible collections.' % target_z)

        if start_tiles is None:
            start_tiles = ((self._collections[0].name, img) for img in self._collections[0].images)

        for start_tile in start_tiles:
            # recursively drill down to the images at the target zoom
            domain = start_tile[1].bbox()
            if domain.intersects(target_domain):
                if start_tile[0] == target_z:
                        yield start_tile
                else:
                    for tile in self.subtiles(start_tile):
                        for result in self.find_images(target_domain, target_z, start_tiles=[tile]):
                            yield result

    def subtiles(self, collection_image):
        return iter(self._ancestry.get(collection_image, []))

    desired_tile_form = 'RGB'
    def get_image(self, collection_image):
        img = collection_image[1]
        img_data = Image.open(img.filename)
        img_data = img_data.convert(self.desired_tile_form)
        return img_data, img.extent, img.origin

    @classmethod
    def from_configuration(cls, name, crs, name_dir_pairs,
                           glob_pattern='*.tif', img_collection_cls=ImageCollection,
                           ):
        """
        Creates a NestedImageCollection given the [collection name, directory] pairs.
        This is very convenient functionality for simple configuration level creation
        of this complex object.

        For example, to produce a nested collection of OS map tiles::

            r = NestedImageCollection.from_configuration('os',
                                                 ccrs.OSGB(),
                                                 [['OS 1:1,000,000', '/directory/to/1_to_1m'],
                                                  ['OS 1:250,000', '/directory/to/1_to_250k'],
                                                  ['OS 1:50,000', '/directory/to/1_to_50k'],
                                                  ],
                                                 )

        """
        collections = []
        for collection_name, collection_dir in name_dir_pairs:
            collection = img_collection_cls(collection_name, crs)
            collection.scan_dir_for_imgs(collection_dir, glob_pattern=glob_pattern)
            collections.append(collection)
        return cls(name, crs, collections)