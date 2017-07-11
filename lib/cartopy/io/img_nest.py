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

import collections
import glob
import os.path

import numpy as np
from PIL import Image
import shapely.geometry as sgeom
from six.moves import zip


_img_class_attrs = ['filename', 'extent', 'origin', 'pixel_size']


class Img(collections.namedtuple('Img', _img_class_attrs)):
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
        """
        Represents a simple geo-located image.

        Args:

        * filename:
            Filename of the image tile.

        * extent:
            The (x_lower, x_upper, y_lower, y_upper) extent of the image
            in units of the native projection.

        * origin:
            Name of the origin.

        * pixel_size:
            The (x_scale, y_scale) pixel width, in units of the native
            projection per pixel.

        .. note::
            API is likely to change in the future to include a CRS.

        """
        self._bbox = None

    def __getstate__(self):
        """
        Override the default to ensure when pickling that any new attributes
        introduced are included in the pickled object.

        """
        return self.__dict__

    def bbox(self):
        """
        Return a :class:`~shapely.geometry.polygon.Polygon` instance for
        this image's extents.

        """
        if self._bbox is None:
            x0, x1, y0, y1 = self.extent
            self._bbox = sgeom.box(x0, y0, x1, y1)
        return self._bbox

    @staticmethod
    def world_files(fname):
        """
        Determine potential world filename combinations, without checking
        their existence.

        For example, a '*.tif' file may have one of the following
        popular conventions for world file extensions '*.tifw',
        '*.tfw', '*.TIFW' or '*.TFW'.

        Given the possible world file extensions, the upper case basename
        combinations are also generated. For example, the file 'map.tif'
        will generate the following world file variations, 'map.tifw',
        'map.tfw', 'map.TIFW', 'map.TFW', 'MAP.tifw', 'MAP.tfw', 'MAP.TIFW'
        and 'MAP.TFW'.

        Args:

        * fname:
            Name of the file for which to get all the possible world
            filename combinations.

        Returns:
            A list of possible world filename combinations.

        Examples:

        >>> from cartopy.io.img_nest import Img
        >>> Img.world_files('img.png')[:6]
        ['img.pngw', 'img.pgw', 'img.PNGW', 'img.PGW', 'IMG.pngw', 'IMG.pgw']
        >>> Img.world_files('/path/to/img.TIF')[:2]
        ['/path/to/img.tifw', '/path/to/img.tfw']
        >>> Img.world_files('/path/to/img/with_no_extension')[0]
        '/path/to/img/with_no_extension.w'

        """
        froot, fext = os.path.splitext(fname)
        # If there was no extension to the filename.
        if froot == fname:
            result = ['{}.{}'.format(fname, 'w'),
                      '{}.{}'.format(fname, 'W')]
        else:
            fext = fext[1::].lower()
            if len(fext) < 3:
                result = ['{}.{}'.format(fname, 'w'),
                          '{}.{}'.format(fname, 'W')]
            else:
                fext_types = [fext + 'w', fext[0] + fext[-1] + 'w']
                fext_types.extend([ext.upper() for ext in fext_types])
                result = ['{}.{}'.format(froot, ext) for ext in fext_types]

        def _convert_basename(name):
            dirname, basename = os.path.dirname(name), os.path.basename(name)
            base, ext = os.path.splitext(basename)
            if base == base.upper():
                result = base.lower() + ext
            else:
                result = base.upper() + ext
            if dirname:
                result = os.path.join(dirname, result)
            return result

        result += [_convert_basename(r) for r in result]
        return result

    def __array__(self):
        return np.array(Image.open(self.filename))

    @classmethod
    def from_world_file(cls, img_fname, world_fname):
        """
        Return an Img instance from the given image filename and
        worldfile filename.

        """
        im = Image.open(img_fname)
        with open(world_fname) as world_fh:
            extent, pix_size = cls.world_file_extent(world_fh, im.size)
        if hasattr(im, 'close'):
            im.close()
        return cls(img_fname, extent, 'lower', pix_size)

    @staticmethod
    def world_file_extent(worldfile_handle, im_shape):
        """
        Return the extent ``(x0, x1, y0, y1)`` and pixel size
        ``(x_width, y_width)`` as defined in the given worldfile file handle
        and associated image shape ``(x, y)``.

        """
        lines = worldfile_handle.readlines()
        if len(lines) != 6:
            raise ValueError('Only world files with 6 lines are supported.')

        pix_size = (float(lines[0]), float(lines[3]))
        pix_rotation = (float(lines[1]), float(lines[2]))
        if pix_rotation != (0., 0.):
            raise ValueError('Rotated pixels in world files is not currently '
                             'supported.')
        ul_corner = (float(lines[4]), float(lines[5]))

        min_x, max_x = (ul_corner[0] - pix_size[0]/2.,
                        ul_corner[0] + pix_size[0]*im_shape[0] -
                        pix_size[0]/2.)
        min_y, max_y = (ul_corner[1] - pix_size[1]/2.,
                        ul_corner[1] + pix_size[1]*im_shape[1] -
                        pix_size[1]/2.)
        return (min_x, max_x, min_y, max_y), pix_size


class ImageCollection(object):
    def __init__(self, name, crs, images=None):
        """
        Represents a collection of images at the same logical level.

        Typically these are images at the same zoom level or resolution.

        Args:

        * name:
            The name of the image collection.

        * crs:
            The :class:`~cartopy.crs.Projection` instance.

        Kwargs:

        * images:
            A list of one or more :class:`~cartopy.io.img_nest.Img` instances.

        """
        self.name = name
        self.crs = crs
        self.images = images or []

    def scan_dir_for_imgs(self, directory, glob_pattern='*.tif',
                          img_class=Img):
        """
        Search the given directory for the associated world files
        of the image files.

        Args:

        * directory:
            The directory path to search for image files.

        Kwargs:

        * glob_pattern:
            The image filename glob pattern to search with.
            Defaults to '*.tif'.

        * img_class
            The class used to construct each image in the Collection.

        .. note::
            Does not recursively search sub-directories.

        """
        imgs = glob.glob(os.path.join(directory, glob_pattern))

        for img in imgs:
            dirname, fname = os.path.split(img)
            worlds = img_class.world_files(fname)
            for fworld in worlds:
                fworld = os.path.join(dirname, fworld)
                if os.path.exists(fworld):
                    break
            else:
                msg = 'Image file {!r} has no associated world file'
                raise ValueError(msg.format(img))

            self.images.append(img_class.from_world_file(img, fworld))


class NestedImageCollection(object):
    def __init__(self, name, crs, collections, _ancestry=None):
        """
        Represents a complex nest of ImageCollections.

        On construction, the image collections are scanned for ancestry,
        leading to fast image finding capabilities.

        A complex (and time consuming to create) NestedImageCollection instance
        can be saved as a pickle file and subsequently be (quickly) restored.

        There is a simplified creation interface for NestedImageCollection
        ``from_configuration`` for more detail.

        Args:

        * name:
            The name of the nested image collection.

        * crs:
            The native :class:`~cartopy.crs.Projection` of all the image
            collections.

        * collections:
            A list of one or more :class:`~cartopy.io.img_nest.ImageCollection`
            instances.

        """
        # NOTE: all collections must have the same crs.
        _names = set([collection.name for collection in collections])
        assert len(_names) == len(collections), \
            'The collections must have unique names.'

        self.name = name
        self.crs = crs
        self._collections_by_name = {collection.name: collection
                                     for collection in collections}

        def sort_func(c):
            return np.max([image.bbox().area for image in c.images])

        self._collections = sorted(collections, key=sort_func, reverse=True)
        self._ancestry = {}
        """
        maps (collection name, image) to a list of children
        (collection name, image).
        """
        if _ancestry is not None:
            self._ancestry = _ancestry
        else:
            parent_wth_children = zip(self._collections,
                                      self._collections[1:])
            for parent_collection, collection in parent_wth_children:
                for parent_image in parent_collection.images:
                    for image in collection.images:
                        if self._is_parent(parent_image, image):
                            # append the child image to the parent's ancestry
                            key = (parent_collection.name, parent_image)
                            self._ancestry.setdefault(key, []).append(
                                (collection.name, image))

            # TODO check that the ancestry is in a good state (i.e. that each
            # collection has child images)

    @staticmethod
    def _is_parent(parent, child):
        """
        Returns whether the given Image is the parent of image.
        Used by __init__.

        """
        result = False
        pbox = parent.bbox()
        cbox = child.bbox()
        if pbox.area > cbox.area:
            result = pbox.intersects(cbox) and not pbox.touches(cbox)
        return result

    def image_for_domain(self, target_domain, target_z):
        """
        Determine the image that provides complete coverage of target location.

        The composed image is merged from one or more image tiles that overlay
        the target location and provide complete image coverage of the target
        location.

        Args:

        * target_domain:
            A :class:`~shapely.geometry.linestring.LineString` instance that
            specifies the target location requiring image coverage.

        * target_z:
            The name of the target
            :class`~cartopy.io.img_nest.ImageCollection` which specifies the
            target zoom level (resolution) of the required images.

        Returns:
            A tuple containing three items, consisting of the target
            location :class:`numpy.ndarray` image data, the
            (x-lower, x-upper, y-lower, y-upper) extent of the image, and the
            origin for the target location.

        """
        # XXX Copied from cartopy.io.img_tiles
        if target_z not in self._collections_by_name:
            # TODO: Handle integer depths also?
            msg = '{!r} is not one of the possible collections.'
            raise ValueError(msg.format(target_z))

        tiles = []
        for tile in self.find_images(target_domain, target_z):
            try:
                img, extent, origin = self.get_image(tile)
            except IOError:
                continue

            img = np.array(img)

            x = np.linspace(extent[0], extent[1], img.shape[1],
                            endpoint=False)
            y = np.linspace(extent[2], extent[3], img.shape[0],
                            endpoint=False)
            tiles.append([np.array(img), x, y, origin])

        from cartopy.io.img_tiles import _merge_tiles
        img, extent, origin = _merge_tiles(tiles)
        return img, extent, origin

    def find_images(self, target_domain, target_z, start_tiles=None):
        """
        A generator that finds all images that overlap the bounded
        target location.

        Args:

        * target_domain:
            A :class:`~shapely.geometry.linestring.LineString` instance that
            specifies the target location requiring image coverage.

        * target_z:
            The name of the target
            :class:`~cartopy.io.img_nest.ImageCollection` which specifies
            the target zoom level (resolution) of the required images.

        Kwargs:

        * start_tiles:
            A list of one or more tuple pairs, composed of a
            :class:`~cartopy.io.img_nest.ImageCollection` name and an
            :class:`~cartopy.io.img_nest.Img` instance, from which to search
            for the target images.

        Returns:
            A generator tuple pair composed of a
            :class:`~cartopy.io.img_nest.ImageCollection` name and an
            :class:`~cartopy.io.img_nest.Img` instance.

        """
        # XXX Copied from cartopy.io.img_tiles
        if target_z not in self._collections_by_name:
            # TODO: Handle integer depths also?
            msg = '{!r} is not one of the possible collections.'
            raise ValueError(msg.format(target_z))

        if start_tiles is None:
            start_tiles = ((self._collections[0].name, img)
                           for img in self._collections[0].images)

        for start_tile in start_tiles:
            # recursively drill down to the images at the target zoom
            domain = start_tile[1].bbox()
            if target_domain.intersects(domain) and \
                    not target_domain.touches(domain):
                if start_tile[0] == target_z:
                        yield start_tile
                else:
                    for tile in self.subtiles(start_tile):
                        for result in self.find_images(target_domain,
                                                       target_z,
                                                       start_tiles=[tile]):
                            yield result

    def subtiles(self, collection_image):
        """
        Find the higher resolution image tiles that compose this parent
        image tile.

        Args:

        * collection_image:
            A tuple pair containing the parent
            :class:`~cartopy.io.img_nest.ImageCollection` name and
            :class:`~cartopy.io.img_nest.Img` instance.

        Returns:
            An iterator of tuple pairs containing the higher resolution child
            :class:`~cartopy.io.img_nest.ImageCollection` name and
            :class:`~cartopy.io.img_nest.Img` instance that compose the parent.

        """
        return iter(self._ancestry.get(collection_image, []))

    desired_tile_form = 'RGB'

    def get_image(self, collection_image):
        """
        Retrieve the data of the target image from file.

        .. note::
          The format of the retrieved image file data is controlled by
          :attr:`~cartopy.io.img_nest.NestedImageCollection.desired_tile_form`,
          which defaults to 'RGB' format.

        Args:

        * collection_image:
            A tuple pair containing the target
            :class:`~cartopy.io.img_nest.ImageCollection` name and
            :class:`~cartopy.io.img_nest.Img` instance.

        Returns:
            A tuple containing three items, consisting of the associated image
            file data, the (x_lower, x_upper, y_lower, y_upper) extent of the
            image, and the image origin.

        """
        img = collection_image[1]
        img_data = Image.open(img.filename)
        img_data = img_data.convert(self.desired_tile_form)
        return img_data, img.extent, img.origin

    @classmethod
    def from_configuration(cls, name, crs, name_dir_pairs,
                           glob_pattern='*.tif',
                           img_class=Img):
        """
        Creates a :class:`~cartopy.io.img_nest.NestedImageCollection` instance
        given the list of image collection name and directory path pairs.

        This is very convenient functionality for simple configuration level
        creation of this complex object.

        For example, to produce a nested collection of OS map tiles::

            files = [['OS 1:1,000,000', '/directory/to/1_to_1m'],
                     ['OS 1:250,000', '/directory/to/1_to_250k'],
                     ['OS 1:50,000', '/directory/to/1_to_50k'],
                    ]
            r = NestedImageCollection.from_configuration('os',
                                                         ccrs.OSGB(),
                                                         files)

        .. important::
            The list of image collection name and directory path pairs must be
            given in increasing resolution order i.e. from low resolution to
            high resolution.

        Args:

        * name:
            The name for the
            :class:`~cartopy.io.img_nest.NestedImageCollection` instance.

        * crs:
            The :class:`~cartopy.crs.Projection` of the image collection.

        * name_dir_pairs:
            A list of image collection name and directory path pairs.

        Kwargs:

        * glob_pattern:
            The image collection filename glob pattern.
            Defaults to '*.tif'.

        * img_class:
            The class of images created in the image collection.

        Returns:
            A :class:`~cartopy.io.img_nest.NestedImageCollection` instance.

        """
        collections = []
        for collection_name, collection_dir in name_dir_pairs:
            collection = ImageCollection(collection_name, crs)
            collection.scan_dir_for_imgs(collection_dir,
                                         glob_pattern=glob_pattern,
                                         img_class=img_class)
            collections.append(collection)
        return cls(name, crs, collections)
