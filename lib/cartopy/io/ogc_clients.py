# (C) British Crown Copyright 2014 - 2017, Met Office
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
"""
Implements RasterSource classes which can retrieve imagery from web services
such as WMS and WMTS.

The matplotlib interface can make use of RasterSources via the
:meth:`cartopy.mpl.geoaxes.GeoAxes.add_raster` method,
with additional specific methods which make use of this for WMS and WMTS
(:meth:`~cartopy.mpl.geoaxes.GeoAxes.add_wms` and
:meth:`~cartopy.mpl.geoaxes.GeoAxes.add_wmts`). An example of using WMTS in
this way can be found at :ref:`examples-wmts`.


"""

from __future__ import (absolute_import, division, print_function)

import six

import collections
import io
import math
import warnings
import weakref
from xml.etree import ElementTree

from PIL import Image
import numpy as np
import shapely.geometry as sgeom

try:
    from owslib.wms import WebMapService
    from owslib.wfs import WebFeatureService
    import owslib.util
    import owslib.wmts
    _OWSLIB_AVAILABLE = True
except ImportError:
    WebMapService = None
    WebFeatureService = None
    _OWSLIB_AVAILABLE = False

import cartopy.crs as ccrs
from cartopy.io import LocatedImage, RasterSource
from cartopy.img_transform import warp_array


_OWSLIB_REQUIRED = 'OWSLib is required to use OGC web services.'

# Hardcode some known EPSG codes for now.
# The order given here determines the preferred SRS for WMS retrievals.
_CRS_TO_OGC_SRS = collections.OrderedDict(
    [(ccrs.PlateCarree(), 'EPSG:4326'),
     (ccrs.Mercator.GOOGLE, 'EPSG:900913'),
     (ccrs.OSGB(), 'EPSG:27700')])

# Standard pixel size of 0.28 mm as defined by WMTS.
METERS_PER_PIXEL = 0.28e-3

_WGS84_METERS_PER_UNIT = 2 * math.pi * 6378137 / 360

METERS_PER_UNIT = {
    'urn:ogc:def:crs:EPSG::27700': 1,
    'urn:ogc:def:crs:EPSG::900913': 1,
    'urn:ogc:def:crs:OGC:1.3:CRS84': _WGS84_METERS_PER_UNIT,
    'urn:ogc:def:crs:EPSG::3031': 1,
    'urn:ogc:def:crs:EPSG::3413': 1
}

_URN_TO_CRS = collections.OrderedDict([
    ('urn:ogc:def:crs:OGC:1.3:CRS84', ccrs.PlateCarree()),
    ('urn:ogc:def:crs:EPSG::4326', ccrs.PlateCarree()),
    ('urn:ogc:def:crs:EPSG::900913', ccrs.GOOGLE_MERCATOR),
    ('urn:ogc:def:crs:EPSG::27700', ccrs.OSGB()),
    ('urn:ogc:def:crs:EPSG::3031', ccrs.Stereographic(
        central_latitude=-90,
        true_scale_latitude=-71)),
    ('urn:ogc:def:crs:EPSG::3413', ccrs.Stereographic(
        central_longitude=-45,
        central_latitude=90,
        true_scale_latitude=70))
])

# XML namespace definitions
_MAP_SERVER_NS = '{http://mapserver.gis.umn.edu/mapserver}'
_GML_NS = '{http://www.opengis.net/gml}'


def _warped_located_image(image, source_projection, source_extent,
                          output_projection, output_extent, target_resolution):
    """
    Reproject an Image from one source-projection and extent to another.

    Returns:
        A reprojected LocatedImage, the extent of which is >= the requested
        'output_extent'.

    """
    if source_projection == output_projection:
        extent = output_extent
    else:
        # Convert Image to numpy array (flipping so that origin
        # is 'lower').
        img, extent = warp_array(np.asanyarray(image)[::-1],
                                 source_proj=source_projection,
                                 source_extent=source_extent,
                                 target_proj=output_projection,
                                 target_res=target_resolution,
                                 target_extent=output_extent,
                                 mask_extrapolated=True)

        # Convert arrays with masked RGB(A) values to non-masked RGBA
        # arrays, setting the alpha channel to zero for masked values.
        # This avoids unsightly grey boundaries appearing when the
        # extent is limited (i.e. not global).
        if np.ma.is_masked(img):
            if img.shape[2:3] == (3,):
                # RGB
                old_img = img
                img = np.zeros(img.shape[:2] + (4,), dtype=img.dtype)
                img[:, :, 0:3] = old_img
                img[:, :, 3] = ~ np.any(old_img.mask, axis=2)
                if img.dtype.kind == 'u':
                    img[:, :, 3] *= 255
            elif img.shape[2:3] == (4,):
                # RGBA
                img[:, :, 3] = np.where(np.any(img.mask, axis=2), 0,
                                        img[:, :, 3])
                img = img.data

        # Convert warped image array back to an Image, undoing the
        # earlier flip.
        image = Image.fromarray(img[::-1])

    return LocatedImage(image, extent)


def _target_extents(extent, requested_projection, available_projection):
    """
    Translate the requested extent in the display projection into a list of
    extents in the projection available from the service (multiple if it
    crosses seams).

    The extents are represented as (min_x, max_x, min_y, max_y).

    """
    # Start with the requested area.
    min_x, max_x, min_y, max_y = extent
    target_box = sgeom.box(min_x, min_y, max_x, max_y)

    # If the requested area (i.e. target_box) is bigger (or nearly bigger) than
    # the entire output requested_projection domain, then we erode the request
    # area to avoid re-projection instabilities near the projection boundary.
    buffered_target_box = target_box.buffer(requested_projection.threshold,
                                            resolution=1)
    fudge_mode = buffered_target_box.contains(requested_projection.domain)
    if fudge_mode:
        target_box = requested_projection.domain.buffer(
            -requested_projection.threshold)

    # Translate the requested area into the server projection.
    polys = available_projection.project_geometry(target_box,
                                                  requested_projection)

    # Return the polygons' rectangular bounds as extent tuples.
    target_extents = []
    for poly in polys:
        min_x, min_y, max_x, max_y = poly.bounds
        if fudge_mode:
            # If we shrunk the request area before, then here we
            # need to re-inflate.
            radius = min(max_x - min_x, max_y - min_y) / 5.0
            radius = min(radius, available_projection.threshold * 15)
            poly = poly.buffer(radius, resolution=1)
            # Prevent the expanded request going beyond the
            # limits of the requested_projection.
            poly = available_projection.domain.intersection(poly)
            min_x, min_y, max_x, max_y = poly.bounds
        target_extents.append((min_x, max_x, min_y, max_y))

    return target_extents


class WMSRasterSource(RasterSource):
    """
    A WMS imagery retriever which can be added to a map.

    .. note:: Requires owslib and Pillow to work.

    .. note::

        No caching of retrieved maps is done with this WMSRasterSource.

        To reduce load on the WMS server it is encouraged to tile
        map requests and subsequently stitch them together to recreate
        a single raster, thus allowing for a more aggressive caching scheme,
        but this WMSRasterSource does not currently implement WMS tile
        fetching.

        Whilst not the same service, there is also a WMTSRasterSource which
        makes use of tiles and comes with built-in caching for fast repeated
        map retrievals.

    """

    def __init__(self, service, layers, getmap_extra_kwargs=None):
        """
        Parameters
        ----------
        service : string or WebMapService instance
            The WebMapService instance, or URL of a WMS service, from whence
            to retrieve the image.
        layers : string or list of strings
            The name(s) of layers to use from the WMS service.
        getmap_extra_kwargs : dict or None
            Extra keywords to pass through to the service's getmap method.
            If None, a dictionary with ``{'transparent': True}`` will
            be defined.

        """
        if WebMapService is None:
            raise ImportError(_OWSLIB_REQUIRED)

        if isinstance(service, six.string_types):
            service = WebMapService(service)

        if isinstance(layers, six.string_types):
            layers = [layers]

        if getmap_extra_kwargs is None:
            getmap_extra_kwargs = {'transparent': True}

        if len(layers) == 0:
            raise ValueError('One or more layers must be defined.')
        for layer in layers:
            if layer not in service.contents:
                raise ValueError('The {!r} layer does not exist in '
                                 'this service.'.format(layer))

        #: The OWSLib WebMapService instance.
        self.service = service

        #: The names of the layers to fetch.
        self.layers = layers

        #: Extra kwargs passed through to the service's getmap request.
        self.getmap_extra_kwargs = getmap_extra_kwargs

    def _native_srs(self, projection):
        # Return the SRS which corresponds to the given projection when
        # known, otherwise return None.
        return _CRS_TO_OGC_SRS.get(projection)

    def _fallback_proj_and_srs(self):
        """
        Return a :class:`cartopy.crs.Projection` and corresponding
        SRS string in which the WMS service can supply the requested
        layers.

        """
        contents = self.service.contents
        for proj, srs in six.iteritems(_CRS_TO_OGC_SRS):
            missing = any(srs not in contents[layer].crsOptions for
                          layer in self.layers)
            if not missing:
                break
        if missing:
            raise ValueError('The requested layers are not available in a '
                             'known SRS.')
        return proj, srs

    def validate_projection(self, projection):
        if self._native_srs(projection) is None:
            self._fallback_proj_and_srs()

    def _image_and_extent(self, wms_proj, wms_srs, wms_extent, output_proj,
                          output_extent, target_resolution):
        min_x, max_x, min_y, max_y = wms_extent
        wms_image = self.service.getmap(layers=self.layers,
                                        srs=wms_srs,
                                        bbox=(min_x, min_y, max_x, max_y),
                                        size=target_resolution,
                                        format='image/png',
                                        **self.getmap_extra_kwargs)
        wms_image = Image.open(io.BytesIO(wms_image.read()))

        return _warped_located_image(wms_image, wms_proj, wms_extent,
                                     output_proj, output_extent,
                                     target_resolution)

    def fetch_raster(self, projection, extent, target_resolution):
        target_resolution = [int(np.ceil(val)) for val in target_resolution]
        wms_srs = self._native_srs(projection)
        if wms_srs is not None:
            wms_proj = projection
            wms_extents = [extent]
        else:
            # The SRS for the requested projection is not known, so
            # attempt to use the fallback and perform the necessary
            # transformations.
            wms_proj, wms_srs = self._fallback_proj_and_srs()

            # Calculate the bounding box(es) in WMS projection.
            wms_extents = _target_extents(extent, projection, wms_proj)

        located_images = []
        for wms_extent in wms_extents:
            located_images.append(self._image_and_extent(wms_proj, wms_srs,
                                                         wms_extent,
                                                         projection, extent,
                                                         target_resolution))
        return located_images


class WMTSRasterSource(RasterSource):
    """
    A WMTS imagery retriever which can be added to a map.

    Uses tile caching for fast repeated map retrievals.

    .. note:: Requires owslib and Pillow to work.

    """

    _shared_image_cache = weakref.WeakKeyDictionary()
    """
    A nested mapping from WMTS, layer name, tile matrix name, tile row
    and tile column to the resulting PIL image::

        {wmts: {(layer_name, tile_matrix_name): {(row, column): Image}}}

    This provides a significant boost when producing multiple maps of the
    same projection or with an interactive figure.

    """

    def __init__(self, wmts, layer_name, gettile_extra_kwargs=None):
        """
        Args:

            * wmts - The URL of the WMTS, or an
                     owslib.wmts.WebMapTileService instance.
            * layer_name - The name of the layer to use.

        Kwargs:

            * gettile_extra_kwargs : dict or None
                Extra keywords (e.g. time) to pass through to the
                service's gettile method.

        """
        if WebMapService is None:
            raise ImportError(_OWSLIB_REQUIRED)

        if not (hasattr(wmts, 'tilematrixsets') and
                hasattr(wmts, 'contents') and
                hasattr(wmts, 'gettile')):
            wmts = owslib.wmts.WebMapTileService(wmts)

        try:
            layer = wmts.contents[layer_name]
        except KeyError:
            raise ValueError('Invalid layer name {!r} for WMTS at {!r}'.format(
                layer_name, wmts.url))

        #: The OWSLib WebMapTileService instance.
        self.wmts = wmts

        #: The layer to fetch.
        self.layer = layer

        #: Extra kwargs passed through to the service's gettile request.
        if gettile_extra_kwargs is None:
            gettile_extra_kwargs = {}
        self.gettile_extra_kwargs = gettile_extra_kwargs

        self._matrix_set_name_map = {}

    def _matrix_set_name(self, target_projection):
        key = id(target_projection)
        matrix_set_name = self._matrix_set_name_map.get(key)
        if matrix_set_name is None:
            if hasattr(self.layer, 'tilematrixsetlinks'):
                matrix_set_names = self.layer.tilematrixsetlinks.keys()
            else:
                matrix_set_names = self.layer.tilematrixsets

            def find_projection(match_projection):
                result = None
                for tile_matrix_set_name in matrix_set_names:
                    matrix_sets = self.wmts.tilematrixsets
                    tile_matrix_set = matrix_sets[tile_matrix_set_name]
                    crs_urn = tile_matrix_set.crs
                    tms_crs = _URN_TO_CRS.get(crs_urn)
                    if tms_crs == match_projection:
                        result = tile_matrix_set_name
                        break
                return result

            # First search for a matrix set in the target projection.
            matrix_set_name = find_projection(target_projection)
            if matrix_set_name is None:
                # Search instead for a set in _any_ projection we can use.
                for possible_projection in _URN_TO_CRS.values():
                    # Look for supported projections (in a preferred order).
                    matrix_set_name = find_projection(possible_projection)
                    if matrix_set_name is not None:
                        break
                if matrix_set_name is None:
                    # Fail completely.
                    available_urns = sorted(set(
                        self.wmts.tilematrixsets[name].crs
                        for name in matrix_set_names))
                    msg = 'Unable to find tile matrix for projection.'
                    msg += '\n    Projection: ' + str(target_projection)
                    msg += '\n    Available tile CRS URNs:'
                    msg += '\n        ' + '\n        '.join(available_urns)
                    raise ValueError(msg)
            self._matrix_set_name_map[key] = matrix_set_name
        return matrix_set_name

    def validate_projection(self, projection):
        self._matrix_set_name(projection)

    def fetch_raster(self, projection, extent, target_resolution):
        matrix_set_name = self._matrix_set_name(projection)
        wmts_projection = _URN_TO_CRS[
            self.wmts.tilematrixsets[matrix_set_name].crs]

        if wmts_projection == projection:
            wmts_extents = [extent]
        else:
            # Calculate (possibly multiple) extents in the given projection.
            wmts_extents = _target_extents(extent, projection, wmts_projection)
            # Bump resolution by a small factor, as a weak alternative to
            # delivering a minimum projected resolution.
            # Generally, the desired area is smaller than the enclosing extent
            # in projection space and may have varying scaling, so the ideal
            # solution is a hard problem !
            resolution_factor = 1.4
            target_resolution = np.array(target_resolution) * resolution_factor

        width, height = target_resolution
        located_images = []
        for wmts_desired_extent in wmts_extents:
            # Calculate target resolution for the actual polygon.  Note that
            # this gives *every* polygon enough pixels for the whole result,
            # which is potentially excessive!
            min_x, max_x, min_y, max_y = wmts_desired_extent
            if wmts_projection == projection:
                max_pixel_span = min((max_x - min_x) / width,
                                     (max_y - min_y) / height)
            else:
                # X/Y orientation is arbitrary, so use a worst-case guess.
                max_pixel_span = (min(max_x - min_x, max_y - min_y) /
                                  max(width, height))

            # Fetch a suitable image and its actual extent.
            wmts_image, wmts_actual_extent = self._wmts_images(
                self.wmts, self.layer, matrix_set_name,
                extent=wmts_desired_extent,
                max_pixel_span=max_pixel_span)

            # Return each (image, extent) as a LocatedImage.
            if wmts_projection == projection:
                located_image = LocatedImage(wmts_image, wmts_actual_extent)
            else:
                # Reproject the image to the desired projection.
                located_image = _warped_located_image(
                    wmts_image,
                    wmts_projection, wmts_actual_extent,
                    output_projection=projection, output_extent=extent,
                    target_resolution=target_resolution)

            located_images.append(located_image)

        return located_images

    def _choose_matrix(self, tile_matrices, meters_per_unit, max_pixel_span):
        # Get the tile matrices in order of increasing resolution.
        tile_matrices = sorted(tile_matrices,
                               key=lambda tm: tm.scaledenominator,
                               reverse=True)

        # Find which tile matrix has the appropriate resolution.
        max_scale = max_pixel_span * meters_per_unit / METERS_PER_PIXEL
        for tm in tile_matrices:
            if tm.scaledenominator <= max_scale:
                return tm
        return tile_matrices[-1]

    def _tile_span(self, tile_matrix, meters_per_unit):
        pixel_span = tile_matrix.scaledenominator * (
            METERS_PER_PIXEL / meters_per_unit)
        tile_span_x = tile_matrix.tilewidth * pixel_span
        tile_span_y = tile_matrix.tileheight * pixel_span
        return tile_span_x, tile_span_y

    def _select_tiles(self, tile_matrix, tile_matrix_limits,
                      tile_span_x, tile_span_y, extent):
        # Convert the requested extent from CRS coordinates to tile
        # indices. See annex H of the WMTS v1.0.0 spec.
        # NB. The epsilons get rid of any tiles which only just
        # (i.e. one part in a million) intrude into the requested
        # extent. Since these wouldn't be visible anyway there's nothing
        # to be gained by spending the time downloading them.
        min_x, max_x, min_y, max_y = extent
        matrix_min_x, matrix_max_y = tile_matrix.topleftcorner
        epsilon = 1e-6
        min_col = int((min_x - matrix_min_x) / tile_span_x + epsilon)
        max_col = int((max_x - matrix_min_x) / tile_span_x - epsilon)
        min_row = int((matrix_max_y - max_y) / tile_span_y + epsilon)
        max_row = int((matrix_max_y - min_y) / tile_span_y - epsilon)
        # Clamp to the limits of the tile matrix.
        min_col = max(min_col, 0)
        max_col = min(max_col, tile_matrix.matrixwidth - 1)
        min_row = max(min_row, 0)
        max_row = min(max_row, tile_matrix.matrixheight - 1)
        # Clamp to any layer-specific limits on the tile matrix.
        if tile_matrix_limits:
            min_col = max(min_col, tile_matrix_limits.mintilecol)
            max_col = min(max_col, tile_matrix_limits.maxtilecol)
            min_row = max(min_row, tile_matrix_limits.mintilerow)
            max_row = min(max_row, tile_matrix_limits.maxtilerow)
        return min_col, max_col, min_row, max_row

    def _wmts_images(self, wmts, layer, matrix_set_name, extent,
                     max_pixel_span):
        """
        Add images from the specified WMTS layer and matrix set to cover
        the specified extent at an appropriate resolution.

        The zoom level (aka. tile matrix) is chosen to give the lowest
        possible resolution which still provides the requested quality.
        If insufficient resolution is available, the highest available
        resolution is used.

        Args:

            * wmts - The owslib.wmts.WebMapTileService providing the tiles.
            * layer - The owslib.wmts.ContentMetadata (aka. layer) to draw.
            * matrix_set_name - The name of the matrix set to use.
            * extent - Tuple of (left, right, bottom, top) in Axes coordinates.
            * max_pixel_span - Preferred maximum pixel width or height
                               in Axes coordinates.

        """

        # Find which tile matrix has the appropriate resolution.
        tile_matrix_set = wmts.tilematrixsets[matrix_set_name]
        tile_matrices = tile_matrix_set.tilematrix.values()
        meters_per_unit = METERS_PER_UNIT[tile_matrix_set.crs]
        tile_matrix = self._choose_matrix(tile_matrices, meters_per_unit,
                                          max_pixel_span)

        # Determine which tiles are required to cover the requested extent.
        tile_span_x, tile_span_y = self._tile_span(tile_matrix,
                                                   meters_per_unit)
        tile_matrix_set_links = getattr(layer, 'tilematrixsetlinks', None)
        if tile_matrix_set_links is None:
            tile_matrix_limits = None
        else:
            tile_matrix_set_link = tile_matrix_set_links[matrix_set_name]
            tile_matrix_limits = tile_matrix_set_link.tilematrixlimits.get(
                tile_matrix.identifier)
        min_col, max_col, min_row, max_row = self._select_tiles(
            tile_matrix, tile_matrix_limits, tile_span_x, tile_span_y, extent)

        # Find the relevant section of the image cache.
        tile_matrix_id = tile_matrix.identifier
        cache_by_wmts = WMTSRasterSource._shared_image_cache
        cache_by_layer_matrix = cache_by_wmts.setdefault(wmts, {})
        image_cache = cache_by_layer_matrix.setdefault((layer.id,
                                                        tile_matrix_id), {})

        # To avoid nasty seams between the individual tiles, we
        # accumulate the tile images into a single image.
        big_img = None
        n_rows = 1 + max_row - min_row
        n_cols = 1 + max_col - min_col
        # Ignore out-of-range errors if the current version of OWSLib
        # doesn't provide the regional information.
        ignore_out_of_range = tile_matrix_set_links is None
        for row in range(min_row, max_row + 1):
            for col in range(min_col, max_col + 1):
                # Get the tile's Image from the cache if possible.
                img_key = (row, col)
                img = image_cache.get(img_key)
                if img is None:
                    try:
                        tile = wmts.gettile(
                            layer=layer.id,
                            tilematrixset=matrix_set_name,
                            tilematrix=tile_matrix_id,
                            row=row, column=col,
                            **self.gettile_extra_kwargs)
                    except owslib.util.ServiceException as exception:
                        if ('TileOutOfRange' in exception.message and
                                ignore_out_of_range):
                            continue
                        raise exception
                    img = Image.open(io.BytesIO(tile.read()))
                    image_cache[img_key] = img
                if big_img is None:
                    size = (img.size[0] * n_cols, img.size[1] * n_rows)
                    big_img = Image.new('RGBA', size, (255, 255, 255, 255))
                top = (row - min_row) * tile_matrix.tileheight
                left = (col - min_col) * tile_matrix.tilewidth
                big_img.paste(img, (left, top))

        if big_img is None:
            img_extent = None
        else:
            matrix_min_x, matrix_max_y = tile_matrix.topleftcorner
            min_img_x = matrix_min_x + tile_span_x * min_col
            max_img_y = matrix_max_y - tile_span_y * min_row
            img_extent = (min_img_x, min_img_x + n_cols * tile_span_x,
                          max_img_y - n_rows * tile_span_y, max_img_y)
        return big_img, img_extent


class WFSGeometrySource(object):
    """Web Feature Service (WFS) retrieval for Cartopy."""

    def __init__(self, service, features, getfeature_extra_kwargs=None):
        """
        Args:

        * service:
            The URL of a WFS, or an instance of
            :class:`owslib.wfs.WebFeatureService`.
        * features:
            The typename(s) of the features from the WFS that
            will be retrieved and made available as geometries.

        Kwargs:

        * getfeature_extra_kwargs:
            Extra keyword args to pass to the service's `getfeature` call.

        """
        if WebFeatureService is None:
            raise ImportError(_OWSLIB_REQUIRED)

        if isinstance(service, six.string_types):
            service = WebFeatureService(service)

        if isinstance(features, six.string_types):
            features = [features]

        if getfeature_extra_kwargs is None:
            getfeature_extra_kwargs = {}

        if not features:
            raise ValueError('One or more features must be specified.')
        for feature in features:
            if feature not in service.contents:
                raise ValueError('The {!r} feature does not exist in this '
                                 'service.'.format(feature))

        self.service = service
        self.features = features
        self.getfeature_extra_kwargs = getfeature_extra_kwargs

        self._default_urn = None

    def default_projection(self):
        """
        Return a :class:`cartopy.crs.Projection` in which the WFS
        service can supply the requested features.

        """
        # Using first element in crsOptions (default).
        if self._default_urn is None:
            default_urn = set(self.service.contents[feature].crsOptions[0] for
                              feature in self.features)
            if len(default_urn) != 1:
                ValueError('Failed to find a single common default SRS '
                           'across all features (typenames).')
            else:
                default_urn = default_urn.pop()
                default_srs = default_urn.id

            if six.text_type(default_urn) not in _URN_TO_CRS:
                raise ValueError('Unknown mapping from SRS/CRS_URN {!r} to '
                                 'cartopy projection.'.format(default_urn))

            self._default_urn = default_urn

        return _URN_TO_CRS[six.text_type(self._default_urn)]

    def fetch_geometries(self, projection, extent):
        """
        Return any Point, Linestring or LinearRing geometries available
        from the WFS that lie within the specified extent.

        Args:

        * projection: :class:`cartopy.crs.Projection`
            The projection in which the extent is specified and in
            which the geometries should be returned. Only the default
            (native) projection is supported.

        * extent: four element tuple
            (min_x, max_x, min_y, max_y) tuple defining the geographic extent
            of the geometries to obtain.

        Returns:
            A list of Shapely geometries.

        """
        if self.default_projection() != projection:
            raise ValueError('Geometries are only available in projection '
                             '{!r}.'.format(self.default_projection()))

        min_x, max_x, min_y, max_y = extent
        response = self.service.getfeature(typename=self.features,
                                           bbox=(min_x, min_y, max_x, max_y),
                                           **self.getfeature_extra_kwargs)
        geoms_by_srs = self._to_shapely_geoms(response)
        if not geoms_by_srs:
            geoms = []
        elif len(geoms_by_srs) > 1:
            raise ValueError('Unexpected response from the WFS server. The '
                             'geometries are in multiple SRSs, when only one '
                             'was expected.')
        else:
            srs, geoms = list(geoms_by_srs.items())[0]
            # Attempt to verify the SRS associated with the geometries (if any)
            # matches the specified projection.
            if srs is not None:
                if srs in _URN_TO_CRS:
                    geom_proj = _URN_TO_CRS[srs]
                    if geom_proj != projection:
                        raise ValueError('The geometries are not in expected '
                                         'projection. Expected {!r}, got '
                                         '{!r}.'.format(projection, geom_proj))
                else:
                    msg = 'Unable to verify matching projections due ' \
                          'to incomplete mappings from SRS identifiers ' \
                          'to cartopy projections. The geometries have ' \
                          'an SRS of {!r}.'.format(srs)
                    warnings.warn(msg)
        return geoms

    def _to_shapely_geoms(self, response):
        """
        Convert polygon coordinate strings in WFS response XML to Shapely
        geometries.

        Args:

        * response: (file-like object)
            WFS response XML data.

        Returns:
            A dictionary containing geometries, with key-value pairs of
            the form {srsname: [geoms]}.

        """
        linear_rings_data = []
        linestrings_data = []
        points_data = []
        tree = ElementTree.parse(response)

        for node in tree.findall('.//{}msGeometry'.format(_MAP_SERVER_NS)):
            # Find LinearRing geometries in our msGeometry node.
            find_str = './/{gml}LinearRing'.format(gml=_GML_NS)
            if self._node_has_child(node, find_str):
                data = self._find_polygon_coords(node, find_str)
                linear_rings_data.extend(data)

            # Find LineString geometries in our msGeometry node.
            find_str = './/{gml}LineString'.format(gml=_GML_NS)
            if self._node_has_child(node, find_str):
                data = self._find_polygon_coords(node, find_str)
                linestrings_data.extend(data)

            # Find Point geometries in our msGeometry node.
            find_str = './/{gml}Point'.format(gml=_GML_NS)
            if self._node_has_child(node, find_str):
                data = self._find_polygon_coords(node, find_str)
                points_data.extend(data)

        geoms_by_srs = {}
        for srs, x, y in linear_rings_data:
            geoms_by_srs.setdefault(srs, []).append(
                sgeom.LinearRing(zip(x, y)))
        for srs, x, y in linestrings_data:
            geoms_by_srs.setdefault(srs, []).append(
                sgeom.LineString(zip(x, y)))
        for srs, x, y in points_data:
            geoms_by_srs.setdefault(srs, []).append(
                sgeom.Point(zip(x, y)))
        return geoms_by_srs

    def _find_polygon_coords(self, node, find_str):
        """
        Return the x, y coordinate values for all the geometries in
        a given`node`.

        Args:

        * node: :class:`xml.etree.ElementTree.Element`
            Node of the parsed XML response.

        * find_str: string
            A search string used to match subelements that contain
            the coordinates of interest, for example:
            './/{http://www.opengis.net/gml}LineString'

        Returns:
            A list of (srsName, x_vals, y_vals) tuples.

        """

        data = []
        for polygon in node.findall(find_str):
            feature_srs = polygon.attrib.get('srsName')
            x, y = [], []

            # We can have nodes called `coordinates` or `coord`.
            coordinates_find_str = '{}coordinates'.format(_GML_NS)
            coords_find_str = '{}coord'.format(_GML_NS)

            if self._node_has_child(polygon, coordinates_find_str):
                points = polygon.findtext(coordinates_find_str)
                coords = points.strip().split(' ')
                for coord in coords:
                    x_val, y_val = coord.split(',')
                    x.append(float(x_val))
                    y.append(float(y_val))
            elif self._node_has_child(polygon, coords_find_str):
                for coord in polygon.findall(coords_find_str):
                    x.append(float(coord.findtext('{}X'.format(_GML_NS))))
                    y.append(float(coord.findtext('{}Y'.format(_GML_NS))))
            else:
                raise ValueError('Unable to find or parse coordinate values '
                                 'from the XML.')

            data.append((feature_srs, x, y))
        return data

    @staticmethod
    def _node_has_child(node, find_str):
        """
        Return whether `node` contains (at any sub-level), a node with name
        equal to `find_str`.

        """
        element = node.find(find_str)
        return element is not None
