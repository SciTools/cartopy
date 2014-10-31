# (C) British Crown Copyright 2014, Met Office
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

import io
import math
import os
import weakref
from xml.etree import ElementTree

from PIL import Image
from shapely.geometry import LinearRing, LineString, Point

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

from cartopy.io import RasterSource
import cartopy.crs as ccrs


_OWSLIB_REQUIRED = 'OWSLib is required to use OGC web services.'

# Hardcode some known EPSG codes for now.
_CRS_TO_OGC_SRS = {ccrs.PlateCarree(): 'EPSG:4326'
                   }

# Standard pixel size of 0.28 mm as defined by WMTS.
METERS_PER_PIXEL = 0.28e-3

_WGS84_METERS_PER_UNIT = 2 * math.pi * 6378137 / 360

METERS_PER_UNIT = {
    'urn:ogc:def:crs:EPSG::900913': 1,
    'urn:ogc:def:crs:OGC:1.3:CRS84': _WGS84_METERS_PER_UNIT,
}

_URN_TO_CRS = {
    'urn:ogc:def:crs:EPSG::4326': ccrs.PlateCarree(),
    'urn:ogc:def:crs:EPSG::900913': ccrs.GOOGLE_MERCATOR,
    'urn:ogc:def:crs:EPSG::32661': ccrs.Mercator(),
    'urn:ogc:def:crs:OGC:1.3:CRS84': ccrs.PlateCarree(),
}


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

        if isinstance(service, basestring):
            service = WebMapService(service)

        if isinstance(layers, basestring):
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

        self._srs_for_projection_id = {}

    def _srs(self, projection):
        key = id(projection)
        srs = self._srs_for_projection_id.get(key)
        if srs is None:
            srs = _CRS_TO_OGC_SRS.get(projection)
            if srs is None:
                raise ValueError('The projection {!r} was not convertible to '
                                 'a suitable WMS SRS.'.format(projection))
            for layer in self.layers:
                if srs not in self.service.contents[layer].crsOptions:
                    raise ValueError('The SRS {} is not a valid SRS for the '
                                     '{!r} WMS layer.'.format(srs, layer))
            self._srs_for_projection_id[key] = srs
        return srs

    def validate_projection(self, projection):
        self._srs(projection)

    def fetch_raster(self, projection, extent, target_resolution):
        service = self.service
        min_x, max_x, min_y, max_y = extent
        wms_image = service.getmap(layers=self.layers,
                                   srs=self._srs(projection),
                                   bbox=(min_x, min_y, max_x, max_y),
                                   size=target_resolution, format='image/png',
                                   **self.getmap_extra_kwargs)
        wms_image = Image.open(io.BytesIO(wms_image.read()))
        return wms_image, extent


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

    def __init__(self, wmts, layer_name):
        """
        Args:

            * wmts - The URL of the WMTS, or an
                     owslib.wmts.WebMapTileService instance.
            * layer_name - The name of the layer to use.

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

        self._matrix_set_name_map = {}

    def _matrix_set_name(self, projection):
        key = id(projection)
        matrix_set_name = self._matrix_set_name_map.get(key)
        if matrix_set_name is None:
            wmts = self.wmts
            if hasattr(self.layer, 'tilematrixsetlinks'):
                matrix_set_names = self.layer.tilematrixsetlinks.keys()
            else:
                matrix_set_names = self.layer.tilematrixsets
            for tile_matrix_set_name in matrix_set_names:
                tile_matrix_set = wmts.tilematrixsets[tile_matrix_set_name]
                crs_urn = tile_matrix_set.crs
                if crs_urn in _URN_TO_CRS:
                    tms_crs = _URN_TO_CRS[crs_urn]
                    if tms_crs == projection:
                        matrix_set_name = tile_matrix_set_name
                        break
            if matrix_set_name is None:
                available_urns = sorted(set(
                    wmts.tilematrixsets[name].crs for name in
                    matrix_set_names))
                msg = 'Unable to find tile matrix for projection.'
                msg += '\n    Projection: ' + str(projection)
                msg += '\n    Available tile CRS URNs:'
                msg += '\n        ' + '\n        '.join(available_urns)
                raise ValueError(msg)
            self._matrix_set_name_map[key] = matrix_set_name
        return matrix_set_name

    def validate_projection(self, projection):
        self._matrix_set_name(projection)

    def fetch_raster(self, projection, extent, target_resolution):
        matrix_set_name = self._matrix_set_name(projection)
        min_x, max_x, min_y, max_y = extent
        width, height = target_resolution
        max_pixel_span = min((max_x - min_x) / width,
                             (max_y - min_y) / height)
        image, extent = self._wmts_images(self.wmts, self.layer,
                                          matrix_set_name, extent,
                                          max_pixel_span)
        return image, extent

    def _choose_matrix(self, tile_matrices, meters_per_unit, max_pixel_span):
        # Get the tile matrices in order of increasing resolution.
        tile_matrices = sorted(tile_matrices,
                               key=lambda tm: tm.scaledenominator,
                               reverse=True)

        # Find which tile matrix has the appropriate resolution.
        max_scale = max_pixel_span * meters_per_unit / METERS_PER_PIXEL
        ok_tile_matrices = filter(lambda tm: tm.scaledenominator <= max_scale,
                                  tile_matrices)
        if ok_tile_matrices:
            tile_matrix = ok_tile_matrices[0]
        else:
            tile_matrix = tile_matrices[-1]
        return tile_matrix

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
                            row=row, column=col)
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


class WFSRasterSource(RasterSource):
    """Web Feature Service (WFS) retrieval for Cartopy."""

    def __init__(self, service, features, getfeature_extra_kwargs=None):
        """
        Args:

        * service:
            The URL of a WFS, or an :class:`owslib.wfs.WebFeatureService`
            instance.
        * features:
            The name(s) of the features from the WFS to use.

        Kwargs:

        * getfeature_extra_kwargs:
            Extra keyword args to pass to the service's `getfeature` call.

        """
        if WebFeatureService is None:
            raise ImportError(_OWSLIB_REQUIRED)

        if isinstance(service, basestring):
            service = WebFeatureService(service)

        if isinstance(features, basestring):
            features = [features]

        # Populate an empty kwargs dict with something harmless.
        if getfeature_extra_kwargs is None:
            getfeature_extra_kwargs = {'propertyname': ['*']}

        if len(features) == 0:
            raise ValueError('One or more features must be specified.')
        for feature in features:
            if feature not in service.contents:
                raise ValueError('The {!r} feature does not exist in this '
                                 'service.'.format(feature))

        self.service = service
        self.features = features
        self.getfeature_extra_kwargs = getfeature_extra_kwargs

        # For parsing GML. The enclosing {} are required as part of this.
        self.ms = "{http://mapserver.gis.umn.edu/mapserver}"
        self.gml = "{http://www.opengis.net/gml}"

        self._crs_urn_for_projection_id = {}

    def _crs_urn(self, projection):
        """
        Confirm that the SRS of the features in the WFS response is recognised
        and that it lines up with the requested `projection`.

        """
        key = id(projection)
        crs_urn = self._crs_urn_for_projection_id.get(key)
        if crs_urn is None:
            # Confirm the CRS URNs are recognised and match `projection`.
            for feature in self.features:
                crs_options = self.service.contents[feature].crsOptions
                for features_crs_urn in crs_options:
                    if features_crs_urn in _URN_TO_CRS:
                        this_crs = _URN_TO_CRS[features_crs_urn]
                        if this_crs == projection:
                            crs_urn = features_crs_urn
                            break
            # Fail informatively if not.
            if crs_urn is None:
                raise ValueError('The projection {!r} was not convertible to '
                                 'a suitable WFS SRS.'.format(projection))
            self._crs_urn_for_projection_id[key] = crs_urn
        return crs_urn

    def validate_projection(self, projection):
        self._crs_urn(projection)

    def fetch_raster(self, projection, extent):
        """
        Return Shapely geometries of any Point, Linestring or LinearRing
        geometries found in the WFS request.

        .. note::
            Kwargs passed to the get feature request are not guaranteed to
            have an impact on the response data as the WFS server will not
            necessarily be able to supply data that precisely matches the
            request. In such cases the response will contain the default
            data provided by the WFS server.

        """
        service = self.service
        min_x, max_x, min_y, max_y = extent
        response = service.getfeature(typename=self.features,
                                      bbox=(min_x, min_y, max_x, max_y),
                                      **self.getfeature_extra_kwargs)
        wfs_features = self._to_shapely_geoms(response)
        return wfs_features, extent

    def _to_shapely_geoms(self, response):
        """
        Convert polygon coordinate strings in WFS response XML to Shapely
        geometries.

        Args:

        * response:
            WFS response XML data.

        """
        linear_rings_data = []
        linestrings_data = []
        points_data = []
        tree = ElementTree.parse(response)

        for node in tree.findall('.//{ms}msGeometry'.format(ms=self.ms)):
            # Find LinearRing geometries in our Polygon node.
            find_str = './/{gml}LinearRing'.format(gml=self.gml)
            if self._node_has_child(node, find_str):
                data = self._find_polygon_coords(node, find_str)
                linear_rings_data.extend(data)

            # Find LineString geometries in our Polygon node.
            find_str = './/{gml}LineString'.format(gml=self.gml)
            if _node_has_child(node, find_str):
                data = self._find_polygon_coords(node, find_str)
                linestrings_data.extend(data)

            # Find Point geometries in our Polygon node.
            find_str = './/{gml}Point'.format(gml=self.gml)
            if _node_has_child(node, find_str):
                data = self._find_polygon_coords(node, find_str)
                points_data.extend(data)

        linear_rings = {}
        for k, x, y in linear_rings_data:
            linear_rings.setdefault(k, []).append(LinearRing(zip(x, y)))
        linestrings = {}
        for k, x, y in linestrings_data:
            linestrings.setdefault(k, []).append(LineString(zip(x, y)))
        points = {}
        for k, x, y in points_data:
            points.setdefault(k, []).append(Point(zip(x, y)))
        return linear_rings, linestrings, points

    def _find_polygon_coords(self, node, find_str):
        """Find all coordinates data for a given Polygon `node`."""

        data = []
        for polygon in node.findall(find_str):
            feature_srs = polygon.attrib.get('srsName')
            # Assume we can either have nodes called `coordinates` or `coords`.
            find_str = '{}coordinates'.format(self.gml)
            if _node_has_child(polygon, find_str): 
                points = polygon.findtext(find_str)
                x, y = self._coordinates_to_points(points)
            else:
                x, y = [], []
                for coord in polygon.findall('{}coord'.format(self.gml)):
                    x.append(float(coord.findtext('{}X'.format(self.gml))))
                    y.append(float(coord.findtext('{}Y'.format(self.gml))))

            data.append([feature_srs, x, y])
        return data

    def _node_has_child(self, node, find_str):
        """
        Determine whether `node` contains (at any sub-level), a node with name
        equal to `find_str`.

        """
        found = list(node.iterfind(find_str))
        return len(found) > 0

    def _coordinates_to_points(self, points):
        """
        If our XML has a coordinates node then we'll get one string of
        whitespace-separated points:

            "x,y x,y ... x,y"

        Parses this string to return lists of x and y points.

        """
        x = []
        y = []
        coords = points.strip().split(' ')
        for coord in coords:
            x_val, y_val = coord.split(',')
            x.append(float(x_val))
            y.append(float(y_val))
        return x, y
