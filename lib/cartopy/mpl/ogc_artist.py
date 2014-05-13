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
This module defines the :class:`WMTSArtist` class, for drawing
WMTS layers with matplotlib.

"""
import io
import math
import weakref

from matplotlib.image import AxesImage
import matplotlib.artist

import cartopy.crs as ccrs


# Standard pixel size of 0.28 mm as defined by WMTS.
METERS_PER_PIXEL = 0.28e-3

_WGS84_METERS_PER_UNIT = 2 * math.pi * 6378137 / 360

METERS_PER_UNIT = {
    'urn:ogc:def:crs:OGC:1.3:CRS84': _WGS84_METERS_PER_UNIT,
}

_URN_TO_CRS = {
    'urn:ogc:def:crs:OGC:1.3:CRS84': ccrs.PlateCarree(),
}


class WMTSArtist(matplotlib.artist.Artist):
    """
    A subclass of :class:`~matplotlib.artist.Artist` capable of
    drawing a WMTS layer.

    Requires owslib and PIL to work.

    """

    _shared_image_cache = weakref.WeakKeyDictionary()
    """
    A nested mapping from WMTS, layer name, tile matrix name, tile row
    and tile column to the resulting PIL image::

        {wmts: {(layer_name, tile_matrix_name): {(row, column): Image}}}

    This provides a significant boost when producing multiple maps of the
    same projection or with an interactive figure.

    """

    def __init__(self, wmts, layer_name, matrix_set_name=None):
        """
        Args:

            * wmts - The URL of the WMTS, or an
                     owslib.wmts.WebMapTileService instance.
            * layer_name - The name of the layer to use.

        Kwargs:

            * matrix_set_name - Optional matrix set name. Defaults to
                                a matrix set which matches the current
                                projection parameters.

        """
        super(WMTSArtist, self).__init__()

        if not (hasattr(wmts, 'tilematrixsets') and
                hasattr(wmts, 'contents') and
                hasattr(wmts, 'gettile')):
            try:
                import owslib.wmts
            except ImportError:
                raise ImportError('OWSLib is required for showing WMTS layers')
            wmts = owslib.wmts.WebMapTileService(wmts)

        try:
            import PIL.Image
        except ImportError:
            raise ImportError('PIL is required for showing WMTS layers')

        self._wmts = wmts
        self._layer_name = layer_name
        self._matrix_set_name = matrix_set_name
        self._pil_image = PIL.Image

    def set_axes(self, axes):
        """
        Set the :class:`~matplotlib.axes.Axes` instance in which the
        artist resides, if any.

        ACCEPTS: an :class:`~matplotlib.axes.Axes` instance

        """
        super(WMTSArtist, self).set_axes(axes)

        if self._matrix_set_name is None:
            wmts = self._wmts
            layer = wmts.contents[self._layer_name]
            for tile_matrix_set_name in layer.tilematrixsets:
                tile_matrix_set = wmts.tilematrixsets[tile_matrix_set_name]
                crs_urn = tile_matrix_set.crs
                if crs_urn in _URN_TO_CRS:
                    tms_crs = _URN_TO_CRS[crs_urn]
                    if tms_crs == self.axes.projection:
                        self._matrix_set_name = tile_matrix_set_name
                        break
            if self._matrix_set_name is None:
                available_urns = sorted(set(
                    wmts.tilematrixsets[name].crs for name in
                    layer.tilematrixsets))
                msg = 'Unable to find tile matrix for current projection.'
                msg += '\n    Projection: ' + str(axes.projection)
                msg += '\n    Available tile CRS URNs:'
                msg += '\n        ' + '\n        '.join(available_urns)
                raise ValueError(msg)

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer, *args, **kwargs):
        """
        Draws the tiles from the WMTS layer that intersect with the
        extent of the :class:`cartopy.mpl.GeoAxes` instance to which
        this object has been added.

        """
        if not self.get_visible():
            return

        ax = self.get_axes()
        window_extent = ax.get_window_extent()
        max_pixel_span = min(ax.viewLim.width / window_extent.width,
                             ax.viewLim.height / window_extent.height)
        self._add_wmts_images(self._wmts, self._layer_name,
                              self._wmts.tilematrixsets[self._matrix_set_name],
                              ax.get_extent(), max_pixel_span, renderer)

    def _add_wmts_images(self, wmts, layer_name, tile_matrix_set, extent,
                         max_pixel_span, renderer):
        """
        Add images from the specified WMTS layer and matrix set to cover
        the specified extent at an appropriate resolution.

        The zoom level (aka. tile matrix) is chosen to give the lowest
        possible resolution which still provides the requested quality.
        If insufficient resolution is available, the highest available
        resolution is used.

        Args:

            * wmts - The owslib.wmts.WebMapTileService providing the tiles.
            * layer_name - The name of the layer to use.
            * tile_matrix_set - A owslib.wmts.TileMatrixSet relevant to
                                the layer.
            * extent - Tuple of (left, right, bottom, top) in Axes coordinates.
            * max_pixel_span - Preferred maximum pixel width or height
                               in Axes coordinates.

        """
        min_x, max_x, min_y, max_y = extent

        # Get the tile matrices in order of increasing resolution.
        tile_matrices = sorted(tile_matrix_set.tilematrix.values(),
                               key=lambda tm: tm.scaledenominator,
                               reverse=True)

        # Find which tile matrix has the appropriate resolution.
        meters_per_unit = METERS_PER_UNIT[tile_matrix_set.crs]
        max_scale = max_pixel_span * meters_per_unit / METERS_PER_PIXEL
        ok_tile_matrices = filter(lambda tm: tm.scaledenominator <= max_scale,
                                  tile_matrices)
        if ok_tile_matrices:
            tile_matrix = ok_tile_matrices[0]
        else:
            tile_matrix = tile_matrices[-1]

        # Determine which tiles are required to cover the requested extent.
        pixel_span = tile_matrix.scaledenominator * (
            METERS_PER_PIXEL / meters_per_unit)
        tile_span_x = tile_matrix.tilewidth * pixel_span
        tile_span_y = tile_matrix.tileheight * pixel_span

        # Convert the requested extent from CRS coordinates to tile
        # indices. See annex H of the WMTS v1.0.0 spec.
        # NB. The epsilons get rid of any tiles which only just
        # (i.e. one part in a million) intrude into the requested
        # extent. Since these wouldn't be visible anyway there's nothing
        # to be gained by spending the time downloading them.
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

        tile_matrix_id = tile_matrix.identifier
        image_cache = self._shared_image_cache.setdefault(wmts, {}) \
            .setdefault((layer_name, tile_matrix_id), {})

        # To avoid nasty seams between the individual tiles, we
        # accumulate the tile images into a single image.
        big_img = None
        n_rows = 1 + max_row - min_row
        n_cols = 1 + max_col - min_col
        for row in range(min_row, max_row + 1):
            for col in range(min_col, max_col + 1):
                # Get the tile's Image from the cache if possible.
                img_key = (row, col)
                img = image_cache.get(img_key)
                if img is None:
                    tile = wmts.gettile(layer=layer_name,
                                        tilematrix=tile_matrix_id,
                                        row=row, column=col)
                    img = self._pil_image.open(io.BytesIO(tile.read()))
                    image_cache[img_key] = img
                if big_img is None:
                    size = (img.size[0] * n_cols, img.size[1] * n_rows)
                    big_img = self._pil_image.new('RGB', size, None)
                top = (row - min_row) * tile_matrix.tileheight
                left = (col - min_col) * tile_matrix.tilewidth
                big_img.paste(img, (left, top))

        # Draw the combined image.
        ax = self.get_axes()
        min_img_x = matrix_min_x + tile_span_x * min_col
        max_img_y = matrix_max_y - tile_span_y * min_row
        img_extent = (min_img_x, min_img_x + n_cols * tile_span_x,
                      max_img_y - n_rows * tile_span_y, max_img_y)
        img_artist = AxesImage(ax, extent=img_extent, origin='upper')
        img_artist.set_data(big_img)
        img_artist.set_clip_path(ax.outline_patch)
        img_artist.draw(renderer)
