# (C) British Crown Copyright 2011 - 2018, Met Office
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
Implements image tile identification and fetching from various sources.


The Matplotlib interface can make use of tile objects (defined below) via the
:meth:`cartopy.mpl.geoaxes.GeoAxes.add_image` method. For example, to add a
:class:`MapQuest Open Aerial tileset <MapQuestOpenAerial>` to an existing axes
at zoom level 2, do ``ax.add_image(MapQuestOpenAerial(), 2)``. An example of
using tiles in this way can be found at the
:ref:`sphx_glr_gallery_eyja_volcano.py` example.

"""

from __future__ import (absolute_import, division, print_function)

from abc import ABCMeta, abstractmethod
import warnings

from PIL import Image
import shapely.geometry as sgeom
import numpy as np
import six

import cartopy.crs as ccrs


class GoogleWTS(six.with_metaclass(ABCMeta, object)):
    """
    Implement web tile retrieval using the Google WTS coordinate system.

    A "tile" in this class refers to the coordinates (x, y, z).

    """
    def __init__(self, desired_tile_form='RGB'):
        self.imgs = []
        self.crs = ccrs.Mercator.GOOGLE
        self.desired_tile_form = desired_tile_form

    def image_for_domain(self, target_domain, target_z):
        tiles = []
        for tile in self.find_images(target_domain, target_z):
            try:
                img, extent, origin = self.get_image(tile)
            except IOError:
                continue
            img = np.array(img)
            x = np.linspace(extent[0], extent[1], img.shape[1])
            y = np.linspace(extent[2], extent[3], img.shape[0])
            tiles.append([img, x, y, origin])

        img, extent, origin = _merge_tiles(tiles)
        return img, extent, origin

    def _find_images(self, target_domain, target_z, start_tile=(0, 0, 0)):
        """Target domain is a shapely polygon in native coordinates."""

        assert isinstance(target_z, int) and target_z >= 0, ('target_z must '
                                                             'be an integer '
                                                             '>=0.')

        # Recursively drill down to the images at the target zoom.
        x0, x1, y0, y1 = self._tileextent(start_tile)
        domain = sgeom.box(x0, y0, x1, y1)
        if domain.intersects(target_domain):
            if start_tile[2] == target_z:
                    yield start_tile
            else:
                for tile in self._subtiles(start_tile):
                    for result in self._find_images(target_domain, target_z,
                                                    start_tile=tile):
                        yield result

    find_images = _find_images

    def subtiles(self, x_y_z):
        x, y, z = x_y_z
        # Google tile specific (i.e. up->down).
        for xi in range(0, 2):
            for yi in range(0, 2):
                yield x * 2 + xi, y * 2 + yi, z + 1

    _subtiles = subtiles

    def tile_bbox(self, x, y, z, y0_at_north_pole=True):
        """
        Return the ``(x0, x1), (y0, y1)`` bounding box for the given x, y, z
        tile position.

        Parameters
        ----------
        x
            The x tile coordinate in the Google tile numbering system.
        y
            The y tile coordinate in the Google tile numbering system.
        z
            The z tile coordinate in the Google tile numbering system.

        y0_at_north_pole: optional
            Boolean representing whether the numbering of the y coordinate
            starts at the north pole (as is the convention for Google tiles)
            or not (in which case it will start at the south pole, as is the
            convention for TMS). Defaults to True.


        """
        n = 2 ** z
        assert 0 <= x <= (n - 1), ("Tile's x index is out of range. Upper "
                                   "limit %s. Got %s" % (n, x))
        assert 0 <= y <= (n - 1), ("Tile's y index is out of range. Upper "
                                   "limit %s. Got %s" % (n, y))

        x0, x1 = self.crs.x_limits
        y0, y1 = self.crs.y_limits

        # Compute the box height and width in native coordinates
        # for this zoom level.
        box_h = (y1 - y0) / n
        box_w = (x1 - x0) / n

        # Compute the native x & y extents of the tile.
        n_xs = x0 + (x + np.arange(0, 2, dtype=np.float64)) * box_w
        n_ys = y0 + (y + np.arange(0, 2, dtype=np.float64)) * box_h

        if y0_at_north_pole:
            n_ys = -1 * n_ys[::-1]

        return n_xs, n_ys

    def tileextent(self, x_y_z):
        """Return extent tuple ``(x0,x1,y0,y1)`` in Mercator coordinates."""
        x, y, z = x_y_z
        x_lim, y_lim = self.tile_bbox(x, y, z, y0_at_north_pole=True)
        return tuple(x_lim) + tuple(y_lim)

    _tileextent = tileextent

    @abstractmethod
    def _image_url(self, tile):
        pass

    def get_image(self, tile):
        if six.PY3:
            from urllib.request import urlopen
        else:
            from urllib2 import urlopen

        url = self._image_url(tile)

        fh = urlopen(url)
        im_data = six.BytesIO(fh.read())
        fh.close()
        img = Image.open(im_data)

        img = img.convert(self.desired_tile_form)

        return img, self.tileextent(tile), 'lower'


class GoogleTiles(GoogleWTS):
    def __init__(self, desired_tile_form='RGB', style="street",
                 url=('https://mts0.google.com/vt/lyrs={style}'
                      '@177000000&hl=en&src=api&x={x}&y={y}&z={z}&s=G')):
        """
        Parameters
        ----------
        desired_tile_form: optional
            Defaults to 'RGB'.
        style: optional
            The style for the Google Maps tiles.  One of 'street',
            'satellite', 'terrain', and 'only_streets'.  Defaults to 'street'.
        url: optional
            URL pointing to a tile source and containing {x}, {y}, and {z}.
            Such as: ``'https://server.arcgisonline.com/ArcGIS/rest/services/\
World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg'``

        """
        styles = ["street", "satellite", "terrain", "only_streets"]
        style = style.lower()
        self.url = url
        if style not in styles:
            msg = "Invalid style '%s'. Valid styles: %s" % \
                (style, ", ".join(styles))
            raise ValueError(msg)
        self.style = style

        # The 'satellite' and 'terrain' styles require pillow with a jpeg
        # decoder.
        if self.style in ["satellite", "terrain"] and \
                not hasattr(Image.core, "jpeg_decoder") or \
                not Image.core.jpeg_decoder:
            msg = "The '%s' style requires pillow with jpeg decoding support."
            raise ValueError(msg % self.style)
        return super(GoogleTiles, self).__init__(
            desired_tile_form=desired_tile_form)

    def _image_url(self, tile):
        style_dict = {
            "street": "m",
            "satellite": "s",
            "terrain": "t",
            "only_streets": "h"}
        url = self.url.format(
            style=style_dict[self.style],
            x=tile[0], X=tile[0],
            y=tile[1], Y=tile[1],
            z=tile[2], Z=tile[2])
        return url


class MapQuestOSM(GoogleWTS):
    # http://developer.mapquest.com/web/products/open/map for terms of use
    # http://devblog.mapquest.com/2016/06/15/
    # modernization-of-mapquest-results-in-changes-to-open-tile-access/
    # this now requires a sign up to a plan
    def _image_url(self, tile):
        x, y, z = tile
        url = 'http://otile1.mqcdn.com/tiles/1.0.0/osm/%s/%s/%s.jpg' % (
            z, x, y)
        mqdevurl = ('http://devblog.mapquest.com/2016/06/15/'
                    'modernization-of-mapquest-results-in-changes'
                    '-to-open-tile-access/')
        warnings.warn('{} will require a log in and and will likely'
                      ' fail. see {} for more details.'.format(url, mqdevurl))
        return url


class MapQuestOpenAerial(GoogleWTS):
    # http://developer.mapquest.com/web/products/open/map for terms of use
    # The following attribution should be included in the resulting image:
    # "Portions Courtesy NASA/JPL-Caltech and U.S. Depart. of Agriculture,
    #  Farm Service Agency"
    def _image_url(self, tile):
        x, y, z = tile
        url = 'http://oatile1.mqcdn.com/tiles/1.0.0/sat/%s/%s/%s.jpg' % (
            z, x, y)
        return url


class OSM(GoogleWTS):
    # http://developer.mapquest.com/web/products/open/map for terms of use
    def _image_url(self, tile):
        x, y, z = tile
        url = 'https://a.tile.openstreetmap.org/%s/%s/%s.png' % (z, x, y)
        return url


class Stamen(GoogleWTS):
    """
    Retrieves tiles from maps.stamen.com. Styles include
    ``terrain-background``, ``terrain``, ``toner`` and ``watercolor``.

    For a full reference on the styles available please see
    http://maps.stamen.com. Of particular note are the sub-styles
    that are made available (e.g. ``terrain`` and ``terrain-background``).
    To determine the name of the particular [sub-]style you want,
    follow the link on http://maps.stamen.com to your desired style and
    observe the style name in the URL. Your style name will be in the
    form of: ``http://maps.stamen.com/{STYLE_NAME}/#9/37/-122``.

    Except otherwise noted, the Stamen map tile sets are copyright Stamen
    Design, under a Creative Commons Attribution (CC BY 3.0) license.

    Please see the attribution notice at http://maps.stamen.com on how to
    attribute this imagery.

    """
    def __init__(self, style='toner'):
        super(Stamen, self).__init__()
        self.style = style

    def _image_url(self, tile):
        return ('http://tile.stamen.com/{self.style}/{z}/{x}/{y}.png'
                .format(self=self, x=tile[0], y=tile[1], z=tile[2]))


class StamenTerrain(Stamen):
    """
    **DEPRECATED:** This class is deprecated. Please use
    ``Stamen('terrain-background')`` instead.

    Terrain tiles defined for the continental United States, and include land
    color and shaded hills. The land colors are a custom palette developed by
    Gem Spear for the National Atlas 1km land cover data set, which defines
    twenty-four land classifications including five kinds of forest,
    combinations of shrubs, grasses and crops, and a few tundras and wetlands.
    The colors are at their highest contrast when fully zoomed-out to the
    whole U.S., and they slowly fade out to pale off-white as you zoom in to
    leave room for foreground data and break up the weirdness of large areas
    of flat, dark green.

    References
    ----------

     * http://mike.teczno.com/notes/osm-us-terrain-layer/background.html
     * http://maps.stamen.com/
     * https://wiki.openstreetmap.org/wiki/List_of_OSM_based_Services
     * https://github.com/migurski/DEM-Tools


    """
    def __init__(self):
        warnings.warn(
            "The StamenTerrain class was deprecated in v0.17. "
            "Please use Stamen('terrain-background') instead.")

        # NOTE: This subclass of Stamen exists for legacy reasons.
        # No further Stamen subclasses will be accepted as
        # they can easily be created in user code with Stamen(style_name).
        return super(StamenTerrain, self).__init__(style='terrain-background')


class MapboxTiles(GoogleWTS):
    """
    Implement web tile retrieval from Mapbox.

    For terms of service, see https://www.mapbox.com/tos/.

    """
    def __init__(self, access_token, map_id):
        """
        Set up a new Mapbox tiles instance.

        Access to Mapbox web services requires an access token and a map ID.
        See https://www.mapbox.com/api-documentation/ for details.

        Parameters
        ----------
        access_token
            A valid Mapbox API access token.
        map_id
            An ID for a publicly accessible map (provided by Mapbox).
            This is the map whose tiles will be retrieved through this process.

        """
        self.access_token = access_token
        self.map_id = map_id
        super(MapboxTiles, self).__init__()

    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://api.mapbox.com/v4/mapbox.{id}/{z}/{x}/{y}.png'
               '?access_token={token}'.format(z=z, y=y, x=x,
                                              id=self.map_id,
                                              token=self.access_token))
        return url


class MapboxStyleTiles(GoogleWTS):
    """
    Implement web tile retrieval from a user-defined Mapbox style. For more
    details on Mapbox styles, see
    https://www.mapbox.com/studio-manual/overview/map-styling/.

    For terms of service, see https://www.mapbox.com/tos/.

    """
    def __init__(self, access_token, username, map_id):
        """
        Set up a new instance to retrieve tiles from a Mapbox style.

        Access to Mapbox web services requires an access token and a map ID.
        See https://www.mapbox.com/api-documentation/ for details.

        Parameters
        ----------
        access_token
            A valid Mapbox API access token.
        username
            The username for the Mapbox user who defined the Mapbox style.
        map_id
            A map ID for a map defined by a Mapbox style. This is the map whose
            tiles will be retrieved through this process. Note that this style
            may be private and if your access token does not have permissions
            to view this style, then map tile retrieval will fail.

        """
        self.access_token = access_token
        self.username = username
        self.map_id = map_id
        super(MapboxStyleTiles, self).__init__()

    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://api.mapbox.com/styles/v1/'
               '{user}/{mapid}/tiles/256/{z}/{x}/{y}'
               '?access_token={token}'.format(z=z, y=y, x=x,
                                              user=self.username,
                                              mapid=self.map_id,
                                              token=self.access_token))
        return url


class QuadtreeTiles(GoogleWTS):
    """
    Implement web tile retrieval using the Microsoft WTS quadkey coordinate
    system.

    A "tile" in this class refers to a quadkey such as "1", "14" or "141"
    where the length of the quatree is the zoom level in Google Tile terms.

    """
    def _image_url(self, tile):
        url = ('http://ecn.dynamic.t1.tiles.virtualearth.net/comp/'
               'CompositionHandler/{tile}?mkt=en-'
               'gb&it=A,G,L&shading=hill&n=z'.format(tile=tile))
        return url

    def tms_to_quadkey(self, tms, google=False):
        quadKey = ""
        x, y, z = tms
        # this algorithm works with google tiles, rather than tms, so convert
        # to those first.
        if not google:
            y = (2 ** z - 1) - y
        for i in range(z, 0, -1):
            digit = 0
            mask = 1 << (i - 1)
            if (x & mask) != 0:
                digit += 1
            if (y & mask) != 0:
                digit += 2
            quadKey += str(digit)
        return quadKey

    def quadkey_to_tms(self, quadkey, google=False):
        # algorithm ported from
        # https://msdn.microsoft.com/en-us/library/bb259689.aspx
        assert isinstance(quadkey, six.string_types), \
            'quadkey must be a string'

        x = y = 0
        z = len(quadkey)
        for i in range(z, 0, -1):
            mask = 1 << (i - 1)
            if quadkey[z - i] == '0':
                pass
            elif quadkey[z - i] == '1':
                x |= mask
            elif quadkey[z - i] == '2':
                y |= mask
            elif quadkey[z - i] == '3':
                x |= mask
                y |= mask
            else:
                raise ValueError('Invalid QuadKey digit '
                                 'sequence.' + str(quadkey))
        # the algorithm works to google tiles, so convert to tms
        if not google:
            y = (2 ** z - 1) - y
        return (x, y, z)

    def subtiles(self, quadkey):
        for i in range(4):
            yield quadkey + str(i)

    def tileextent(self, quadkey):
        x_y_z = self.quadkey_to_tms(quadkey, google=True)
        return GoogleWTS.tileextent(self, x_y_z)

    def find_images(self, target_domain, target_z, start_tile=None):
        """
        Find all the quadtrees at the given target zoom, in the given
        target domain.

        target_z must be a value >= 1.

        """
        if target_z == 0:
            raise ValueError('The empty quadtree cannot be returned.')

        if start_tile is None:
            start_tiles = ['0', '1', '2', '3']
        else:
            start_tiles = [start_tile]

        for start_tile in start_tiles:
            start_tile = self.quadkey_to_tms(start_tile, google=True)
            for tile in GoogleWTS.find_images(self, target_domain, target_z,
                                              start_tile=start_tile):
                yield self.tms_to_quadkey(tile, google=True)


def _merge_tiles(tiles):
    """Return a single image, merging the given images."""
    if not tiles:
        raise ValueError('A non-empty list of tiles should '
                         'be provided to merge.')
    xset = [set(x) for i, x, y, _ in tiles]
    yset = [set(y) for i, x, y, _ in tiles]

    xs = xset[0]
    xs.update(*xset[1:])
    ys = yset[0]
    ys.update(*yset[1:])
    xs = sorted(xs)
    ys = sorted(ys)

    other_len = tiles[0][0].shape[2:]
    img = np.zeros((len(ys), len(xs)) + other_len, dtype=np.uint8) - 1

    for tile_img, x, y, origin in tiles:
        y_first, y_last = y[0], y[-1]
        yi0, yi1 = np.where((y_first == ys) | (y_last == ys))[0]
        if origin == 'upper':
            yi0 = tile_img.shape[0] - yi0 - 1
            yi1 = tile_img.shape[0] - yi1 - 1
        start, stop, step = yi0, yi1, 1 if yi0 < yi1 else -1
        if step == 1 and stop == img.shape[0] - 1:
            stop = None
        elif step == -1 and stop == 0:
            stop = None
        else:
            stop += step
        y_slice = slice(start, stop, step)

        xi0, xi1 = np.where((x[0] == xs) | (x[-1] == xs))[0]

        start, stop, step = xi0, xi1, 1 if xi0 < xi1 else -1

        if step == 1 and stop == img.shape[1] - 1:
            stop = None
        elif step == -1 and stop == 0:
            stop = None
        else:
            stop += step

        x_slice = slice(start, stop, step)

        img_slice = (y_slice, x_slice, Ellipsis)

        if origin == 'lower':
            tile_img = tile_img[::-1, ::]

        img[img_slice] = tile_img

    return img, [min(xs), max(xs), min(ys), max(ys)], 'lower'
