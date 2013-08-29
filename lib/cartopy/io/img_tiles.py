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


"""
Implements image tile identification and fetching from various sources.

Tile generation is explicitly not yet implemented.

"""

import math

import PIL.Image as Image
import shapely.geometry as sgeom
import numpy as np

import cartopy.crs as ccrs


class GoogleTiles(object):
    """
    Implements web tile retrieval using the Google WTS coordinate system.

    A "tile" in this class refers to the coordinates (x, y, z).

    """
    def __init__(self, desired_tile_form='RGB'):
        # XXX consider fixing the CRS???
        self.imgs = []
        self.crs = ccrs.Mercator(
            min_latitude=-85., max_latitude=85.,
            globe=ccrs.Globe(semimajor_axis=math.degrees(1)))
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

        # recursively drill down to the images at the target zoom
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
        # google tile specific (i.e. up->down)
        for xi in range(0, 2):
            for yi in range(0, 2):
                yield x * 2 + xi, y * 2 + yi, z + 1

    _subtiles = subtiles

    @staticmethod
    def tile_bbox(prj, x, y, z, lat_extent_at_z0=(-85., 85.),
                  bottom_up=True):
        """
        Returns the ``(x0, x1), (y0, y1)`` bounding box for the given x, y, z
        tile position.

        NOTE: This interface is highly liable to change in the future.
        """
        n = 2 ** z
        assert 0 <= x <= (n - 1), ("Tile's x index is out of range. Upper "
                                   "limit %s. Got %s" % (n, x))
        assert 0 <= y <= (n - 1), ("Tile's y index is out of range. Upper "
                                   "limit %s. Got %s" % (n, y))

        x0 = -180.
        # compute the box height in native coordinates for this zoom level
        box_w = 360. / n

        result = prj.transform_points(ccrs.PlateCarree(),
                                      np.array([0., 0]),
                                      np.array(lat_extent_at_z0))
        y1, y0 = result[:, 1]
        # compute the box height in native coordinates for this zoom level
        box_h = (y1 - y0) / n

        # compute the native x & native y extents of the box
    #    x_lower, x_upper = x0 + x*box_w, x0 + (x+1)*box_w
    #    y_lower, y_upper = y0 + y*box_h, y0 + (y+1)*box_h
        n_xs = x0 + (x + np.arange(0, 2, dtype=np.float64)) * box_w
        n_ys = y0 + (y + np.arange(0, 2, dtype=np.float64)) * box_h

        if not bottom_up:
            # n.b. assumes that the projection is symmetric
            n_ys = -1 * n_ys[::-1]

        return n_xs, n_ys

    def tileextent(self, x_y_z):
        """Returns extent tuple in MercatorCoords ``(x0, x1, y0, y1)``."""
        x, y, z = x_y_z
        x_lim, y_lim = self.tile_bbox(self.crs, x, y, z, bottom_up=True)
        return tuple(x_lim) + tuple(y_lim)

    _tileextent = tileextent

    def _image_url(self, tile):
        url = ('http://chart.apis.google.com/chart?chst=d_text_outline&'
               'chs=256x256&chf=bg,s,00000055&chld=FFFFFF|16|h|000000|b||||'
               'Google:%20%20(' + str(tile[0]) + ',' + str(tile[1]) + ')'
               '|Zoom%20' + str(tile[2]) + '||||||______________________'
               '______')
#        print url
#        url = ('http://mts0.google.com/vt/lyrs=m@177000000&hl=en&src=api&'
#               'x=%s&y=%s&z=%s&s=G' % tile)
        return url

    def get_image(self, tile):
        import cStringIO  # *much* faster than StringIO
        import urllib

        url = self._image_url(tile)

        fh = urllib.urlopen(url)
        im_data = cStringIO.StringIO(fh.read())
        fh.close()
        img = Image.open(im_data)

        img = img.convert(self.desired_tile_form)

        return img, self.tileextent(tile), 'lower'


class MapQuestOSM(GoogleTiles):
    # http://developer.mapquest.com/web/products/open/map for terms of use
    def _image_url(self, tile):
        x, y, z = tile
        url = 'http://otile1.mqcdn.com/tiles/1.0.0/osm/%s/%s/%s.jpg' % (
            z, x, y)
        return url


class MapQuestOpenAerial(GoogleTiles):
    # http://developer.mapquest.com/web/products/open/map for terms of use
    # The following attribution should be included in the resulting image:
    # "Portions Courtesy NASA/JPL-Caltech and U.S. Depart. of Agriculture,
    #  Farm Service Agency"
    def _image_url(self, tile):
        x, y, z = tile
        url = 'http://oatile1.mqcdn.com/tiles/1.0.0/sat/%s/%s/%s.jpg' % (
            z, x, y)
        return url


class OSM(GoogleTiles):
    # http://developer.mapquest.com/web/products/open/map for terms of use
    def _image_url(self, tile):
        x, y, z = tile
        url = 'http://a.tile.openstreetmap.org/%s/%s/%s.png' % (z, x, y)
        return url


class QuadtreeTiles(GoogleTiles):
    """
    Implements web tile retrieval using the Microsoft WTS quadkey coordinate
    system.

    A "tile" in this class refers to a quadkey such as "1", "14" or "141"
    where the length of the quatree is the zoom level in Google Tile terms.

    """
    def _image_url(self, tile):
        url = ('http://ecn.dynamic.t1.tiles.virtualearth.net/comp/'
               'CompositionHandler/' + str(tile) + '?mkt=en-'
               'gb&it=A,G,L&shading=hill&n=z')
        url = ('http://chart.apis.google.com/chart?chst=d_text_outline&'
               'chs=256x256&chf=bg,s,00000055&chld=FFFFFF|16|h|000000|b'
               '||||Quadkey:%20%20(' + str(tile) + ')||||||____________'
               '________________')

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
        # http://msdn.microsoft.com/en-us/library/bb259689.aspx
        assert isinstance(quadkey, basestring), 'quadkey must be a string'

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
        return GoogleTiles.tileextent(self, x_y_z)

    def find_images(self, target_domain, target_z, start_tile=None):
        """
        Find all the quadtree's at the given target zoom, in the given
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
            for tile in GoogleTiles.find_images(self, target_domain, target_z,
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
