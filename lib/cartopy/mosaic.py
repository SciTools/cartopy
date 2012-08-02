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

Please note: This is a proof of concept implementation of retrieval of the Google and 
Bing web map tiles. Storage of the images using tiles from the Google and Bing services
contravenes the respective service policies and is not condoned. 

"""
# TODO: IMPLEMENT A FILE INSPECTION VERSION WHICH CACHES THE ZOOM STRUCTURE...


import numpy

import cartopy.crs as ccrs

import PIL.Image as Image

import shapely.geometry


class GoogleTiles(object):
    def __init__(self, desired_tile_form='RGB'):
        # XXX consider fixing the CRS???
        self.imgs = []
        self.desired_tile_form = desired_tile_form
    
    def image_for_domain(self, target_domain, target_z):
        tiles = []
        for tile in self.find_images(target_domain, target_z):
            try:
                img, extent, origin = self.get_image(tile)
            except IOError:
                continue
#            img, extent, orientation = self.get_image(tile)
#            print 'img: ', tile, extent
#            
            img = numpy.array(img)
            import numpy as np
            x = np.linspace(extent[0], extent[1], img.shape[1], endpoint=False)
            y = np.linspace(extent[2], extent[3], img.shape[0], endpoint=False)
            tiles.append([numpy.array(img), x, y, origin])

        img, extent, origin = merge_tiles(tiles)
        return img, extent, origin
    
    def _find_images(self, target_domain, target_z, start_tile=(0, 0, 0)):
        """Target domain is a shapely polygon in native coordinates."""
        
        assert isinstance(target_z, int) and target_z >= 0, 'target_z must be an integer >=0.'
        
        # recursively drill down to the images at the target zoom
        domain = self.tiledomain(start_tile)
        if domain.intersects(target_domain):
            if start_tile[2] == target_z:
                    yield start_tile
            else:
                for tile in self._subtiles(start_tile):
                    for result in self._find_images(target_domain, target_z, start_tile=tile):
                        yield result
    
    find_images = _find_images
        
    def draw_all_domains(self, ax, **kwargs):
        for _, _, domain, crs in self.imgs:
            self.draw_domain(ax, domain, crs, **kwargs)
                
    def draw_found_domains(self, ax, target_domain, **kwargs):
        for _, _, domain, crs in self.find_images(target_domain):
            self.draw_domain(ax, domain, crs, **kwargs)
                
    def draw_domain(self, ax, domain, crs, **kwargs):
        import matplotlib.patches as mpatches
        import cartopy.mpl_integration.patch as patch
        
        for path in patch.geos_to_path(domain):
            p = mpatches.PathPatch(path, transform=crs._as_mpl_transform(ax),
                                   **kwargs)
            ax.add_patch(p)

    def subtiles(self, x_y_z):
        x, y, z = x_y_z
        # google tile specific (i.e. up->down)
        
#        yield x * 2, y * 2, z + 1 
#        return
        for xi in range(0, 2):
            for yi in range(0, 2):
                yield x * 2 + xi, y * 2 + yi, z + 1

    _subtiles = subtiles 

    def tileextent(self, x_y_z):
        x, y, z = x_y_z
        # this was a copy from tiledomain
        import cartopy.examples.gmaptiles as gmt
        prj = ccrs.Mercator()
        x_lim, y_lim = gmt.tile_bbox(prj, x, y, z, bottom_up=True)
        
        return tuple(x_lim) + tuple(y_lim)
#        result = ccrs.PlateCarree().transform_points(prj, x_lim.astype(numpy.float64), y_lim.astype(numpy.float64))
#        x_lim = result[:, 0]
#        y_lim = result[:, 1]
#        return tuple(x_lim) + tuple(y_lim)
        
    def tiledomain(self, x_y_z):
#        print x_y_z
        x, y, z = x_y_z
        import cartopy.examples.gmaptiles as gmt
        prj = ccrs.Mercator()
        x_lim, y_lim = gmt.tile_bbox(prj, x, y, z, bottom_up=True)
        
        result = ccrs.PlateCarree().transform_points(prj, x_lim.astype(numpy.float64), y_lim.astype(numpy.float64))
        x_lim = result[:, 0]
        y_lim = result[:, 1]
        
        domain = shapely.geometry.Polygon([[x_lim[0], y_lim[0]], 
                                  [x_lim[1], y_lim[0]], 
                                  [x_lim[1], y_lim[1]], 
                                  [x_lim[0], y_lim[1]], 
                                  [x_lim[0], y_lim[0]]])
#        return x_lim, y_lim
        return domain
    
    def _image_url(self, tile):
        url = 'http://chart.apis.google.com/chart?chst=d_text_outline&chs=256x256&chf=bg,' + \
              's,00000055&chld=FFFFFF|16|h|000000|b||||Google:%20%20('+str(tile[0])+','+str(tile[1])+')' + \
              '|Zoom%20'+str(tile[2])+'||||||____________________________'
              
        url = 'http://mts0.google.com/vt/lyrs=m@177000000&hl=en&src=api&x=%s&y=%s&z=%s&s=G' % tile
        return url
    
    def get_image(self, tile):
        # get the tile from google, and turn it into RGB
        import cStringIO
        import urllib
        
        url = self._image_url(tile)      
                
        fh = urllib.urlopen(url)
        im_data = cStringIO.StringIO(fh.read())
        fh.close()
        img = Image.open(im_data)
        
        img = img.convert(self.desired_tile_form)
        
        return img, self.tileextent(tile), 'lower'
        
        
class BingTiles(GoogleTiles):
    def _image_url(self, tile):
        url = 'http://ecn.dynamic.t1.tiles.virtualearth.net/comp/CompositionHandler/' + str(tile) +'?mkt=en-gb&it=A,G,L&shading=hill&n=z'
#        url = 'http://chart.apis.google.com/chart?chst=d_text_outline&chs=256x256&chf=bg,' + \
#              's,00000055&chld=FFFFFF|16|h|000000|b||||Bing:%20%20('+str(tile)+')' + \
#              '||||||____________________________'
        
        return url
    
    def tms_to_quadkey(self, tms, google=False):
        quadKey = ""
        x, y, z = tms
        # this algorithm works with google tiles, rather than tms, so convert to those first.
        if not google:
            y = (2**z - 1) - y
        for i in range(z, 0, -1):
            digit = 0
            mask = 1 << (i-1)
            if (x & mask) != 0:
                digit += 1
            if (y & mask) != 0:
                digit += 2
            quadKey += str(digit)
        return quadKey
    
    def quadkey_to_tms(self, quadkey, google=False):
        # algorithm ported from http://msdn.microsoft.com/en-us/library/bb259689.aspx
        assert isinstance(quadkey, basestring), 'quadkey must be a string'
        
        x = y = 0;
        z = len(quadkey)
        for i in range(z, 0, -1):
            mask = 1 << (i -1)
            if quadkey[z-i] == '0':
                pass
            elif quadkey[z-i] == '1':
                x |= mask
            elif quadkey[z-i] == '2':
                y |= mask
            elif quadkey[z-i] == '3':
                x |= mask
                y |= mask
            else:
                raise ValueError('Invalid QuadKey digit sequence.' + str(quadkey))
        # the algorithm works to google tiles, so convert to tms
        if not google:
            y = (2**z - 1) - y
        return (x, y, z)
    
    def subtiles(self, quadkey):
        for i in range(4):
            yield quadkey + str(i)
            
    def tileextent(self, quadkey):
        x_y_z = self.quadkey_to_tms(quadkey, google=True)
#        print 'tile extent: ', quadkey, x_y_z, GoogleTiles.tileextent(self, x_y_z)
        return GoogleTiles.tileextent(self, x_y_z)
        
    def find_images(self, target_domain, target_z, start_tile=None):
        if start_tile is None:
            start_tiles = ['0', '1', '2', '3']
        else:
            start_tiles = [start_tile]
            
        for start_tile in start_tiles:
#            st = start_tile
            start_tile = self.quadkey_to_tms(start_tile, google=True)
#            print 'start tile', start_tile, st
            for tile in GoogleTiles.find_images(self, target_domain, target_z, start_tile=start_tile):
#                print 'yielded ', tile, self.tms_to_quadkey(tile, google=True) 
                yield self.tms_to_quadkey(tile, google=True)
                

def draw_child_domains(tile_system):

    gt = tile_system

    import matplotlib.pyplot as plt
    ax = plt.axes(projection=ccrs.OSGB())
    ax.set_global()

    start_tile = ('SO', '1:250,000')
    target_domain = gt.tiledomain(start_tile)
    gt.draw_domain(ax, target_domain, ccrs.OSGB(), facecolor='coral', alpha='0.3')

    for img in gt.subtiles(start_tile):
        gt.draw_domain(ax, gt.tiledomain(img), ccrs.OSGB(), facecolor='green', alpha='0.3')
    ax.coastlines()
    plt.show()


def merge_tiles(tiles):
    """return a single image when the given images are merged."""
    xset = [set(x) for i, x, y, _ in tiles]
    yset = [set(y) for i, x, y, _ in tiles]
    
    xs = xset[0]
    xs.update(*xset[1:])
    ys = yset[0]
    ys.update(*yset[1:]) 
    xs = sorted(xs)
    ys = sorted(ys)

    import numpy as np
    other_len = tiles[0][0].shape[2:]
    img = np.zeros((len(ys), len(xs)) + other_len, dtype=np.uint8) - 1
    
    for tile_img, x, y, origin in tiles:
        y_first, y_last = y[0], y[-1]
        yi0, yi1 = np.where((y_first == ys) | (y_last == ys))[0]
        if origin == 'upper':
            yi0 = tile_img.shape[0] - yi0 - 1
            yi1 = tile_img.shape[0] - yi1 - 1
        start, stop, step = yi0, yi1, 1 if yi0 < yi1 else -1
        if step == 1 and stop == img.shape[0]-1:
            stop = None
        elif step == -1 and stop == 0:
            stop = None
        else:
            stop += step
        y_slice = slice(start, stop, step)
        
        xi0, xi1 = np.where((x[0] == xs) | (x[-1] == xs))[0]
#        if origin == 'lower':
#            xi0 = tile_img.shape[1] - xi0 - 1
#            xi1 = tile_img.shape[1] - xi1 - 1
        start, stop, step = xi0, xi1, 1 if xi0 < xi1 else -1
        
        if step == 1 and stop == img.shape[1]-1:
            stop = None
        elif step == -1 and stop == 0:
            stop = None
        else:
            stop += step
            
        x_slice = slice(start, stop, step)
        
        img_slice = (y_slice, x_slice, Ellipsis)  
#        print img_slice
#        # XXX probably also depends on origin...
#        if origin == 'upper':
#            direction *= -1
#            
#            yi = img.shape[0] - yi
#            img_slice = list(img_slice)
#            y2 = yi - len(y) -1
#            if y2 <= 0:
#                y2 = None
#            img_slice[0] = slice(yi-1, y2, -1)
#            img_slice = tuple(img_slice)
#            print img_slice
#        else:
#            tile_img = tile_img[::direction, ...] 
        
        if origin == 'lower':
            tile_img = tile_img[::-1, ::]
        
        img[img_slice] = tile_img

    return img, [min(xs), max(xs), min(ys), max(ys)], 'lower'
        

def draw_tiles_in_domain(tile_thing, target_domain, target_z):
    import matplotlib.pyplot as plt
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_global()
    ax.coastlines()
    tile_thing.draw_domain(ax, target_domain, ccrs.PlateCarree(), facecolor='coral', alpha=0.3)
    
    tiles = []
    for tile in tile_thing.find_images(target_domain, target_z):
        img, extent, origin = tile_thing.get_image(tile)
        img = numpy.array(img)
        import numpy as np
        x = np.linspace(extent[0], extent[1], img.shape[1], endpoint=False)
        y = np.linspace(extent[2], extent[3], img.shape[0], endpoint=False)
        tiles.append([numpy.array(img), x, y, origin])
#        ax.imshow(numpy.array(img), extent=extent, transform=ccrs.Mercator())
#        gt.draw_domain(ax, gt.tiledomain(*tile), ccrs.PlateCarree(), facecolor='green', alpha='0.3')

    img, extent, origin = merge_tiles(tiles)
    ax.imshow(numpy.array(img), extent=extent, transform=ccrs.Mercator(), origin=origin,
              interpolation='bilinear', regrid_shape=(1500, 750))
    plt.show()


def google_run():
    gt = GoogleTiles()

    extent = [-15, 00, 50, 60]
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]], 
                                              [extent[2], extent[1]], 
                                              [extent[2], extent[3]], 
                                              [extent[0], extent[3]], 
                                              [extent[0], extent[1]]])

    # XXX implement as tests

#    for img in gt.find_images(target_domain, 2):
#        print img
#
#    for img in gt.subtiles(0, 0, 1):
#        print img
#    
#    print
#    for img in gt.subtiles(1, 0, 1):
#        print img
#        
#    for img in gt.subtiles(2, 2, 2):
#        print img
#
#    gt.tiledomain(0, 0, 0)
#    print gt.tiledomain(0, 0, 0)
#    print gt.tiledomain(2**1-1, 0, 1)
#    print gt.tiledomain(2**2-1, 0, 2)
#    print 'f', gt.tiledomain((2**3)-1, 0, 3)
#    print gt.tiledomain(0, 0, 3)
#    print gt.tiledomain(0, 0, 4)
#    assert '11', gt.tiledomain(2, 1, 4)
#    
#    print gt.tiledomain(2, 2, 2)
#    print gt.tiledomain(8, 9, 4)

#    draw_child_domains(gt)
    draw_tiles_in_domain(gt, target_domain, 3)
    
def bing_run():
    bt = BingTiles()

    extent = [-15, 00, 50, 60]
    
#    extent = [-3.133617, 50.1218, -3.933617, 51.2218,]
    
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]], 
                                              [extent[2], extent[1]], 
                                              [extent[2], extent[3]], 
                                              [extent[0], extent[3]], 
                                              [extent[0], extent[1]]])

    # XXX Implement as tests    
#    print bt.tms_to_quadkey((1, 1, 1)), bt.quadkey_to_tms(bt.tms_to_quadkey((1, 1, 1)))
#    print bt.tms_to_quadkey((1, 0, 1))
#    print bt.quadkey_to_tms('1')
#    print bt.quadkey_to_tms('3'), bt.tms_to_quadkey(bt.quadkey_to_tms('3'))

#    for img in bt.find_images(target_domain, 2):
#        print 'foo', img
#
#    for img in bt.subtiles('1'):
#        print img

    draw_tiles_in_domain(bt, target_domain, 3)


if __name__ == '__main__':
    google_run()
    bing_run()