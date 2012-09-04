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


import matplotlib.axes
from matplotlib.image import imread
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.collections as mcollections
import numpy
import shapely.geometry

#import cartopy.cartesian
import cartopy
import cartopy.crs as ccrs
import cartopy.img_transform
import cartopy.mpl_integration.patch as patch



import matplotlib
assert matplotlib.__version__ >= '1.2', 'Cartopy can only work with matplotlib 1.2 or greater.'



colors = {
          'land': numpy.array((240, 240, 220))/256.,
          'sea': numpy.array((152, 183, 226))/256.,
          }

# XXX call this InterCRSTransform
class InterProjectionTransform(mtransforms.Transform):
    """Transforms coordinates from the source_projection to the target_projection."""
    input_dims = 2
    output_dims = 2
    is_separable = False
    has_inverse = True
    
    def __init__(self, source_projection, target_projection):
        # assert target_projection is cartopy.crs.Projection
        # assert source_projection is cartopy.crs.CRS
        self.source_projection = source_projection
        self.target_projection = target_projection
        mtransforms.Transform.__init__(self)

    def __repr__(self):
        return '< InterProjectionTransform %s -> %s >' % (self.source_projection, self.target_projection)

    def transform_non_affine(self, xy):
        """Transforms from native coordinates to lat lons."""
        if isinstance(xy, numpy.ndarray):
#            print xy.shape
#            x, y = xy[:, 0:1], xy[:, 1:2]
#            r = []
#            for xpt, ypt in zip(x, y):
#                print xpt, ypt
            return self.target_projection.transform_points(self.source_projection, xy[:, 0], xy[:, 1])[:, 0:2]
        else:
            x, y = xy
            x, y = self.target_projection.transform_point(x, y, self.source_projection)
            return x, y

    
    def transform_point(self, point):
        # XXX Needs testing
        return numpy.array([[point[0], point[1]]])
    
    def transform_path_non_affine(self, path):
        if path.vertices.shape == (1, 2):
            return mpath.Path(self.transform(path.vertices))
        
        transformed_geoms = []
        for geom in patch.path_to_geos(path):
            transformed_geoms.append(self.target_projection.project_geometry(geom, self.source_projection))
        
        if not transformed_geoms:
            return mpath.Path(numpy.empty([0, 2]))
        else:
            paths = patch.geos_to_path(transformed_geoms)
            if not paths:
                return mpath.Path(numpy.empty([0, 2]))
            points, codes = zip(*[patch.path_segments(path, curves=False, simplify=False) for path in paths])
            return mpath.Path(numpy.concatenate(points, 0), numpy.concatenate(codes))

    def inverted(self):
        return InterProjectionTransform(self.target_projection, self.source_projection)


class GenericProjectionAxes(matplotlib.axes.Axes):
    def __init__(self, *args, **kwargs):
        self.projection = kwargs.pop('map_projection')        
        super(GenericProjectionAxes, self).__init__(*args, **kwargs)
        
    def __str__(self):
        return '< GenericProjectionAxes: %s >' % self.projection

    def cla(self):
        result = matplotlib.axes.Axes.cla(self)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.autoscale_view(tight=True)
        self.set_aspect('equal')
#        self.patch.set_transform(self.transData)
        self._boundary()
        
        # gives a better display of data. However, does need to be made clear in the docs.
        self._xmargin = 0.15
        self._ymargin = 0.15
        
        return result
    
    def format_coord(self, x, y):
        """Return a format string formatting the coordinate value for GUI purposes only."""
        lon, lat = self.projection.transform_point(x, y, ccrs.Geodetic())
        
        ns = 'N' if lat >= 0.0 else 'S'
        ew = 'E' if lon >= 0.0 else 'W'

        return u'%.4g, %.4g (%f\u00b0%s, %f\u00b0%s)' % (x, y, abs(lat), ns, abs(lon), ew)
    
    def coastlines(self):
        import cartopy.io.shapereader as shapereader
        
        coastline_path = shapereader.natural_earth(resolution='110m', 
                                                   category='physical', 
                                                   name='coastline')
    
        shapereader.mpl_axes_plot(self, shapereader.Reader(coastline_path).geometries())
        
    def coastlines_land(self, facecolor=colors['land'], **kwargs):
        import cartopy.io.shapereader as shapereader
        
        land_path = shapereader.natural_earth(resolution='110m', 
                                               category='physical', 
                                               name='land')
        
        paths = []
        for geom in shapereader.Reader(land_path).geometries():
            
            paths.extend(patch.geos_to_path(self.projection.project_geometry(geom)))
        self.add_collection(mcollections.PathCollection(paths, facecolor=facecolor, **kwargs), autolim=False)
        
    def gshhs(self, outline_color='k', land_fill='green', ocean_fill='None', resolution='coarse', domain=None):
        import cartopy.gshhs as gshhs
        from matplotlib.collections import PatchCollection
        
        if ocean_fill is not None and land_fill is None:
            land_fill = self.outline_patch.get_facecolor()
        
        if domain is None:
            domain = self.ll_boundary_poly()
            
        paths = []
        for points in gshhs.read_gshhc(gshhs.fnames[resolution], poly=True, domain=domain):
            # XXX Sometimes we only want to do lines...
            poly = shapely.geometry.Polygon(points[::-1, :])
            projected = self.projection.project_polygon(poly)
            
            paths.extend(patch.geos_to_path(projected))
            
        if ocean_fill is not None:
            self.outline_patch.set_facecolor(ocean_fill) 
            
        collection = PatchCollection([mpatches.PathPatch(pth) for pth in paths], 
                                     edgecolor=outline_color,
                                     facecolor=land_fill, 
                                     zorder=2,
                                     )
        # XXX Should it update the limits??? (AND HOW???)
        self.add_collection(collection)
    
    def native_extents(self):
        west, south, east, north = self.viewLim.get_points().flatten()
        if west == 0 and south == 0 and east == 1 and north == 1:
            import itertools
            west, east, south, north = itertools.chain(self.projection.x_limits, self.projection.y_limits)
        return west, east, south, north
        
    def boundary_poly(self):
        # returns the boundary of the axes at its current state
        x1, x2, y1, y2 = self.native_extents()
        boundary = shapely.geometry.Polygon([[x1, y1], [x1, y2], [x2, y2], [x2, y1], [x1, y1]])
        return boundary
    
    def map_domain(self, crs):
        x1, x2, y1, y2 = self.native_extents()
        native_domain = shapely.geometry.LineString([[x1, y1], 
                                                  [x2, y1], 
                                                  [x2, y2], 
                                                  [x1, y2], 
                                                  [x1, y1]])
        return crs.project_geometry(native_domain, self.projection)
    
    def ll_boundary_poly(self):
        native_boundary = self.boundary_poly()
        foo = native_boundary
        return cartopy.prj.PlateCarree().project_geometry(foo, self.projection)
    
    def ll_boundary_poly_draw(self):
        for path in patch.geos_to_path(self.ll_boundary_poly()):
        # XXX Seems like great circle interpolation is making the plot strange...
            pp = mpatches.PathPatch(path, color='red', transform=cartopy.prj.PlateCarree(), alpha=0.5)
            dl = self.viewLim.get_points().copy()
            self.add_patch(pp)
            self.viewLim.set_points(dl)
    
    def set_global(self):
        self.set_xlim(self.projection.x_limits)
        self.set_ylim(self.projection.y_limits)
    
    def geod_circle_meters(self, lon_0, lat_0, radius, npts=80, **kwargs):
        # radius is in meters
        geod = self.projection.as_geodetic()
        
        az = numpy.linspace(0, 360, npts)
        lats = numpy.zeros(npts) + lat_0
        lons = numpy.zeros(npts) + lon_0
        distances = numpy.zeros(npts) + radius
        
        lons, lats, _reverse_az = geod.fwd(lons, lats, az, distances, radians=False)
        ll = numpy.concatenate([lons[:, None], lats[:, None]], 1)
        from matplotlib.patches import Polygon
        poly = Polygon(ll, transform=cartopy.prj.PlateCarree(), **kwargs)
        self.add_patch(poly)
        return poly
    
    def gshhs_line(self, outline_color='k', domain=None, resolution='low', **kwargs):
        # domain is a shapely geometry (Polygon or MultiPolygon)
        import cartopy.gshhs as gshhs
#        import cartopy.spherical as spherical
        from matplotlib.collections import PatchCollection, LineCollection
        
        paths = []
            
        projection = self.projection
        
        if domain is None:
            domain = self.map_domain(ccrs.PlateCarree())
            
        for points in gshhs.read_gshhc(gshhs.fnames[resolution], poly=False, domain=domain):
            paths.extend(patch.geos_to_path(shapely.geometry.LineString(points)))
            
#            slinestring = shapely.geometry.LineString(points)
#            projected = projection.project_geometry(slinestring)            
#            paths.extend(patch.geos_to_path(projected))

        collection = PatchCollection([mpatches.PathPatch(pth) for pth in paths], 
                             edgecolor=outline_color, facecolor='none',
                             transform=ccrs.PlateCarree(),
                             **kwargs
                             )
        
        self.add_collection(collection, autolim=False)
    
    def bluemarble(self):
        # XXX remove this method.
        return self.stock_img('bluemarble')

    def stock_img(self, img_name):
        # XXX Turn into a dictionary (inside the method)?
        if img_name == 'bluemarble':
            source_proj = ccrs.PlateCarree()
            fname = '/data/local/dataZoo/cartography/raster/blue_marble_720_360.png'
#            fname = '/data/local/dataZoo/cartography/raster/blue_marble_2000_1000.jpg'            
            img_origin = 'lower'
            img = imread(fname)
            img = img[::-1]
            return self.imshow(img, origin=img_origin, transform=source_proj, extent=[-180, 180, -90, 90])
        elif img_name == 'bm_high':
            source_proj = ccrs.PlateCarree()
            fname = '/data/local/dataZoo/cartography/raster/blue_marble_2000_1000.jpg'            
            img_origin = 'lower'
            img = imread(fname)
            return self.imshow(img, origin=img_origin, transform=source_proj, extent=[-180, 180, -90, 90])
        elif img_name == 'ne_shaded':
            source_proj = ccrs.PlateCarree()
            fname = '/data/local/dataZoo/cartography/raster/NE1_50M_SR_W/NE1_50M_SR_W_720_360.png'
            img_origin = 'lower'
            img = imread(fname)
            img = img[::-1]
            return self.imshow(img, origin=img_origin, transform=source_proj, extent=[-180, 180, -90, 90])
        else:
            raise ValueError('Unknown stock image.')
    
#    def margins(self, *args, **kw):
#        result = matplotlib.axes.Axes.margins(self, *args, **kw)
#        if result is not None:
    
    def imshow(self, img, *args, **kwargs):
        """
        Add the "transform" keyword to imshow.
        
        Extra kwarg:
        transform - is actually a PROJECTION NOT a transform
        regrid_shape - default is (750, 375). But may be changed to "auto" in the future...
        extent = (left, right, bottom, top) - transform coordinates for the extent of the source image.
        target_extent = (left, right, bottom, top) - native coordinates for the extent of the desired image.
        origin - default is changed to 'lower'
        update_datalim - flag whether the image should affect the data limits (default: True)
        
        """
        transform = kwargs.pop('transform', None)
        regrid_shape = kwargs.pop('regrid_shape', (750, 375))
        update_datalim = kwargs.pop('update_datalim', True)
               
        kwargs.setdefault('origin', 'lower') 
        
        same_projection = isinstance(transform, ccrs.Projection) and self.projection == transform
        
        if not update_datalim:
            data_lim = self.dataLim.frozen().get_points()
            view_lim = self.viewLim.frozen().get_points()
        
        if transform is None or transform == self.transData or same_projection:
            if isinstance(transform, ccrs.Projection):
                transform = transform._as_mpl_transform(self)
            result = matplotlib.axes.Axes.imshow(self, img, *args, **kwargs)
        else:
            extent = kwargs.pop('extent', None)
            
            if not isinstance(transform, ccrs.Projection):
                raise ValueError('Expected a projection. Cannot handle a %s in imshow.' % type(transform))
                    
            # XXX adaptive resolution depending on incoming img?
            img, extent = cartopy.img_transform.warp_array(img, 
                                                           source_proj=transform,
                                                           source_extent=extent,
                                                           target_proj=self.projection,
                                                           target_res=regrid_shape,
                                                           target_extent=self.native_extents(),
                                                           )
            result = matplotlib.axes.Axes.imshow(self, img, *args, extent=extent, **kwargs)
        
        # clip the image. This does not work as the patch moves with mouse movement, but the clip path doesn't
        # This could definitely be fixed in matplotlib 
#        if result.get_clip_path() in [None, self.patch]:
#            # image does not already have clipping set, clip to axes patch
#            result.set_clip_path(self.outline_patch)
        
        if not update_datalim:
            self.dataLim.set_points(data_lim)
            self.viewLim.set_points(view_lim)

        return result
        
    def gridlines(self, res):
        if isinstance(self.projection, ccrs.Orthographic):
            import warnings
            warnings.warn('Gridlines and orthographic projections are causing a problem, skipping this operation.')
            return
        # XXX needs to use proper mpl axes
        step = 15
        lons = range(0, 360, step)
        lats = [-75, 0, 75]
        for lon in lons:
            self.plot([lon] * len(lats), lats, linestyle=':', color='k', 
                      transform=ccrs.Geodetic(), 
                      scalex=False, scaley=False
                      )
            
        lons = lons + [lons[0] + 360]
        lats = range(-90 + step, 90, step)
        for lat in lats:
            self.plot(lons, [lat] * len(lons), linestyle=':', color='k', 
                      transform=ccrs.Geodetic(), 
                      scalex=False, scaley=False
                      )

    def _gen_axes_spines(self, locations=None, offset=0.0, units='inches'):
        # generate some axes spines, as some Axes super class machinery requires them. Just make them invisible
        spines = matplotlib.axes.Axes._gen_axes_spines(self, locations=locations, offset=offset, units=units)
        for spine in spines.itervalues():
            spine.set_visible(False)
        return spines

    def _boundary(self):
        """
        Adds the map's boundary. Note, the boundary is not the axes.patch, which provides rectilinear 
        clipping for all of the map's artists.
        The axes.patch will have its visibility set to False inside GeoAxes.gca()
        """
        import cartopy.mpl_integration.patch as p
        path, = p.geos_to_path(self.projection.boundary)
        
#        from matplotlib.collections import PatchCollection
            
        self.sct = sct = SimpleClippedTransform(self.transScale + self.transLimits, self.transAxes)
        self.smt = SimpleMaskingTransform(self.transScale + self.transLimits, self.transAxes)
    
        # XXX Should be exactly one path...
        collection = mpatches.PathPatch(path,
                                        facecolor='none', edgecolor='k', zorder=1000, 
#                                        transform=self.transData,
                                        transform=sct, clip_on=False,
                                        )
        self.outline_patch = collection
        # XXX autolim = False
        self.add_patch(collection)
        
        # put a color patch for background color
        # XXX Should be exactly one path...
        collection = mpatches.PathPatch(path,
                                        facecolor='w', edgecolor='none', zorder=-1, 
                                        transform=sct, clip_on=False,
                                        )
        self.background_patch = collection
        # XXX autolim = False
        self.add_patch(collection)
        
        
        self.patch.set_facecolor((1, 1, 1, 0))
        self.patch.set_edgecolor((0.5, 0.5, 0.5))
        self.patch.set_linewidth(0.0)


class SimpleClippedTransform(mtransforms.Transform):
    """
    Transforms the values using a pre transform, clips them, then post transforms them.
    
    This transform should not be widely used, but is useful for transforming a background patch 
    and clipping the patch to a desired extent.
    
    """
    input_dims = 2
    output_dims = 2
    def __init__(self, pre_clip_transform, post_clip_transform, xclip=(0, 1), yclip=(0, 1)):
        mtransforms.Transform.__init__(self)
        self.pre_clip_transform = pre_clip_transform
        self.post_clip_transform = post_clip_transform
    
        self.x_clips = xclip
        self.y_clips = yclip
                
    def transform_non_affine(self, values):
        new_vals = self.pre_clip_transform.transform(values)
        x, y = new_vals[:, 0:1], new_vals[:, 1:2]
        numpy.clip(x, self.x_clips[0], self.x_clips[1], x)
        numpy.clip(y, self.y_clips[0], self.y_clips[1], y)
        # XXX support ma's?
        return self.post_clip_transform.transform(new_vals)

class SimpleMaskingTransform(SimpleClippedTransform):
    def transform_non_affine(self, values):
        values = numpy.array(values)
        single_point = False
        if values.ndim == 1:
            single_point = True
            values.shape = (1, 2)
        
        new_vals = self.pre_clip_transform.transform(values)
        x, y = new_vals[:, 0:1], new_vals[:, 1:2]
        new_vals = numpy.ma.masked_array(new_vals, mask=None)
        new_vals.mask[:, 1:2] = (x < self.x_clips[0]) | (x > self.x_clips[1])
#        new_vals.mask[:, 0:1] = (y < self.y_clips[0]) | (y > self.y_clips[1])
        numpy.clip(y, self.y_clips[0], self.y_clips[1], y)
        result =  self.post_clip_transform.transform(new_vals)
        
        if single_point:
            return result[0, :]
        else:
            return result


def add_holey_poly(transform):
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection
    from matplotlib.path import Path
    
    
    poly = mpatches.RegularPolygon( (140, 10), 4, 81.0) 
    # XXX internal rings are not yet supported...
    pth = Path([[0, 45], [60, 45], [60, -45], [0, -45], [0, -45], [10, 20], [10, -20], [40, -20], [40, 20], [10, -20]], [1, 2, 2, 2, 79, 1, 2, 2 ,2, 79]) 
    pth = Path([[0, 45], [60, 45], [60, -45], [0, -45], [0, -45]], [1, 2, 2, 2, 79])
#    pth = Path([[0, 45], [0, -45], [60, -45], [60, 45], [0, -45]], [1, 2, 2, 2, 79])
    poly = mpatches.PathPatch(pth) 
    collection = PatchCollection([poly], cmap=matplotlib.cm.jet, alpha=0.4, 
                                 transform=transform
                                 ) 
    plt.gca().add_collection(collection)
   
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
      
#    proj = ccrs.Robinson()  
#    proj = ccrs.PlateCarree(central_longitude=90)
    proj = ccrs.PlateCarree()
#    proj = ccrs.NorthPolarStereo()
#    proj = ccrs.Orthographic()
#    proj = ccrs.RotatedPole(pole_longitude=177, pole_latitude=37)    
    
    # XXX Interrupted isn't working so well yet.
#    proj = ccrs.InterruptedGoodeHomolosine()
    

    # Initialise the projection
#    ax = plt.subplot(211, projection=proj)
    ax = plt.axes(projection=proj)
    
    # define a trans
    ll = ccrs.PlateCarree()
#    ll_transform = InterProjectionTransform(ll, proj) + ax.transData

#    plt.axis(xrange=[0, 50], yrange=[0, 50])
    # ================

#    ax.bluemarble()
    ax.coastlines()
    ax.gridlines()
#    ax.gshhs()
#    ax.gshhs(ocean_fill='blue', land_fill='green')
    
#    ax.plot(0, 0, 'bx')
#    ax.plot(-14000000, 0, 'yx')
#    add_holey_poly(ll)

    # cardiff in NorthPolar...
#    ax.plot(-2.568e5, -4.399e6, 'bo', transform=ccrs.NorthPolarStereo(), scalex=False, scaley=False)
    
    # cape town
#    ax.plot(18.25, -33.5, 'ro', transform=ll)
#    plt.ion()
    coords = [(-0.08, 51.53), (132.00, 43.17)] # London to Vladivostock
    r = ax.plot(*zip(*coords), transform=ll)
    ax.set_global()
#    import cartopy.test
#    cs = plt.contourf(*cartopy.test.wave_data(), transform=ll)
#    cs = plt.contour(*cartopy.test.wave_data(), transform=ll)
    
#    x, y, s = cartopy.test.wave_data()
#    cs = plt.scatter(x.flat, y.flat, transform=ll)
#    print x.flat[101:102], y.flat[201:202]
       
#    plt.draw()
#    print plt.gca().outline_patch.get_transform()
    plt.show()
    
       
#    # ================
#    plt.show()
#    print r[0].get_transform()._invalid
#    plt.draw()
#    print r[0].get_transform()._invalid
#    import time
#    time.sleep(2)
#    plt.gcf().set_size_inches([10, 10])
#    print r[0].get_transform()._invalid
#    plt.draw()
#    print '------------\n' * 4
##    print r[0].get_transform()
#    print r[0].get_transform()._invalid
#    time.sleep(6)
##    plt.show()
#    
