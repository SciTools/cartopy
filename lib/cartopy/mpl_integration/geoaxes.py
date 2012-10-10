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


_colors = {
          'land': numpy.array((240, 240, 220)) / 256.,
          'sea': numpy.array((152, 183, 226)) / 256.,
          }

_PATH_TRANSFORM_CACHE = {}
"""
Stores a mapping between (id(path), hash(self.source_projection), hash(self.target_projection))
and the resulting transformed path.

Provides a significant performance boost for contours which, at matplotlib 1.2.0 called
transform_path_non_affine twice unnecessarily.

"""

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
            return self.target_projection.transform_points(self.source_projection, xy[:, 0], xy[:, 1])[:, 0:2]
        else:
            x, y = xy
            x, y = self.target_projection.transform_point(x, y, self.source_projection)
            return x, y

    def transform_path_non_affine(self, path):
        orig_id = (id(path), hash(self.source_projection), hash(self.target_projection))
        result = _PATH_TRANSFORM_CACHE.get(orig_id, None)
        if result is not None:
            return result
        
        bypass = self.source_projection == self.target_projection
        if bypass:
            projection = self.source_projection
            if isinstance(projection, ccrs._CylindricalProjection):
                x = path.vertices[:, 0]
                x_limits = projection.x_limits
                bypass = x.min() >= x_limits[0] and x.max() <= x_limits[1]
        if bypass:
            return path

        if path.vertices.shape == (1, 2):
            return mpath.Path(self.transform(path.vertices))

        transformed_geoms = []
        for geom in patch.path_to_geos(path):
            transformed_geoms.append(self.target_projection.project_geometry(geom, self.source_projection))

        if not transformed_geoms:
            result = mpath.Path(numpy.empty([0, 2]))
        else:
            paths = patch.geos_to_path(transformed_geoms)
            if not paths:
                return mpath.Path(numpy.empty([0, 2]))
            points, codes = zip(*[patch.path_segments(path, curves=False, simplify=False) for path in paths])
            result = mpath.Path(numpy.concatenate(points, 0), numpy.concatenate(codes))
        
        # store the result in the cache for future performance boosts    
        _PATH_TRANSFORM_CACHE[orig_id] = result
        
        return result

    def inverted(self):
        return InterProjectionTransform(self.target_projection, self.source_projection)


class GeoAxes(matplotlib.axes.Axes):
    def __init__(self, *args, **kwargs):
        self.projection = kwargs.pop('map_projection')
        super(GeoAxes, self).__init__(*args, **kwargs)
        self._gridliners = []
        self.img_factories = []
        self._done_img_factory = False

    def add_image(self, factory, *args, **kwargs):
        """
        Adds an image "factory" to the Axes.

        Any image "factory" added, will be asked to retrieve an image
        with associated metadata for a given bounding box at draw time.
        The advantage of this approach is that the limits of the map
        do not need to be known when adding the image factory, but can
        be deferred until everything which can effect the limits has been
        added.

        Currently an image "factory" is just an object with
        a ``image_for_domain`` method. Examples of image factories
        are :class:`cartopy.io.img_nest.NestedImageCollection` and
        :class:`cartopy.io.image_tiles.GoogleTiles`.

        """
        # XXX TODO: Needs working on
        self.img_factories.append([factory, args, kwargs])

    @matplotlib.axes.allow_rasterization
    def draw(self, renderer=None, inframe=False):
        # if no data has been added, and no extents set, then make the map global
        if self.ignore_existing_data_limits and self._autoscaleXon and self._autoscaleYon:
            self.set_global()
            self.ignore_existing_data_limits = True


        for gl in self._gridliners:
            gl.do_gridlines(background_patch=self.background_patch)
        self._gridliners = []

        # XXX This interface needs a tidy up: 
        #       image drawing on pan/zoom; 
        #       caching the resulting image;
        #       buffering the result by 10%...; 
        if not self._done_img_factory:
            for factory, args, kwargs in self.img_factories:
                img, extent, origin = factory.image_for_domain(self._get_extent_geom(factory.crs), args[0])
                self.imshow(img, extent=extent, origin=origin, transform=factory.crs, *args[1:], **kwargs)
        self._done_img_factory = True
        
        
        return matplotlib.axes.Axes.draw(self, renderer=renderer, inframe=inframe)
    
    def __str__(self):
        return '< GeoAxes: %s >' % self.projection

    def cla(self):
        result = matplotlib.axes.Axes.cla(self)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.autoscale_view(tight=True)
        self.set_aspect('equal')

        pre_bounary = self.ignore_existing_data_limits
        self._boundary()
        self.ignore_existing_data_limits = pre_bounary

        # XXX consider a margin - but only when the map is not global... 
        # self._xmargin = 0.15
        # self._ymargin = 0.15

        return result

    def format_coord(self, x, y):
        """Return a string formatted for the matplotlib GUI status bar."""
        lon, lat = ccrs.Geodetic().transform_point(x, y, self.projection)

        ns = 'N' if lat >= 0.0 else 'S'
        ew = 'E' if lon >= 0.0 else 'W'

        return u'%.4g, %.4g (%f\u00b0%s, %f\u00b0%s)' % (x, y, abs(lat), ns, abs(lon), ew)

    _coastlines_cache = {}
    """
    Caches a mapping between "resolution" and a tuple of the resulting geometries.
    Provides a significant performance benefit (when combined with object id caching 
    in shapereader.mpl_axes_plot) when producing multiple maps of the same projection.
    """
    
    def coastlines(self, resolution='110m', **kwargs):
        """
        Adds coastal **outlines** to the current axes from the Natural Earth
        "coastline" shapefile collection.

        Kwargs:

            * resolution - a named resolution to use from the Natural Earth
                           dataset. Currently can be one of "110m", "50m", and
                           "10m".

        .. note::

            Currently no clipping is done on the coastlines before adding them to the axes.
            This means, if very high resolution coastlines are being used, performance is
            likely to be severely effected. This should be resolved transparently by v0.5.

        """
        import cartopy.io.shapereader as shapereader

        if resolution not in self._coastlines_cache: 
            coastline_path = shapereader.natural_earth(resolution=resolution,
                                                       category='physical',
                                                       name='coastline')
            geoms = tuple(shapereader.Reader(coastline_path).geometries())
            # put the geoms in the cache
            self._coastlines_cache[resolution] = geoms 
        else:
            geoms = self._coastlines_cache[resolution]
        
        shapereader.mpl_axes_plot(self, geoms, **kwargs)

    # TODO: expose an interface similar to ax.add_image for shapely things, 
    # and another for adding paths/patches (consider the land shapefile and gshhs as the primary usecases).
#    def _coastlines_land(self, facecolor=colors['land'], **kwargs):
#        import cartopy.io.shapereader as shapereader
#
#        land_path = shapereader.natural_earth(resolution='110m',
#                                               category='physical',
#                                               name='land')
#
#        paths = []
#        for geom in shapereader.Reader(land_path).geometries():
#
#            paths.extend(patch.geos_to_path(self.projection.project_geometry(geom)))
#        self.add_collection(mcollections.PathCollection(paths, facecolor=facecolor, **kwargs), autolim=False)
#
#    def gshhs(self, outline_color='k', land_fill='green', ocean_fill='None', resolution='coarse', domain=None):
#        import cartopy.gshhs as gshhs
#        from matplotlib.collections import PatchCollection
#
#        if ocean_fill is not None and land_fill is None:
#            land_fill = self.outline_patch.get_facecolor()
#
#        if domain is None:
#            domain = self.ll_boundary_poly()
#
#        paths = []
#        for points in gshhs.read_gshhc(gshhs.fnames[resolution], poly=True, domain=domain):
#            # XXX Sometimes we only want to do lines...
#            poly = shapely.geometry.Polygon(points[::-1, :])
#            projected = self.projection.project_polygon(poly)
#
#            paths.extend(patch.geos_to_path(projected))
#
#        if ocean_fill is not None:
#            self.outline_patch.set_facecolor(ocean_fill)
#
#        collection = PatchCollection([mpatches.PathPatch(pth) for pth in paths],
#                                     edgecolor=outline_color,
#                                     facecolor=land_fill,
#                                     zorder=2,
#                                     )
#        # XXX Should it update the limits??? (AND HOW???)
#        self.add_collection(collection)

    def get_extent(self, crs=None):
        """
        Get the extent (x0, x1, y0, y1) of the map in the given coordinate system.

        If no crs is given, the returned extents' coordinate system will be assumed 
        to be the Geodetic version of this axes' projection.
        
        """
        p = self._get_extent_geom(crs)
        r = p.bounds
        x1, y1, x2, y2 = r
        return x1, x2, y1, y2
    
    def _get_extent_geom(self, crs=None):
        x1, x2 = self.get_xlim()
        y1, y2 = self.get_ylim()
        
        if x1 == 0 and y1 == 0 and x2 == 1 and y2 == 1:
            x1, x2 = self.projection.x_limits
            y1, y2 = self.projection.y_limits
        
        domain_in_crs = shapely.geometry.LineString([[x1, y1], [x2, y1],
                                                     [x2, y2], [x1, y2], [x1, y1]])
        
        if self.projection != crs:
            domain_in_crs = self.projection.project_geometry(domain_in_crs, crs)
        
        return domain_in_crs
        
    def set_extent(self, extents, crs=None):
        """
        Set the extent (x0, x1, y0, y1) of the map in the given coordinate system.
        
        If no crs is given, the extents' coordinate system will be assumed to be the
        Geodetic version of this axes' projection.
        
        """
        # TODO: Implement the same semantics as plt.xlim and 
        # plt.ylim - allowing users to set None for a minimum and/or maximum value
        x1, x2, y1, y2 = extents
        domain_in_crs = shapely.geometry.LineString([[x1, y1], [x2, y1],
                                                     [x2, y2], [x1, y2], [x1, y1]])
    
        r = self.projection.project_geometry(domain_in_crs, crs)
        x1, y1, x2, y2 = r.bounds
        self.set_xlim([x1, x2])
        self.set_ylim([y1, y2])
        
    def set_global(self):
        """
        Set the extent of the Axes to the limits of the projection.

        .. note::

            In some cases where the projection has a limited sensible range
            the ``set_global`` method does not actually make the whole globe
            visible. Instead, the most appropriate extents will be used (e.g.
            Ordnance Survey UK will set the extents to be around the British
            Isles.

        """
        self.set_xlim(self.projection.x_limits)
        self.set_ylim(self.projection.y_limits)

#    def geod_circle_meters(self, lon_0, lat_0, radius, npts=80, **kwargs):
#        # radius is in meters
#        geod = self.projection.as_geodetic()
#
#        az = numpy.linspace(0, 360, npts)
#        lats = numpy.zeros(npts) + lat_0
#        lons = numpy.zeros(npts) + lon_0
#        distances = numpy.zeros(npts) + radius
#
#        lons, lats, _reverse_az = geod.fwd(lons, lats, az, distances, radians=False)
#        ll = numpy.concatenate([lons[:, None], lats[:, None]], 1)
#        from matplotlib.patches import Polygon
#        poly = Polygon(ll, transform=cartopy.prj.PlateCarree(), **kwargs)
#        self.add_patch(poly)
#        return poly
#
#    def gshhs_line(self, outline_color='k', domain=None, resolution='low', **kwargs):
#        # domain is a shapely geometry (Polygon or MultiPolygon)
#        import cartopy.gshhs as gshhs
##        import cartopy.spherical as spherical
#        from matplotlib.collections import PatchCollection, LineCollection
#
#        paths = []
#
#        projection = self.projection
#
#        if domain is None:
#            domain = self.map_domain(ccrs.PlateCarree())
#
#        for points in gshhs.read_gshhc(gshhs.fnames[resolution], poly=False, domain=domain):
#            paths.extend(patch.geos_to_path(shapely.geometry.LineString(points)))
#
##            slinestring = shapely.geometry.LineString(points)
##            projected = projection.project_geometry(slinestring)            
##            paths.extend(patch.geos_to_path(projected))
#
#        collection = PatchCollection([mpatches.PathPatch(pth) for pth in paths],
#                             edgecolor=outline_color, facecolor='none',
#                             transform=ccrs.PlateCarree(),
#                             **kwargs
#                             )
#
#        self.add_collection(collection, autolim=False)

    def stock_img(self, name='ne_shaded'):
        """
        Add a standard image to the map. 
        
        Currently, the only (and default) option is a downsampled version of the Natural Earth
        shaded relief raster. 
        
        """
        # XXX Turn into a dictionary (inside the method)?
#        if img_name == 'bluemarble':
#            source_proj = ccrs.PlateCarree()
#            fname = '/data/local/dataZoo/cartography/raster/blue_marble_720_360.png'
##            fname = '/data/local/dataZoo/cartography/raster/blue_marble_2000_1000.jpg'            
#            img_origin = 'lower'
#            img = imread(fname)
#            img = img[::-1]
#            return self.imshow(img, origin=img_origin, transform=source_proj, extent=[-180, 180, -90, 90])
#        elif img_name == 'bm_high':
#            source_proj = ccrs.PlateCarree()
#            fname = '/data/local/dataZoo/cartography/raster/blue_marble_2000_1000.jpg'
#            img_origin = 'lower'
#            img = imread(fname)
#            return self.imshow(img, origin=img_origin, transform=source_proj, extent=[-180, 180, -90, 90])
        if name == 'ne_shaded':
            import os
            source_proj = ccrs.PlateCarree()
            fname = os.path.join(os.path.dirname(os.path.dirname(__file__)), 
                                 'data', 'raster', 'natural_earth',
                                 '50-natural-earth-1-downsampled.png')
            img_origin = 'lower'
            img = imread(fname)
            img = img[::-1]
            return self.imshow(img, origin=img_origin, transform=source_proj, extent=[-180, 180, -90, 90])
        else:
            raise ValueError('Unknown stock image %r.' % name)

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
                raise ValueError('Expected a projection subclass. Cannot handle a %s in imshow.' % type(transform))

            # XXX adaptive resolution depending on incoming img?
            img, extent = cartopy.img_transform.warp_array(img,
                                                           source_proj=transform,
                                                           source_extent=extent,
                                                           target_proj=self.projection,
                                                           target_res=regrid_shape,
                                                           target_extent=self.get_extent(self.projection),
                                                           )
            # as a workaround to a matplotlib limitation, turn any images which are RGB with a mask into 
            # RGBA images with an alpha channel.
            if isinstance(img, numpy.ma.MaskedArray) and img.shape[2:3] == (3, ):
                old_img = img
                img = numpy.zeros(img.shape[:2] + (4, ))
                img[:, :, 0:3] = old_img
                # put an alpha channel in if the image was masked
                img[:, :, 3] = ~ numpy.any(old_img.mask, axis=2) 
                
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

    def gridlines(self, crs=None, **kwargs):
        """
        Automatically adds gridlines to the axes, in the given coordinate system, at draw time.

        ``**kwargs`` - are passed through to the created :class:`matplotlib.collections.Collection`
                       allowing control of colors and linewidths etc.

        """
        if crs is None:
            crs = ccrs.PlateCarree()
        from cartopy.mpl_integration.gridliner import Gridliner
        gl = Gridliner(self, crs=crs, collection_kwargs=kwargs)
        self._gridliners.append(gl)
        return gl
    
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
                                        facecolor='w', edgecolor='none', zorder= -1,
                                        transform=sct, clip_on=False,
                                        )
        self.background_patch = collection
        # XXX autolim = False
        self.add_patch(collection)


        self.patch.set_facecolor((1, 1, 1, 0))
        self.patch.set_edgecolor((0.5, 0.5, 0.5))
        self.patch.set_linewidth(0.0)


    # mpl 1.2.0rc2 compatibility. To be removed once 1.2 is released
    def contour(self, *args, **kwargs):
        t = kwargs.get('transform', None)
        # Keep this bit - even at mpl v1.2
        if t is None:
            t = self.projection
        if hasattr(t, '_as_mpl_transform'):
            kwargs['transform'] = t._as_mpl_transform(self)
        return matplotlib.axes.Axes.contour(self, *args, **kwargs)
    
    # mpl 1.2.0rc2 compatibility. To be removed once 1.2 is released
    def contourf(self, *args, **kwargs):
        t = kwargs.get('transform', None)
        # Keep this bit - even at mpl v1.2
        if t is None:
            t = self.projection
        if hasattr(t, '_as_mpl_transform'):
            kwargs['transform'] = t._as_mpl_transform(self)
        return matplotlib.axes.Axes.contourf(self, *args, **kwargs)
    
    # mpl 1.2.0rc2 compatibility. To be removed once 1.2 is released
    def scatter(self, *args, **kwargs):
        t = kwargs.get('transform', None)
        # Keep this bit - even at mpl v1.2
        if t is None:
            t = self.projection
        if hasattr(t, '_as_mpl_transform'):
            kwargs['transform'] = t._as_mpl_transform(self)
        return matplotlib.axes.Axes.scatter(self, *args, **kwargs)

    # mpl 1.2.0rc2 compatibility. To be removed once 1.2 is released
    def pcolormesh(self, *args, **kwargs):
        import warnings
        import numpy as np
        import numpy.ma as ma
        import matplotlib as mpl
        import matplotlib.cbook as cbook
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm
        from matplotlib import docstring
        import matplotlib.transforms as transforms
        import matplotlib.artist as artist
        from matplotlib.artist import allow_rasterization
        import matplotlib.backend_bases as backend_bases
        import matplotlib.path as mpath
        import matplotlib.mlab as mlab
        import matplotlib.collections as mcoll
        
        if not self._hold: self.cla()

        alpha = kwargs.pop('alpha', None)
        norm = kwargs.pop('norm', None)
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        shading = kwargs.pop('shading', 'flat').lower()
        antialiased = kwargs.pop('antialiased', False)
        kwargs.setdefault('edgecolors', 'None')

        X, Y, C = self._pcolorargs('pcolormesh', *args)
        Ny, Nx = X.shape

        # convert to one dimensional arrays
        if shading != 'gouraud':
            C = ma.ravel(C[0:Ny-1, 0:Nx-1]) # data point in each cell is value at
                                            # lower left corner
        else:
            C = C.ravel()
        X = X.ravel()
        Y = Y.ravel()

        coords = np.zeros(((Nx * Ny), 2), dtype=float)
        coords[:, 0] = X
        coords[:, 1] = Y

        collection = mcoll.QuadMesh(
            Nx - 1, Ny - 1, coords,
            antialiased=antialiased, shading=shading, **kwargs)
        collection.set_alpha(alpha)
        collection.set_array(C)
        if norm is not None: assert(isinstance(norm, mcolors.Normalize))
        collection.set_cmap(cmap)
        collection.set_norm(norm)
        collection.set_clim(vmin, vmax)
        collection.autoscale_None()

        self.grid(False)

        # Transform from native to data coordinates?
        t = collection._transform
        if (not isinstance(t, mtransforms.Transform)
            and hasattr(t, '_as_mpl_transform')):
            t = t._as_mpl_transform(self.axes)

        if t and any(t.contains_branch_seperately(self.transData)):
            trans_to_data = t - self.transData
            pts = np.vstack([X, Y]).T.astype(np.float)
            transformed_pts = trans_to_data.transform(pts)
            X = transformed_pts[..., 0]
            Y = transformed_pts[..., 1]

            # XXX Not a mpl 1.2 thing...
            no_inf = (X != np.inf) & (Y != np.inf)
            X = X[no_inf]
            Y = Y[no_inf]

        minx = np.amin(X)
        maxx = np.amax(X)
        miny = np.amin(Y)
        maxy = np.amax(Y)
        
        corners = (minx, miny), (maxx, maxy)
        self.update_datalim( corners)
        self.autoscale_view()
        self.add_collection(collection)
        return collection
    
    # mpl 1.2.0rc2 compatibility. To be removed once 1.2 is released
    def pcolor(self, *args, **kwargs):
        t = kwargs.get('transform', None)
        if t is None:
            kwargs['transform'] = self.projection

        import warnings
        import numpy as np
        import numpy.ma as ma
        import matplotlib as mpl
        import matplotlib.cbook as cbook
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm
        from matplotlib import docstring
        import matplotlib.transforms as transforms
        import matplotlib.artist as artist
        from matplotlib.artist import allow_rasterization
        import matplotlib.backend_bases as backend_bases
        import matplotlib.path as mpath
        import matplotlib.mlab as mlab
        import matplotlib.collections as mcoll

        if not self._hold: self.cla()

        alpha = kwargs.pop('alpha', None)
        norm = kwargs.pop('norm', None)
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        shading = kwargs.pop('shading', 'flat')

        X, Y, C = self._pcolorargs('pcolor', *args)
        Ny, Nx = X.shape

        # convert to MA, if necessary.
        C = ma.asarray(C)
        X = ma.asarray(X)
        Y = ma.asarray(Y)
        mask = ma.getmaskarray(X)+ma.getmaskarray(Y)
        xymask = mask[0:-1,0:-1]+mask[1:,1:]+mask[0:-1,1:]+mask[1:,0:-1]
        # don't plot if C or any of the surrounding vertices are masked.
        mask = ma.getmaskarray(C)[0:Ny-1,0:Nx-1]+xymask

        newaxis = np.newaxis
        compress = np.compress

        ravelmask = (mask==0).ravel()
        X1 = compress(ravelmask, ma.filled(X[0:-1,0:-1]).ravel())
        Y1 = compress(ravelmask, ma.filled(Y[0:-1,0:-1]).ravel())
        X2 = compress(ravelmask, ma.filled(X[1:,0:-1]).ravel())
        Y2 = compress(ravelmask, ma.filled(Y[1:,0:-1]).ravel())
        X3 = compress(ravelmask, ma.filled(X[1:,1:]).ravel())
        Y3 = compress(ravelmask, ma.filled(Y[1:,1:]).ravel())
        X4 = compress(ravelmask, ma.filled(X[0:-1,1:]).ravel())
        Y4 = compress(ravelmask, ma.filled(Y[0:-1,1:]).ravel())
        npoly = len(X1)

        xy = np.concatenate((X1[:,newaxis], Y1[:,newaxis],
                             X2[:,newaxis], Y2[:,newaxis],
                             X3[:,newaxis], Y3[:,newaxis],
                             X4[:,newaxis], Y4[:,newaxis],
                             X1[:,newaxis], Y1[:,newaxis]),
                             axis=1)
        verts = xy.reshape((npoly, 5, 2))

        C = compress(ravelmask, ma.filled(C[0:Ny-1,0:Nx-1]).ravel())

        linewidths = (0.25,)
        if 'linewidth' in kwargs:
            kwargs['linewidths'] = kwargs.pop('linewidth')
        kwargs.setdefault('linewidths', linewidths)

        if shading == 'faceted':
            edgecolors = 'k',
        else:
            edgecolors = 'none'
        if 'edgecolor' in kwargs:
            kwargs['edgecolors'] = kwargs.pop('edgecolor')
        ec = kwargs.setdefault('edgecolors', edgecolors)

        # aa setting will default via collections to patch.antialiased
        # unless the boundary is not stroked, in which case the
        # default will be False; with unstroked boundaries, aa
        # makes artifacts that are often disturbing.
        if 'antialiased' in kwargs:
            kwargs['antialiaseds'] = kwargs.pop('antialiased')
        if 'antialiaseds' not in kwargs and ec.lower() == "none":
                kwargs['antialiaseds'] = False


        collection = mcoll.PolyCollection(verts, **kwargs)

        collection.set_alpha(alpha)
        collection.set_array(C)
        if norm is not None: assert(isinstance(norm, mcolors.Normalize))
        collection.set_cmap(cmap)
        collection.set_norm(norm)
        collection.set_clim(vmin, vmax)
        collection.autoscale_None()
        self.grid(False)

        x = X.compressed()
        y = Y.compressed()
        
        # Transform from native to data coordinates?
        t = collection._transform
        if (not isinstance(t, mtransforms.Transform)
            and hasattr(t, '_as_mpl_transform')):
            t = t._as_mpl_transform(self.axes)

        if t and any(t.contains_branch_seperately(self.transData)):
            trans_to_data = t - self.transData
            pts = np.vstack([x, y]).T.astype(np.float)
            transformed_pts = trans_to_data.transform(pts)
            x = transformed_pts[..., 0]
            y = transformed_pts[..., 1]
            
            # XXX Not a mpl 1.2 thing...
            no_inf = (x != np.inf) & (y != np.inf)
            x = x[no_inf]
            y = y[no_inf]
        
        minx = np.amin(x)
        maxx = np.amax(x)
        miny = np.amin(y)
        maxy = np.amax(y)

        corners = (minx, miny), (maxx, maxy)
        self.update_datalim( corners)
        self.autoscale_view()
        self.add_collection(collection)
        return collection


# alias GeoAxes - NOTE: THIS WAS NOT IN v0.4.0rc1
GenericProjectionAxes = GeoAxes


class SimpleClippedTransform(mtransforms.Transform):
    """
    Transforms the values using a pre transform, clips them, then post transforms them.
    
    This transform should not be widely used, but is useful for transforming a background patch 
    and clipping the patch to a desired extent.
    
    """
    input_dims = 2
    output_dims = 2
    has_inverse = True
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

    def inverted(self):
        return (self.pre_clip_transform + self.post_clip_transform).inverted()


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
        result = self.post_clip_transform.transform(new_vals)

        if single_point:
            return result[0, :]
        else:
            return result
