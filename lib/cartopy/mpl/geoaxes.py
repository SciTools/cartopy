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
"""
This module defines the :class:`GeoAxes` class, for use with matplotlib.

When a matplotlib figure contains a GeoAxes the plotting commands can transform
plot results from source coordinates to the GeoAxes' target projection.

"""

from __future__ import (absolute_import, division, print_function)

import collections
import contextlib
import warnings
import weakref

import matplotlib as mpl
import matplotlib.artist
import matplotlib.axes
from matplotlib.image import imread
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import numpy as np
import numpy.ma as ma
import shapely.geometry as sgeom

from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature
import cartopy.img_transform
from cartopy.mpl.clip_path import clip_path
import cartopy.mpl.feature_artist as feature_artist
import cartopy.mpl.patch as cpatch
from cartopy.mpl.slippy_image_artist import SlippyImageArtist
from cartopy.vector_transform import vector_scalar_to_grid


assert matplotlib.__version__ >= '1.3', ('Cartopy is only supported with '
                                         'matplotlib 1.3 or greater.')


_PATH_TRANSFORM_CACHE = weakref.WeakKeyDictionary()
"""
A nested mapping from path, source CRS, and target projection to the
resulting transformed paths::

    {path: {(source_crs, target_projection): list_of_paths}}

Provides a significant performance boost for contours which, at
matplotlib 1.2.0 called transform_path_non_affine twice unnecessarily.

"""

_BACKG_IMG_CACHE = {}
"""
A dictionary of pre-loaded images for large background images, kept as a
dictionary so that large images are loaded only once.
"""

_USER_BG_IMGS = {}
"""
A dictionary of background images in the directory specified by the
CARTOPY_USER_BACKGROUNDS environment variable.
"""


# XXX call this InterCRSTransform
class InterProjectionTransform(mtransforms.Transform):
    """
    Transforms coordinates from the source_projection to
    the ``target_projection``.

    """
    input_dims = 2
    output_dims = 2
    is_separable = False
    has_inverse = True

    def __init__(self, source_projection, target_projection):
        """
        Create the transform object from the given projections.

        Args:

            * source_projection - A :class:`~cartopy.crs.CRS`.
            * target_projection - A :class:`~cartopy.crs.CRS`.

        """
        # assert target_projection is cartopy.crs.Projection
        # assert source_projection is cartopy.crs.CRS
        self.source_projection = source_projection
        self.target_projection = target_projection
        mtransforms.Transform.__init__(self)

    def __repr__(self):
        return ('< {!s} {!s} -> {!s} >'.format(self.__class__.__name__,
                                               self.source_projection,
                                               self.target_projection))

    def transform_non_affine(self, xy):
        """
        Transforms from source to target coordinates.

        Args:

            * xy - An (n,2) array of points in source coordinates.

        Returns:

            * An (n,2) array of transformed points in target coordinates.

        """
        prj = self.target_projection
        if isinstance(xy, np.ndarray):
            return prj.transform_points(self.source_projection,
                                        xy[:, 0], xy[:, 1])[:, 0:2]
        else:
            x, y = xy
            x, y = prj.transform_point(x, y, self.source_projection)
            return x, y

    def transform_path_non_affine(self, src_path):
        """
        Transforms from source to target coordinates.

        Caches results, so subsequent calls with the same *src_path* argument
        (and the same source and target projections) are faster.

        Args:

            * src_path - A matplotlib :class:`~matplotlib.path.Path` object
                         with vertices in source coordinates.

        Returns

            * A matplotlib :class:`~matplotlib.path.Path` with vertices
              in target coordinates.

        """
        mapping = _PATH_TRANSFORM_CACHE.get(src_path)
        if mapping is not None:
            key = (self.source_projection, self.target_projection)
            result = mapping.get(key)
            if result is not None:
                return result

        # Allow the vertices to be quickly transformed, if
        # quick_vertices_transform allows it.
        new_vertices = self.target_projection.quick_vertices_transform(
            src_path.vertices, self.source_projection)
        if new_vertices is not None:
            if new_vertices is src_path.vertices:
                return src_path
            else:
                return mpath.Path(new_vertices, src_path.codes)

        if src_path.vertices.shape == (1, 2):
            return mpath.Path(self.transform(src_path.vertices))

        transformed_geoms = []
        # Check whether this transform has the "force_path_ccw" attribute set.
        # This is a cartopy extension to the Transform API to allow finer
        # control of Path orientation handling (Path ordering is not important
        # in matplotlib, but is in Cartopy).
        geoms = cpatch.path_to_geos(src_path,
                                    getattr(self, 'force_path_ccw', False))

        for geom in geoms:
            proj_geom = self.target_projection.project_geometry(
                geom, self.source_projection)
            transformed_geoms.append(proj_geom)

        if not transformed_geoms:
            result = mpath.Path(np.empty([0, 2]))
        else:
            paths = cpatch.geos_to_path(transformed_geoms)
            if not paths:
                return mpath.Path(np.empty([0, 2]))
            points, codes = list(zip(*[cpatch.path_segments(path,
                                                            curves=False,
                                                            simplify=False)
                                       for path in paths]))
            result = mpath.Path(np.concatenate(points, 0),
                                np.concatenate(codes))

        # store the result in the cache for future performance boosts
        key = (self.source_projection, self.target_projection)
        if mapping is None:
            _PATH_TRANSFORM_CACHE[src_path] = {key: result}
        else:
            mapping[key] = result

        return result

    def inverted(self):
        """
        Return a matplotlib :class:`~matplotlib.transforms.Transform`
        from target to source coordinates.

        """
        return InterProjectionTransform(self.target_projection,
                                        self.source_projection)


class GeoAxes(matplotlib.axes.Axes):
    """
    A subclass of :class:`matplotlib.axes.Axes` which represents a
    map :class:`~cartopy.crs.Projection`.

    This class replaces the matplotlib :class:`~matplotlib.axes.Axes` class
    when created with the *projection* keyword. For example::

        # Set up a standard map for latlon data.
        geo_axes = pyplot.axes(projection=cartopy.crs.PlateCarree())

        # Set up an OSGB map.
        geo_axes = pyplot.subplot(2, 2, 1, projection=cartopy.crs.OSGB())

    When a source projection is provided to one of it's plotting methods,
    using the *transform* keyword, the standard matplotlib plot result is
    transformed from source coordinates to the target projection. For example::

        # Plot latlon data on an OSGB map.
        pyplot.axes(projection=cartopy.crs.OSGB())
        pyplot.contourf(x, y, data, transform=cartopy.crs.PlateCarree())

    """
    def __init__(self, *args, **kwargs):
        """
        Create a GeoAxes object using standard matplotlib
        :class:`~matplotlib.axes.Axes` args and kwargs.

        Kwargs:

            * map_projection - The target :class:`~cartopy.crs.Projection` of
                               this Axes object.

        All other args and keywords are passed through to
        :class:`matplotlib.axes.Axes`.

        """
        self.projection = kwargs.pop('map_projection')
        """The :class:`cartopy.crs.Projection` of this GeoAxes."""

        self.outline_patch = None
        """The patch that provides the line bordering the projection."""

        self.background_patch = None
        """The patch that provides the filled background of the projection."""

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
        if hasattr(factory, 'image_for_domain'):
            # XXX TODO: Needs deprecating.
            self.img_factories.append([factory, args, kwargs])
        else:
            # Args and kwargs not allowed.
            assert not bool(args) and not bool(kwargs)
            image = factory
            try:
                super(GeoAxes, self).add_image(image)
            except AttributeError:
                # If add_image method doesn't exist (only available from
                # v1.4 onwards) we implement it ourselves.
                self._set_artist_props(image)
                self.images.append(image)
                image._remove_method = lambda h: self.images.remove(h)
            return image

    @contextlib.contextmanager
    def hold_limits(self, hold=True):
        """
        Keep track of the original view and data limits for the life of this
        context manager, optionally reverting any changes back to the original
        values after the manager exits.

        Parameters
        ----------
        hold : bool (default True)
            Whether to revert the data and view limits after the context
            manager exits.

        """
        data_lim = self.dataLim.frozen().get_points()
        view_lim = self.viewLim.frozen().get_points()
        other = (self.ignore_existing_data_limits,
                 self._autoscaleXon, self._autoscaleYon)
        try:
            yield
        finally:
            if hold:
                self.dataLim.set_points(data_lim)
                self.viewLim.set_points(view_lim)
                (self.ignore_existing_data_limits,
                    self._autoscaleXon, self._autoscaleYon) = other

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer=None, inframe=False):
        """
        Extends the standard behaviour of :func:`matplotlib.axes.Axes.draw`.

        Draws grid lines and image factory results before invoking standard
        matplotlib drawing. A global range is used if no limits have yet
        been set.

        """
        # If data has been added (i.e. autoscale hasn't been turned off)
        # then we should autoscale the view.
        if self.get_autoscale_on() and self.ignore_existing_data_limits:
            self.autoscale_view()

        if self.outline_patch.reclip or self.background_patch.reclip:
            clipped_path = clip_path(self.outline_patch.orig_path,
                                     self.viewLim)
            self.outline_patch._path = clipped_path
            self.background_patch._path = clipped_path

        for gl in self._gridliners:
            gl._draw_gridliner(background_patch=self.background_patch)
        self._gridliners = []

        # XXX This interface needs a tidy up:
        #       image drawing on pan/zoom;
        #       caching the resulting image;
        #       buffering the result by 10%...;
        if not self._done_img_factory:
            for factory, args, kwargs in self.img_factories:
                img, extent, origin = factory.image_for_domain(
                    self._get_extent_geom(factory.crs), args[0])
                self.imshow(img, extent=extent, origin=origin,
                            transform=factory.crs, *args[1:], **kwargs)
        self._done_img_factory = True

        return matplotlib.axes.Axes.draw(self, renderer=renderer,
                                         inframe=inframe)

    def __str__(self):
        return '< GeoAxes: %s >' % self.projection

    def cla(self):
        """Clears the current axes and adds boundary lines."""
        result = matplotlib.axes.Axes.cla(self)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        # Enable tight autoscaling.
        self._tight = True
        self.set_aspect('equal')

        with self.hold_limits():
            self._boundary()

        # XXX consider a margin - but only when the map is not global...
        # self._xmargin = 0.15
        # self._ymargin = 0.15

        self.dataLim.intervalx = self.projection.x_limits
        self.dataLim.intervaly = self.projection.y_limits

        return result

    def format_coord(self, x, y):
        """Return a string formatted for the matplotlib GUI status bar."""
        lon, lat = ccrs.Geodetic().transform_point(x, y, self.projection)

        ns = 'N' if lat >= 0.0 else 'S'
        ew = 'E' if lon >= 0.0 else 'W'

        return u'%.4g, %.4g (%f\u00b0%s, %f\u00b0%s)' % (x, y, abs(lat),
                                                         ns, abs(lon), ew)

    def coastlines(self, resolution='110m', color='black', **kwargs):
        """
        Adds coastal **outlines** to the current axes from the Natural Earth
        "coastline" shapefile collection.

        Kwargs:

            * resolution - a named resolution to use from the Natural Earth
                           dataset. Currently can be one of "110m", "50m", and
                           "10m".

        .. note::

            Currently no clipping is done on the coastlines before adding
            them to the axes. This means, if very high resolution coastlines
            are being used, performance is likely to be severely effected.
            This should be resolved transparently by v0.5.

        """
        kwargs['edgecolor'] = color
        kwargs['facecolor'] = 'none'
        feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline',
                                                      resolution, **kwargs)
        return self.add_feature(feature)

    def tissot(self, rad_km=5e5, lons=None, lats=None, n_samples=80, **kwargs):
        """
        Adds Tissot's indicatrices to the axes.

        Kwargs:

            * rad_km - The radius in km of the the circles to be drawn.

            * lons - A numpy.ndarray, list or tuple of longitude values that
                     locate the centre of each circle. Specifying more than one
                     dimension allows individual points to be drawn whereas a
                     1D array produces a grid of points.

            * lats - A numpy.ndarray, list or tuple of latitude values that
                     that locate the centre of each circle. See lons.

            * n_samples - Integer number of points sampled around the
                          circumference of each circle.

        **kwargs are passed through to `class:ShapelyFeature`.

        """
        from cartopy import geodesic

        geod = geodesic.Geodesic()
        geoms = []

        if lons is None:
            lons = np.linspace(-180, 180, 6, endpoint=False)
        else:
            lons = np.asarray(lons)
        if lats is None:
            lats = np.linspace(-80, 80, 6)
        else:
            lats = np.asarray(lats)

        if lons.ndim == 1 or lats.ndim == 1:
            lons, lats = np.meshgrid(lons, lats)
        lons, lats = lons.flatten(), lats.flatten()

        if lons.shape != lats.shape:
            raise ValueError('lons and lats must have the same shape.')

        for i in range(len(lons)):
                circle = geod.circle(lons[i], lats[i], rad_km,
                                     n_samples=n_samples)
                geoms.append(sgeom.Polygon(circle))

        feature = cartopy.feature.ShapelyFeature(geoms, ccrs.Geodetic(),
                                                 **kwargs)
        return self.add_feature(feature)

    def natural_earth_shp(self, name='land', resolution='110m',
                          category='physical', **kwargs):
        """
        Adds the geometries from the specified Natural Earth shapefile to the
        Axes as a :class:`~matplotlib.collections.PathCollection`.

        ``**kwargs`` are passed through to the
        :class:`~matplotlib.collections.PathCollection` constructor.

        Returns the created :class:`~matplotlib.collections.PathCollection`.

        .. note::

            Currently no clipping is done on the geometries before adding them
            to the axes. This means, if very high resolution geometries are
            being used, performance is likely to be severely effected. This
            should be resolved transparently by v0.5.

        """
        warnings.warn('This method has been deprecated.'
                      ' Please use `add_feature` instead.')
        kwargs.setdefault('edgecolor', 'face')
        kwargs.setdefault('facecolor', cartopy.feature.COLORS['land'])
        feature = cartopy.feature.NaturalEarthFeature(category, name,
                                                      resolution, **kwargs)
        return self.add_feature(feature)

    def add_feature(self, feature, **kwargs):
        """
        Adds the given :class:`~cartopy.feature.Feature` instance to the axes.

        Args:

        * feature:
            An instance of :class:`~cartopy.feature.Feature`.

        Kwargs:
            Keyword arguments to be used when drawing the feature. This allows
            standard matplotlib control over aspects such as 'facecolor',
            'alpha', etc.

        Returns:
            * A :class:`cartopy.mpl.feature_artist.FeatureArtist`
              instance responsible for drawing the feature.

        """
        # Instantiate an artist to draw the feature and add it to the axes.
        artist = feature_artist.FeatureArtist(feature, **kwargs)
        return self.add_artist(artist)

    def add_geometries(self, geoms, crs, **kwargs):
        """
        Add the given shapely geometries (in the given crs) to the axes.

        Args:

        * geoms:
            A collection of shapely geometries.
        * crs:
            The cartopy CRS in which the provided geometries are defined.

        Kwargs:
            Keyword arguments to be used when drawing this feature.

        Returns:
             A :class:`cartopy.mpl.feature_artist.FeatureArtist`
             instance responsible for drawing the geometries.

        """
        feature = cartopy.feature.ShapelyFeature(geoms, crs, **kwargs)
        return self.add_feature(feature)

    def get_extent(self, crs=None):
        """
        Get the extent (x0, x1, y0, y1) of the map in the given coordinate
        system.

        If no crs is given, the returned extents' coordinate system will be
        the CRS of this Axes.

        """
        p = self._get_extent_geom(crs)
        r = p.bounds
        x1, y1, x2, y2 = r
        return x1, x2, y1, y2

    def _get_extent_geom(self, crs=None):
        # Perform the calculations for get_extent(), which just repackages it.
        with self.hold_limits():
            if self.get_autoscale_on():
                self.autoscale_view()
            [x1, y1], [x2, y2] = self.viewLim.get_points()

        domain_in_src_proj = sgeom.Polygon([[x1, y1], [x2, y1],
                                            [x2, y2], [x1, y2],
                                            [x1, y1]])

        # Determine target projection based on requested CRS.
        if crs is None:
            proj = self.projection
        elif isinstance(crs, ccrs.Projection):
            proj = crs
        else:
            # Attempt to select suitable projection for
            # non-projection CRS.
            if isinstance(crs, ccrs.RotatedGeodetic):
                proj = ccrs.RotatedPole(crs.proj4_params['lon_0'] - 180,
                                        crs.proj4_params['o_lat_p'])
                warnings.warn('Approximating coordinate system {!r} with a '
                              'RotatedPole projection.'.format(crs))
            elif hasattr(crs, 'is_geodetic') and crs.is_geodetic():
                proj = ccrs.PlateCarree(crs.globe)
                warnings.warn('Approximating coordinate system {!r} with the '
                              'PlateCarree projection.'.format(crs))
            else:
                raise ValueError('Cannot determine extent in'
                                 ' coordinate system {!r}'.format(crs))

        # Calculate intersection with boundary and project if necesary.
        boundary_poly = sgeom.Polygon(self.projection.boundary)
        if proj != self.projection:
            # Erode boundary by threshold to avoid transform issues.
            # This is a workaround for numerical issues at the boundary.
            eroded_boundary = boundary_poly.buffer(-self.projection.threshold)
            geom_in_src_proj = eroded_boundary.intersection(
                domain_in_src_proj)
            geom_in_crs = proj.project_geometry(geom_in_src_proj,
                                                self.projection)
        else:
            geom_in_crs = boundary_poly.intersection(domain_in_src_proj)

        return geom_in_crs

    def set_extent(self, extents, crs=None):
        """
        Set the extent (x0, x1, y0, y1) of the map in the given
        coordinate system.

        If no crs is given, the extents' coordinate system will be assumed
        to be the Geodetic version of this axes' projection.

        """
        # TODO: Implement the same semantics as plt.xlim and
        # plt.ylim - allowing users to set None for a minimum and/or
        # maximum value
        x1, x2, y1, y2 = extents
        domain_in_crs = sgeom.polygon.LineString([[x1, y1], [x2, y1],
                                                  [x2, y2], [x1, y2],
                                                  [x1, y1]])

        projected = None

        # Sometimes numerical issues cause the projected vertices of the
        # requested extents to appear outside the projection domain.
        # This results in an empty geometry, which has an empty `bounds`
        # tuple, which causes an unpack error.
        # This workaround avoids using the projection when the requested
        # extents are obviously the same as the projection domain.
        try_workaround = ((crs is None and
                           isinstance(self.projection, ccrs.PlateCarree)) or
                          crs == self.projection)
        if try_workaround:
            boundary = self.projection.boundary
            if boundary.equals(domain_in_crs):
                projected = boundary

        if projected is None:
            projected = self.projection.project_geometry(domain_in_crs, crs)
        try:
            # This might fail with an unhelpful error message ('need more
            # than 0 values to unpack') if the specified extents fall outside
            # the projection extents, so try and give a better error message.
            x1, y1, x2, y2 = projected.bounds
        except ValueError:
            msg = ('Failed to determine the required bounds in projection '
                   'coordinates. Check that the values provided are within '
                   'the valid range (x_limits=[{xlim[0]}, {xlim[1]}], '
                   'y_limits=[{ylim[0]}, {ylim[1]}]).')
            raise ValueError(msg.format(xlim=self.projection.x_limits,
                                        ylim=self.projection.y_limits))
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

    def set_xticks(self, ticks, minor=False, crs=None):
        """
        Set the x ticks.

        Args:

            * ticks - list of floats denoting the desired position of x ticks.

        Kwargs:

            * minor - boolean flag indicating whether the ticks should be minor
                      ticks i.e. small and unlabelled (default is False).

            * crs - An instance of :class:`~cartopy.crs.CRS` indicating the
                    coordinate system of the provided tick values. If no
                    coordinate system is specified then the values are assumed
                    to be in the coordinate system of the projection.
                    Only transformations from one rectangular coordinate system
                    to another rectangular coordinate system are supported.

        .. note::

            This interface is subject to change whilst functionality is added
            to support other map projections.

        """
        # Project ticks if crs differs from axes' projection
        if crs is not None and crs != self.projection:
            if not isinstance(crs, (ccrs._RectangularProjection,
                                    ccrs.Mercator)) or \
                    not isinstance(self.projection,
                                   (ccrs._RectangularProjection,
                                    ccrs.Mercator)):
                raise RuntimeError('Cannot handle non-rectangular coordinate '
                                   'systems.')
            proj_xyz = self.projection.transform_points(crs,
                                                        np.asarray(ticks),
                                                        np.zeros(len(ticks)))
            xticks = proj_xyz[..., 0]
        else:
            xticks = ticks

        # Switch on drawing of x axis
        self.xaxis.set_visible(True)

        return super(GeoAxes, self).set_xticks(xticks, minor)

    def set_yticks(self, ticks, minor=False, crs=None):
        """
        Set the y ticks.

        Args:

            * ticks - list of floats denoting the desired position of y ticks.

        Kwargs:

            * minor - boolean flag indicating whether the ticks should be minor
                      ticks i.e. small and unlabelled (default is False).

            * crs - An instance of :class:`~cartopy.crs.CRS` indicating the
                    coordinate system of the provided tick values. If no
                    coordinate system is specified then the values are assumed
                    to be in the coordinate system of the projection.
                    Only transformations from one rectangular coordinate system
                    to another rectangular coordinate system are supported.

        .. note::

            This interface is subject to change whilst functionality is added
            to support other map projections.

        """
        # Project ticks if crs differs from axes' projection
        if crs is not None and crs != self.projection:
            if not isinstance(crs, (ccrs._RectangularProjection,
                                    ccrs.Mercator)) or \
                    not isinstance(self.projection,
                                   (ccrs._RectangularProjection,
                                    ccrs.Mercator)):
                raise RuntimeError('Cannot handle non-rectangular coordinate '
                                   'systems.')
            proj_xyz = self.projection.transform_points(crs,
                                                        np.zeros(len(ticks)),
                                                        np.asarray(ticks))
            yticks = proj_xyz[..., 1]
        else:
            yticks = ticks

        # Switch on drawing of y axis
        self.yaxis.set_visible(True)

        return super(GeoAxes, self).set_yticks(yticks, minor)

    def stock_img(self, name='ne_shaded'):
        """
        Add a standard image to the map.

        Currently, the only (and default) option is a downsampled version of
        the Natural Earth shaded relief raster.

        """
        if name == 'ne_shaded':
            import os
            source_proj = ccrs.PlateCarree()
            fname = os.path.join(config["repo_data_dir"],
                                 'raster', 'natural_earth',
                                 '50-natural-earth-1-downsampled.png')

            return self.imshow(imread(fname), origin='upper',
                               transform=source_proj,
                               extent=[-180, 180, -90, 90])
        else:
            raise ValueError('Unknown stock image %r.' % name)

    def background_img(self, name='ne_shaded', resolution='low', extent=None,
                       cache=False):
        """
        Adds a background image to the map, from a selection of pre-prepared
        images held in a directory specified by the CARTOPY_USER_BACKGROUNDS
        environment variable. That directory is checked with
        func:`self.read_user_background_images` and needs to contain a JSON
        file which defines for the image metadata.

        Kwargs:

            * name - the name of the image to read according to the contents
                     of the JSON file. A typical file might have, for instance:
                     'ne_shaded' : Natural Earth Shaded Relief
                     'ne_grey' : Natural Earth Grey Earth

            * resolution - the resolution of the image to read, according to
                           the contents of the JSON file. A typical file might
                           have the following for each name of the image:
                           'low', 'med', 'high', 'vhigh', 'full'.

            * extent - using a high resolution background image, zoomed into
                       a small area, will take a very long time to render as
                       the image is prepared globally, even though only a small
                       area is used. Adding the extent will only render a
                       particular geographic region. Specified as
                       [longitude start, longitude end,
                        latitude start, latitude end].

                       e.g. [-11, 3, 48, 60] for the UK
                       or [167.0, 193.0, 47.0, 68.0] to cross the date line.

            * cache - logical flag as to whether or not to cache the loaded
                      images into memory. The images are stored before the
                      extent is used.
        """
        # read in the user's background image directory:
        if len(_USER_BG_IMGS) == 0:
            self.read_user_background_images()
        import os
        bgdir = os.getenv('CARTOPY_USER_BACKGROUNDS')
        if bgdir is None:
            bgdir = os.path.join(config["repo_data_dir"],
                                 'raster', 'natural_earth')
        # now get the filename we want to use:
        try:
            fname = _USER_BG_IMGS[name][resolution]
        except KeyError:
            msg = ('Image "{}" and resolution "{}" are not present in '
                   'the user background image metadata in directory "{}"')
            raise ValueError(msg.format(name, resolution, bgdir))
        # Now obtain the image data from file or cache:
        fpath = os.path.join(bgdir, fname)
        if cache:
            if fname in _BACKG_IMG_CACHE:
                img = _BACKG_IMG_CACHE[fname]
            else:
                img = imread(fpath)
                _BACKG_IMG_CACHE[fname] = img
        else:
            img = imread(fpath)
        if len(img.shape) == 2:
            # greyscale images are only 2-dimensional, so need replicating
            # to 3 colour channels:
            img = np.repeat(img[:, :, np.newaxis], 3, axis=2)
        # now get the projection from the metadata:
        if _USER_BG_IMGS[name]['__projection__'] == 'PlateCarree':
            # currently only PlateCarree is defined:
            source_proj = ccrs.PlateCarree()
        else:
            raise NotImplementedError('Background image projection undefined')

        if extent is None:
            # not specifying an extent, so return all of it:
            return self.imshow(img, origin='upper',
                               transform=source_proj,
                               extent=[-180, 180, -90, 90])
        else:
            # return only a subset of the image:
            # set up coordinate arrays:
            d_lat = 180.0 / img.shape[0]
            d_lon = 360.0 / img.shape[1]
            # latitude starts at 90N for this image:
            lat_pts = (np.arange(img.shape[0]) * -d_lat - (d_lat / 2.0)) + 90.0
            lon_pts = (np.arange(img.shape[1]) * d_lon + (d_lon / 2.0)) - 180.0

            # which points are in range:
            lat_in_range = np.logical_and(lat_pts >= extent[2],
                                          lat_pts <= extent[3])
            if extent[0] < 180 and extent[1] > 180:
                # we have a region crossing the dateline
                # this is the westerly side of the input image:
                lon_in_range1 = np.logical_and(lon_pts >= extent[0],
                                               lon_pts <= 180.0)
                img_subset1 = img[lat_in_range, :, :][:, lon_in_range1, :]
                # and the eastward half:
                lon_in_range2 = lon_pts + 360. <= extent[1]
                img_subset2 = img[lat_in_range, :, :][:, lon_in_range2, :]
                # now join them up:
                img_subset = np.concatenate((img_subset1, img_subset2), axis=1)
                # now define the extent for output that matches those points:
                ret_extent = [lon_pts[lon_in_range1][0] - d_lon / 2.0,
                              lon_pts[lon_in_range2][-1] + d_lon / 2.0 + 360,
                              lat_pts[lat_in_range][-1] - d_lat / 2.0,
                              lat_pts[lat_in_range][0] + d_lat / 2.0]
            else:
                # not crossing the dateline, so just find the region:
                lon_in_range = np.logical_and(lon_pts >= extent[0],
                                              lon_pts <= extent[1])
                img_subset = img[lat_in_range, :, :][:, lon_in_range, :]
                # now define the extent for output that matches those points:
                ret_extent = [lon_pts[lon_in_range][0] - d_lon / 2.0,
                              lon_pts[lon_in_range][-1] + d_lon / 2.0,
                              lat_pts[lat_in_range][-1] - d_lat / 2.0,
                              lat_pts[lat_in_range][0] + d_lat / 2.0]

            return self.imshow(img_subset, origin='upper',
                               transform=source_proj,
                               extent=ret_extent)

    def read_user_background_images(self, verify=True):
        """
        Reads the metadata in the specified CARTOPY_USER_BACKGROUNDS
        environment variable to populate the dictionaries for background_img.

        If CARTOPY_USER_BACKGROUNDS is not set then by default the image in
        lib/cartopy/data/raster/natural_earth/ will be made available.

        The metadata should be a standard JSON file which specifies a two
        level dictionary. The first level is the image type.
        For each image type there must be the fields:
        __comment__, __source__ and __projection__
        and then an element giving the filename for each resolution.

        An example JSON file can be found at:
        lib/cartopy/data/raster/natural_earth/images.json

        """
        import os
        import json

        bgdir = os.getenv('CARTOPY_USER_BACKGROUNDS')
        if bgdir is None:
            bgdir = os.path.join(config["repo_data_dir"],
                                 'raster', 'natural_earth')
        json_file = os.path.join(bgdir, 'images.json')

        with open(json_file, 'r') as js_obj:
            dict_in = json.load(js_obj)
        for img_type in dict_in:
            _USER_BG_IMGS[img_type] = dict_in[img_type]

        if verify:
            required_info = ['__comment__', '__source__', '__projection__']
            for img_type in _USER_BG_IMGS:
                if img_type == '__comment__':
                    # the top level comment doesn't need verifying:
                    pass
                else:
                    # check that this image type has the required info:
                    for required in required_info:
                        if required not in _USER_BG_IMGS[img_type]:
                            msg = ('User background metadata file "{}", '
                                   'image type "{}", does not specify '
                                   'metadata item "{}"')
                            raise ValueError(msg.format(json_file, img_type,
                                                        required))
                    for resln in _USER_BG_IMGS[img_type]:
                        # the required_info items are not resolutions:
                        if resln not in required_info:
                            img_it_r = _USER_BG_IMGS[img_type][resln]
                            test_file = os.path.join(bgdir, img_it_r)
                            if not os.path.isfile(test_file):
                                msg = 'File "{}" not found'
                                raise ValueError(msg.format(test_file))

    def add_raster(self, raster_source, **slippy_image_kwargs):
        """
        Add the given raster source to the GeoAxes.

        Parameters
        ----------
        raster_source : :class:`cartopy.io.RasterSource` like instance
            ``raster_source`` may be any object which implements the
            RasterSource interface, including instances of objects such as
            :class:`~cartopy.io.ogc_clients.WMSRasterSource` and
            :class:`~cartopy.io.ogc_clients.WMTSRasterSource`. Note that image
            retrievals are done at draw time, not at creation time.

        """
        # Allow a fail-fast error if the raster source cannot provide
        # images in the current projection.
        raster_source.validate_projection(self.projection)
        img = SlippyImageArtist(self, raster_source, **slippy_image_kwargs)
        with self.hold_limits():
            self.add_image(img)
        return img

    def _regrid_shape_aspect(self, regrid_shape, target_extent):
        """
        Helper for setting regridding shape which is used in several
        plotting methods.

        """
        if not isinstance(regrid_shape, collections.Sequence):
            target_size = int(regrid_shape)
            x_range, y_range = np.diff(target_extent)[::2]
            desired_aspect = x_range / y_range
            if x_range >= y_range:
                regrid_shape = (target_size * desired_aspect, target_size)
            else:
                regrid_shape = (target_size, target_size / desired_aspect)
        return regrid_shape

    def imshow(self, img, *args, **kwargs):
        """
        Add the "transform" keyword to :func:`~matplotlib.pyplot.imshow'.

        Parameters
        ----------

        transform : :class:`~cartopy.crs.Projection` or matplotlib transform
            The coordinate system in which the given image is rectangular.
        regrid_shape : int or pair of ints
            The shape of the desired image if it needs to be transformed.
            If a single integer is given then that will be used as the minimum
            length dimension, while the other dimension will be scaled up
            according to the target extent's aspect ratio. The default is for
            the minimum dimension of a transformed image to have length 750,
            so for an image being transformed into a global PlateCarree
            projection the resulting transformed image would have a shape of
            ``(750, 1500)``.
        extent : tuple
            The corner coordinates of the image in the form
            ``(left, right, bottom, top)``. The coordinates should be in the
            coordinate system passed to the transform keyword.
        origin : {'lower', 'upper'}
            The origin of the vertical pixels. See
            :func:`matplotlib.pyplot.imshow` for further details. Default
            is ``'lower'``.

        """
        transform = kwargs.pop('transform', None)
        if 'update_datalim' in kwargs:
            raise ValueError('The update_datalim keyword has been removed in '
                             'imshow. To hold the data and view limits see '
                             'GeoAxes.hold_limits.')

        kwargs.setdefault('origin', 'lower')

        same_projection = (isinstance(transform, ccrs.Projection) and
                           self.projection == transform)

        if transform is None or transform == self.transData or same_projection:
            if isinstance(transform, ccrs.Projection):
                transform = transform._as_mpl_transform(self)
            result = matplotlib.axes.Axes.imshow(self, img, *args, **kwargs)
        else:
            extent = kwargs.pop('extent', None)
            img = np.asanyarray(img)
            if kwargs['origin'] == 'upper':
                # It is implicitly assumed by the regridding operation that the
                # origin of the image is 'lower', so simply adjust for that
                # here.
                img = img[::-1]
                kwargs['origin'] = 'lower'

            if not isinstance(transform, ccrs.Projection):
                raise ValueError('Expected a projection subclass. Cannot '
                                 'handle a %s in imshow.' % type(transform))

            target_extent = self.get_extent(self.projection)
            regrid_shape = kwargs.pop('regrid_shape', 750)
            regrid_shape = self._regrid_shape_aspect(regrid_shape,
                                                     target_extent)
            warp_array = cartopy.img_transform.warp_array
            img, extent = warp_array(img,
                                     source_proj=transform,
                                     source_extent=extent,
                                     target_proj=self.projection,
                                     target_res=regrid_shape,
                                     target_extent=target_extent,
                                     mask_extrapolated=True,
                                     )

            # As a workaround to a matplotlib limitation, turn any images
            # which are RGB with a mask into RGBA images with an alpha
            # channel.
            if (isinstance(img, np.ma.MaskedArray) and
                    img.shape[2:3] == (3, ) and
                    img.mask is not False):
                old_img = img
                img = np.zeros(img.shape[:2] + (4, ), dtype=img.dtype)
                img[:, :, 0:3] = old_img
                # Put an alpha channel in if the image was masked.
                img[:, :, 3] = ~ np.any(old_img.mask, axis=2)
                if img.dtype.kind == 'u':
                    img[:, :, 3] *= 255

            result = matplotlib.axes.Axes.imshow(self, img, *args,
                                                 extent=extent, **kwargs)

        # clip the image. This does not work as the patch moves with mouse
        # movement, but the clip path doesn't
        # This could definitely be fixed in matplotlib
#        if result.get_clip_path() in [None, self.patch]:
#            # image does not already have clipping set, clip to axes patch
#            result.set_clip_path(self.outline_patch)
        return result

    def gridlines(self, crs=None, draw_labels=False, xlocs=None,
                  ylocs=None, **kwargs):
        """
        Automatically adds gridlines to the axes, in the given coordinate
        system, at draw time.

        Kwargs:

        * crs
            The :class:`cartopy._crs.CRS` defining the coordinate system in
            which gridlines are drawn.
            Default is :class:`cartopy.crs.PlateCarree`.

        * draw_labels
            Label gridlines like axis ticks, around the edge.

        * xlocs
            An iterable of gridline locations or a
            :class:`matplotlib.ticker.Locator` instance which will be used to
            determine the locations of the gridlines in the x-coordinate of
            the given CRS. Defaults to None, which implies automatic locating
            of the gridlines.

        * ylocs
            An iterable of gridline locations or a
            :class:`matplotlib.ticker.Locator` instance which will be used to
            determine the locations of the gridlines in the y-coordinate of
            the given CRS. Defaults to None, which implies automatic locating
            of the gridlines.

        Returns:

            A :class:`cartopy.mpl.gridliner.Gridliner` instance.

        All other keywords control line properties.  These are passed through
        to :class:`matplotlib.collections.Collection`.

        """
        if crs is None:
            crs = ccrs.PlateCarree()
        from cartopy.mpl.gridliner import Gridliner
        if xlocs is not None and not isinstance(xlocs, mticker.Locator):
            xlocs = mticker.FixedLocator(xlocs)
        if ylocs is not None and not isinstance(ylocs, mticker.Locator):
            ylocs = mticker.FixedLocator(ylocs)
        gl = Gridliner(
            self, crs=crs, draw_labels=draw_labels, xlocator=xlocs,
            ylocator=ylocs, collection_kwargs=kwargs)
        self._gridliners.append(gl)
        return gl

    def _gen_axes_spines(self, locations=None, offset=0.0, units='inches'):
        # generate some axes spines, as some Axes super class machinery
        # requires them. Just make them invisible
        spines = matplotlib.axes.Axes._gen_axes_spines(self,
                                                       locations=locations,
                                                       offset=offset,
                                                       units=units)
        for spine in spines.values():
            spine.set_visible(False)
        return spines

    def _boundary(self):
        """
        Adds the map's boundary to this GeoAxes, attaching the appropriate
        artists to :data:`.outline_patch` and :data:`.background_patch`.

        .. note::

            The boundary is not the ``axes.patch``. ``axes.patch``
            is made invisible by this method - its only remaining
            purpose is to provide a rectilinear clip patch for
            all Axes artists.

        """
        # Hide the old "background" patch used by matplotlib - it is not
        # used by cartopy's GeoAxes.
        self.patch.set_facecolor((1, 1, 1, 0))
        self.patch.set_edgecolor((0.5, 0.5, 0.5))
        self.patch.set_visible(False)
        self.background_patch = None
        self.outline_patch = None

        path, = cpatch.geos_to_path(self.projection.boundary)

        # Get the outline path in terms of self.transData
        proj_to_data = self.projection._as_mpl_transform(self) - self.transData
        trans_path = proj_to_data.transform_path(path)

        # Set the boundary - we can make use of the rectangular clipping.
        self.set_boundary(trans_path, use_as_clip_path=False)

        # Attach callback events for when the xlim or ylim are changed. This
        # is what triggers the patches to be re-clipped at draw time.
        self.callbacks.connect('xlim_changed', _trigger_patch_reclip)
        self.callbacks.connect('ylim_changed', _trigger_patch_reclip)

    def set_boundary(self, path, transform=None, use_as_clip_path=True):
        """
        Given a path, update the :data:`.outline_patch` and
        :data:`.background_patch` to take its shape.

        Parameters
        ----------

        path : :class:`matplotlib.path.Path`
            The path of the desired boundary.
        transform : None or :class:`matplotlib.transforms.Transform`
            The coordinate system of the given path. Currently this must be
            convertible to data coordinates, and therefore cannot extend beyond
            the limits of the axes' projection.
        use_as_clip_path : bool
            Whether axes.patch should be updated. Updating axes.patch means
            that any artists subsequently created will inherit clipping from
            this path, rather than the standard unit square in axes
            coordinates.

        """
        if transform is None:
            transform = self.transData

        if isinstance(transform, cartopy.crs.CRS):
            transform = transform._as_mpl_transform(self)

        if self.background_patch is None:
            background = matplotlib.patches.PathPatch(path, edgecolor='none',
                                                      facecolor='white',
                                                      zorder=-1, clip_on=False,
                                                      transform=transform)
        else:
            background = matplotlib.patches.PathPatch(path, zorder=-1,
                                                      clip_on=False)
            background.update_from(self.background_patch)
            self.background_patch.remove()
            background.set_transform(transform)

        if self.outline_patch is None:
            outline = matplotlib.patches.PathPatch(path, edgecolor='black',
                                                   facecolor='none',
                                                   zorder=2.5, clip_on=False,
                                                   transform=transform)
        else:
            outline = matplotlib.patches.PathPatch(path, zorder=2.5,
                                                   clip_on=False)
            outline.update_from(self.outline_patch)
            self.outline_patch.remove()
            outline.set_transform(transform)

        # Attach the original path to the patches. This will be used each time
        # a new clipped path is calculated.
        outline.orig_path = path
        background.orig_path = path

        # Attach a "reclip" attribute, which determines if the patch's path is
        # reclipped before drawing. A callback is used to change the "reclip"
        # state.
        outline.reclip = True
        background.reclip = True

        # Add the patches to the axes, and also make them available as
        # attributes.
        self.background_patch = background
        self.outline_patch = outline

        if use_as_clip_path:
            self.patch = background

        with self.hold_limits():
            self.add_patch(outline)
            self.add_patch(background)

    def contour(self, *args, **kwargs):
        """
        Add the "transform" keyword to :func:`~matplotlib.pyplot.contour'.

        Extra kwargs:

            transform - a :class:`~cartopy.crs.Projection`.

        """
        t = kwargs.get('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical contouring is not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        if isinstance(t, ccrs.Projection):
            kwargs['transform'] = t._as_mpl_transform(self)
        else:
            kwargs['transform'] = t
        result = matplotlib.axes.Axes.contour(self, *args, **kwargs)

        self.autoscale_view()
        return result

    def contourf(self, *args, **kwargs):
        """
        Add the "transform" keyword to :func:`~matplotlib.pyplot.contourf'.

        Extra kwargs:

            transform - a :class:`~cartopy.crs.Projection`.

        """
        t = kwargs.get('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical contouring is not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        if isinstance(t, ccrs.Projection):
            kwargs['transform'] = t = t._as_mpl_transform(self)
        else:
            kwargs['transform'] = t

        # Set flag to indicate correcting orientation of paths if not ccw
        if isinstance(t, mtransforms.Transform):
            for sub_trans, _ in t._iter_break_from_left_to_right():
                if isinstance(sub_trans, InterProjectionTransform):
                    if not hasattr(sub_trans, 'force_path_ccw'):
                        sub_trans.force_path_ccw = True

        result = matplotlib.axes.Axes.contourf(self, *args, **kwargs)

        # We need to compute the dataLim correctly for contours.
        if matplotlib.__version__ >= '1.4':
            extent = mtransforms.Bbox.union([col.get_datalim(self.transData)
                                             for col in result.collections
                                             if col.get_paths()])
            self.dataLim.update_from_data_xy(extent.get_points())

        self.autoscale_view()
        return result

    def scatter(self, *args, **kwargs):
        """
        Add the "transform" keyword to :func:`~matplotlib.pyplot.scatter'.

        Extra kwargs:

            transform - a :class:`~cartopy.crs.Projection`.

        """
        t = kwargs.get('transform', None)
        # Keep this bit - even at mpl v1.2
        if t is None:
            t = self.projection
        if hasattr(t, '_as_mpl_transform'):
            kwargs['transform'] = t._as_mpl_transform(self)

        # exclude Geodetic as a vaild source CS
        if (isinstance(kwargs.get('transform', None),
                       InterProjectionTransform) and
                kwargs['transform'].source_projection.is_geodetic()):
            raise ValueError('Cartopy cannot currently do spherical '
                             'contouring. The source CRS cannot be a '
                             'geodetic, consider using the cyllindrical form '
                             '(PlateCarree or RotatedPole).')

        result = matplotlib.axes.Axes.scatter(self, *args, **kwargs)
        self.autoscale_view()
        return result

    def pcolormesh(self, *args, **kwargs):
        """
        Add the "transform" keyword to :func:`~matplotlib.pyplot.pcolormesh'.

        Extra kwargs:

            transform - a :class:`~cartopy.crs.Projection`.

        """
        t = kwargs.get('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical pcolormesh is not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        kwargs.setdefault('transform', t)
        result = self._pcolormesh_patched(*args, **kwargs)
        self.autoscale_view()
        return result

    def _pcolormesh_patched(self, *args, **kwargs):
        """
        A temporary, modified duplicate of
        :func:`~matplotlib.pyplot.pcolormesh'.

        This function contains a workaround for a matplotlib issue
        and will be removed once the issue has been resolved.
        https://github.com/matplotlib/matplotlib/pull/1314

        """
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

        if not self._hold:
            self.cla()

        alpha = kwargs.pop('alpha', None)
        norm = kwargs.pop('norm', None)
        cmap = kwargs.pop('cmap', None)
        vmin = kwargs.pop('vmin', None)
        vmax = kwargs.pop('vmax', None)
        shading = kwargs.pop('shading', 'flat').lower()
        antialiased = kwargs.pop('antialiased', False)
        kwargs.setdefault('edgecolors', 'None')

        allmatch = (shading == 'gouraud')

        X, Y, C = self._pcolorargs('pcolormesh', *args, allmatch=allmatch)
        Ny, Nx = X.shape

        # convert to one dimensional arrays
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
        if norm is not None:
            assert(isinstance(norm, mcolors.Normalize))
        collection.set_cmap(cmap)
        collection.set_norm(norm)
        collection.set_clim(vmin, vmax)
        collection.autoscale_None()

        self.grid(False)

        # Transform from native to data coordinates?
        t = collection._transform
        if (not isinstance(t, mtransforms.Transform) and
                hasattr(t, '_as_mpl_transform')):
            t = t._as_mpl_transform(self.axes)

        if t and any(t.contains_branch_seperately(self.transData)):
            trans_to_data = t - self.transData
            pts = np.vstack([X, Y]).T.astype(np.float)
            transformed_pts = trans_to_data.transform(pts)
            X = transformed_pts[..., 0]
            Y = transformed_pts[..., 1]

            ########################
            # PATCH
            # XXX Non-standard matplotlib thing.
            no_inf = (X != np.inf) & (Y != np.inf)
            X = X[no_inf]
            Y = Y[no_inf]
            # END OF PATCH
            ##############

        ########################
        # PATCH
        # XXX Non-standard matplotlib thing (length check).
        if len(X):
            minx = np.amin(X)
            maxx = np.amax(X)
        else:
            minx = maxx = np.nan
        if len(Y):
            miny = np.amin(Y)
            maxy = np.amax(Y)
        else:
            miny = maxy = np.nan
        # END OF PATCH
        ##############

        corners = (minx, miny), (maxx, maxy)
        ########################
        # PATCH
        # XXX Non-standard matplotlib thing.
        collection._corners = mtransforms.Bbox(corners)
        collection.get_datalim = lambda transData: collection._corners
        # END OF PATCH
        ##############

        self.update_datalim(corners)
        self.add_collection(collection)
        self.autoscale_view()

        ########################
        # PATCH
        # XXX Non-standard matplotlib thing.
        # Handle a possible wrap around for rectangular projections.
        t = kwargs.get('transform', None)
        if isinstance(t, ccrs.CRS):
            wrap_proj_types = (ccrs._RectangularProjection,
                               ccrs._WarpedRectangularProjection,
                               ccrs.InterruptedGoodeHomolosine,
                               ccrs.Mercator)
            if isinstance(t, wrap_proj_types) and \
                    isinstance(self.projection, wrap_proj_types):

                C = C.reshape((Ny - 1, Nx - 1))
                transformed_pts = transformed_pts.reshape((Ny, Nx, 2))

                # compute the vertical line angles of the pcolor in
                # transformed coordinates
                with np.errstate(invalid='ignore'):
                    horizontal_vert_angles = np.arctan2(
                        np.diff(transformed_pts[..., 0], axis=1),
                        np.diff(transformed_pts[..., 1], axis=1)
                    )

                # if the change in angle is greater than 90 degrees (absolute),
                # then mark it for masking later on.
                dx_horizontal = np.diff(horizontal_vert_angles)
                to_mask = ((np.abs(dx_horizontal) > np.pi / 2) |
                           np.isnan(dx_horizontal))

                if np.any(to_mask):
                    if collection.get_cmap()._rgba_bad[3] != 0.0:
                        warnings.warn("The colormap's 'bad' has been set, but "
                                      "in order to wrap pcolormesh across the "
                                      "map it must be fully transparent.")

                    # at this point C has a shape of (Ny-1, Nx-1), to_mask has
                    # a shape of (Ny, Nx-2) and pts has a shape of (Ny*Nx, 2)

                    mask = np.zeros(C.shape, dtype=np.bool)

                    # mask out the neighbouring cells if there was a cell
                    # found with an angle change of more than pi/2 . NB.
                    # Masking too much only has a detrimental impact on
                    # performance.
                    to_mask_y_shift = to_mask[:-1, :]
                    mask[:, :-1][to_mask_y_shift] = True
                    mask[:, 1:][to_mask_y_shift] = True

                    to_mask_x_shift = to_mask[1:, :]
                    mask[:, :-1][to_mask_x_shift] = True
                    mask[:, 1:][to_mask_x_shift] = True

                    C_mask = getattr(C, 'mask', None)
                    if C_mask is not None:
                        dmask = mask | C_mask
                    else:
                        dmask = mask

                    # create the masked array to be used with this pcolormesh
                    pcolormesh_data = np.ma.array(C, mask=mask)

                    collection.set_array(pcolormesh_data.ravel())

                    # now that the pcolormesh has masked the bad values,
                    # create a pcolor with just those values that were masked
                    pcolor_data = pcolormesh_data.copy()
                    # invert the mask
                    pcolor_data.mask = ~pcolor_data.mask

                    # remember to re-apply the original data mask to the array
                    if C_mask is not None:
                        pcolor_data.mask = pcolor_data.mask | C_mask

                    pts = pts.reshape((Ny, Nx, 2))
                    if np.any(~pcolor_data.mask):
                        # plot with slightly lower zorder to avoid odd issue
                        # where the main plot is obscured
                        zorder = collection.zorder - .1
                        kwargs.pop('zorder', None)
                        kwargs.setdefault('snap', False)
                        pcolor_col = self.pcolor(pts[..., 0], pts[..., 1],
                                                 pcolor_data, zorder=zorder,
                                                 **kwargs)

                        pcolor_col.set_cmap(cmap)
                        pcolor_col.set_norm(norm)
                        pcolor_col.set_clim(vmin, vmax)
                        # scale the data according to the *original* data
                        pcolor_col.norm.autoscale_None(C)

                        # put the pcolor_col on the pcolormesh collection so
                        # that if really necessary, users can do things post
                        # this method
                        collection._wrapped_collection_fix = pcolor_col

            # Clip the QuadMesh to the projection boundary, which is required
            # to keep the shading inside the projection bounds.
            collection.set_clip_path(self.outline_patch)

        # END OF PATCH
        ##############

        return collection

    def pcolor(self, *args, **kwargs):
        """
        Add the "transform" keyword to :func:`~matplotlib.pyplot.pcolor'.

        Extra kwargs:

            transform - a :class:`~cartopy.crs.Projection`.

        """
        t = kwargs.get('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical pcolor is not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        kwargs.setdefault('transform', t)
        result = matplotlib.axes.Axes.pcolor(self, *args, **kwargs)

        # Update the datalim for this pcolor.
        limits = result.get_datalim(self.axes.transData)
        self.axes.update_datalim(limits)

        self.autoscale_view()
        return result

    def quiver(self, x, y, u, v, *args, **kwargs):
        """
        Plot a field of arrows.

        Extra Kwargs:

        * transform: :class:`cartopy.crs.Projection` or matplotlib transform
            The coordinate system in which the vectors are defined.

        * regrid_shape: int or 2-tuple of ints
            If given, specifies that the points where the arrows are
            located will be interpolated onto a regular grid in
            projection space. If a single integer is given then that
            will be used as the minimum grid length dimension, while the
            other dimension will be scaled up according to the target
            extent's aspect ratio. If a pair of ints are given they
            determine the grid length in the x and y directions
            respectively.

        * target_extent: 4-tuple
            If given, specifies the extent in the target CRS that the
            regular grid defined by *regrid_shape* will have. Defaults
            to the current extent of the map projection.

        See :func:`matplotlib.pyplot.quiver` for details on arguments
        and other keyword arguments.

        .. note::

           The vector components must be defined as grid eastward and
           grid northward.

        """
        t = kwargs.get('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical quiver is not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        if isinstance(t, ccrs.Projection):
            kwargs['transform'] = t._as_mpl_transform(self)
        else:
            kwargs['transform'] = t
        regrid_shape = kwargs.pop('regrid_shape', None)
        target_extent = kwargs.pop('target_extent',
                                   self.get_extent(self.projection))
        if regrid_shape is not None:
            # If regridding is required then we'll be handling transforms
            # manually and plotting in native coordinates.
            regrid_shape = self._regrid_shape_aspect(regrid_shape,
                                                     target_extent)
            if args:
                # Interpolate color array as well as vector components.
                x, y, u, v, c = vector_scalar_to_grid(
                    t, self.projection, regrid_shape, x, y, u, v, args[0],
                    target_extent=target_extent)
                args = (c,) + args[1:]
            else:
                x, y, u, v = vector_scalar_to_grid(
                    t, self.projection, regrid_shape, x, y, u, v,
                    target_extent=target_extent)
            kwargs.pop('transform', None)
        elif t != self.projection:
            # Transform the vectors if the projection is not the same as the
            # data transform.
            if (x.ndim == 1 and y.ndim == 1) and (x.shape != u.shape):
                x, y = np.meshgrid(x, y)
            u, v = self.projection.transform_vectors(t, x, y, u, v)
        return matplotlib.axes.Axes.quiver(self, x, y, u, v, *args, **kwargs)

    def barbs(self, x, y, u, v, *args, **kwargs):
        """
        Plot a field of barbs.

        Extra Kwargs:

        * transform: :class:`cartopy.crs.Projection` or matplotlib transform
            The coordinate system in which the vectors are defined.

        * regrid_shape: int or 2-tuple of ints
            If given, specifies that the points where the arrows are
            located will be interpolated onto a regular grid in
            projection space. If a single integer is given then that
            will be used as the minimum grid length dimension, while the
            other dimension will be scaled up according to the target
            extent's aspect ratio. If a pair of ints are given they
            determine the grid length in the x and y directions
            respectively.

        * target_extent: 4-tuple
            If given, specifies the extent in the target CRS that the
            regular grid defined by *regrid_shape* will have. Defaults
            to the current extent of the map projection.

        See :func:`matplotlib.pyplot.barbs` for details on arguments
        and keyword arguments.

        .. note::

           The vector components must be defined as grid eastward and
           grid northward.

        """
        t = kwargs.get('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical barbs are not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        if isinstance(t, ccrs.Projection):
            kwargs['transform'] = t._as_mpl_transform(self)
        else:
            kwargs['transform'] = t
        regrid_shape = kwargs.pop('regrid_shape', None)
        target_extent = kwargs.pop('target_extent',
                                   self.get_extent(self.projection))
        if regrid_shape is not None:
            # If regridding is required then we'll be handling transforms
            # manually and plotting in native coordinates.
            regrid_shape = self._regrid_shape_aspect(regrid_shape,
                                                     target_extent)
            if args:
                # Interpolate color array as well as vector components.
                x, y, u, v, c = vector_scalar_to_grid(
                    t, self.projection, regrid_shape, x, y, u, v, args[0],
                    target_extent=target_extent)
                args = (c,) + args[1:]
            else:
                x, y, u, v = vector_scalar_to_grid(
                    t, self.projection, regrid_shape, x, y, u, v,
                    target_extent=target_extent)
            kwargs.pop('transform', None)
        elif t != self.projection:
            # Transform the vectors if the projection is not the same as the
            # data transform.
            if (x.ndim == 1 and y.ndim == 1) and (x.shape != u.shape):
                x, y = np.meshgrid(x, y)
            u, v = self.projection.transform_vectors(t, x, y, u, v)
        return matplotlib.axes.Axes.barbs(self, x, y, u, v, *args, **kwargs)

    def streamplot(self, x, y, u, v, **kwargs):
        """
        Draws streamlines of a vector flow.

        Extra Kwargs:

        * transform: :class:`cartopy.crs.Projection` or matplotlib transform
            The coordinate system in which the vector field is defined.

        See :func:`matplotlib.pyplot.streamplot` for details on arguments
        and keyword arguments.

        .. note::

           The vector components must be defined as grid eastward and
           grid northward.

        """
        t = kwargs.pop('transform', None)
        if t is None:
            t = self.projection
        if isinstance(t, ccrs.CRS) and not isinstance(t, ccrs.Projection):
            raise ValueError('invalid transform:'
                             ' Spherical streamplot is not supported - '
                             ' consider using PlateCarree/RotatedPole.')
        # Regridding is required for streamplot, it must have an evenly spaced
        # grid to work correctly. Choose our destination grid based on the
        # density keyword. The grid need not be bigger than the grid used by
        # the streamplot integrator.
        density = kwargs.get('density', 1)
        if np.isscalar(density):
            regrid_shape = [int(30 * density)] * 2
        else:
            regrid_shape = [int(25 * d) for d in density]
        # The color and linewidth keyword arguments can be arrays so they will
        # need to be gridded also.
        c = kwargs.get('color', None)
        l = kwargs.get('linewidth', None)
        scalars = []
        color_array = isinstance(c, np.ndarray)
        linewidth_array = isinstance(l, np.ndarray)
        if color_array:
            scalars.append(c)
        if linewidth_array:
            scalars.append(l)
        # Do the regridding including any scalar fields.
        target_extent = self.get_extent(self.projection)
        gridded = vector_scalar_to_grid(t, self.projection, regrid_shape,
                                        x, y, u, v, *scalars,
                                        target_extent=target_extent)
        x, y, u, v = gridded[:4]
        # If scalar fields were regridded then replace the appropriate keyword
        # arguments with the gridded arrays.
        scalars = list(gridded[4:])
        if linewidth_array:
            kwargs['linewidth'] = scalars.pop()
        if color_array:
            kwargs['color'] = ma.masked_invalid(scalars.pop())
        with warnings.catch_warnings():
            # The workaround for nan values in streamplot colors gives rise to
            # a warning which is not at all important so it is hidden from the
            # user to avoid confusion.
            message = 'Warning: converting a masked element to nan.'
            warnings.filterwarnings('ignore', message=message,
                                    category=UserWarning)
            sp = matplotlib.axes.Axes.streamplot(self, x, y, u, v, **kwargs)
        return sp

    def add_wmts(self, wmts, layer_name, wmts_kwargs=None, **kwargs):
        """
        Add the specified WMTS layer to the axes.

        This function requires owslib and PIL to work.

        Args:

            * wmts - The URL of the WMTS, or an
                     owslib.wmts.WebMapTileService instance.
            * layer_name - The name of the layer to use.

        Kwargs:

            * wmts_kwargs - dict or None. Passed through to the
                :class:`~cartopy.io.ogc_clients.WMTSRasterSource`
                constructor's ``gettile_extra_kwargs`` (e.g. time).

        All other keywords are passed through to the construction of the
        image artist. See :meth:`~matplotlib.axes.Axes.imshow()` for
        more details.

        """
        from cartopy.io.ogc_clients import WMTSRasterSource
        wmts = WMTSRasterSource(wmts, layer_name,
                                gettile_extra_kwargs=wmts_kwargs)
        return self.add_raster(wmts, **kwargs)

    def add_wms(self, wms, layers, wms_kwargs=None, **kwargs):
        """
        Add the specified WMS layer to the axes.

        This function requires owslib and PIL to work.

        Parameters
        ----------
        wms : string or :class:`owslib.wms.WebMapService` instance
            The web map service URL or owslib WMS instance to use.
        layers : string or iterable of string
            The name of the layer(s) to use.
        wms_kwargs : dict or None
            Passed through to the
            :class:`~cartopy.io.ogc_clients.WMSRasterSource`
            constructor's ``getmap_extra_kwargs`` for defining getmap time
            keyword arguments.

        All other keywords are passed through to the construction of the
        image artist. See :meth:`~matplotlib.axes.Axes.imshow()` for
        more details.

        """
        from cartopy.io.ogc_clients import WMSRasterSource
        wms = WMSRasterSource(wms, layers, getmap_extra_kwargs=wms_kwargs)
        return self.add_raster(wms, **kwargs)


# Define the GeoAxesSubplot class, so that a type(ax) will emanate from
# cartopy.mpl.geoaxes, not matplotlib.axes.
class GeoAxesSubplot(matplotlib.axes.SubplotBase, GeoAxes):
    _axes_class = GeoAxes


try:
    matplotlib.axes._subplots._subplot_classes[GeoAxes] = GeoAxesSubplot
except AttributeError:
    matplotlib.axes._subplot_classes[GeoAxes] = GeoAxesSubplot


def _trigger_patch_reclip(event):
    """
    Defines an event callback for a GeoAxes which forces the outline and
    background patches to be re-clipped next time they are drawn.

    """
    axes = event.axes
    # trigger the outline and background patches to be re-clipped
    axes.outline_patch.reclip = True
    axes.background_patch.reclip = True
