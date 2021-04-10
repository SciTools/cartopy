# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

"""
This module defines the :class:`FeatureArtist` class, for drawing
:class:`Feature` instances with matplotlib.

"""

from collections import OrderedDict
import warnings
import weakref

import numpy as np
import matplotlib.artist
import matplotlib.collections
import matplotlib.figure as figure


import cartopy.mpl.patch as cpatch
from cartopy.crs import PlateCarree
from .style import merge as style_merge, finalize as style_finalize

import shapely.ops as ops
import shapely.geometry as sgeom


class _GeomKey:
    """
    Provide id() based equality and hashing for geometries.

    Instances of this class must be treated as immutable for the caching
    to operate correctly.

    A workaround for Shapely polygons no longer being hashable as of 1.5.13.

    """
    def __init__(self, geom):
        self._id = id(geom)

    def __eq__(self, other):
        return self._id == other._id

    def __hash__(self):
        return hash(self._id)


def _freeze(obj):
    """
    Recursively freeze the given object so that it might be suitable for
    use as a hashable.

    """
    if isinstance(obj, dict):
        obj = frozenset(((k, _freeze(v)) for k, v in obj.items()))
    elif isinstance(obj, list):
        obj = tuple(_freeze(item) for item in obj)
    elif isinstance(obj, np.ndarray):
        obj = tuple(obj)
    return obj


class FeatureArtist(matplotlib.artist.Artist):
    """
    A subclass of :class:`~matplotlib.artist.Artist` capable of
    drawing a :class:`cartopy.feature.Feature`.

    """

    _geom_key_to_geometry_cache = weakref.WeakValueDictionary()
    """
    A mapping from _GeomKey to geometry to assist with the caching of
    transformed Matplotlib paths.

    """
    _geom_key_to_path_cache = weakref.WeakKeyDictionary()
    """
    A nested mapping from geometry (converted to a _GeomKey) and target
    projection to the resulting transformed Matplotlib paths::

        {geom: {target_projection: list_of_paths}}

    This provides a significant boost when producing multiple maps of the
    same projection.

    """

    def __init__(self, feature, **kwargs):
        """
        Parameters
        ----------
        feature
            An instance of :class:`cartopy.feature.Feature` to draw.
        styler
            A callable that given a gemometry, returns matplotlib styling
            parameters.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments to be used when drawing the feature. These
            will override those shared with the feature.

        """
        super().__init__()

        if kwargs is None:
            kwargs = {}
        self._styler = kwargs.pop('styler', None)
        self._kwargs = dict(kwargs)

        if 'color' in self._kwargs:
            # We want the user to be able to override both face and edge
            # colours if the original feature already supplied it.
            color = self._kwargs.pop('color')
            self._kwargs['facecolor'] = self._kwargs['edgecolor'] = color

        # Set default zorder so that features are drawn before
        # lines e.g. contours but after images.
        # Note that the zorder of Patch, PatchCollection and PathCollection
        # are all 1 by default. Assuming equal zorder drawing takes place in
        # the following order: collections, patches, lines (default zorder=2),
        # text (default zorder=3), then other artists e.g. FeatureArtist.
        if self._kwargs.get('zorder') is not None:
            self.set_zorder(self._kwargs['zorder'])
        elif feature.kwargs.get('zorder') is not None:
            self.set_zorder(feature.kwargs['zorder'])
        else:
            # The class attribute matplotlib.collections.PathCollection.zorder
            # was removed after mpl v1.2.0, so the hard-coded value of 1 is
            # used instead.
            self.set_zorder(1)

        self._feature = feature

    def xy_samples(self, lon_sample, lat_sample):
        feature_crs = self._feature.crs
        if not isinstance(feature_crs, PlateCarree):
            x_lims = feature_crs.x_limits
            y_lims = feature_crs.y_limits

            x_samples, x_sep = np.linspace(*x_lims, 400,
                                           endpoint=True, retstep=True)

            if 0. not in x_samples:
                x_samples = np.insert(x_samples,
                                      np.searchsorted(x_samples, 0.), 0.)

            y_samples, y_sep = np.linspace(*y_lims, 400,
                                           endpoint=True, retstep=True)

            if 0. not in y_samples:
                y_samples = np.insert(y_samples,
                                      np.searchsorted(y_samples, 0.), 0.)

            return x_samples, x_sep, y_samples, y_sep
        else:
            # We need to sample the central longitude and latitude, the
            # antipode, and -1 * (central latitude)
            x_sample, y_sample = lon_sample, lat_sample

            x_sep = y_sep = 1

            x_offset = x_sample - np.floor(x_sample)
            x_samples = np.arange(-180. + x_offset, 180. + x_offset, x_sep)

            y_offset = y_sample - np.floor(y_sample)
            y_samples = np.arange(-90. + y_offset, 90. + y_offset, y_sep)

            if x_offset == 0.:
                x_samples = np.append(x_samples, 180.)
            if y_offset == 0.:
                y_samples = np.append(y_samples, 90.)

            if -y_sample not in y_samples:
                ind = np.searchsorted(y_samples, -y_sample)
                y_samples = np.insert(y_samples, ind, -y_sample)

            return x_samples, x_sep, y_samples, y_sep

    def build_invalid_geom(self):
        ax = self.axes
        # the following may work after improvements to multipolygon creation
        # from projection results
#         feature_crs = self._feature.crs
#         projection_crs = ax.projection
#         projection_domain = projection_crs.domain
#         feature_proj_domain = feature_crs.project_geometry(projection_domain,
#         projection_crs)
#         feature_domain = feature_crs.domain
#         invalid_geom = feature_domain.difference(feature_proj_domain)
        lon_sample = ax.projection._lon_sample
        lat_sample = ax.projection._lat_sample

        x_samples, x_sep, y_samples, y_sep = self.xy_samples(
            lon_sample, lat_sample)
        x_grid, y_grid = np.meshgrid(x_samples, y_samples)
        xyz = ax.projection.transform_points(self._feature.crs, x_grid, y_grid)
        x_proj_grid, y_proj_grid = xyz[:, :, 0], xyz[:, :, 1]
        inds_x = ~np.isfinite(x_proj_grid)
        inds_y = ~np.isfinite(y_proj_grid)

        inds = inds_x | inds_y
        inds_bin = inds.astype(int)

        fig_ = figure.Figure()
        ax_ = fig_.add_subplot()
        # Method: create contours around the points that project to infinity
        # Then contract the exterior rings and expand the interior rings to
        # buffer the contours
        # Create Polygons from the rings
        fcontour = ax_.contourf(x_grid, y_grid,
                                inds_bin, levels=[0.5, 1.5])

        paths = fcontour.collections[0].get_paths()
        invalid_geoms = []
        for path in paths:
            poly = path.to_polygons()
            if poly:
                # 0th path is polygon exterior
                # Subsequent paths are interior rings
                exterior = poly[0]
                exterior = sgeom.LinearRing(exterior)
                offset = max(x_sep, y_sep)
                exterior = exterior.parallel_offset(offset, "right",
                                                    join_style=3)
                # Point of discussion:
                # parallel_offset can output a MultiLineString from LinearRing
                # input. Do we want to catch this behavior?
                if isinstance(exterior, sgeom.MultiLineString):
                    raise NotImplementedError
                exterior.coords = list(exterior.coords)[::-1]
                interiors = poly[1:]
                interiors_shrunk = []
                for interior in interiors:
                    interior = sgeom.LinearRing(interior)
                    interior_shrunk = interior.parallel_offset(offset, "right",
                                                               join_style=3)
                    # Again the possibility of a MultiLineString
                    if not interior_shrunk.is_empty:
                        interior_shrunk.coords = list(
                                                interior_shrunk.coords)[::-1]
                        interiors_shrunk.append(interior_shrunk)
                invalid_geom = sgeom.Polygon(exterior, holes=interiors_shrunk)
                invalid_geoms.append(invalid_geom)

        invalid_geom = ops.unary_union(invalid_geoms)
        return invalid_geom

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer, *args, **kwargs):
        """
        Draw the geometries of the feature that intersect with the extent of
        the :class:`cartopy.mpl.GeoAxes` instance to which this
        object has been added.

        """
        if not self.get_visible():
            return

        ax = self.axes
        feature_crs = self._feature.crs

        # Get geometries that we need to draw.
        extent = None
        try:
            extent = ax.get_extent(feature_crs)
        except ValueError:
            warnings.warn('Unable to determine extent. Defaulting to global.')
        geoms = self._feature.intersecting_geometries(extent)

        # Combine all the keyword args in priority order.
        prepared_kwargs = style_merge(self._feature.kwargs,
                                      self._kwargs,
                                      kwargs)

        # Freeze the kwargs so that we can use them as a dict key. We will
        # need to unfreeze this with dict(frozen) before passing to mpl.
        prepared_kwargs = _freeze(prepared_kwargs)

        # Project (if necessary) and convert geometries to matplotlib paths.
        stylised_paths = OrderedDict()
        key = ax.projection

        skip_geoms = ["LineString", "LinearRing", "MultiLineString"]

        for geom in geoms:
            # As Shapely geometries cannot be relied upon to be
            # hashable, we have to use a WeakValueDictionary to manage
            # their weak references. The key can then be a simple,
            # "disposable", hashable geom-key object that just uses the
            # id() of a geometry to determine equality and hash value.
            # The only persistent, strong reference to the geom-key is
            # in the WeakValueDictionary, so when the geometry is
            # garbage collected so is the geom-key.
            # The geom-key is also used to access the WeakKeyDictionary
            # cache of transformed geometries. So when the geom-key is
            # garbage collected so are the transformed geometries.
            geom_key = _GeomKey(geom)
            FeatureArtist._geom_key_to_geometry_cache.setdefault(
                geom_key, geom)
            mapping = FeatureArtist._geom_key_to_path_cache.setdefault(
                geom_key, {})
            geom_paths = mapping.get(key)
            if not self._styler:
                style = prepared_kwargs
            else:
                # Unfreeze, then add the computed style,
                # and then re-freeze.
                style = style_merge(dict(prepared_kwargs),
                                    self._styler(geom))
                style = _freeze(style)
            if geom_paths is None:
                # Shapely can't represent the difference between a line and a
                # polygon. If all geoms are lines, bypass building the polygon
                # containing points projecting to infinity and use different
                # masking method.
                if (ax.projection != feature_crs) & \
                   (geom.geom_type not in skip_geoms):
                    if feature_crs not in ax.projection.invalid_geoms.keys():
                        invalid_geom = self.build_invalid_geom()
                        ax.projection.invalid_geoms[feature_crs] = invalid_geom
                    else:
                        invalid_geom = ax.projection.invalid_geoms[feature_crs]
                    if not geom.is_valid:
                        geom = geom.buffer(0)
                    cleaned_geom = geom.difference(invalid_geom)
                    if not cleaned_geom.is_empty:
                        projected_geom = ax.projection.project_geometry(
                                                         cleaned_geom,
                                                         feature_crs)
                        geom_paths = cpatch.geos_to_path(projected_geom)
                        mapping[key] = geom_paths

                elif (ax.projection != feature_crs) & \
                     (geom.geom_type in skip_geoms):
                    if isinstance(geom, sgeom.LineString):
                        geom = sgeom.MultiLineString([geom])
                    elif isinstance(geom, sgeom.LinearRing):
                        geom = [geom]
                    for subgeom in geom:
                        x, y = np.array(subgeom.xy)
                        xyz = ax.projection.transform_points(feature_crs, x, y)
                        x_proj, y_proj = xyz[:, 0], xyz[:, 1]

                        inds_x = ~np.isfinite(x_proj)
                        inds_y = ~np.isfinite(y_proj)
                        inds = inds_x | inds_y

                        x_proj = np.ma.masked_array(x_proj, mask=inds)
                        clump_slices = np.ma.clump_unmasked(x_proj)
                        geom_paths = []
                        for clump_slice in clump_slices:
                            if clump_slice.stop - clump_slice.start > 1:
                                xt = x[clump_slice]
                                yt = y[clump_slice]
                                cleaned_geom = sgeom.LineString(zip(xt, yt))
                                projected_geom = \
                                    ax.projection.project_geometry(
                                            cleaned_geom, feature_crs
                                    )
                                geom_paths.extend(
                                    cpatch.geos_to_path(projected_geom))
                        if key not in mapping.keys():
                            mapping[key] = geom_paths
                        else:
                            mapping[key].extend(geom_paths)
                else:
                    geom_paths = cpatch.geos_to_path(geom)
                    mapping[key] = geom_paths
            if geom_paths:
                stylised_paths.setdefault(style, []).extend(geom_paths)

        transform = ax.projection._as_mpl_transform(ax)

        # Draw one PathCollection per style. We could instead pass an array
        # of style items through to a single PathCollection, but that
        # complexity does not yet justify the effort.
        for style, paths in stylised_paths.items():
            style = style_finalize(dict(style))
            # Build path collection and draw it.
            c = matplotlib.collections.PathCollection(paths,
                                                      transform=transform,
                                                      **style)
            c.set_clip_path(ax.patch)
            c.set_figure(ax.figure)
            c.draw(renderer)
        # n.b. matplotlib.collection.Collection.draw returns None
        return None
