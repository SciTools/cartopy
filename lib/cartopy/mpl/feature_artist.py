# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

"""
This module defines the :class:`FeatureArtist` class, for drawing
:class:`Feature` instances through an extension of the Matplotlib Artist interfaces.

"""

import warnings
import weakref

import matplotlib.artist
import matplotlib.collections
from matplotlib.path import Path
import numpy as np
import shapely

import cartopy.feature as cfeature
import cartopy.mpl.path as cpath


# Key in the per-geometry mapping of FeatureArtist._geom_key_to_path_cache
# under which the parts of a multi-part geometry are stored.
_PARTS_KEY = object()


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


def _is_partially_outside(geom, extent):
    """
    Return whether the given geometry consists of multiple parts, and does not
    lie completely within the given extent (x0, x1, y0, y1).

    """
    if not isinstance(geom, (shapely.MultiPolygon, shapely.MultiLineString)):
        return False
    x0, y0, x1, y1 = geom.bounds
    return x0 < extent[0] or x1 > extent[1] or y0 < extent[2] or y1 > extent[3]


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


class FeatureArtist(matplotlib.collections.Collection):
    """
    A subclass of :class:`~matplotlib.collections.Collection` capable of
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

    For multi-part geometries, the mapping additionally stores the parts of
    the geometry under ``_PARTS_KEY``, so that the paths of the individual
    parts can be cached in the same way.

    """

    def __init__(self, feature, **kwargs):
        """
        Parameters
        ----------
        feature
            An instance of :class:`cartopy.feature.Feature` to draw.
        styler
            A callable that given a geometry, returns matplotlib styling
            parameters.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments to be used when drawing the feature. These
            will override those shared with the feature.

        """
        super().__init__()

        self._styler = kwargs.pop('styler', None)
        self._kwargs = dict(kwargs)

        if 'color' in self._kwargs:
            # We want the user to be able to override both face and edge
            # colours if the original feature already supplied it.
            color = self._kwargs.pop('color')
            self._kwargs['facecolor'] = self._kwargs['edgecolor'] = color

        # Paths are worked out at draw, but add_collection fails if paths is
        # left to the default of None.
        self.set_paths([])

        # Set default zorder so that features are drawn under
        # lines e.g. contours but over images and filled patches.
        # Note that the zorder of Patch, PatchCollection and PathCollection
        # are all 1 by default. Assuming default zorder, drawing takes place in
        # the following order: collections, patches, FeatureArtist, lines,
        # text.
        self.set_zorder(1.5)

        # Update drawing styles from the feature and **kwargs.
        self.set(**feature.kwargs)
        self.set(**self._kwargs)

        self._feature = feature

    def set_facecolor(self, c):
        """
        Set the facecolor(s) of the `.FeatureArtist`.  If set to 'never' then
        subsequent calls will have no effect.  Otherwise works the same as
        `matplotlib.collections.Collection.set_facecolor`.
        """
        if isinstance(c, str) and c == 'never':
            self._never_fc = True
            super().set_facecolor('none')

        elif (getattr(self, '_never_fc', False) and
                (not isinstance(c, str) or c != 'none')):
            warnings.warn('facecolor will have no effect as it has been '
                          'defined as "never".')
        else:
            super().set_facecolor(c)

    def _get_geoms_paths(self):
        ax = self.axes
        feature_crs = self._feature.crs

        # Get geometries that we need to draw.
        extent = None
        try:
            extent = ax.get_extent(feature_crs)
        except ValueError:
            warnings.warn('Unable to determine extent. Defaulting to global.')

        if isinstance(self._feature, cfeature.ShapelyFeature):
            # User passed a specific list of geometries.  If they also passed
            # `array` or a list of facecolors then we should keep the colours
            # consistent after pan/zoom.  Do this by creating a Path for every
            # geometry regardless of whether they are currently in view.
            geoms = self._feature.geometries()
        else:
            # For efficiency on local maps with high resolution features (e.g
            # from Natural Earth), only create paths for geometries that are
            # in view.
            geoms = self._feature.intersecting_geometries(extent)

        extent_geom = None
        # shapely 2.0 returns tuple of NaNs instead of None for empty geometry
        # -> check for both
        if extent is not None and not np.isnan(extent[0]):
            extent_geom = shapely.box(extent[0], extent[2], extent[1], extent[3])
            shapely.prepare(extent_geom)

        # Use the CRS of the feature from before the geometries are read above
        # (which is deferred until iterating over ``geoms``), as reading them
        # can replace it with an equivalent CRS, e.g. from the .prj file of a
        # Natural Earth shapefile, that no longer compares equal to the
        # projection of the axes.
        src_crs = feature_crs if ax.projection != feature_crs else None

        # Project (if necessary) and convert geometries to matplotlib paths.
        for geom in geoms:
            mapping = self._get_path_mapping(geom)
            if extent_geom is not None and _is_partially_outside(geom, extent):
                # Multi-part geometries, e.g. from Natural Earth, often span
                # the whole globe even if only a few of their parts are in
                # view. Only project the parts that are in view and combine
                # their paths, which gives the same path as projecting the
                # whole geometry, minus the parts that are not visible.
                parts = mapping.get(_PARTS_KEY)
                if parts is None:
                    parts = mapping[_PARTS_KEY] = shapely.get_parts(geom)
                visible_parts = parts[shapely.intersects(extent_geom, parts)]
                geom_path = Path.make_compound_path(
                    *[self._get_path(part, self._get_path_mapping(part), src_crs)
                      for part in visible_parts])
            else:
                geom_path = self._get_path(geom, mapping, src_crs)

            yield geom, geom_path

    @staticmethod
    def _get_path_mapping(geom):
        """
        Return the cached mapping from target projection to the path of the
        given geometry.

        """
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
        FeatureArtist._geom_key_to_geometry_cache.setdefault(geom_key, geom)
        return FeatureArtist._geom_key_to_path_cache.setdefault(geom_key, {})

    def _get_path(self, geom, mapping, src_crs):
        """
        Return the path of the given geometry in the projection of the axes,
        using and updating the cache in ``mapping``.

        The geometry is projected from ``src_crs``, unless that is None.

        """
        key = self.axes.projection
        geom_path = mapping.get(key)
        if geom_path is None:
            if src_crs is not None:
                projected_geom = key.project_geometry(geom, src_crs)
            else:
                projected_geom = geom

            geom_path = cpath.shapely_to_path(projected_geom)
            mapping[key] = geom_path
        return geom_path

    def get_paths(self):
        paths = super().get_paths()
        if paths:
            # When we are drawing, there is an explicit list of paths set.
            # Return these for the renderer.
            return paths

        # When not drawing, the path list is empty.  Find all the relevant paths for
        # the current axes extent.
        return [path for _, path in self._get_geoms_paths()]

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer):
        """
        Draw the geometries of the feature that intersect with the extent of
        the :class:`cartopy.mpl.geoaxes.GeoAxes` instance to which this
        object has been added.

        """
        if not self.get_visible():
            return

        stylised_paths = {}
        # Make an empty placeholder style dictionary for when styler is not
        # used.  Freeze it so that we can use it as a dict key.  We will need
        # to unfreeze all style dicts with dict(frozen) before passing to mpl.
        no_style = _freeze({})
        for geom, geom_path in self._get_geoms_paths():
            if self._styler is None:
                stylised_paths.setdefault(no_style, []).append(geom_path)
            else:
                style = _freeze(self._styler(geom))
                stylised_paths.setdefault(style, []).append(geom_path)

        self.set_clip_path(self.axes.patch)

        # Draw each style individually.  Note that there will only be multiple
        # styles if styler was used.
        for style, paths in stylised_paths.items():
            style = dict(style)

            # Temporarily replace properties.
            orig_style = {k: getattr(self, f"get_{k}")() for k in style}
            self.set(paths=paths, **style)

            super().draw(renderer)

            self.set(paths=[], **orig_style)
