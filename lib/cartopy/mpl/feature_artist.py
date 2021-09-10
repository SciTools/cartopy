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

import cartopy.mpl.patch as cpatch
from .style import merge as style_merge, finalize as style_finalize


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

        # Set default zorder so that features are drawn under
        # lines e.g. contours but over images and filled patches.
        # Note that the zorder of Patch, PatchCollection and PathCollection
        # are all 1 by default. Assuming default zorder, drawing takes place in
        # the following order: collections, patches, FeatureArtist, lines,
        # text.
        if self._kwargs.get('zorder') is not None:
            self.set_zorder(self._kwargs['zorder'])
        elif feature.kwargs.get('zorder') is not None:
            self.set_zorder(feature.kwargs['zorder'])
        else:
            self.set_zorder(1.5)

        self._feature = feature

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
            if geom_paths is None:
                if ax.projection != feature_crs:
                    projected_geom = ax.projection.project_geometry(
                        geom, feature_crs)
                else:
                    projected_geom = geom
                geom_paths = cpatch.geos_to_path(projected_geom)
                mapping[key] = geom_paths

            if not self._styler:
                style = prepared_kwargs
            else:
                # Unfreeze, then add the computed style, and then re-freeze.
                style = style_merge(dict(prepared_kwargs), self._styler(geom))
                style = _freeze(style)

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
