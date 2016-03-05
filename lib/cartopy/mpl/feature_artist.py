# (C) British Crown Copyright 2011 - 2016, Met Office
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
This module defines the :class:`FeatureArtist` class, for drawing
:class:`Feature` instances with matplotlib.

"""

from __future__ import (absolute_import, division, print_function)

import warnings
import weakref

import matplotlib.artist
import matplotlib.collections

import cartopy.mpl.patch as cpatch


class _GeomKey(object):
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


class FeatureArtist(matplotlib.artist.Artist):
    """
    A subclass of :class:`~matplotlib.artist.Artist` capable of
    drawing a :class:`cartopy.feature.Feature`.

    """

    _geom_key_to_geometry_cache = weakref.WeakValueDictionary()
    """
    A mapping from _GeomKey to geometry to assist with the caching of
    transformed matplotlib paths.

    """
    _geom_key_to_path_cache = weakref.WeakKeyDictionary()
    """
    A nested mapping from geometry (converted to a _GeomKey) and target
    projection to the resulting transformed matplotlib paths::

        {geom: {target_projection: list_of_paths}}

    This provides a significant boost when producing multiple maps of the
    same projection.

    """

    def __init__(self, feature, **kwargs):
        """
        Args:

        * feature:
            an instance of :class:`cartopy.feature.Feature` to draw.
        * kwargs:
            keyword arguments to be used when drawing the feature. These
            will override those shared with the feature.

        """
        super(FeatureArtist, self).__init__()

        if kwargs is None:
            kwargs = {}
        self._kwargs = dict(kwargs)

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

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer, *args, **kwargs):
        """
        Draws the geometries of the feature that intersect with the extent of
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

        # Project (if necessary) and convert geometries to matplotlib paths.
        paths = []
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
            paths.extend(geom_paths)

        # Build path collection and draw it.
        transform = ax.projection._as_mpl_transform(ax)
        # Combine all the keyword args in priority order
        final_kwargs = dict(self._feature.kwargs)
        final_kwargs.update(self._kwargs)
        final_kwargs.update(kwargs)
        c = matplotlib.collections.PathCollection(paths,
                                                  transform=transform,
                                                  **final_kwargs)
        c.set_clip_path(ax.patch)
        c.set_figure(ax.figure)
        return c.draw(renderer)
