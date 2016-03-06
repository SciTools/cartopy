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
Provides shapely geometry <-> matplotlib path support.


See also `Shapely Geometric Objects <see_also_shapely>`_
and `Matplotlib Path API <http://matplotlib.org/api/path_api.html>`_.

.. see_also_shapely:
   http://toblerity.org/shapely/manual.html#geometric-objects

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
import matplotlib.path
from matplotlib.path import Path
import shapely.geometry as sgeom


def geos_to_path(shape):
    """
    Creates a list of :class:`matplotlib.path.Path` objects that describe
    a shape.

    Args:

    * shape
        A list, tuple or single instance of any of the following
        types: :class:`shapely.geometry.point.Point`,
        :class:`shapely.geometry.linestring.LineString`,
        :class:`shapely.geometry.polygon.Polygon`,
        :class:`shapely.geometry.multipoint.MultiPoint`,
        :class:`shapely.geometry.multipolygon.MultiPolygon`,
        :class:`shapely.geometry.multilinestring.MultiLineString`,
        :class:`shapely.geometry.collection.GeometryCollection`,
        or any type with a _as_mpl_path() method.

    Returns:
        A list of :class:`matplotlib.path.Path` objects.

    """
    if isinstance(shape, (list, tuple)):
        paths = []
        for shp in shape:
            paths.extend(geos_to_path(shp))
        return paths

    if isinstance(shape, (sgeom.LineString, sgeom.Point)):
        return [Path(np.vstack(shape.xy).T)]
    elif isinstance(shape, sgeom.Polygon):
        def poly_codes(poly):
            codes = np.ones(len(poly.xy[0])) * Path.LINETO
            codes[0] = Path.MOVETO
            return codes
        if shape.is_empty:
            return []
        vertices = np.concatenate([np.array(shape.exterior.xy)] +
                                  [np.array(ring.xy) for ring in
                                   shape.interiors], 1).T
        codes = np.concatenate([poly_codes(shape.exterior)] +
                               [poly_codes(ring) for ring in shape.interiors])
        return [Path(vertices, codes)]
    elif isinstance(shape, (sgeom.MultiPolygon, sgeom.GeometryCollection,
                            sgeom.MultiLineString, sgeom.MultiPoint)):
        paths = []
        for geom in shape.geoms:
            paths.extend(geos_to_path(geom))
        return paths
    elif hasattr(shape, '_as_mpl_path'):
        vertices, codes = shape._as_mpl_path()
        return [Path(vertices, codes)]
    else:
        raise ValueError('Unsupported shape type {}.'.format(type(shape)))


def path_segments(path, transform=None, remove_nans=False, clip=None,
                  quantize=False, simplify=False, curves=False,
                  stroke_width=1.0, snap=False):
    """
    Creates an array of vertices and a corresponding array of codes from a
    :class:`matplotlib.path.Path`.

    Args:

    * path
        A :class:`matplotlib.path.Path` instance.

    Kwargs:
        See :func:`matplotlib.path.iter_segments` for details of the keyword
        arguments.

    Returns:
        A (vertices, codes) tuple, where vertices is a numpy array of
        coordinates, and codes is a numpy array of matplotlib path codes.
        See :class:`matplotlib.path.Path` for information on the types of
        codes and their meanings.

    """
    # XXX assigned to avoid a ValueError inside the mpl C code...
    a = transform, remove_nans, clip, quantize, simplify, curves

    # Series of cleanups and conversions to the path e.g. it
    # can convert curved segments to line segments.
    vertices, codes = matplotlib.path.cleanup_path(path, transform,
                                                   remove_nans, clip,
                                                   snap, stroke_width,
                                                   simplify, curves)

    # Remove the final vertex (with code 0)
    return vertices[:-1, :], codes[:-1]


# Matplotlib v1.3+ deprecates the use of matplotlib.path.cleanup_path. Instead
# there is a method on a Path instance to simplify this.
if hasattr(matplotlib.path.Path, 'cleaned'):
    _path_segments_doc = path_segments.__doc__

    def path_segments(path, **kwargs):
        pth = path.cleaned(**kwargs)
        return pth.vertices[:-1, :], pth.codes[:-1]

    path_segments.__doc__ = _path_segments_doc


def path_to_geos(path, force_ccw=False):
    """
    Creates a list of Shapely geometric objects from a
    :class:`matplotlib.path.Path`.

    Args:

    * path
        A :class:`matplotlib.path.Path` instance.

    Kwargs:

    * force_ccw
        Boolean flag determining whether the path can be inverted to enforce
        ccw.

    Returns:
        A list of :class:`shapely.geometry.polygon.Polygon`,
        :class:`shapely.geometry.linestring.LineString` and/or
        :class:`shapely.geometry.multilinestring.MultiLineString` instances.

    """
    # Convert path into numpy array of vertices (and associated codes)
    path_verts, path_codes = path_segments(path, curves=False)

    # Split into subarrays such that each subarray consists of connected
    # line segments based on the start of each one being marked by a
    # matplotlib MOVETO code.
    verts_split_inds = np.where(path_codes == Path.MOVETO)[0]
    verts_split = np.split(path_verts, verts_split_inds)
    codes_split = np.split(path_codes, verts_split_inds)

    # Iterate through the vertices generating a list of
    # (external_geom, [internal_polygons]) tuples.
    other_result_geoms = []
    collection = []
    for path_verts, path_codes in zip(verts_split, codes_split):
        if len(path_verts) == 0:
            continue

        # XXX A path can be given which does not end with close poly, in that
        # situation, we have to guess?
        verts_same_as_first = np.all(path_verts[0, :] == path_verts[1:, :],
                                     axis=1)
        if all(verts_same_as_first):
            geom = sgeom.Point(path_verts[0, :])
        elif (path_verts.shape[0] > 2 and
                (path_codes[-1] == Path.CLOSEPOLY or
                 verts_same_as_first[-1])):
            if path_codes[-1] == Path.CLOSEPOLY:
                geom = sgeom.Polygon(path_verts[:-1, :])
            else:
                geom = sgeom.Polygon(path_verts)
        else:
            geom = sgeom.LineString(path_verts)

        # If geom is a Polygon and is contained within the last geom in
        # collection, add it to its list of internal polygons, otherwise
        # simply append it as a new external geom.
        if geom.is_empty:
            pass
        elif (len(collection) > 0 and
                isinstance(collection[-1][0], sgeom.Polygon) and
                isinstance(geom, sgeom.Polygon) and
                collection[-1][0].contains(geom.exterior)):
            collection[-1][1].append(geom.exterior)
        elif isinstance(geom, sgeom.Point):
            other_result_geoms.append(geom)
        else:
            collection.append((geom, []))

    # Convert each (external_geom, [internal_polygons]) pair into a
    # a shapely Polygon that encapsulates the internal polygons, if the
    # external geom is a LineString leave it alone.
    geom_collection = []
    for external_geom, internal_polys in collection:
        if internal_polys:
            # XXX worry about islands within lakes
            geom = sgeom.Polygon(external_geom.exterior, internal_polys)
        else:
            geom = external_geom

        # Correctly orientate the polygon (ccw)
        if isinstance(geom, sgeom.Polygon):
            if force_ccw and not geom.exterior.is_ccw:
                geom = sgeom.polygon.orient(geom)

        geom_collection.append(geom)

    # If the geom_collection only contains LineStrings combine them
    # into a single MultiLinestring.
    if geom_collection and all(isinstance(geom, sgeom.LineString) for
                               geom in geom_collection):
        geom_collection = [sgeom.MultiLineString(geom_collection)]

    # Remove any zero area Polygons
    def not_zero_poly(geom):
        return ((isinstance(geom, sgeom.Polygon) and not geom._is_empty and
                 geom.area != 0) or
                not isinstance(geom, sgeom.Polygon))

    result = list(filter(not_zero_poly, geom_collection))

    return result + other_result_geoms
