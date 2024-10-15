# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Extra functionality that is primarily intended for developers, providing support for
transforming between Shapely geometries and Matplotlib paths.

See also `Shapely Geometric Objects
<https://shapely.readthedocs.io/en/latest/manual.html#geometric-objects>`_
and `Matplotlib Path API <https://matplotlib.org/stable/api/path_api.html>`_.

"""

import warnings

from matplotlib.path import Path
import numpy as np
import shapely.geometry as sgeom

from cartopy.mpl import _MPL_38


def _ensure_path_closed(path):
    """
    Method to ensure that a path contains only closed sub-paths.

    Parameters
    ----------
    path
        A :class:`matplotlib.path.Path` instance.

    Returns
    -------
    path
        A :class:`matplotlib.path.Path` instance with only closed polygons.

    """
    # Split path into potential sub-paths and close all polygons
    # (explicitly disable path simplification applied in to_polygons)
    should_simplify = path.should_simplify
    try:
        path.should_simplify = False
        polygons = path.to_polygons()
    finally:
        path.should_simplify = should_simplify

    codes, vertices = [], []
    for poly in polygons:
        vertices.extend([poly[0], *poly])
        codes.extend([Path.MOVETO, *[Path.LINETO]*(len(poly) - 1), Path.CLOSEPOLY])

    return Path(vertices, codes)

def geos_to_path(shape):
    """
    Create a list of :class:`matplotlib.path.Path` objects that describe
    a shape.

    .. deprecated:: 0.25
       Use `shapely_to_path` instead.

    Parameters
    ----------
    shape
        A list, tuple or single instance of any of the following
        types: :class:`shapely.geometry.point.Point`,
        :class:`shapely.geometry.linestring.LineString`,
        :class:`shapely.geometry.polygon.LinearRing`,
        :class:`shapely.geometry.polygon.Polygon`,
        :class:`shapely.geometry.multipoint.MultiPoint`,
        :class:`shapely.geometry.multipolygon.MultiPolygon`,
        :class:`shapely.geometry.multilinestring.MultiLineString`,
        :class:`shapely.geometry.collection.GeometryCollection`,
        or any type with a _as_mpl_path() method.

    Returns
    -------
    paths
        A list of :class:`matplotlib.path.Path` objects.

    """
    warnings.warn("geos_to_path is deprecated and will be removed in a future release."
                  "  Use shapely_to_path instead.",
                  DeprecationWarning, stacklevel=2)
    if isinstance(shape, (list, tuple)):
        paths = []
        for shp in shape:
            paths.extend(geos_to_path(shp))
        return paths

    if isinstance(shape, sgeom.LinearRing):
        return [Path(np.column_stack(shape.xy), closed=True)]
    elif isinstance(shape, (sgeom.LineString, sgeom.Point)):
        return [Path(np.column_stack(shape.xy))]
    elif isinstance(shape, sgeom.Polygon):
        def poly_codes(poly):
            codes = np.ones(len(poly.xy[0])) * Path.LINETO
            codes[0] = Path.MOVETO
            codes[-1] = Path.CLOSEPOLY
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
        raise ValueError(f'Unsupported shape type {type(shape)}.')


def path_segments(path, **kwargs):
    """
    Create an array of vertices and a corresponding array of codes from a
    :class:`matplotlib.path.Path`.

    Parameters
    ----------
    path
        A :class:`matplotlib.path.Path` instance.

    Other Parameters
    ----------------
    kwargs
        See :func:`matplotlib.path.iter_segments` for details of the keyword
        arguments.

    Returns
    -------
    vertices, codes
        A (vertices, codes) tuple, where vertices is a numpy array of
        coordinates, and codes is a numpy array of matplotlib path codes.
        See :class:`matplotlib.path.Path` for information on the types of
        codes and their meanings.

    """
    pth = path.cleaned(**kwargs)
    return pth.vertices[:-1, :], pth.codes[:-1]


def path_to_geos(path, force_ccw=False):
    """
    Create a list of Shapely geometric objects from a
    :class:`matplotlib.path.Path`.

    .. deprecated:: 0.25
       Use `path_to_shapely` instead.

    Parameters
    ----------
    path
        A :class:`matplotlib.path.Path` instance.

    Other Parameters
    ----------------
    force_ccw
        Boolean flag determining whether the path can be inverted to enforce
        ccw. Defaults to False.

    Returns
    -------
    A list of instances of the following type(s):
        :class:`shapely.geometry.polygon.Polygon`,
        :class:`shapely.geometry.linestring.LineString` and/or
        :class:`shapely.geometry.multilinestring.MultiLineString`.

    """
    warnings.warn("path_to_geos is deprecated and will be removed in a future release."
                  "  Use path_to_shapely instead.",
                  DeprecationWarning, stacklevel=2)
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

        if path_codes[-1] == Path.CLOSEPOLY:
            path_verts[-1, :] = path_verts[0, :]

        verts_same_as_first = np.isclose(path_verts[0, :], path_verts[1:, :],
                                         rtol=1e-10, atol=1e-13)
        verts_same_as_first = np.logical_and.reduce(verts_same_as_first,
                                                    axis=1)

        if all(verts_same_as_first):
            geom = sgeom.Point(path_verts[0, :])
        elif path_verts.shape[0] > 4 and path_codes[-1] == Path.CLOSEPOLY:
            geom = sgeom.Polygon(path_verts[:-1, :])
        else:
            geom = sgeom.LineString(path_verts)

        # If geom is a Polygon and is contained within the last geom in
        # collection, it usually needs to be an interior to that geom (e.g. a
        # lake within a land mass).  Sometimes there is a further geom within
        # this interior (e.g. an island in a lake, or some instances of
        # contours).  This needs to be a new external geom in the collection.
        if geom.is_empty:
            pass
        elif (len(collection) > 0 and
                isinstance(collection[-1][0], sgeom.Polygon) and
                isinstance(geom, sgeom.Polygon) and
                collection[-1][0].contains(geom.exterior)):
            if any(internal.contains(geom) for internal in collection[-1][1]):
                collection.append((geom, []))
            else:
                collection[-1][1].append(geom)
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
            exteriors = [geom.exterior for geom in internal_polys]
            geom = sgeom.Polygon(external_geom.exterior, exteriors)
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
        return ((isinstance(geom, sgeom.Polygon) and not geom.is_empty and
                 geom.area != 0) or
                not isinstance(geom, sgeom.Polygon))

    result = list(filter(not_zero_poly, geom_collection))

    return result + other_result_geoms


def path_to_shapely(path, force_ccw=False):
    """
    Create a Shapely geometric objects from a :class:`matplotlib.path.Path`.

    Parameters
    ----------
    path
        A :class:`matplotlib.path.Path` instance.

    Other Parameters
    ----------------
    force_ccw
        Boolean flag determining whether the path can be inverted to enforce
        ccw. Defaults to False.

    Returns
    -------
    One of the following Shapely objects:

        :class:`shapely.geometry.polygon.Polygon`,
        :class:`shapely.geometry.linestring.LineString`
        :class:`shapely.geometry.point.Point`,
        :class:`shapely.geometry.multipolygon.MultiPolygon`
        :class:`shapely.geometry.multilinestring.MultiLineString`
        :class:`shapely.geometry.multipoint.MultiPoint`
        :class:`shapely.geometry.collection.GeometryCollection`.

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
    # (external_poly, [internal_polygons]) tuples for the polygons and separate
    # lists for linestrings and points.
    points = []
    linestrings = []
    polygon_bits = []
    for path_verts, path_codes in zip(verts_split, codes_split):
        if len(path_verts) == 0:
            continue

        if path_codes[-1] == Path.CLOSEPOLY:
            path_verts[-1, :] = path_verts[0, :]

        verts_same_as_first = np.isclose(path_verts[0, :], path_verts[1:, :],
                                         rtol=1e-10, atol=1e-13)
        verts_same_as_first = np.logical_and.reduce(verts_same_as_first,
                                                    axis=1)

        if all(verts_same_as_first):
            points.append(sgeom.Point(path_verts[0, :]))
        elif not(path_verts.shape[0] > 4 and path_codes[-1] == Path.CLOSEPOLY):
            linestrings.append(sgeom.LineString(path_verts))
        else:
            geom = sgeom.Polygon(path_verts[:-1, :])
            # If geom is a Polygon and is contained within the last geom in
            # polygon_bits, it usually needs to be an interior to that geom (e.g. a
            # lake within a land mass).  Sometimes there is a further geom within
            # this interior (e.g. an island in a lake, or some instances of
            # contours).  This needs to be a new external geom in polygon_bits.
            if (len(polygon_bits) > 0 and polygon_bits[-1][0].contains(geom.exterior)):
                if any(internal.contains(geom) for internal in polygon_bits[-1][1]):
                    polygon_bits.append((geom, []))
                else:
                    polygon_bits[-1][1].append(geom)
            else:
                polygon_bits.append((geom, []))

    # Convert each (external_polygon, [internal_polygons]) pair into a
    # a shapely Polygon that encapsulates the internal polygons.
    polygons = []
    for external_poly, internal_polys in polygon_bits:
        if internal_polys:
            exteriors = [geom.exterior for geom in internal_polys]
            geom = sgeom.Polygon(external_poly.exterior, exteriors)
        else:
            geom = external_poly

        # Correctly orientate the polygon (ccw)
        if force_ccw and not geom.exterior.is_ccw:
            geom = sgeom.polygon.orient(geom)

        polygons.append(geom)

    # Remove any zero area Polygons
    def not_zero_poly(geom):
        return (not geom.is_empty and geom.area != 0)

    polygons = list(filter(not_zero_poly, polygons))

    # Figure out what type of object to return
    if not polygons:
        if not linestrings:
            if not points:
                # No geometries.  Return an empty point
                return sgeom.Point()
            elif len(points) > 1:
                return sgeom.MultiPoint(points)
            else:
                return points[0]
        elif not points:
            if len(linestrings) > 1:
                return sgeom.MultiLineString(linestrings)
            else:
                return linestrings[0]
    else:
        if not linestrings and not points:
            if len(polygons) > 1:
                return sgeom.MultiPolygon(polygons)
            else:
                return polygons[0]

    # If we got to here, we have at least two types of geometry, so return
    # a geometry collection.
    return sgeom.GeometryCollection(polygons + linestrings + points)


def shapely_to_path(shape):
    """
    Create a :class:`matplotlib.path.Path` object that describes a shape.

    Parameters
    ----------
    shape
        :class:`shapely.geometry.point.Point`,
        :class:`shapely.geometry.linestring.LineString`,
        :class:`shapely.geometry.polygon.LinearRing`,
        :class:`shapely.geometry.polygon.Polygon`,
        :class:`shapely.geometry.multipoint.MultiPoint`,
        :class:`shapely.geometry.multipolygon.MultiPolygon`,
        :class:`shapely.geometry.multilinestring.MultiLineString`,
        :class:`shapely.geometry.collection.GeometryCollection`,
        or any type with a _as_mpl_path() method.

    Returns
    -------
    path
        :class:`matplotlib.path.Path`

    """
    if shape.is_empty:
        return Path(np.empty([0, 2]))
    elif isinstance(shape, sgeom.LinearRing):
        return Path(np.column_stack(shape.xy), closed=True)
    elif isinstance(shape, (sgeom.LineString, sgeom.Point)):
        return Path(np.column_stack(shape.xy))
    elif isinstance(shape, sgeom.Polygon):
        def poly_codes(poly):
            codes = np.ones(len(poly.xy[0])) * Path.LINETO
            codes[0] = Path.MOVETO
            codes[-1] = Path.CLOSEPOLY
            return codes
        vertices = np.concatenate([np.array(shape.exterior.xy)] +
                                  [np.array(ring.xy) for ring in
                                   shape.interiors], 1).T
        codes = np.concatenate([poly_codes(shape.exterior)] +
                               [poly_codes(ring) for ring in shape.interiors])
        return Path(vertices, codes)
    elif isinstance(shape, (sgeom.MultiPolygon, sgeom.GeometryCollection,
                            sgeom.MultiLineString, sgeom.MultiPoint)):
        paths = []
        for geom in shape.geoms:
            path = shapely_to_path(geom)
            if _MPL_38 or path.vertices.size > 0:
                paths.append(path)
        return Path.make_compound_path(*paths)
    elif hasattr(shape, '_as_mpl_path'):
        vertices, codes = shape._as_mpl_path()
        return Path(vertices, codes)
    else:
        raise ValueError(f'Unsupported shape type {type(shape)}.')
