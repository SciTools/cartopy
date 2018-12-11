# (C) British Crown Copyright 2011 - 2018, Met Office
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
#
# cython: embedsignature=True

"""
This module pulls together ``_trace.cpp``, proj, GEOS and ``_crs.pyx`` to
implement a function to project a `~shapely.geometry.LinearRing` /
`~shapely.geometry.LineString`. In general, this should never be called
manually, instead leaving the processing to be done by the
:class:`cartopy.crs.Projection` subclasses.
"""

from libc.stdint cimport uintptr_t as ptr


cdef extern from "geos_c.h":
    ctypedef void *GEOSContextHandle_t
    ctypedef struct GEOSGeometry:
        pass

from cartopy._crs cimport CRS


import shapely.geometry as sgeom
from shapely.geos import lgeos


cdef extern from "proj_api.h":
    ctypedef void *projPJ


cdef extern from "_trace.h":
    cdef cppclass Interpolator:
        pass

    cdef cppclass SphericalInterpolator:
        SphericalInterpolator(projPJ src_proj, projPJ dest_proj)

    cdef cppclass CartesianInterpolator:
        CartesianInterpolator(projPJ src_proj, projPJ dest_proj)

    # XXX Rename? It handles LinearRings too.
    GEOSGeometry *_project_line_string(GEOSContextHandle_t handle,
                                       GEOSGeometry *g_line_string,
                                       Interpolator *interpolator,
                                       GEOSGeometry *g_domain,
                                       double threshold)


cdef GEOSContextHandle_t get_geos_context_handle():
    cdef ptr handle = lgeos.geos_handle
    return <GEOSContextHandle_t>handle


cdef GEOSGeometry *geos_from_shapely(shapely_geom) except *:
    """Get the GEOS pointer from the given shapely geometry."""
    cdef ptr geos_geom = shapely_geom._geom
    return <GEOSGeometry *>geos_geom


cdef shapely_from_geos(GEOSGeometry *geom):
    """Turn the given GEOS geometry pointer into a shapely geometry."""
    # This is the "correct" way to do it...
    #   return geom_factory(<ptr>geom)
    # ... but it's quite slow, so we do it by hand.
    multi_line_string = sgeom.base.BaseGeometry()
    multi_line_string.__class__ = sgeom.MultiLineString
    multi_line_string.__geom__ = <ptr>geom
    multi_line_string.__parent__ = None
    multi_line_string._ndim = 2
    return multi_line_string


def project_linear(geometry not None, CRS src_crs not None,
                   dest_projection not None):
    """
    Project a geometry from one projection to another.

    Parameters
    ----------
    geometry : `shapely.geometry.LineString` or `shapely.geometry.LinearRing`
        A geometry to be projected.
    src_crs : cartopy.crs.CRS
        The coordinate system of the line to be projected.
    dest_projection : cartopy.crs.Projection
        The projection for the resulting projected line.

    Returns
    -------
    `shapely.geometry.MultiLineString`
        The result of projecting the given geometry from the source projection
        into the destination projection.

    """
    cdef:
        double threshold = dest_projection.threshold
        GEOSContextHandle_t handle = get_geos_context_handle()
        GEOSGeometry *g_linear = geos_from_shapely(geometry)
        Interpolator *interpolator
        GEOSGeometry *g_domain
        GEOSGeometry *g_multi_line_string

    g_domain = geos_from_shapely(dest_projection.domain)

    if src_crs.is_geodetic():
        interpolator = <Interpolator *>new SphericalInterpolator(
                src_crs.proj4, (<CRS>dest_projection).proj4)
    else:
        interpolator = <Interpolator *>new CartesianInterpolator(
                src_crs.proj4, (<CRS>dest_projection).proj4)

    g_multi_line_string = _project_line_string(handle, g_linear,
                                               interpolator, g_domain, threshold)
    del interpolator
    multi_line_string = shapely_from_geos(g_multi_line_string)
    return multi_line_string
