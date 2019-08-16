# (C) British Crown Copyright 2011 - 2019, Met Office
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
This module pulls together proj, GEOS and ``_crs.pyx`` to implement a function
to project a `~shapely.geometry.LinearRing` / `~shapely.geometry.LineString`.
In general, this should never be called manually, instead leaving the
processing to be done by the :class:`cartopy.crs.Projection` subclasses.
"""

cimport cython
from libc.math cimport HUGE_VAL, sqrt
from numpy.math cimport isfinite, isnan
from libc.stdint cimport uintptr_t as ptr
from libcpp cimport bool
from libcpp.list cimport list
from libcpp.vector cimport vector

cdef bool DEBUG = False

cdef extern from "geos_c.h":
    ctypedef void *GEOSContextHandle_t
    ctypedef struct GEOSGeometry:
        pass
    ctypedef struct GEOSCoordSequence
    ctypedef struct GEOSPreparedGeometry
    GEOSCoordSequence *GEOSCoordSeq_create_r(GEOSContextHandle_t, unsigned int, unsigned int) nogil
    GEOSGeometry *GEOSGeom_createPoint_r(GEOSContextHandle_t, GEOSCoordSequence *) nogil
    GEOSGeometry *GEOSGeom_createLineString_r(GEOSContextHandle_t, GEOSCoordSequence *) nogil
    GEOSGeometry *GEOSGeom_createCollection_r(GEOSContextHandle_t, int, GEOSGeometry **, unsigned int) nogil
    GEOSGeometry *GEOSGeom_createEmptyCollection_r(GEOSContextHandle_t, int) nogil
    void GEOSGeom_destroy_r(GEOSContextHandle_t, GEOSGeometry *) nogil
    GEOSCoordSequence *GEOSGeom_getCoordSeq_r(GEOSContextHandle_t, GEOSGeometry *) nogil
    int GEOSCoordSeq_getSize_r(GEOSContextHandle_t handle, const GEOSCoordSequence* s, unsigned int *size) nogil
    int GEOSCoordSeq_getX_r(GEOSContextHandle_t, GEOSCoordSequence *, int, double *) nogil
    int GEOSCoordSeq_getY_r(GEOSContextHandle_t, GEOSCoordSequence *, int, double *) nogil
    int GEOSCoordSeq_setX_r(GEOSContextHandle_t, GEOSCoordSequence *, int, double) nogil
    int GEOSCoordSeq_setY_r(GEOSContextHandle_t, GEOSCoordSequence *, int, double) nogil
    const GEOSPreparedGeometry *GEOSPrepare_r(GEOSContextHandle_t handle, const GEOSGeometry* g) nogil
    char GEOSPreparedCovers_r(GEOSContextHandle_t, const GEOSPreparedGeometry*, const GEOSGeometry*) nogil
    char GEOSPreparedDisjoint_r(GEOSContextHandle_t, const GEOSPreparedGeometry*, const GEOSGeometry*) nogil
    void GEOSPreparedGeom_destroy_r(GEOSContextHandle_t handle, const GEOSPreparedGeometry* g) nogil
    cdef int GEOS_MULTILINESTRING

from cartopy._crs cimport CRS
from ._proj4 cimport (projPJ, projLP, pj_get_spheroid_defn, pj_transform,
                      pj_strerrno, DEG_TO_RAD)
from .geodesic._geodesic cimport (geod_geodesic, geod_geodesicline,
                                  geod_init, geod_geninverse,
                                  geod_lineinit, geod_genposition,
                                  GEOD_ARCMODE, GEOD_LATITUDE, GEOD_LONGITUDE)


import shapely.geometry as sgeom
from shapely.geos import lgeos


cdef GEOSContextHandle_t get_geos_context_handle():
    cdef ptr handle = lgeos.geos_handle
    return <GEOSContextHandle_t>handle


cdef GEOSGeometry *geos_from_shapely(shapely_geom) except *:
    """Get the GEOS pointer from the given shapely geometry."""
    cdef ptr geos_geom = shapely_geom._geom
    return <GEOSGeometry *>geos_geom


cdef shapely_from_geos(GEOSGeometry *geom):
    """Turn the given GEOS geometry pointer into a shapely geometry."""
    return sgeom.base.geom_factory(<ptr>geom)


ctypedef struct Point:
    double x
    double y

ctypedef list[Point] Line


cdef bool degenerate_line(const Line &value):
    return value.size() < 2


cdef bool close(double a, double b):
    return abs(a - b) <= (1e-8 + 1e-5 * abs(b))


@cython.final
cdef class LineAccumulator:
    cdef list[Line] lines

    def __init__(self):
        self.new_line()

    cdef void new_line(self):
        cdef Line line
        self.lines.push_back(line)

    cdef void add_point(self, const Point &point):
        self.lines.back().push_back(point)

    cdef void add_point_if_empty(self, const Point &point):
        if self.lines.back().empty():
            self.add_point(point)

    cdef GEOSGeometry *as_geom(self, GEOSContextHandle_t handle):
        from cython.operator cimport dereference, preincrement
        # self.lines.remove_if(degenerate_line) is not available in Cython.
        cdef list[Line].iterator it = self.lines.begin()
        while it != self.lines.end():
            if degenerate_line(dereference(it)):
                it = self.lines.erase(it)
            else:
                preincrement(it)

        cdef Point first, last
        if self.lines.size() > 1:
            first = self.lines.front().front()
            last = self.lines.back().back()
            if close(first.x, last.x) and close(first.y, last.y):
                self.lines.front().pop_front()
                self.lines.back().splice(self.lines.back().end(),
                                         self.lines.front())
                self.lines.pop_front()

        cdef Line ilines
        cdef Point ipoints
        cdef vector[GEOSGeometry *] geoms
        cdef int i
        cdef GEOSCoordSequence *coords
        for ilines in self.lines:
            coords = GEOSCoordSeq_create_r(handle, ilines.size(), 2)
            for i, ipoints in enumerate(ilines):
                GEOSCoordSeq_setX_r(handle, coords, i, ipoints.x)
                GEOSCoordSeq_setY_r(handle, coords, i, ipoints.y)

            geoms.push_back(GEOSGeom_createLineString_r(handle, coords))

        cdef GEOSGeometry *geom
        if geoms.empty():
            geom = GEOSGeom_createEmptyCollection_r(handle,
                                                    GEOS_MULTILINESTRING)
        else:
            geom = GEOSGeom_createCollection_r(handle, GEOS_MULTILINESTRING,
                                               &geoms[0], geoms.size())
        return geom

    cdef size_t size(self):
        return self.lines.size()


cdef class Interpolator:
    cdef Point start
    cdef Point end
    cdef projPJ src_proj
    cdef projPJ dest_proj

    cdef void init(self, projPJ src_proj, projPJ dest_proj):
        self.src_proj = src_proj
        self.dest_proj = dest_proj

    cdef void set_line(self, const Point &start, const Point &end):
        self.start = start
        self.end = end

    cdef Point interpolate(self, double t):
        raise NotImplementedError

    cdef Point project(self, const Point &point):
        raise NotImplementedError


cdef class CartesianInterpolator(Interpolator):
    cdef Point interpolate(self, double t):
        cdef Point xy
        xy.x = self.start.x + (self.end.x - self.start.x) * t
        xy.y = self.start.y + (self.end.y - self.start.y) * t
        return self.project(xy)

    cdef Point project(self, const Point &src_xy):
        cdef Point dest_xy
        cdef projLP xy

        xy.u = src_xy.x
        xy.v = src_xy.y

        cdef int status = pj_transform(self.src_proj, self.dest_proj,
                                       1, 1, &xy.u, &xy.v, NULL)
        if status in (-14, -20):
            # -14 => "latitude or longitude exceeded limits"
            # -20 => "tolerance condition error"
            xy.u = xy.v = HUGE_VAL
        elif status != 0:
            raise Exception('pj_transform failed: %d\n%s' % (
                status,
                pj_strerrno(status)))

        dest_xy.x = xy.u
        dest_xy.y = xy.v
        return dest_xy


cdef class SphericalInterpolator(Interpolator):
    cdef geod_geodesic geod
    cdef geod_geodesicline geod_line
    cdef double a13

    cdef void init(self, projPJ src_proj, projPJ dest_proj):
        self.src_proj = src_proj
        self.dest_proj = dest_proj

        cdef double major_axis
        cdef double eccentricity_squared
        pj_get_spheroid_defn(self.src_proj, &major_axis, &eccentricity_squared)
        geod_init(&self.geod, major_axis, 1 - sqrt(1 - eccentricity_squared))

    cdef void set_line(self, const Point &start, const Point &end):
        cdef double azi1
        self.a13 = geod_geninverse(&self.geod,
                                   start.y, start.x, end.y, end.x,
                                   NULL, &azi1, NULL, NULL, NULL, NULL, NULL)
        geod_lineinit(&self.geod_line, &self.geod, start.y, start.x, azi1,
                      GEOD_LATITUDE | GEOD_LONGITUDE);

    cdef Point interpolate(self, double t):
        cdef Point lonlat

        geod_genposition(&self.geod_line, GEOD_ARCMODE, self.a13 * t,
                         &lonlat.y, &lonlat.x, NULL, NULL, NULL, NULL, NULL,
                         NULL)

        return self.project(lonlat)

    cdef Point project(self, const Point &lonlat):
        cdef Point xy
        cdef projLP dest

        dest.u = lonlat.x * DEG_TO_RAD
        dest.v = lonlat.y * DEG_TO_RAD

        cdef int status = pj_transform(self.src_proj, self.dest_proj,
                                       1, 1, &dest.u, &dest.v, NULL)
        if status in (-14, -20):
            # -14 => "latitude or longitude exceeded limits"
            # -20 => "tolerance condition error"
            dest.u = dest.v = HUGE_VAL
        elif status != 0:
            raise Exception('pj_transform failed: %d\n%s' % (
                status,
                pj_strerrno(status)))

        xy.x = dest.u
        xy.y = dest.v
        return xy


cdef enum State:
    POINT_IN = 1,
    POINT_OUT,
    POINT_NAN


cdef State get_state(const Point &point, const GEOSPreparedGeometry *gp_domain,
                     GEOSContextHandle_t handle):
    cdef State state
    cdef GEOSCoordSequence *coords
    cdef GEOSGeometry *g_point

    if isfinite(point.x) and isfinite(point.y):
        # TODO: Avoid create-destroy
        coords = GEOSCoordSeq_create_r(handle, 1, 2)
        GEOSCoordSeq_setX_r(handle, coords, 0, point.x)
        GEOSCoordSeq_setY_r(handle, coords, 0, point.y)
        g_point = GEOSGeom_createPoint_r(handle, coords)
        state = (POINT_IN
                 if GEOSPreparedCovers_r(handle, gp_domain, g_point)
                 else POINT_OUT)
        GEOSGeom_destroy_r(handle, g_point)
    else:
        state = POINT_NAN
    return state


@cython.cdivision(True)  # Want divide-by-zero to produce NaN.
cdef bool straightAndDomain(double t_start, const Point &p_start,
                            double t_end, const Point &p_end,
                            Interpolator interpolator, double threshold,
                            GEOSContextHandle_t handle,
                            const GEOSPreparedGeometry *gp_domain,
                            bool inside):
    """
    Return whether the given line segment is suitable as an
    approximation of the projection of the source line.

    t_start: Interpolation parameter for the start point.
    p_start: Projected start point.
    t_end: Interpolation parameter for the end point.
    p_start: Projected end point.
    interpolator: Interpolator for current source line.
    threshold: Lateral tolerance in target projection coordinates.
    handle: Thread-local context handle for GEOS.
    gp_domain: Prepared polygon of target map domain.
    inside: Whether the start point is within the map domain.

    """
    # Straight and in-domain (de9im[7] == 'F')
    cdef bool valid
    cdef double t_mid
    cdef Point p_mid
    cdef double seg_dx, seg_dy
    cdef double mid_dx, mid_dy
    cdef double seg_hypot_sq
    cdef double along
    cdef double separation
    cdef double hypot
    cdef GEOSCoordSequence *coords
    cdef GEOSGeometry *g_segment

    # This could be optimised out of the loop.
    if not (isfinite(p_start.x) and isfinite(p_start.y)):
        valid = False
    elif not (isfinite(p_end.x) and isfinite(p_end.y)):
        valid = False
    else:
        # Find the projected mid-point
        t_mid = (t_start + t_end) * 0.5
        p_mid = interpolator.interpolate(t_mid)

        # Determine the closest point on the segment to the midpoint, in
        # normalized coordinates. We could use GEOSProjectNormalized_r
        # here, but since it's a single line segment, it's much easier to
        # just do the math ourselves:
        #     ○̩ (x1, y1) (assume that this is not necessarily vertical)
        #     │
        #     │   D
        #    ╭├───────○ (x, y)
        #    ┊│┘     ╱
        #    ┊│     ╱
        #    ┊│    ╱
        #    L│   ╱
        #    ┊│  ╱
        #    ┊│θ╱
        #    ┊│╱
        #    ╰̍○̍
        #  (x0, y0)
        # The angle θ can be found by arctan2:
        #     θ = arctan2(y1 - y0, x1 - x0) - arctan2(y - y0, x - x0)
        # and the projection onto the line is simply:
        #     L = hypot(x - x0, y - y0) * cos(θ)
        # with the normalized form being:
        #     along = L / hypot(x1 - x0, y1 - y0)
        #
        # Plugging those into SymPy and .expand().simplify(), we get the
        # following equations (with a slight refactoring to reuse some
        # intermediate values):
        seg_dx = p_end.x - p_start.x
        seg_dy = p_end.y - p_start.y
        mid_dx = p_mid.x - p_start.x
        mid_dy = p_mid.y - p_start.y
        seg_hypot_sq = seg_dx*seg_dx + seg_dy*seg_dy

        along = (seg_dx*mid_dx + seg_dy*mid_dy) / seg_hypot_sq

        if isnan(along):
            valid = True
        else:
            valid = 0.0 < along < 1.0
            if valid:
                # For the distance of the point from the line segment, using
                # the same geometry above, use sin instead of cos:
                #     D = hypot(x - x0, y - y0) * sin(θ)
                # and then simplify with SymPy again:
                separation = (abs(mid_dx*seg_dy - mid_dy*seg_dx) /
                              sqrt(seg_hypot_sq))
                if inside:
                    # Scale the lateral threshold by the distance from
                    # the nearest end. I.e. Near the ends the lateral
                    # threshold is much smaller; it only has its full
                    # value in the middle.
                    valid = (separation <=
                             threshold * 2.0 * (0.5 - abs(0.5 - along)))
                else:
                    # Check if the mid-point makes less than ~11 degree
                    # angle with the straight line.
                    # sin(11') => 0.2
                    # To save the square-root we just use the square of
                    # the lengths, hence:
                    # 0.2 ^ 2 => 0.04
                    hypot = mid_dx*mid_dx + mid_dy*mid_dy
                    valid = ((separation * separation) / hypot) < 0.04

        if valid:
            # TODO: Re-use geometries, instead of create-destroy!

            # Create a LineString for the current end-point.
            coords = GEOSCoordSeq_create_r(handle, 2, 2)
            GEOSCoordSeq_setX_r(handle, coords, 0, p_start.x)
            GEOSCoordSeq_setY_r(handle, coords, 0, p_start.y)
            GEOSCoordSeq_setX_r(handle, coords, 1, p_end.x)
            GEOSCoordSeq_setY_r(handle, coords, 1, p_end.y)
            g_segment = GEOSGeom_createLineString_r(handle, coords)

            if inside:
                valid = GEOSPreparedCovers_r(handle, gp_domain, g_segment)
            else:
                valid = GEOSPreparedDisjoint_r(handle, gp_domain, g_segment)

            GEOSGeom_destroy_r(handle, g_segment)

    return valid


cdef void bisect(double t_start, const Point &p_start, const Point &p_end,
                 GEOSContextHandle_t handle,
                 const GEOSPreparedGeometry *gp_domain, const State &state,
                 Interpolator interpolator, double threshold,
                 double &t_min, Point &p_min, double &t_max, Point &p_max):
    cdef double t_current
    cdef Point p_current
    cdef bool valid

    # Initialise our bisection range to the start and end points.
    (&t_min)[0] = t_start
    (&p_min)[0] = p_start
    (&t_max)[0] = 1.0
    (&p_max)[0] = p_end

    # Start the search at the end.
    t_current = t_max
    p_current = p_max

    # TODO: See if we can convert the 't' threshold into one based on the
    # projected coordinates - e.g. the resulting line length.

    while abs(t_max - t_min) > 1.0e-6:
        if DEBUG:
            print("t: ", t_current)

        if state == POINT_IN:
            # Straight and entirely-inside-domain
            valid = straightAndDomain(t_start, p_start, t_current, p_current,
                                      interpolator, threshold,
                                      handle, gp_domain, True)

        elif state == POINT_OUT:
            # Straight and entirely-outside-domain
            valid = straightAndDomain(t_start, p_start, t_current, p_current,
                                      interpolator, threshold,
                                      handle, gp_domain, False)
        else:
            valid = not isfinite(p_current.x) or not isfinite(p_current.y)

        if DEBUG:
            print("   => valid: ", valid)

        if valid:
            (&t_min)[0] = t_current
            (&p_min)[0] = p_current
        else:
            (&t_max)[0] = t_current
            (&p_max)[0] = p_current

        t_current = (t_min + t_max) * 0.5
        p_current = interpolator.interpolate(t_current)


cdef void _project_segment(GEOSContextHandle_t handle,
                           const GEOSCoordSequence *src_coords,
                           unsigned int src_idx_from, unsigned int src_idx_to,
                           Interpolator interpolator,
                           const GEOSPreparedGeometry *gp_domain,
                           double threshold, LineAccumulator lines):
    cdef Point p_current, p_min, p_max, p_end
    cdef double t_current, t_min, t_max
    cdef State state

    GEOSCoordSeq_getX_r(handle, src_coords, src_idx_from, &p_current.x)
    GEOSCoordSeq_getY_r(handle, src_coords, src_idx_from, &p_current.y)
    GEOSCoordSeq_getX_r(handle, src_coords, src_idx_to, &p_end.x)
    GEOSCoordSeq_getY_r(handle, src_coords, src_idx_to, &p_end.y)
    if DEBUG:
        print("Setting line:")
        print("   ", p_current.x, ", ", p_current.y)
        print("   ", p_end.x, ", ", p_end.y)

    interpolator.set_line(p_current, p_end)
    p_current = interpolator.project(p_current)
    p_end = interpolator.project(p_end)
    if DEBUG:
        print("Projected as:")
        print("   ", p_current.x, ", ", p_current.y)
        print("   ", p_end.x, ", ", p_end.y)

    t_current = 0.0
    state = get_state(p_current, gp_domain, handle)

    cdef size_t old_lines_size = lines.size()
    while t_current < 1.0 and (lines.size() - old_lines_size) < 100:
        if DEBUG:
            print("Bisecting from: ", t_current, " (")
            if state == POINT_IN:
                print("IN")
            elif state == POINT_OUT:
                print("OUT")
            else:
                print("NAN")
            print(")")
            print("   ", p_current.x, ", ", p_current.y)
            print("   ", p_end.x, ", ", p_end.y)

        bisect(t_current, p_current, p_end, handle, gp_domain, state,
               interpolator, threshold,
               t_min, p_min, t_max, p_max)
        if DEBUG:
            print("   => ", t_min, "to", t_max)
            print("   => (", p_min.x, ", ", p_min.y, ") to (",
                  p_max.x, ", ", p_max.y, ")")

        if state == POINT_IN:
            lines.add_point_if_empty(p_current)
            if t_min != t_current:
                lines.add_point(p_min)
                t_current = t_min
                p_current = p_min
            else:
                t_current = t_max
                p_current = p_max
                state = get_state(p_current, gp_domain, handle)
                if state == POINT_IN:
                    lines.new_line()

        elif state == POINT_OUT:
            if t_min != t_current:
                t_current = t_min
                p_current = p_min
            else:
                t_current = t_max
                p_current = p_max
                state = get_state(p_current, gp_domain, handle)
                if state == POINT_IN:
                    lines.new_line()

        else:
            t_current = t_max
            p_current = p_max
            state = get_state(p_current, gp_domain, handle)
            if state == POINT_IN:
                lines.new_line()


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
        Interpolator interpolator
        GEOSGeometry *g_domain
        const GEOSCoordSequence *src_coords
        unsigned int src_size, src_idx
        const GEOSPreparedGeometry *gp_domain
        LineAccumulator lines
        GEOSGeometry *g_multi_line_string

    g_domain = geos_from_shapely(dest_projection.domain)

    if src_crs.is_geodetic():
        interpolator = SphericalInterpolator()
    else:
        interpolator = CartesianInterpolator()
    interpolator.init(src_crs.proj4, (<CRS>dest_projection).proj4)

    src_coords = GEOSGeom_getCoordSeq_r(handle, g_linear)
    gp_domain = GEOSPrepare_r(handle, g_domain)

    GEOSCoordSeq_getSize_r(handle, src_coords, &src_size)  # check exceptions

    lines = LineAccumulator()
    for src_idx in range(1, src_size):
        _project_segment(handle, src_coords, src_idx - 1, src_idx,
                         interpolator, gp_domain, threshold, lines);

    GEOSPreparedGeom_destroy_r(handle, gp_domain)

    g_multi_line_string = lines.as_geom(handle)

    del lines, interpolator
    multi_line_string = shapely_from_geos(g_multi_line_string)
    return multi_line_string
