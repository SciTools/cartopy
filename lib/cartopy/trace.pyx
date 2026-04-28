# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
#
# cython: embedsignature=True

"""
Trace pulls together proj, GEOS and ``_crs.pyx`` to implement a function
to project a `~shapely.geometry.LinearRing` / `~shapely.geometry.LineString`.
In general, this should never be called manually, instead leaving the
processing to be done by the :class:`cartopy.crs.Projection` subclasses.
"""
from functools import lru_cache

cimport cython
from libc.math cimport HUGE_VAL, sqrt, isfinite, isnan
from libcpp cimport bool
from libcpp.list cimport list

cdef bool DEBUG = False

from cython.operator cimport dereference, preincrement
import numpy as np
import shapely
import shapely.geometry as sgeom
import shapely.lib as _shapely_lib
import shapely.prepared as sprep
from pyproj import Geod, Transformer, proj_version_str

# Module-level reusable buffers for fast segment construction in straightAndDomain.
# Using shapely.lib.set_coordinates on a pre-built 0-d array avoids the overhead.
# NOTE: these buffers are written on every call to straightAndDomain and are
# therefore not thread-safe.
_seg_geom_arr = np.asarray(sgeom.LineString([(0.0, 0.0), (1.0, 1.0)]))
_seg_coords_buf = np.empty((2, 2), dtype=np.float64)


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

    cdef object as_geom(self):
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
        cdef unsigned int n, j
        geoms = []
        for ilines in self.lines:
            # Fill a pre-allocated numpy array then hand it straight to
            # shapely.linestrings() so no Python tuple objects are created.
            n = ilines.size()
            coords_arr = np.empty((n, 2), dtype=np.float64)
            j = 0
            for ipoints in ilines:
                coords_arr[j, 0] = ipoints.x
                coords_arr[j, 1] = ipoints.y
                j += 1
            geoms.append(shapely.linestrings(coords_arr))

        geom = sgeom.MultiLineString(geoms)
        return geom

    cdef object as_ring_fragments(self):
        """Return an ordered list of LineStrings

        Representing projected ring fragments in ring-traversal order.
        The caller is responsible for closing each fragment pair with a
        boundary arc. Unlike as_geom(), this method does not attempt to join
        the first and last fragments.
        """
        cdef list[Line].iterator it = self.lines.begin()
        while it != self.lines.end():
            if degenerate_line(dereference(it)):
                it = self.lines.erase(it)
            else:
                preincrement(it)

        cdef Line ilines
        cdef Point ipoints
        cdef unsigned int n, j
        result = []
        for ilines in self.lines:
            n = ilines.size()
            coords_arr = np.empty((n, 2), dtype=np.float64)
            j = 0
            for ipoints in ilines:
                coords_arr[j, 0] = ipoints.x
                coords_arr[j, 1] = ipoints.y
                j += 1
            result.append(shapely.linestrings(coords_arr))
        return result

    cdef size_t size(self):
        return self.lines.size()

cdef class Interpolator:
    cdef Point start
    cdef Point end
    cdef readonly transformer
    cdef double src_scale
    cdef double dest_scale
    cdef bint to_180
    # Minimum number of source vertices for batch midpoint precomputation
    cdef int batch_threshold

    def __cinit__(self):
        self.src_scale = 1
        self.dest_scale = 1
        self.to_180 = False
        # emprically found, may change in the future
        self.batch_threshold = 22

    cdef void init(self, src_crs, dest_crs) except *:
        self.transformer = Transformer.from_crs(src_crs, dest_crs, always_xy=True)
        self.to_180 = (
            self.transformer.name == "noop" and
            src_crs.__class__.__name__ in ("PlateCarree", "RotatedPole")
        )

    cdef void set_line(self, const Point &start, const Point &end):
        self.start = start
        self.end = end

    cdef Point project(self, const Point &src_xy) except *:
        cdef Point dest_xy
        cdef double xx, yy

        # errcheck=False: PROJ returns inf (== HUGE_VAL) for out-of-domain
        # points without raising a Python exception. get_state() classifies
        # non-finite values as POINT_NAN, so we don't need the extra error handling
        if self.src_scale != 1.0:
            xx, yy = self.transformer.transform(
                src_xy.x * self.src_scale,
                src_xy.y * self.src_scale,
                errcheck=False
            )
        else:
            xx, yy = self.transformer.transform(src_xy.x, src_xy.y,
                                                errcheck=False)

        # Only wrap finite values; inf/nan from failed projections must not
        # be passed to the modulo operator (inf % 360 == nan)
        if self.to_180 and isfinite(xx) and (xx > 180 or xx < -180):
            xx = (((xx + 180) % 360) - 180)

        if self.dest_scale != 1.0:
            dest_xy.x = xx * self.dest_scale
            dest_xy.y = yy * self.dest_scale
        else:
            dest_xy.x = xx
            dest_xy.y = yy
        return dest_xy

    cdef double[:, :] project_points(self, double[:, :] src_xy) except *:
        # Batch-project an array of source coordinates.
        # errcheck=False: out-of-domain points silently return inf instead of
        # raising ProjError for performance
        src = np.asarray(src_xy)
        if self.src_scale != 1.0:
            src = src * self.src_scale

        xx, yy = self.transformer.transform(src[:, 0], src[:, 1],
                                            errcheck=False)

        result = np.empty((len(xx), 2), dtype=np.float64)
        result[:, 0] = xx
        result[:, 1] = yy

        if self.to_180:
            # Only wrap finite values; inf % 360 == nan and must be excluded.
            wrap_locs = np.isfinite(result[:, 0]) & (
                (result[:, 0] > 180) | (result[:, 0] < -180)
            )
            result[wrap_locs, 0] = (
                ((result[wrap_locs, 0] + 180) % 360) - 180
            )

        if self.dest_scale != 1.0:
            result *= self.dest_scale
        return result

    cdef Point interpolate(self, double t) except *:
        raise NotImplementedError

    cdef object batch_midpoints(self, double[:, :] src_coords):
        """Project all midpoints at once for performance

        Return a (N-1, 2) float64 array of projected midpoints for all
        adjacent source-coordinate pairs, or None when the batch optimization
        should not be applied (too few segments, or not implemented)."""
        return None


cdef class CartesianInterpolator(Interpolator):
    cdef Point interpolate(self, double t) except *:
        cdef Point xy
        xy.x = self.start.x + (self.end.x - self.start.x) * t
        xy.y = self.start.y + (self.end.y - self.start.y) * t
        return self.project(xy)

    cdef object batch_midpoints(self, double[:, :] src_coords):
        # Arithmetic midpoint is exact for linear interpolation.
        if len(src_coords) < self.batch_threshold:
            return None
        src = np.asarray(src_coords)
        return np.asarray(self.project_points((src[:-1] + src[1:]) * 0.5))


cdef class SphericalInterpolator(Interpolator):
    cdef object geod
    cdef double azim
    cdef double s12
    # Cache populated by batch_midpoints() so that the subsequent
    # per-segment set_line() calls can consume pre-computed geod.inv
    # results instead of re-running the scalar geod.inv each time.
    cdef object _batch_az12   # 1-D float64 ndarray or None
    cdef object _batch_s12    # 1-D float64 ndarray or None
    cdef int    _batch_idx    # next slot to consume; -1 means no live cache

    def __cinit__(self):
        self._batch_idx = -1

    cdef void init(self, src_crs, dest_crs) except *:
        self.transformer = Transformer.from_crs(src_crs, dest_crs, always_xy=True)

        cdef double major_axis = src_crs.ellipsoid.semi_major_metre
        cdef double flattening = 0
        if src_crs.ellipsoid.inverse_flattening > 0:
            flattening = 1 / src_crs.ellipsoid.inverse_flattening
        self.geod = Geod(a=major_axis, f=flattening)

    cdef void set_line(self, const Point &start, const Point &end):
        Interpolator.set_line(self, start, end)
        if self._batch_idx >= 0:
            # Reuse the geod.inv result pre-computed in batch_midpoints().
            self.azim = self._batch_az12[self._batch_idx]
            self.s12  = self._batch_s12[self._batch_idx]
            self._batch_idx += 1
        else:
            self.azim, _, self.s12 = self.geod.inv(start.x, start.y, end.x, end.y)

    cdef Point interpolate(self, double t) except *:
        cdef Point lonlat

        lonlat.x, lonlat.y, _ = self.geod.fwd(self.start.x, self.start.y, self.azim, self.s12 * t)
        return self.project(lonlat)

    cdef object batch_midpoints(self, double[:, :] src_coords):
        # Walk the geodesic to its midpoint for each adjacent pair:
        #   az12, s12  = geod.inv(lon1, lat1, lon2, lat2)
        #   lon_m, lat_m = geod.fwd(lon1, lat1, az12, s12 / 2)
        #
        # Crucially, the az12/s12 arrays are also cached on self so that the
        # subsequent per-segment set_line() calls (one per segment in
        # project_linear's loop) can read the pre-computed value instead of
        # calling geod.inv() a second time for the same pair.

        # Reset cache regardless of whether we proceed; ensures that a short
        # follow-up call does not accidentally consume stale cached values.
        self._batch_idx = -1

        if len(src_coords) < self.batch_threshold:
            return None

        src = np.asarray(src_coords)
        lons1 = src[:-1, 0]
        lats1 = src[:-1, 1]
        lons2 = src[1:, 0]
        lats2 = src[1:, 1]
        az12, _, s12 = self.geod.inv(lons1, lats1, lons2, lats2)

        # Store for set_line() to consume, then compute midpoints.
        self._batch_az12 = az12
        self._batch_s12  = s12
        self._batch_idx  = 0

        mid_lons, mid_lats, _ = self.geod.fwd(lons1, lats1, az12, s12 * 0.5)
        return np.asarray(
            self.project_points(np.stack([mid_lons, mid_lats], axis=-1)))


cdef enum State:
    POINT_IN = 1,
    POINT_OUT,
    POINT_NAN


cdef State get_state(const Point &point, object gp_domain, bool geom_fully_inside=False):
    cdef State state
    if geom_fully_inside:
        # Fast-path return because the geometry is fully inside
        return POINT_IN
    if isfinite(point.x) and isfinite(point.y):
        # Use shapely.lib directly to skip the shapely.predicates decorator
        state = POINT_IN if _shapely_lib.intersects_xy(
            gp_domain.context, point.x, point.y) else POINT_OUT
    else:
        state = POINT_NAN
    return state


@cython.cdivision(True)  # Want divide-by-zero to produce NaN.
cdef bool _check_straight(const Point &p_start, const Point &p_end,
                          const Point &p_mid,
                          double threshold, bool inside):
    """
    Pure geometry check: is the projected midpoint p_mid close enough to the
    straight line from p_start to p_end?

    Assumes p_start and p_end are finite. p_mid may be non-finite (returns
    True for degenerate/zero-length segments where along is NaN).

    Determine the closest point on the segment to the midpoint, in
    normalized coordinates.
        ○̩ (x1, y1) (assume that this is not necessarily vertical)
        │
        │   D
       ╭├───────○ (x, y)
       ┊│┘     ╱
       ┊│     ╱
       ┊│    ╱
       L│   ╱
       ┊│  ╱
       ┊│θ╱
       ┊│╱
       ╰̍○̍
     (x0, y0)
    The angle θ can be found by arctan2:
        θ = arctan2(y1 - y0, x1 - x0) - arctan2(y - y0, x - x0)
    and the projection onto the line is simply:
        L = hypot(x - x0, y - y0) * cos(θ)
    with the normalized form being:
        along = L / hypot(x1 - x0, y1 - y0)

    Plugging those into SymPy and .expand().simplify(), we get the
    following equations (with a slight refactoring to reuse some
    intermediate values):
    """
    cdef double seg_dx, seg_dy
    cdef double mid_dx, mid_dy
    cdef double seg_hypot_sq
    cdef double along, separation, hypot

    seg_dx = p_end.x - p_start.x
    seg_dy = p_end.y - p_start.y
    mid_dx = p_mid.x - p_start.x
    mid_dy = p_mid.y - p_start.y
    seg_hypot_sq = seg_dx*seg_dx + seg_dy*seg_dy

    along = (seg_dx*mid_dx + seg_dy*mid_dy) / seg_hypot_sq

    if isnan(along):
        return True
    if not (0.0 < along < 1.0):
        return False

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
        return separation <= threshold * 2.0 * (0.5 - abs(0.5 - along))
    else:
        # Check if the mid-point makes less than ~11 degree
        # angle with the straight line.
        # sin(11') => 0.2
        # To save the square-root we just use the square of
        # the lengths, hence:
        # 0.2 ^ 2 => 0.04
        hypot = mid_dx*mid_dx + mid_dy*mid_dy
        return ((separation * separation) / hypot) < 0.04


cdef bool straightAndDomain(const Point &p_start, const Point &p_end,
                            const Point &p_mid,
                            double threshold,
                            object gp_domain,
                            bool inside,
                            bool geom_fully_inside=False) except *:
    """
    Return whether the given line segment is suitable as an
    approximation of the projection of the source line.

    p_start: Projected start point.
    p_end:   Projected end point.
    p_mid:   Pre-computed projected midpoint (at t = (t_start + t_end) / 2).
    threshold: Lateral tolerance in target projection coordinates.
    gp_domain: Prepared polygon of target map domain.
    inside: Whether the start point is within the map domain.
    geom_fully_inside: Whether all points are within the map domain.
    """
    cdef bool valid

    if not (isfinite(p_start.x) and isfinite(p_start.y)):
        return False
    if not (isfinite(p_end.x) and isfinite(p_end.y)):
        return False

    valid = _check_straight(p_start, p_end, p_mid, threshold, inside)

    if valid and not geom_fully_inside:
        # Build a 2-point LineString via set_coordinates on a cached
        # geometry array
        _seg_coords_buf[0, 0] = p_start.x
        _seg_coords_buf[0, 1] = p_start.y
        _seg_coords_buf[1, 0] = p_end.x
        _seg_coords_buf[1, 1] = p_end.y
        g_segment = _shapely_lib.set_coordinates(_seg_geom_arr, _seg_coords_buf)
        if inside:
            valid = _shapely_lib.covers(gp_domain.context, g_segment)
        else:
            valid = _shapely_lib.disjoint(gp_domain.context, g_segment)

    return valid


# recursive subdivision limit; matching d3-geo's default
cdef int MAX_DEPTH = 16


cdef void _find_crossing(
    double t_a, const Point &p_a, State state_a,
    double t_b, const Point &p_b,
    Interpolator interpolator,
    object gp_domain,
    double *t_in, Point *p_in,
    double *t_out, Point *p_out,
    bool geom_fully_inside=False
) except *:
    """
    Binary-search for the domain boundary between p_a (state_a) and p_b
    (a different state). On return, t_in/p_in hold the last point still in
    state_a and t_out/p_out hold the first point on the other side.
    """
    cdef double t_lo, t_hi, t_mid
    cdef Point p_lo, p_hi, p_mid
    t_lo = t_a
    p_lo = p_a
    t_hi = t_b
    p_hi = p_b
    while t_hi - t_lo > 1e-6:
        t_mid = (t_lo + t_hi) * 0.5
        p_mid = interpolator.interpolate(t_mid)
        if get_state(p_mid, gp_domain, geom_fully_inside) == state_a:
            t_lo = t_mid
            p_lo = p_mid
        else:
            t_hi = t_mid
            p_hi = p_mid
    t_in[0] = t_lo
    p_in[0] = p_lo
    t_out[0] = t_hi
    p_out[0] = p_hi


cdef void _find_valid_boundary(
    double t_a, const Point &p_a,
    double t_lo, const Point &p_lo,
    double t_hi, const Point &p_hi,
    Interpolator interpolator,
    object gp_domain,
    double threshold,
    double *t_valid, Point *p_valid,
    double *t_invalid, Point *p_invalid,
    bool geom_fully_inside=False
) except *:
    """
    Binary-search for the last t in [t_lo, t_hi] where the projected
    segment from (t_a, p_a) to (t, interpolate(t)) satisfies
    straightAndDomain with inside=True. Used to locate the wrap-around
    boundary (e.g. antimeridian) when both segment endpoints are inside the
    domain but the geodesic midpoint falls outside the [0, 1] along-range.
    On return, t_valid/p_valid is the last valid endpoint and
    t_invalid/p_invalid is the first invalid one.
    """
    cdef double t_lo_, t_hi_, t_test
    cdef Point p_lo_, p_hi_, p_test
    t_lo_ = t_lo
    p_lo_ = p_lo
    t_hi_ = t_hi
    p_hi_ = p_hi
    cdef double t_mid_test
    cdef Point p_mid_test
    while t_hi_ - t_lo_ > 1e-6:
        t_test = (t_lo_ + t_hi_) * 0.5
        p_test = interpolator.interpolate(t_test)
        # Use only the pure-geometry curvature check (no shapely): we already
        # know both endpoints are inside the domain.
        t_mid_test = (t_a + t_test) * 0.5
        p_mid_test = interpolator.interpolate(t_mid_test)
        if _check_straight(p_a, p_test, p_mid_test, threshold, True):
            t_lo_ = t_test
            p_lo_ = p_test
        else:
            t_hi_ = t_test
            p_hi_ = p_test
    t_valid[0] = t_lo_
    p_valid[0] = p_lo_
    t_invalid[0] = t_hi_
    p_invalid[0] = p_hi_


cdef void _resample_recursive(
    double t0, const Point &p0, const State &state0,
    double t1, const Point &p1,
    Interpolator interpolator,
    object gp_domain,
    double threshold,
    LineAccumulator lines,
    int depth,
    bool geom_fully_inside=False,
    double precomp_pmid_x=HUGE_VAL,
    double precomp_pmid_y=HUGE_VAL,
) except *:
    """
    Recursively resample a projected line segment using adaptive midpoint
    subdivision. The algorithm is symmetric: swapping (t0, p0) and (t1, p1)
    produces identical sample points because the midpoint is always computed
    at t = (t0 + t1) / 2 in source-coordinate space.

    On acceptance the endpoint p1 is emitted based on state transitions.
    The caller is responsible for ensuring that, for a POINT_IN segment,
    p0 has already been added to the accumulator (or will be added via
    add_point_if_empty in the acceptance case).

    precomp_pmid_x/y: Pre-computed projected midpoint at t = (t0+t1)/2,
        supplied by project_linear when batch-projecting first-level midpoints
        avoids repeated scalar pyproj calls. HUGE_VAL sentinel means
        "not supplied; compute via interpolator".
    """
    cdef bool inside = (state0 == POINT_IN)
    cdef bool ok
    cdef double t_mid, t_in, t_out
    cdef double seg_dx_w, seg_dy_w, seg_hq_w, mid_dx_w, mid_dy_w, along_w
    cdef Point p_mid, p_in, p_out
    cdef State state1, state_mid, state_out

    # Always compute t_mid (needed for recursive sub-calls even when p_mid is
    # precomputed at the top level).
    t_mid = (t0 + t1) * 0.5
    # Use the pre-computed midpoint if provided; otherwise project via the
    # interpolator. Precomputation is only valid at the top level of a segment
    # (t0=0, t1=1 from project_linear); all recursive sub-calls pass HUGE_VAL.
    if isfinite(precomp_pmid_x):
        p_mid.x = precomp_pmid_x
        p_mid.y = precomp_pmid_y
    else:
        p_mid = interpolator.interpolate(t_mid)

    # Check whether the direct segment is geometrically acceptable and
    # domain-consistent.
    ok = straightAndDomain(p0, p1, p_mid, threshold,
                           gp_domain, inside, geom_fully_inside)

    if ok:
        # Segment is smooth and domain-consistent: emit it directly.
        if state0 == POINT_IN:
            lines.add_point_if_empty(p0)
            lines.add_point(p1)
        else:
            # POINT_OUT or POINT_NAN: only emit p1 if we are entering the domain
            state1 = get_state(p1, gp_domain, geom_fully_inside)
            if state1 == POINT_IN:
                lines.new_line()
                lines.add_point(p1)
        return

    if depth == 0:
        # Reached maximum subdivision depth without finding a smooth segment.
        # This means there is a genuine discontinuity (e.g. antimeridian
        # crossing, domain boundary) that cannot be resolved further.
        # Close the current line and start a new one
        state1 = get_state(p1, gp_domain, geom_fully_inside)
        if state0 == POINT_IN:
            lines.add_point_if_empty(p0)
            # Start a new segment for whatever is on the other side.
            if state1 == POINT_IN:
                lines.new_line()
                lines.add_point(p1)
        else:
            if state1 == POINT_IN:
                lines.new_line()
                lines.add_point(p1)
        return

    # p_mid was already computed before the straightAndDomain check above.
    state_mid = get_state(p_mid, gp_domain, geom_fully_inside)

    if state_mid != state0:
        # Boundary crossing in the left half [t0, t_mid].
        # Locate it precisely, then recurse on clean sub-intervals so that
        # curvature resampling is not confused by the domain boundary.
        _find_crossing(t0, p0, state0, t_mid, p_mid,
                       interpolator, gp_domain,
                       &t_in, &p_in, &t_out, &p_out, geom_fully_inside)
        if t_in > t0 + 1e-6:
            # valid t_in segment, recurse on it
            _resample_recursive(t0, p0, state0, t_in, p_in,
                                interpolator, gp_domain, threshold,
                                lines, depth - 1, geom_fully_inside)
        elif state0 == POINT_IN:
            # p0 is itself on the boundary: just ensure it is recorded.
            lines.add_point_if_empty(p0)
        state_out = get_state(p_out, gp_domain, geom_fully_inside)
        if state0 != POINT_IN and state_out == POINT_IN:
            lines.new_line()
        _resample_recursive(t_out, p_out, state_out, t1, p1,
                            interpolator, gp_domain, threshold,
                            lines, depth - 1, geom_fully_inside)
    else:
        state1 = get_state(p1, gp_domain, geom_fully_inside)
        if state1 != state0:
            # Boundary crossing in the right half [t_mid, t1].
            _find_crossing(t_mid, p_mid, state0, t1, p1,
                           interpolator, gp_domain,
                           &t_in, &p_in, &t_out, &p_out, geom_fully_inside)
            _resample_recursive(t0, p0, state0, t_in, p_in,
                                interpolator, gp_domain, threshold,
                                lines, depth - 1, geom_fully_inside)
            state_out = get_state(p_out, gp_domain, geom_fully_inside)
            if state0 != POINT_IN and state_out == POINT_IN:
                lines.new_line()
            _resample_recursive(t_out, p_out, state_out, t1, p1,
                                interpolator, gp_domain, threshold,
                                lines, depth - 1, geom_fully_inside)
        else:
            # Same state throughout.
            # Check for a wrap-around crossing: the geodesic midpoint lies
            # outside the [0, 1] along-range of the projected segment. This
            # signals that the geodesic crosses the domain boundary (e.g.
            # the antimeridian in PlateCarree) without a detectable state
            # transition, requiring a line split.
            if state0 == POINT_IN:
                seg_dx_w = p1.x - p0.x
                seg_dy_w = p1.y - p0.y
                seg_hq_w = seg_dx_w*seg_dx_w + seg_dy_w*seg_dy_w
                if seg_hq_w > 0.0:
                    mid_dx_w = p_mid.x - p0.x
                    mid_dy_w = p_mid.y - p0.y
                    along_w = (seg_dx_w*mid_dx_w + seg_dy_w*mid_dy_w) / seg_hq_w
                    if along_w < 0.0 or along_w > 1.0:
                        # Wrap detected: binary-search for the last valid
                        # point before the domain boundary, split, then
                        # continue from just past the boundary.
                        # Search the full [t0, t1] interval.
                        _find_valid_boundary(t0, p0, t0, p0, t1, p1,
                                             interpolator, gp_domain, threshold,
                                             &t_in, &p_in, &t_out, &p_out,
                                             geom_fully_inside)
                        if t_in > t0 + 1e-6:
                            _resample_recursive(t0, p0, state0, t_in, p_in,
                                                interpolator, gp_domain, threshold,
                                                lines, depth - 1, geom_fully_inside)
                        else:
                            # Start is right at the boundary: record p0 then split.
                            lines.add_point_if_empty(p0)
                        lines.new_line()
                        _resample_recursive(t_out, p_out, state0, t1, p1,
                                            interpolator, gp_domain, threshold,
                                            lines, depth - 1, geom_fully_inside)
                        return

            # No wrap detected: normal adaptive subdivision for curvature.
            # Only recurse for visible (POINT_IN) segments. Non-visible
            # segments (POINT_NAN, POINT_OUT) have no points to emit, so
            # curvature refinement is wasteful and would cause O(2^depth)
            # recursive calls per segment, slowing things down.
            if state0 != POINT_IN:
                return
            _resample_recursive(t0, p0, state0, t_mid, p_mid,
                                interpolator, gp_domain, threshold,
                                lines, depth - 1, geom_fully_inside)
            _resample_recursive(t_mid, p_mid, state_mid, t1, p1,
                                interpolator, gp_domain, threshold,
                                lines, depth - 1, geom_fully_inside)


cdef void _project_segment(double[:] src_from, double[:] src_to,
                           double[:] dest_from, double[:] dest_to,
                           Interpolator interpolator,
                           object gp_domain,
                           double threshold, LineAccumulator lines,
                           bool geom_fully_inside=False,
                           double precomp_pmid_x=HUGE_VAL,
                           double precomp_pmid_y=HUGE_VAL) except *:
    cdef Point p_src_from, p_src_to
    cdef Point p0, p1
    cdef State state0

    # Set up the interpolator with the source (un-projected) endpoints.
    p_src_from.x = src_from[0]
    p_src_from.y = src_from[1]
    p_src_to.x = src_to[0]
    p_src_to.y = src_to[1]
    interpolator.set_line(p_src_from, p_src_to)

    # Work in destination (projected) coordinates from here on.
    p0.x = dest_from[0]
    p0.y = dest_from[1]
    p1.x = dest_to[0]
    p1.y = dest_to[1]

    state0 = get_state(p0, gp_domain, geom_fully_inside)

    _resample_recursive(0.0, p0, state0, 1.0, p1,
                        interpolator, gp_domain, threshold,
                        lines, MAX_DEPTH, geom_fully_inside,
                        precomp_pmid_x, precomp_pmid_y)


@lru_cache(maxsize=4)
def _interpolator(src_crs, dest_projection):
    # Get an Interpolator from the given CRS and projection.
    # Callers must hold a reference to these systems for the lifetime
    # of the interpolator. If they get garbage-collected while interpolator
    # exists you *will* segfault.

    cdef Interpolator interpolator
    if src_crs.is_geodetic():
        interpolator = SphericalInterpolator()
    else:
        interpolator = CartesianInterpolator()
    interpolator.init(src_crs, dest_projection)
    return interpolator


def project_linear(geometry not None, src_crs not None,
                   dest_projection not None, bint is_ring=False):
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
    is_ring : bool, optional
        Set to ``True`` when *geometry* is a closed ring.  Controls the
        return type: ``True`` returns an ordered list of
        `~shapely.geometry.LineString` fragments (ring-traversal order);
        ``False`` returns a `~shapely.geometry.MultiLineString`.
        Defaults to ``False``.

    Returns
    -------
    `shapely.geometry.MultiLineString` or list of `shapely.geometry.LineString`
        When *is_ring* is ``False``, returns a MultiLineString of projected
        segments.

        When *is_ring* is ``True``, returns an ordered list of LineString
        fragments in ring-traversal order.  Fragment ``i`` ends on the
        projection boundary and fragment ``(i+1) % N`` starts on the same
        boundary; the caller is responsible for inserting the closing
        boundary arc between them.  If the ring projects entirely inside the
        domain with no cuts, a single-element list containing a closed
        LineString is returned.

    """
    cdef:
        double threshold = dest_projection.threshold
        Interpolator interpolator
        object g_domain
        double[:, :] src_coords, dest_coords, mid_dest_coords
        unsigned int src_size, src_idx
        object gp_domain
        LineAccumulator lines
        bool has_mid_precomp

    g_domain = dest_projection.domain

    interpolator = _interpolator(src_crs, dest_projection)

    src_coords = np.asarray(geometry.coords)
    dest_coords = interpolator.project_points(src_coords)
    # Use the cached prepared domain from the projection – avoids rebuilding the
    # prepared geometry (sprep.prep()) on every ring projection call.
    gp_domain = dest_projection._prepared_domain

    src_size = len(src_coords)  # check exceptions

    # Test the entire geometry to see if there are any domain crossings.
    # If there are none, we can skip expensive per-segment shapely checks.
    cdef bool geom_fully_inside = False
    cdef bool all_finite = True
    cdef unsigned int f_idx
    for f_idx in range(src_size):
        if not (isfinite(dest_coords[f_idx, 0]) and isfinite(dest_coords[f_idx, 1])):
            all_finite = False
            break
    if all_finite:
        dest_arr = np.asarray(dest_coords)
        dest_line = shapely.linestrings(dest_arr)
        geom_fully_inside = gp_domain.covers(dest_line)

    # Batch-project all first-level midpoints in one vectorised call.
    mid_dest_np = interpolator.batch_midpoints(src_coords)
    has_mid_precomp = mid_dest_np is not None
    if has_mid_precomp:
        mid_dest_coords = mid_dest_np

    lines = LineAccumulator()
    for src_idx in range(1, src_size):
        if has_mid_precomp:
            _project_segment(src_coords[src_idx - 1, :2], src_coords[src_idx, :2],
                             dest_coords[src_idx - 1, :2], dest_coords[src_idx, :2],
                             interpolator, gp_domain, threshold, lines,
                             geom_fully_inside,
                             mid_dest_coords[src_idx - 1, 0],
                             mid_dest_coords[src_idx - 1, 1])
        else:
            _project_segment(src_coords[src_idx - 1, :2], src_coords[src_idx, :2],
                             dest_coords[src_idx - 1, :2], dest_coords[src_idx, :2],
                             interpolator, gp_domain, threshold, lines,
                             geom_fully_inside);

    if is_ring:
        result = lines.as_ring_fragments()
    else:
        result = lines.as_geom()

    del lines, interpolator
    return result
