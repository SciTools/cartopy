# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
#
# cython: embedsignature=True

cdef extern from "geodesic.h":
    # External imports of Proj4.9 functions
    cdef struct geod_geodesic:
        pass
    cdef struct geod_geodesicline:
        pass

    void geod_init(geod_geodesic*, double, double) nogil
    void geod_direct(geod_geodesic*, double, double, double, double,
                     double*, double*, double*) nogil
    void geod_inverse(geod_geodesic*, double, double, double, double,
                      double*, double*, double*) nogil
    double geod_geninverse(geod_geodesic*, double, double, double, double,
                           double*, double*, double*, double*, double*,
                           double*, double*) nogil
    void geod_lineinit(geod_geodesicline*, geod_geodesic*, double, double,
                       double, int) nogil
    void geod_genposition(geod_geodesicline*, int, double, double*,
                          double*, double*, double*, double*, double*,
                          double*, double*) nogil

    cdef int GEOD_ARCMODE
    cdef int GEOD_LATITUDE
    cdef int GEOD_LONGITUDE
