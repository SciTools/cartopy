# (C) British Crown Copyright 2018, Met Office
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
