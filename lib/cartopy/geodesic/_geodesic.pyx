# (C) British Crown Copyright 2015 - 2018, Met Office
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
This module defines the Geodesic class which can interface with the Proj.4.
geodesic functions.

"""
from cpython.mem cimport PyMem_Malloc
from cpython.mem cimport PyMem_Free
import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport prange

cdef extern from "geodesic.h":
    # External imports of Proj4.9 functions
    cdef struct geod_geodesic:
        pass

    ctypedef geod_geodesic* geodesic_t

    void geod_init(geodesic_t, double, double)
    void geod_direct(geodesic_t, double, double, double, double,
                     double*, double*, double*) nogil
    void geod_inverse(geodesic_t, double, double, double, double,
                      double*, double*, double*) nogil

cdef class Geodesic:
    """
    Defines an ellipsoid on which to solve geodesic problems.

    """
    cdef geod_geodesic* geod
    cdef double radius
    cdef double flattening

    def __cinit__(self, radius=6378137.0, flattening=1/298.257223563):
        """
        Creates an ellipsoid with a given radius and flattening.

        Parameters
        ----------
        radius: optional
            Equatorial radius (metres). Defaults to the WGS84 semimajor axis
            (6378137.0 metres).

        flattening: optional
            Flattening of ellipsoid.  Setting flattening = 0 gives a sphere.
            Negative flattening gives a prolate ellipsoid. If flattening > 1,
            set flattening to 1/flattening.
            Defaults to the WGS84 flattening (1/298.257223563).

        """
        # allocate some memory (filled with random data)
        self.geod = <geod_geodesic*> PyMem_Malloc(sizeof(geod_geodesic))
        if not self.geod:
            raise MemoryError()
        geod_init(self.geod, radius, flattening)
        self.radius = radius
        self.flattening = flattening

    def __dealloc__(self):
        # Free allocated memory.
        PyMem_Free(self.geod)

    def __str__(self):
        fmt = self.radius, 1/self.flattening
        return '<Geodesic: radius=%0.3f, flattening=1/%0.3f>' % fmt

    @cython.boundscheck(False)
    def direct(self, points, azimuths, distances):
        """
        Solve the direct geodesic problem where the length of the geodesic is
        specified in terms of distance.

        Can accept and broadcast length 1 arguments. For example, given a
        single start point and distance, an array of different azimuths can be
        supplied to locate multiple endpoints.

        Parameters
        ----------
        points
            An n (or 1) by 2 numpy.ndarray, list or tuple of lon-lat
            points.  The starting point(s) from which to travel.

        azimuths
            A length n (or 1) numpy.ndarray or list of azimuth
            values (degrees).

        distances
            A length n (or 1) numpy.ndarray or list of distances
            values (metres).

        Returns
        -------
        An n by 3 np.ndarray of lons, lats, and forward azimuths of the
        located endpoint(s).

        """

        cdef int n_points, i
        cdef double[:, :] pts, orig_pts
        cdef double[:] azims, dists

        # Create numpy arrays from inputs, and ensure correct shape. Note:
        # reshape(-1) returns a 1D array from a 0 dimensional array as required
        # for broadcasting.
        pts = np.array(points, dtype=np.float64).reshape((-1, 2))
        azims = np.array(azimuths, dtype=np.float64).reshape(-1)
        dists = np.array(distances, dtype=np.float64).reshape(-1)

        sizes = [pts.shape[0], azims.size, dists.size]
        n_points = max(sizes)
        if not all(size in [1, n_points] for size in sizes):
            raise ValueError("Inputs must have common length n or length one.")

        # Broadcast any length 1 arrays to the correct size.
        if pts.shape[0] == 1:
            orig_pts = pts
            pts = np.empty([n_points, 2], dtype=np.float64)
            pts[:, :] = orig_pts

        if azims.size == 1:
            azims = np.repeat(azims, n_points)

        if dists.size == 1:
            dists = np.repeat(dists, n_points)

        cdef double[:, :] return_pts = np.empty((n_points, 3),
                                                dtype=np.float64)

        with nogil:
            for i in prange(n_points):
                geod_direct(self.geod,
                            pts[i, 1], pts[i, 0], azims[i], dists[i],
                            &return_pts[i, 1], &return_pts[i, 0],
                            &return_pts[i, 2])

        return return_pts

    @cython.boundscheck(False)
    def inverse(self, points, endpoints):
        """
        Solves the inverse geodesic problem.

        Can accept and broadcast length 1 arguments. For example, given a
        single start point, an array of different endpoints can be supplied to
        find multiple distances.

        Parameters
        ----------
        points
            An n (or 1) by 2 numpy.ndarray, list or tuple of lon-lat points.
            The starting point(s) from which to travel.

        endpoints:
            An n (or 1) by 2 numpy.ndarray, list or tuple of lon-lat points.
            The point(s) to travel to.

        Returns
        -------
        An n by 3 np.ndarray of distances, and the (forward) azimuths of
        the start and end points.

        """

        cdef int n_points, i
        cdef double[:, :] pts, epts, orig_pts

        # Create numpy arrays from inputs, and ensure correct shape. Note:
        # reshape(-1) returns a 1D array from a 0 dimensional array as required
        # for broadcasting.
        pts = np.array(points, dtype=np.float64).reshape((-1, 2))
        epts = np.array(endpoints, dtype=np.float64).reshape((-1, 2))

        sizes = [pts.shape[0], epts.shape[0]]
        n_points = max(sizes)
        if not all(size in [1, n_points] for size in sizes):
            raise ValueError("Inputs must have common length n or length one.")

        # Broadcast any length 1 arrays to the correct size.
        if pts.shape[0] == 1:
            orig_pts = pts
            pts = np.empty([n_points, 2], dtype=np.float64)
            pts[:, :] = orig_pts

        if epts.shape[0] == 1:
            orig_pts = epts
            epts = np.empty([n_points, 2], dtype=np.float64)
            epts[:, :] = orig_pts

        cdef double[:, :] results = np.empty((n_points, 3), dtype=np.float64)

        with nogil:
            for i in prange(n_points):
                geod_inverse(self.geod, pts[i, 1], pts[i, 0], epts[i, 1],
                             epts[i, 0], &results[i, 0], &results[i, 1],
                             &results[i, 2])

        return results

    def circle(self, double lon, double lat, double radius, int n_samples=180,
               endpoint=False):
        """
        Finds a geodesic circle of given radius at a given point.

        Parameters
        ----------
        lon
            Longitude coordinate of the centre.
        lat
            Latitude coordinate of the centre.
        radius
            The radius of the circle (metres).
        n_samples: optional
            Integer number of sample points of circle.
        endpoint: optional
            Boolean for whether to repeat endpoint at the end of returned
            array.

        Returns
        -------
        An n_samples by 2 np.ndarray of evenly spaced lon-lat points on the
        circle.

        """

        cdef int i

        # Put the input arguments into c-typed values.
        cdef double[:, :] center = np.array([lon, lat]).reshape((1, 2))
        cdef double[:] radius_m = np.asarray(radius).reshape(1)

        azimuths = np.linspace(360., 0., n_samples,
                               endpoint=endpoint).astype(np.double)

        return self.direct(center, azimuths, radius_m)[:, 0:2]
