# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
#
# cython: embedsignature=True

"""
This module defines the Geodesic class which can interface with the Proj
geodesic functions.

"""
from cpython.mem cimport PyMem_Malloc
from cpython.mem cimport PyMem_Free
import numpy as np
cimport numpy as np
cimport cython
from cython.parallel cimport prange
import shapely.geometry as sgeom

from ._geodesic cimport geod_geodesic, geod_init, geod_direct, geod_inverse


cdef class Geodesic:
    """
    Define an ellipsoid on which to solve geodesic problems.

    """
    cdef geod_geodesic* geod
    cdef double radius
    cdef double flattening

    def __init__(self, radius=6378137.0, flattening=1/298.257223563):
        """
        Parameters
        ----------
        radius: float, optional
            Equatorial radius (metres). Defaults to the WGS84 semimajor axis
            (6378137.0 metres).

        flattening: float, optional
            Flattening of ellipsoid.  Setting flattening = 0 gives a sphere.
            Negative flattening gives a prolate ellipsoid. If flattening > 1,
            set flattening to 1/flattening.
            Defaults to the WGS84 flattening (1/298.257223563).

        """
        # This method exists solely for docstrings. The __cinit__ method is
        # where the real work happens. Make sure to keep the signatures synced.

    def __cinit__(self, radius=6378137.0, flattening=1/298.257223563):
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
        points: array_like, shape=(n *or* 1, 2)
            The starting longitude-latitude point(s) from which to travel.

        azimuths: float or array_like with shape=(n, )
            List of azimuth values (degrees).

        distances : float or array_like with shape(n, )
            List of distances values (metres).

        Returns
        -------
        `numpy.ndarray`, shape=(n, 3)
            The longitudes, latitudes, and forward azimuths of the located
            endpoint(s).

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
        Solve the inverse geodesic problem.

        Can accept and broadcast length 1 arguments. For example, given a
        single start point, an array of different endpoints can be supplied to
        find multiple distances.

        Parameters
        ----------
        points: array_like, shape=(n *or* 1, 2)
            The starting longitude-latitude point(s) from which to travel.

        endpoints: array_like, shape=(n *or* 1, 2)
            The longitude-latitude point(s) to travel to.

        Returns
        -------
        `numpy.ndarray`, shape=(n, 3)
            The distances, and the (forward) azimuths of the start and end
            points.

        """

        cdef int n_points, i
        cdef double[:, :] pts, epts, orig_pts

        # Create numpy arrays from inputs, and ensure correct shape.
        points = np.array(points, dtype=np.float64)
        endpoints = np.array(endpoints, dtype=np.float64)

        if points.ndim > 2 or (points.ndim == 2 and points.shape[1] != 2):
            raise ValueError(
                'Expecting input points to be (N, 2), got {}'
                ''.format(points.shape))

        pts = points.reshape((-1, 2))
        epts = endpoints.reshape((-1, 2))

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
        Find a geodesic circle of given radius at a given point.

        Parameters
        ----------
        lon : float
            Longitude coordinate of the centre.
        lat : float
            Latitude coordinate of the centre.
        radius : float
            The radius of the circle (metres).
        n_samples: int, optional
            Integer number of sample points of circle.
        endpoint: bool, optional
            Whether to repeat endpoint at the end of returned array.

        Returns
        -------
        `numpy.ndarray`, shape=(n_samples, 2)
            The evenly spaced longitude-latitude points on the circle.

        """

        cdef int i

        # Put the input arguments into c-typed values.
        cdef double[:, :] center = np.array([lon, lat]).reshape((1, 2))
        cdef double[:] radius_m = np.asarray(radius).reshape(1)

        azimuths = np.linspace(360., 0., n_samples,
                               endpoint=endpoint).astype(np.double)

        return self.direct(center, azimuths, radius_m)[:, 0:2]

    def geometry_length(self, geometry):
        """
        Return the distance (in physical meters) of the given Shapely geometry.

        The geometry is assumed to be in spherical (lon, lat) coordinates.

        Parameters
        ----------
        geometry : `shapely.geometry.BaseGeometry`
            The Shapely geometry to compute the length of. For polygons, the
            exterior length will be calculated. For multi-part geometries, the
            sum of the parts will be computed.

        """
        result = None
        if hasattr(geometry, 'geoms'):
            # Multi-geometry.
            result = sum(self.geometry_length(geom) for geom in geometry.geoms)

        elif hasattr(geometry, 'exterior'):
            # Polygon.
            result = self.geometry_length(geometry.exterior)

        elif hasattr(geometry, 'coords') and not isinstance(geometry, sgeom.Point):
            coords = np.array(geometry.coords)
            result = self.geometry_length(coords)

        elif isinstance(geometry, np.ndarray):
            coords = geometry
            distances, _, _ = np.array(
                self.inverse(coords[:-1, :], coords[1:, :]).T)
            result = distances.sum()

        else:
            raise TypeError('Unhandled type {}'.format(geometry.__class__))

        return result
