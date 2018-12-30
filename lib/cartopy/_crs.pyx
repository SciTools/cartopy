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
This module defines the core CRS class which can interface with Proj.
The CRS class is the base-class for all projections defined in :mod:`cartopy.crs`.

"""
from collections import OrderedDict
import sys
import warnings

import numpy as np
import six

cimport numpy as np

cdef double NAN = float('nan')


PROJ4_VERSION = (PROJ_VERSION_MAJOR, PROJ_VERSION_MINOR, PROJ_VERSION_PATCH)

WGS84_SEMIMAJOR_AXIS = 6378137.0
WGS84_SEMIMINOR_AXIS = 6356752.3142


def pystrdecode(cstr):
    """
    Decode a string to a python string.
    """
    if sys.version_info[0] > 2:
        try:
            return cstr.decode("utf-8")
        except AttributeError:
            pass
    return cstr


class ProjError(Exception):
    """
    Raised when there has been an exception calling PROJ.

    Add a ``status`` attribute to the exception which has the
    PROJ error reference.

    """
    internal_proj_error = None
    def __init__(self, int status=0):
        internal_proj_error = "unknown"
        cdef char * errno_str = NULL
        if self.internal_proj_error is not None:
            internal_proj_error = self.internal_proj_error
            ProjError.clear()
        elif status:
            errno_str = proj_errno_string(status)
            if errno_str is not NULL:
                internal_proj_error = pystrdecode(errno_str)
        error_message = "Internal PROJ error: {internal_proj_error})".format(
            internal_proj_error=internal_proj_error,
        )
        super(ProjError, self).__init__(error_message)

    @staticmethod
    def clear():
        """
        This will clear the internal PROJ erro message.
        """
        ProjError.internal_proj_error = None


cdef void cartopy_log_function(void *user_data, int level, const char *error_msg):
    """
    Log function for catching PROJ errors.
    """
    if level == PJ_LOG_ERROR:
        ProjError.internal_proj_error = pystrdecode(error_msg)


proj_log_func(NULL, NULL, cartopy_log_function)


class Globe(object):
    """
    Define an ellipsoid and, optionally, how to relate it to the real world.

    """
    def __init__(self, datum=None, ellipse='WGS84',
                 semimajor_axis=None, semiminor_axis=None,
                 flattening=None, inverse_flattening=None,
                 towgs84=None, nadgrids=None):
        """
        Parameters
        ----------
        datum
            Proj "datum" definition. Defaults to None.
        ellipse
            Proj "ellps" definition. Defaults to 'WGS84'.
        semimajor_axis
            Semimajor axis of the spheroid / ellipsoid.  Defaults to None.
        semiminor_axis
            Semiminor axis of the ellipsoid.  Defaults to None.
        flattening
            Flattening of the ellipsoid.  Defaults to None.
        inverse_flattening
            Inverse flattening of the ellipsoid.  Defaults to None.
        towgs84
            Passed through to the Proj definition.  Defaults to None.
        nadgrids
            Passed through to the Proj definition.  Defaults to None.

        """
        self.datum = datum
        self.ellipse = ellipse
        self.semimajor_axis = semimajor_axis
        self.semiminor_axis = semiminor_axis
        self.flattening = flattening
        self.inverse_flattening = inverse_flattening
        self.towgs84 = towgs84
        self.nadgrids = nadgrids

    def to_proj4_params(self):
        """
        Create an OrderedDict of key value pairs which represents this globe
        in terms of proj params.

        """
        proj4_params = (['datum', self.datum], ['ellps', self.ellipse],
                        ['a', self.semimajor_axis], ['b', self.semiminor_axis],
                        ['f', self.flattening], ['rf', self.inverse_flattening],
                        ['towgs84', self.towgs84], ['nadgrids', self.nadgrids])
        return OrderedDict((k, v) for k, v in proj4_params if v is not None)

cdef class TransProj:
    def __cinit__(self):
        self.projpj = NULL

    def __init__(self, CRS p1, CRS p2):
        ProjError.clear()
        self.projpj = proj_create_crs_to_crs(
            NULL,
            six.b(p1.proj4_init),
            six.b(p2.proj4_init),
            NULL)
        if self.projpj == NULL:
            raise ProjError()

    def __dealloc__(self):
        """destroy projection definition"""
        if self.projpj != NULL:
            proj_destroy(self.projpj)


cdef class CRS:
    """
    Define a Coordinate Reference System using proj.

    """

    #: Whether this projection can handle ellipses.
    _handles_ellipses = True

    def __cinit__(self):
        self.proj4 = NULL

    def __dealloc__(self):
        if self.proj4 != NULL:
            proj_destroy(self.proj4)

    def __init__(self, proj4_params, globe=None):
        """
        Parameters
        ----------
        proj4_params: iterable of key-value pairs
            The proj4 parameters required to define the
            desired CRS.  The parameters should not describe
            the desired elliptic model, instead create an
            appropriate Globe instance. The ``proj4_params``
            parameters will override any parameters that the
            Globe defines.
        globe: :class:`~cartopy.crs.Globe` instance, optional
            If omitted, the default Globe instance will be created.
            See :class:`~cartopy.crs.Globe` for details.

        """
        if globe is None:
            if self._handles_ellipses:
                globe = Globe()
            else:
                globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS,
                              ellipse=None)
        if not self._handles_ellipses:
            a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS
            b = globe.semiminor_axis or a
            if a != b or globe.ellipse is not None:
                warnings.warn('The "{}" projection does not handle elliptical '
                              'globes.'.format(self.__class__.__name__))

        self.globe = globe
        self.proj4_params = self.globe.to_proj4_params()
        self.proj4_params.update(proj4_params)

        first_item_inserted = False
        init_items = []
        for k, v in self.proj4_params.items():
            if v is not None:
                if isinstance(v, float):
                    init_items.append('+{}={:.16}'.format(k, v))
                elif isinstance(v, np.float32):
                    init_items.append('+{}={:.8}'.format(k, v))
                elif not first_item_inserted and k in ("init", "proj"):
                    init_items.insert(0, '+{}={}'.format(k, v))
                    initial_item_inserted = True
                else:
                    init_items.append('+{}={}'.format(k, v))
            else:
                init_items.append('+{}'.format(k))
        self.proj4_init = ' '.join(init_items) + ' +no_defs'
        proj4_init_bytes = six.b(self.proj4_init)
        if self.proj4 != NULL:
            proj_destroy(self.proj4)
        self.proj4 = proj_create(NULL, proj4_init_bytes)
        if self.proj4 == NULL:
            raise ProjError()
        self.proj4_type = proj_get_type(self.proj4)

    # Cython uses this method instead of the normal rich comparisons.
    def __richcmp__(self, other, op):
        # We're only interested in:
        #   == -> 2
        #   != -> 3
        result = NotImplemented
        if isinstance(other, CRS):
            if op == 2:
                result = self.proj4_init == other.proj4_init
            elif op == 3:
                result = self.proj4_init != other.proj4_init
        return result

    def __hash__(self):
        """Hash the CRS based on its proj4_init string."""
        return hash(self.proj4_init)

    def __reduce__(self):
        """
        Implement the __reduce__ API so that unpickling produces a stateless
        instance of this class (e.g. an empty tuple). The state will then be
        added via __getstate__ and __setstate__.

        """
        return self.__class__, tuple()

    def __getstate__(self):
        """Return the full state of this instance for reconstruction
        in ``__setstate__``.

        """
        return {'proj4_params': self.proj4_params}

    def __setstate__(self, state):
        """
        Take the dictionary created by ``__getstate__`` and passes it
        through to the class's __init__ method.

        """
        self.__init__(self, **state)

    # TODO
    #def _geod(self): # to return the pyproj.Geod

    def _as_mpl_transform(self, axes=None):
        """
        Cast this CRS instance into a :class:`matplotlib.axes.Axes` using
        the Matplotlib ``_as_mpl_transform`` interface.

        """
        # lazy import mpl.geoaxes (and therefore matplotlib) as mpl
        # is only an optional dependency
        import cartopy.mpl.geoaxes as geoaxes
        if not isinstance(axes, geoaxes.GeoAxes):
            raise ValueError('Axes should be an instance of GeoAxes, got %s' % type(axes))
        return geoaxes.InterProjectionTransform(self, axes.projection) + axes.transData

    property proj4_params:
        def __get__(self):
            return dict(self.proj4_params)

    def as_geocentric(self):
        """
        Return a new Geocentric CRS with the same ellipse/datum as this
        CRS.

        """
        return Geocentric(self.globe)

    def as_geodetic(self):
        """
        Return a new Geodetic CRS with the same ellipse/datum as this
        CRS.

        """
        return Geodetic(self.globe)

    def is_geodetic(self):
        """
        Returns
        -------
        bool: True if projection in geodetic coordinates
        """
        return self.proj4_type is PJ_TYPE_GEODETIC_CRS

    def transform_point(self, double x, double y, CRS src_crs not None, trap=True):
        """
        transform_point(x, y, src_crs)

        Transform the given float64 coordinate pair, in the given source
        coordinate system (``src_crs``), to this coordinate system.

        Parameters
        ----------
        x
            the x coordinate, in ``src_crs`` coordinates, to transform
        y
            the y coordinate, in ``src_crs`` coordinates, to transform
        src_crs
            instance of :class:`CRS` that represents the coordinate
            system of ``x`` and ``y``.
        trap
            Whether proj errors for "latitude or longitude exceeded limits" and
            "tolerance condition error" should be trapped.

        Returns
        -------
        (x, y) in this coordinate system

        """
        pj_trans = TransProj(src_crs, self)

        cdef:
            double cx, cy
        cx = x
        cy = y
        if src_crs.is_geodetic():
            cx = proj_torad(cx)
            cy = proj_torad(cy)

        proj_trans_generic(
            pj_trans.projpj,
            PJ_FWD,
            &cx, sizeof(double), 1,
            &cy, sizeof(double), 1,
            NULL, 0, 0,
            NULL, 0, 0,
        )
        cdef int status = proj_errno(pj_trans.projpj)
        if trap and status == -14 or status == -20:
            # -14 => "latitude or longitude exceeded limits"
            # -20 => "tolerance condition error"
            cx = cy = NAN
        elif trap and status != 0:
            raise ProjError(status)

        if self.is_geodetic():
            cx = proj_todeg(cx)
            cy = proj_todeg(cy)
        return (cx, cy)

    def transform_points(self, CRS src_crs not None,
                                np.ndarray x not None,
                                np.ndarray y not None,
                                np.ndarray z=None):
        """
        transform_points(src_crs, x, y[, z])

        Transform the given coordinates, in the given source
        coordinate system (``src_crs``), to this coordinate system.

        Parameters
        ----------
        src_crs
            instance of :class:`CRS` that represents the
            coordinate system of ``x``, ``y`` and ``z``.
        x
            the x coordinates (array), in ``src_crs`` coordinates,
            to transform.  May be 1 or 2 dimensional.
        y
            the y coordinates (array), in ``src_crs`` coordinates,
            to transform.  Its shape must match that of x.
        z: optional
            the z coordinates (array), in ``src_crs`` coordinates, to
            transform.  Defaults to None.
            If supplied, its shape must match that of x.

        Returns
        -------
            Array of shape ``x.shape + (3, )`` in this coordinate system.

        """
        pj_trans = TransProj(src_crs, self)

        cdef np.ndarray[np.double_t, ndim=2] result

        result_shape = tuple(x.shape[i] for i in range(x.ndim)) + (3, )

        if z is None:
            if x.ndim > 2 or y.ndim > 2:
                raise ValueError('x and y arrays must be 1 or 2 dimensional')
            elif x.ndim != 1 or y.ndim != 1:
                x, y = x.flatten(), y.flatten()

            if x.shape[0] != y.shape[0]:
                raise ValueError('x and y arrays must have the same length')
        else:
            if x.ndim > 2 or y.ndim > 2 or z.ndim > 2:
                raise ValueError('x, y and z arrays must be 1 or 2 '
                                 'dimensional')
            elif x.ndim != 1 or y.ndim != 1 or z.ndim != 1:
                x, y, z = x.flatten(), y.flatten(), z.flatten()

            if not x.shape[0] == y.shape[0] == z.shape[0]:
                raise ValueError('x, y, and z arrays must have the same '
                                 'length')

        npts = x.shape[0]

        result = np.empty([npts, 3], dtype=np.double)
        if src_crs.is_geodetic():
            result[:, 0] = np.deg2rad(x)
            result[:, 1] = np.deg2rad(y)
        else:
            result[:, 0] = x
            result[:, 1] = y
        # if a z has been given, put it in the result array which will be
        # transformed in-place
        if z is None:
            result[:, 2] = 0
        else:
            result[:, 2] = z

        # call proj. The result array is modified in place. This is only
        # safe if npts is not 0.
        cdef int status = 0
        if npts:
            ProjError.clear()
            stride = result.strides[0]
            proj_trans_generic(
                pj_trans.projpj,
                PJ_FWD,
                &result[0, 0], stride, npts,
                &result[0, 1], stride, npts,
                &result[0, 2], stride, npts,
                NULL, 0, 0,
            )

            status = proj_errno(pj_trans.projpj)
            if status:
                raise ProjError(status)
            elif ProjError.internal_proj_error is not None:
                raise ProjError()

        if self.is_geodetic():
            result[:, :2] = np.rad2deg(result[:, :2])

        if len(result_shape) > 2:
            return result.reshape(result_shape)

        return result

    def transform_vectors(self, src_proj, x, y, u, v):
        """
        transform_vectors(src_proj, x, y, u, v)

        Transform the given vector components, with coordinates in the
        given source coordinate system (``src_proj``), to this coordinate
        system. The vector components must be given relative to the
        source projection's coordinate reference system (grid eastward and
        grid northward).

        Parameters
        ----------
        src_proj
            The :class:`CRS.Projection` that represents the coordinate system
            the vectors are defined in.
        x
            The x coordinates of the vectors in the source projection.
        y
            The y coordinates of the vectors in the source projection.
        u
            The grid-eastward components of the vectors.
        v
            The grid-northward components of the vectors.

        Note
        ----
            x, y, u and v may be 1 or 2 dimensional, but must all have matching
            shapes.

        Returns
        -------
            ut, vt: The transformed vector components.

        Note
        ----
           The algorithm used to transform vectors is an approximation
           rather than an exact transform, but the accuracy should be
           good enough for visualization purposes.

        """
        if not (x.shape == y.shape == u.shape == v.shape):
            raise ValueError('x, y, u and v arrays must be the same shape')
        if x.ndim not in (1, 2):
            raise ValueError('x, y, u and v must be 1 or 2 dimensional')
        # Transform the coordinates to the target projection.
        proj_xyz = self.transform_points(src_proj, x, y)
        target_x, target_y = proj_xyz[..., 0], proj_xyz[..., 1]
        # Rotate the input vectors to the projection.
        #
        # 1: Find the magnitude and direction of the input vectors.
        vector_magnitudes = (u**2 + v**2)**0.5
        vector_angles = np.arctan2(v, u)
        # 2: Find a point in the direction of the original vector that is
        #    a small distance away from the base point of the vector (near
        #    the poles the point may have to be in the opposite direction
        #    to be valid).
        factor = 360000.
        delta = (src_proj.x_limits[1] - src_proj.x_limits[0]) / factor
        x_perturbations = delta * np.cos(vector_angles)
        y_perturbations = delta * np.sin(vector_angles)
        # 3: Handle points that are invalid. These come from picking a new
        #    point that is outside the domain of the CRS. The first step is
        #    to apply the native transform to the input coordinates to make
        #    sure they are in the appropriate range. Then detect all the
        #    coordinates where the perturbation takes the point out of the
        #    valid x-domain and fix them. After that do the same for points
        #    that are outside the valid y-domain, which may reintroduce some
        #    points outside of the valid x-domain
        proj_xyz = src_proj.transform_points(src_proj, x, y)
        source_x, source_y = proj_xyz[..., 0], proj_xyz[..., 1]
        #    Detect all the coordinates where the perturbation takes the point
        #    outside of the valid x-domain, and reverse the direction of the
        #    perturbation to fix this.
        eps = 1e-9
        invalid_x = np.logical_or(
            source_x + x_perturbations < src_proj.x_limits[0]-eps,
            source_x + x_perturbations > src_proj.x_limits[1]+eps)
        if invalid_x.any():
            x_perturbations[invalid_x] *= -1
            y_perturbations[invalid_x] *= -1
        #    Do the same for coordinates where the perturbation takes the point
        #    outside of the valid y-domain. This may reintroduce some points
        #    that will be outside the x-domain when the perturbation is
        #    applied.
        invalid_y = np.logical_or(
            source_y + y_perturbations < src_proj.y_limits[0]-eps,
            source_y + y_perturbations > src_proj.y_limits[1]+eps)
        if invalid_y.any():
            x_perturbations[invalid_y] *= -1
            y_perturbations[invalid_y] *= -1
        #    Keep track of the points where the perturbation direction was
        #    reversed.
        reversed_vectors = np.logical_xor(invalid_x, invalid_y)
        #    See if there were any points where we cannot reverse the direction
        #    of the perturbation to get the perturbed point within the valid
        #    domain of the projection, and issue a warning if there are.
        problem_points = np.logical_or(
            source_x + x_perturbations < src_proj.x_limits[0]-eps,
            source_x + x_perturbations > src_proj.x_limits[1]+eps)
        if problem_points.any():
            warnings.warn('Some vectors at source domain corners '
                          'may not have been transformed correctly')
        # 4: Transform this set of points to the projection coordinates and
        #    find the angle between the base point and the perturbed point
        #    in the projection coordinates (reversing the direction at any
        #    points where the original was reversed in step 3).
        proj_xyz = self.transform_points(src_proj,
                                         source_x + x_perturbations,
                                         source_y + y_perturbations)
        target_x_perturbed = proj_xyz[..., 0]
        target_y_perturbed = proj_xyz[..., 1]
        projected_angles = np.arctan2(target_y_perturbed - target_y,
                                      target_x_perturbed - target_x)
        if reversed_vectors.any():
            projected_angles[reversed_vectors] += np.pi
        # 5: Form the projected vector components, preserving the magnitude
        #    of the original vectors.
        projected_u = vector_magnitudes * np.cos(projected_angles)
        projected_v = vector_magnitudes * np.sin(projected_angles)
        return projected_u, projected_v


class Geodetic(CRS):
    """
    Define a latitude/longitude coordinate system with spherical topology,
    geographical distance and coordinates are measured in degrees.

    """
    def __init__(self, globe=None):
        """
        Parameters
        ----------
        globe: A :class:`cartopy.crs.Globe`, optional
            Defaults to a "WGS84" datum.

        """
        proj4_params = [('proj', 'lonlat')]
        globe = globe or Globe(datum='WGS84')
        super(Geodetic, self).__init__(proj4_params, globe)

    # XXX Implement fwd such as Basemap's Geod. Would be used in the tissot example.
    # Could come from http://geographiclib.sourceforge.net


class Geocentric(CRS):
    """
    Define a Geocentric coordinate system, where x, y, z are Cartesian
    coordinates from the center of the Earth.

    """
    def __init__(self, globe=None):
        """
        Parameters
        ----------
        globe: A :class:`cartopy.crs.Globe`, optional
            Defaults to a "WGS84" datum.

        """
        proj4_params = [('proj', 'geocent')]
        globe = globe or Globe(datum='WGS84')
        super(Geocentric, self).__init__(proj4_params, globe)
