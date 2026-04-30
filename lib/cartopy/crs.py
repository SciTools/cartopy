# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

"""
The crs module defines Coordinate Reference Systems and the transformations
between them.

"""

from abc import ABCMeta
import bisect
from collections import Counter, OrderedDict
from functools import lru_cache
import io
import json
import math
import warnings

import numpy as np
import pyproj
from pyproj import Transformer
from pyproj.exceptions import ProjError
import shapely
import shapely.geometry as sgeom
from shapely.prepared import prep

import cartopy.trace


try:
    # https://github.com/pyproj4/pyproj/pull/912
    from pyproj.crs import CustomConstructorCRS as _CRS
except ImportError:
    from pyproj import CRS as _CRS

__document_these__ = ['CRS', 'Geocentric', 'Geodetic', 'Globe']

WGS84_SEMIMAJOR_AXIS = 6378137.0
WGS84_SEMIMINOR_AXIS = 6356752.3142

# cache to avoid recreating in every ring call
_EMPTY_MLS = sgeom.MultiLineString([])


# Cache the transformer creation method
@lru_cache
def _get_transformer_from_crs(src_crs, tgt_crs):
    return Transformer.from_crs(src_crs, tgt_crs, always_xy=True)


def _safe_pj_transform(src_crs, tgt_crs, x, y, z=None, trap=True):
    transformer = _get_transformer_from_crs(src_crs, tgt_crs)

    # if a projection is essentially 2d there
    # should be no harm in setting its z to 0
    if z is None:
        z = np.zeros_like(x)

    with warnings.catch_warnings():
        # pyproj implicitly converts size-1 arrays to scalars, which is
        # deprecated in numpy 1.25, but *also* handles the future error
        # see https://github.com/numpy/numpy/pull/10615
        # and https://github.com/SciTools/cartopy/pull/2194
        warnings.filterwarnings(
            "ignore",
            message="Conversion of an array with ndim > 0"
        )
        return transformer.transform(x, y, z, errcheck=trap)


class Globe:
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
        proj4_params = (
            ['datum', self.datum],
            ['ellps', self.ellipse],
            ['a', self.semimajor_axis],
            ['b', self.semiminor_axis],
            ['f', self.flattening],
            ['rf', self.inverse_flattening],
            ['towgs84', self.towgs84],
            ['nadgrids', self.nadgrids]
        )
        return OrderedDict((k, v) for k, v in proj4_params if v is not None)


class CRS(_CRS):
    """
    Define a Coordinate Reference System using proj. The :class:`cartopy.crs.CRS`
    class is the very core of cartopy, all coordinate reference systems in cartopy
    have :class:`~cartopy.crs.CRS` as a parent class.
    """

    #: Whether this projection can handle ellipses.
    _handles_ellipses = True

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
        self.input = (proj4_params, globe)

        # for compatibility with pyproj.CRS and rasterio.crs.CRS
        try:
            proj4_params = proj4_params.to_wkt()
        except AttributeError:
            pass
        # handle PROJ JSON
        if (
            isinstance(proj4_params, dict) and
            "proj" not in proj4_params and
            "init" not in proj4_params
        ):
            proj4_params = json.dumps(proj4_params)

        if globe is not None and isinstance(proj4_params, str):
            raise ValueError("Cannot have 'globe' with string params.")
        if globe is None and not isinstance(proj4_params, str):
            if self._handles_ellipses:
                globe = Globe()
            else:
                globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS,
                              ellipse=None)
        if globe is not None and not self._handles_ellipses:
            a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS
            b = globe.semiminor_axis or a
            if a != b or globe.ellipse is not None:
                warnings.warn(f'The {self.__class__.__name__!r} projection '
                              'does not handle elliptical globes.')
        self.globe = globe
        if isinstance(proj4_params, str):
            self._proj4_params = {}
            self.proj4_init = proj4_params
        else:
            self._proj4_params = self.globe.to_proj4_params()
            self._proj4_params.update(proj4_params)

            init_items = []
            for k, v in self._proj4_params.items():
                if v is not None:
                    if isinstance(v, float):
                        init_items.append(f'+{k}={v:.16}')
                    elif isinstance(v, np.float32):
                        init_items.append(f'+{k}={v:.8}')
                    else:
                        init_items.append(f'+{k}={v}')
                else:
                    init_items.append(f'+{k}')
            self.proj4_init = ' '.join(init_items) + ' +no_defs'
        super().__init__(self.proj4_init)

    def __eq__(self, other):
        if isinstance(other, CRS) and self.proj4_init == other.proj4_init:
            # Fast path Cartopy's CRS
            return True
        # For everything else, we let pyproj handle the comparison
        return super().__eq__(other)

    def __hash__(self):
        """Hash the CRS based on its proj4_init string."""
        return hash(self.proj4_init)

    def __reduce__(self):
        """
        Implement the __reduce__ method used when pickling or performing deepcopy.
        """
        if type(self) is CRS:
            # State can be reproduced by the proj4_params and globe inputs.
            return self.__class__, self.input
        else:
            # Produces a stateless instance of this class (e.g. an empty tuple).
            # The state will then be added via __getstate__ and __setstate__.
            # We are forced to this approach because a CRS does not store
            # the constructor keyword arguments in its state.
            return self.__class__, (), self.__getstate__()

    def __getstate__(self):
        """Return the full state of this instance for reconstruction
        in ``__setstate__``.
        """
        state = self.__dict__.copy()
        # remove pyproj specific attrs
        state.pop('srs', None)
        state.pop('_local', None)
        # Remove the proj4 instance and the proj4_init string, which can
        # be re-created (in __setstate__) from the other arguments.
        state.pop('proj4', None)
        state.pop('proj4_init', None)
        state['proj4_params'] = self.proj4_params
        return state

    def __setstate__(self, state):
        """
        Take the dictionary created by ``__getstate__`` and passes it
        through to this implementation's __init__ method.
        """
        # Strip out the key state items for a CRS instance.
        CRS_state = {key: state.pop(key) for key in ['proj4_params', 'globe']}
        # Put everything else directly into the dict of the instance.
        self.__dict__.update(state)
        # Call the init of this class to ensure that the projection is
        # properly initialised with proj4.
        CRS.__init__(self, **CRS_state)

    def _as_mpl_transform(self, axes=None):
        """
        Cast this CRS instance into a :class:`matplotlib.axes.Axes` using
        the Matplotlib ``_as_mpl_transform`` interface.

        """
        # lazy import mpl.geoaxes (and therefore matplotlib) as mpl
        # is only an optional dependency
        import cartopy.mpl.geoaxes as geoaxes
        if not isinstance(axes, geoaxes.GeoAxes):
            raise ValueError(
                f'Axes should be an instance of GeoAxes, got {type(axes)}'
            )
        return (
            geoaxes.InterProjectionTransform(self, axes.projection) +
            axes.transData
        )

    @property
    def proj4_params(self):
        return dict(self._proj4_params)

    def as_geocentric(self):
        """
        Return a new Geocentric CRS with the same ellipse/datum as this
        CRS.

        """
        return CRS(
            {
                "$schema": (
                    "https://proj.org/schemas/v0.2/projjson.schema.json"
                ),
                "type": "GeodeticCRS",
                "name": "unknown",
                "datum": self.datum.to_json_dict(),
                "coordinate_system": {
                    "subtype": "Cartesian",
                    "axis": [
                        {
                            "name": "Geocentric X",
                            "abbreviation": "X",
                            "direction": "geocentricX",
                            "unit": "metre"
                        },
                        {
                            "name": "Geocentric Y",
                            "abbreviation": "Y",
                            "direction": "geocentricY",
                            "unit": "metre"
                        },
                        {
                            "name": "Geocentric Z",
                            "abbreviation": "Z",
                            "direction": "geocentricZ",
                            "unit": "metre"
                        }
                    ]
                }
            }
        )

    def as_geodetic(self):
        """
        Return a new Geodetic CRS with the same ellipse/datum as this
        CRS.

        """
        return CRS(self.geodetic_crs.srs)

    def is_geodetic(self):
        return self.is_geographic

    def transform_point(self, x, y, src_crs, trap=True):
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
        result = self.transform_points(
            src_crs, np.array([x]), np.array([y]), trap=trap,
        ).reshape((1, 3))
        return result[0, 0], result[0, 1]

    def transform_points(self, src_crs, x, y, z=None, trap=False):
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
        trap
            Whether proj errors for "latitude or longitude exceeded limits" and
            "tolerance condition error" should be trapped.

        Returns
        -------
            Array of shape ``x.shape + (3, )`` in this coordinate system.

        """
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
        if npts:
            if self == src_crs and (
                    isinstance(src_crs, _CylindricalProjection) or
                    self.is_geodetic()):
                # convert from [0,360] to [-180,180]
                x = np.array(x, copy=True)
                to_180 = (x > 180) | (x < -180)
                x[to_180] = (((x[to_180] + 180) % 360) - 180)
            try:
                result[:, 0], result[:, 1], result[:, 2] = \
                    _safe_pj_transform(src_crs, self, x, y, z, trap=trap)
            except ProjError as err:
                msg = str(err).lower()
                if (
                    "latitude" in msg or
                    "longitude" in msg or
                    "outside of projection domain" in msg or
                    "tolerance condition error" in msg
                ):
                    result[:] = np.nan
                else:
                    raise

            if not trap:
                result[np.isinf(result)] = np.nan

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
        vector_magnitudes = np.hypot(u, v)
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
            source_x + x_perturbations < src_proj.x_limits[0] - eps,
            source_x + x_perturbations > src_proj.x_limits[1] + eps)
        if invalid_x.any():
            x_perturbations[invalid_x] *= -1
            y_perturbations[invalid_x] *= -1
        #    Do the same for coordinates where the perturbation takes the point
        #    outside of the valid y-domain. This may reintroduce some points
        #    that will be outside the x-domain when the perturbation is
        #    applied.
        invalid_y = np.logical_or(
            source_y + y_perturbations < src_proj.y_limits[0] - eps,
            source_y + y_perturbations > src_proj.y_limits[1] + eps)
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
            source_x + x_perturbations < src_proj.x_limits[0] - eps,
            source_x + x_perturbations > src_proj.x_limits[1] + eps)
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
        super().__init__(proj4_params, globe)

    # XXX Implement fwd such as Basemap's Geod.
    # Would be used in the tissot example.
    # Could come from https://geographiclib.sourceforge.io


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
        super().__init__(proj4_params, globe)


class RotatedGeodetic(CRS):
    """
    Define a rotated latitude/longitude coordinate system with spherical
    topology and geographical distance.

    Coordinates are measured in degrees.

    The class uses proj to perform an ob_tran operation, using the
    pole_longitude to set a lon_0 then performing two rotations based on
    pole_latitude and central_rotated_longitude.
    This is equivalent to setting the new pole to a location defined by
    the pole_latitude and pole_longitude values in the GeogCRS defined by
    globe, then rotating this new CRS about it's pole using the
    central_rotated_longitude value.

    """

    def __init__(self, pole_longitude, pole_latitude,
                 central_rotated_longitude=0.0, globe=None):
        """
        Parameters
        ----------
        pole_longitude
            Pole longitude position, in unrotated degrees.
        pole_latitude
            Pole latitude position, in unrotated degrees.
        central_rotated_longitude: optional
            Longitude rotation about the new pole, in degrees.  Defaults to 0.
        globe: optional
            A :class:`cartopy.crs.Globe`.  Defaults to a "WGS84" datum.

        """
        globe = globe or Globe(datum='WGS84')
        proj4_params = [('proj', 'ob_tran'), ('o_proj', 'latlon'),
                        ('o_lon_p', central_rotated_longitude),
                        ('o_lat_p', pole_latitude),
                        ('lon_0', 180 + pole_longitude),
                        ('to_meter', math.radians(1) * (
                            globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS))]

        super().__init__(proj4_params, globe=globe)


class Projection(CRS, metaclass=ABCMeta):
    """
    Define a projected coordinate system with flat topology and Euclidean
    distance.

    """

    _method_map = {
        'Point': '_project_point',
        'LineString': '_project_line_string',
        'LinearRing': '_project_linear_ring',
        'Polygon': '_project_polygon',
        'MultiPoint': '_project_multipoint',
        'MultiLineString': '_project_multiline',
        'MultiPolygon': '_project_multipolygon',
        'GeometryCollection': '_project_geometry_collection'
    }
    # Whether or not this projection can handle wrapped coordinates
    _wrappable = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.bounds = None
        if self.area_of_use:
            # Convert lat/lon bounds to projected bounds.
            # Geographic area of the entire dataset referenced to WGS 84
            # NB. We can't use a polygon transform at this stage because
            # that relies on the existence of the map boundary... the very
            # thing we're trying to work out! ;-)
            x0 = self.area_of_use.west
            x1 = self.area_of_use.east
            y0 = self.area_of_use.south
            y1 = self.area_of_use.north
            lons = np.array([x0, x0, x1, x1])
            lats = np.array([y0, y1, y1, y0])
            points = self.transform_points(
                PlateCarree().as_geodetic(), lons, lats
            )
            x = points[:, 0]
            y = points[:, 1]
            self.bounds = (x.min(), x.max(), y.min(), y.max())
            x0, x1, y0, y1 = self.bounds
            self.threshold = min(x1 - x0, y1 - y0) / 100.
        elif self.is_geographic:
            # If the projection is geographic without an area of use, assume
            # the bounds are the full globe.
            self.bounds = (-180, 180, -90, 90)

    @property
    def boundary(self):
        if self.bounds is None:
            raise NotImplementedError
        x0, x1, y0, y1 = self.bounds
        return sgeom.LineString([(x0, y0), (x0, y1), (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        if self.bounds is None:
            raise NotImplementedError
        x0, x1, y0, y1 = self.bounds
        return (x0, x1)

    @property
    def y_limits(self):
        if self.bounds is None:
            raise NotImplementedError
        x0, x1, y0, y1 = self.bounds
        return (y0, y1)

    @property
    def threshold(self):
        return getattr(self, '_threshold', 0.5)

    @threshold.setter
    def threshold(self, t):
        self._threshold = t

    @property
    def cw_boundary(self):
        try:
            boundary = self._cw_boundary
        except AttributeError:
            boundary = sgeom.LinearRing(self.boundary)
            self._cw_boundary = boundary
        return boundary

    @property
    def ccw_boundary(self):
        try:
            boundary = self._ccw_boundary
        except AttributeError:
            boundary = sgeom.LinearRing(self.boundary.coords[::-1])
            self._ccw_boundary = boundary
        return boundary

    @property
    def domain(self):
        try:
            domain = self._domain
        except AttributeError:
            domain = self._domain = sgeom.Polygon(self.boundary)
        return domain

    @property
    def _ring_threshold(self):
        """Coordinate-match tolerance for ring stitching (1e-5 of domain scale).

        Cached on first access; recomputed only if the projection bounds change.
        """
        try:
            return self.__ring_threshold_cache
        except AttributeError:
            val = max(abs(v) for v in self.x_limits + self.y_limits) * 1e-5
            self.__ring_threshold_cache = val
            return val

    @property
    def _prepared_domain(self):
        """Prepared (indexed) version of :attr:`domain`, cached per projection instance.

        Used by :func:`cartopy.trace.project_linear` for fast point-in-polygon
        tests.  The domain polygon is immutable, so a single prepared instance
        is safe to reuse for the projection's lifetime.
        """
        try:
            return self.__prepared_domain_cache
        except AttributeError:
            self.__prepared_domain_cache = prep(self.domain)
            return self.__prepared_domain_cache

    def is_geodetic(self):
        return False

    def _determine_longitude_bounds(self, central_longitude):
        # In new proj, using exact limits will wrap-around, so subtract a
        # small epsilon:
        epsilon = 1e-10
        minlon = -180 + central_longitude
        maxlon = 180 + central_longitude
        if central_longitude > 0:
            maxlon -= epsilon
        elif central_longitude < 0:
            minlon += epsilon
        return minlon, maxlon

    def _repr_html_(self):
        from html import escape
        try:
            # As matplotlib is not a core cartopy dependency, don't error
            # if it's not available.
            import matplotlib.pyplot as plt
        except ImportError:
            # We can't return an SVG of the CRS, so let Jupyter fall back to
            # a default repr by returning None.
            return None

        # Produce a visual repr of the Projection instance.
        fig, ax = plt.subplots(figsize=(5, 3),
                               subplot_kw={'projection': self})
        ax.set_global()
        ax.coastlines('auto')
        ax.gridlines()
        buf = io.StringIO()
        fig.savefig(buf, format='svg', bbox_inches='tight')
        plt.close(fig)
        # "Rewind" the buffer to the start and return it as an svg string.
        buf.seek(0)
        svg = buf.read()
        return f'{svg}<pre>{escape(object.__repr__(self))}</pre>'

    def _as_mpl_axes(self):
        import cartopy.mpl.geoaxes as geoaxes
        return geoaxes.GeoAxes, {'projection': self}

    def project_geometry(self, geometry, src_crs=None):
        """
        Project the given geometry into this projection.

        Parameters
        ----------
        geometry
            The geometry to (re-)project.
        src_crs: optional
            The source CRS.  Defaults to None.

            If src_crs is None, the source CRS is assumed to be a geodetic
            version of the target CRS.

        Returns
        -------
        geometry
            The projected result (a shapely geometry).

        """
        if src_crs is None:
            src_crs = self.as_geodetic()
        elif not isinstance(src_crs, CRS):
            raise TypeError('Source CRS must be an instance of CRS'
                            ' or one of its subclasses, or None.')
        geom_type = geometry.geom_type
        method_name = self._method_map.get(geom_type)
        if not method_name:
            raise ValueError(f'Unsupported geometry type {geom_type!r}')
        return getattr(self, method_name)(geometry, src_crs)

    def _project_point(self, point, src_crs):
        if point.is_empty:
            return point
        return sgeom.Point(*self.transform_point(point.x, point.y, src_crs))

    def _project_line_string(self, geometry, src_crs):
        return cartopy.trace.project_linear(geometry, src_crs, self, is_ring=False)

    def _project_linear_ring(self, linear_ring, src_crs, is_ccw=None):
        """
        Project the given LinearRing from src_crs into this CRS.

        When called from :meth:`_project_polygon`, ``is_ccw`` is the winding
        of the source polygon exterior and the method returns a
        ``(list[LinearRing], MultiLineString)`` tuple.

        When called via ``project_geometry(linear_ring)`` (i.e. without
        ``is_ccw``), the winding is inferred from the ring itself and the
        method returns a ``GeometryCollection([*rings, MultiLineString])``.
        """
        if is_ccw is None:
            # infer winding and produce GeometryCollection.
            is_ccw = True if src_crs.is_geodetic() else linear_ring.is_ccw
            rings, mls = self._project_linear_ring(linear_ring, src_crs, is_ccw)
            return sgeom.GeometryCollection([*rings, mls])

        # project_linear returns an ordered list of LineString fragments when
        # given a LinearRing. Fragment i exits the domain boundary and
        # fragment (i+1) % N enters it.
        fragments = cartopy.trace.project_linear(
            linear_ring, src_crs, self, is_ring=True
        )

        if not fragments:
            return [], _EMPTY_MLS

        threshold = self._ring_threshold

        if len(fragments) == 1:
            ls = fragments[0]
            if len(ls.coords) < 2:
                return [], _EMPTY_MLS
            c0, c1 = ls.coords[0], ls.coords[-1]
            if abs(c0[0] - c1[0]) <= threshold and abs(c0[1] - c1[1]) <= threshold:
                # Ring projected cleanly inside domain; no boundary stitching.
                return [sgeom.LinearRing(ls.coords[:-1])], _EMPTY_MLS
            # Single open fragment — needs boundary stitching.
            return [], sgeom.MultiLineString([ls])

        # When the ring's first source point was inside the domain,
        # fragments[0].start == fragments[-1].end (both map to the same interior
        # point). Merge last + first[1:] so all fragments have both endpoints
        # on the boundary.
        frags = list(fragments)
        c0, cm1 = frags[0].coords[0], frags[-1].coords[-1]
        if abs(c0[0] - cm1[0]) <= threshold and abs(c0[1] - cm1[1]) <= threshold:
            merged = sgeom.LineString(
                list(frags[-1].coords) + list(frags[0].coords)[1:]
            )
            frags = [merged] + frags[1:-1]

        # Merge consecutive fragments that share an endpoint. This happens
        # when the ring crosses the boundary twice at nearly the same location.
        # Adjacent matches only: the ordered fragment list guarantees that
        # frag[i].end == frag[i+1].start is the only possible adjacent join.
        i = 0
        while i < len(frags) - 1:
            ce, cs = frags[i].coords[-1], frags[i + 1].coords[0]
            if abs(ce[0] - cs[0]) <= threshold and abs(ce[1] - cs[1]) <= threshold:
                stitched = sgeom.LineString(
                    list(frags[i].coords) + list(frags[i + 1].coords)[1:]
                )
                frags = frags[:i] + [stitched] + frags[i + 2:]
            else:
                i += 1

        # A fragment whose endpoints are nearly coincident is a complete ring
        # that projected entirely within the domain — it doesn't cross the
        # boundary and should not go to _attach_lines_to_boundary.
        rings = []
        open_frags = []
        for frag in frags:
            coords = frag.coords
            ca, cb = coords[0], coords[-1]
            if (len(coords) > 3
                    and abs(ca[0] - cb[0]) <= threshold
                    and abs(ca[1] - cb[1]) <= threshold):
                ring = sgeom.LinearRing(coords[:-1])
                # Discard near-zero area rings; they are boundary artefacts
                # produced when a source vertex lies just outside the domain.
                if sgeom.Polygon(ring).area >= threshold ** 2:
                    rings.append(ring)
            else:
                open_frags.append(frag)

        return rings, sgeom.MultiLineString(open_frags)

    def _project_multipoint(self, geometry, src_crs):
        geoms = []
        for geom in geometry.geoms:
            geoms.append(self._project_point(geom, src_crs))
        return sgeom.MultiPoint(geoms)

    def _project_multiline(self, geometry, src_crs):
        geoms = []
        for geom in geometry.geoms:
            r = self._project_line_string(geom, src_crs)
            if r:
                geoms.extend(r.geoms)
        return sgeom.MultiLineString(geoms)

    def _project_multipolygon(self, geometry, src_crs):
        geoms = []
        for geom in geometry.geoms:
            r = self._project_polygon(geom, src_crs)
            if r:
                geoms.extend(r.geoms)
        return sgeom.MultiPolygon(geoms)

    def _project_geometry_collection(self, geometry, src_crs):
        return sgeom.GeometryCollection(
            [self.project_geometry(geom, src_crs) for geom in geometry.geoms])

    def _project_polygon(self, polygon, src_crs):
        """
        Return the projected polygon(s) derived from the given polygon.

        """
        # Determine orientation of polygon.
        if src_crs.is_geodetic():
            is_ccw = True
        else:
            is_ccw = polygon.exterior.is_ccw
        # Project each ring (exterior + interiors).  Self-closed rings are
        # returned immediately; open fragments are pooled for stitching.
        rings = []
        multi_lines = []
        for src_ring in [polygon.exterior] + list(polygon.interiors):
            closed, mls = self._project_linear_ring(src_ring, src_crs, is_ccw)
            rings.extend(closed)
            if len(mls.geoms) > 0:
                multi_lines.append(mls)

        # Convert any open fragments to rings by attaching boundary arcs.
        if multi_lines:
            rings.extend(self._attach_lines_to_boundary(multi_lines, is_ccw))

        # Resolve all inside vs. outside rings and convert to MultiPolygon.
        return self._rings_to_multi_polygon(rings, is_ccw)

    def _attach_lines_to_boundary(self, multi_line_strings, is_ccw):
        """
        Return a list of LinearRings by closing projected line segments with
        boundary arcs.

        Each projected segment has both its start and end point lying on the
        projection boundary. The algorithm connects them into closed rings by
        inserting boundary arcs between consecutive endpoints.

        How the forward walk and shorter arc work together:
        Segments are sorted by their start-point distance along the boundary.
        To close a ring, we walk *forward* (in the boundary's parameterisation
        direction) from each segment's end to find the next segment whose start
        arrives first, or the ring's own start if it comes around sooner.
        This forward walk determines which segment is next and whether to
        close.

        However, the boundary-corner vertices inserted between two consecutive
        endpoints use the *shorter* of the two possible arcs (forward or
        backward). These decisions are independent: the forward walk gives
        correct ring topology, while the shorter-arc rule correctly handles the
        edge case where two crossings straddle a boundary corner on opposite
        sides of the perimeter seam.

        Parameters
        ----------
        multi_line_strings : list of MultiLineString
            One MultiLineString per source ring (exterior or interior).
        is_ccw : bool
            Expected winding of the source polygon exterior ring.
        """
        boundary = self.ccw_boundary if is_ccw else self.cw_boundary
        perimeter = boundary.length
        threshold = self.threshold

        # Pre-compute boundary corner (distance, coord) pairs once, sorted
        # for O(log N) bisect queries inside fwd_arc_corners().
        corner_pairs  = sorted(
            (boundary.project(sgeom.Point(c)), (c[0], c[1]))
            for c in boundary.coords[:-1]
        )
        corner_dists  = [d for d, _ in corner_pairs]
        corner_coords = [c for _, c in corner_pairs]

        def fwd_arc_corners(d_from, arc_len):
            """Boundary corners on the forward arc of *arc_len* from *d_from*."""
            if arc_len == 0:
                return []
            if d_from + arc_len <= perimeter:
                lo = bisect.bisect_right(corner_dists, d_from)
                hi = bisect.bisect_left(corner_dists, d_from + arc_len)
                return corner_coords[lo:hi]
            else:
                wrap_end = d_from + arc_len - perimeter
                lo = bisect.bisect_right(corner_dists, d_from)
                hi = bisect.bisect_left(corner_dists, wrap_end)
                return corner_coords[lo:] + corner_coords[:hi]

        def arc_corners(d_from, d_to):
            """Boundary corners on the shorter arc from *d_from* to *d_to*."""
            d_fwd = (d_to - d_from) % perimeter
            if d_fwd <= perimeter / 2:
                return fwd_arc_corners(d_from, d_fwd)
            else:
                return list(reversed(fwd_arc_corners(d_to, perimeter - d_fwd)))

        def endpoint_info(xy):
            """Return (boundary_distance, is_on_boundary) for a segment endpoint."""
            pt = sgeom.Point(xy)
            d = boundary.project(pt)
            bd_pt = boundary.interpolate(d)
            return d, math.hypot(xy[0] - bd_pt.x, xy[1] - bd_pt.y) < threshold

        # Flatten all source ring segments into one pool, tracking which
        # source MultiLineString (source polygon ring) each segment came from.
        line_strings = []
        mls_idx_of_seg = []  # keep mapping of indices back to original mls
        for mls_i, mls in enumerate(multi_line_strings):
            for ls in mls.geoms:
                line_strings.append(ls)
                mls_idx_of_seg.append(mls_i)
        if not line_strings:
            return []

        completed = []
        remaining_segs = []  # (d_start, d_end, seg_idx, mls_i) for boundary segs

        # Partition: segments with both endpoints on the boundary go through
        # the forward-walk; segments with interior endpoints are already closed
        # rings and go directly to completed.
        for k, ls in enumerate(line_strings):
            d_start, start_on = endpoint_info(ls.coords[0])
            d_end,   end_on   = endpoint_info(ls.coords[-1])
            if start_on and end_on:
                remaining_segs.append((d_start, d_end, k, mls_idx_of_seg[k]))
            else:
                completed.append(list(ls.coords))

        # Sort by d_start so the forward-walk can always find the next segment
        # with a single bisect query.
        remaining_segs.sort(key=lambda t: t[0])
        # Parallel list of just the d_start values for bisect lookups.
        remaining_segs_d = [t[0] for t in remaining_segs]

        # Drop coordinate-wrapping artifact segments: a single segment from
        # its source ring that crosses the perimeter seam (d_end < d_start)
        # and whose closed form has the wrong winding is a PROJ artefact from
        # a source ring that lies entirely outside our domain. Keeping it
        # would cause _rings_to_multi_polygon to invert it into a flood polygon.
        segs_per_mls = Counter(mls_i for _, _, _, mls_i in remaining_segs)
        remaining_segs = [
            (d_s, d_e, k, mls_i) for d_s, d_e, k, mls_i in remaining_segs
            if not (segs_per_mls[mls_i] == 1 and d_e < d_s
                    and _ring_is_ccw(
                        sgeom.LinearRing(
                            list(line_strings[k].coords) + arc_corners(d_e, d_s)
                        )
                    ) != is_ccw)
        ]
        remaining_segs_d = [t[0] for t in remaining_segs]

        def finish_ring(d_end_local, ring_coords_local, d_fwd_local):
            """Close ring_coords with the shorter arc back to d_ring_start."""
            corners = arc_corners(d_end_local, d_ring_start)
            ring_coords_local += corners
            if not corners and d_fwd_local > 0:
                # No corners on the closing arc: insert a boundary midpoint so
                # the ring orientation is unambiguous to _rings_to_multi_polygon.
                if d_fwd_local <= perimeter / 2:
                    d_mid = (d_end_local + d_fwd_local / 2) % perimeter
                else:
                    d_mid = (d_end_local - (perimeter - d_fwd_local) / 2) % perimeter
                mid_pt = boundary.interpolate(d_mid)
                ring_coords_local.append((mid_pt.x, mid_pt.y))
            completed.append(ring_coords_local)

        while remaining_segs:
            # Start a new ring from the segment with the smallest d_start.
            # This choice ensures each ring consumes exactly the segments
            # on its forward arc before reaching any other ring's start.
            d_ring_start, d_end, i, current_mls = remaining_segs.pop(0)
            remaining_segs_d.pop(0)
            ring_coords = list(line_strings[i].coords)
            # d_end tracks the boundary distance of ring_coords[-1] without
            # requiring an extra boundary.project() call each inner iteration.

            while True:
                d_fwd_close = (d_ring_start - d_end) % perimeter

                if not remaining_segs:
                    # No more segments; close this ring.
                    finish_ring(d_end, ring_coords, d_fwd_close)
                    break

                # Find the next segment start going forward from d_end.
                # bisect_right gives the index of the first remaining_segs_d
                # entry > d_end. The modulo handles wrap-around (when all
                # remaining starts are <= d_end we cycle back to the smallest).
                pos = bisect.bisect_right(remaining_segs_d, d_end) % len(remaining_segs)
                d_next, d_next_end, j, j_mls = remaining_segs[pos]

                # Compare forward distance to the ring's own start vs. d_next.
                # On a tie the ring closes without consuming the next segment;
                # that segment will become the start of the following ring.
                d_fwd_next = (d_next - d_end) % perimeter

                if d_fwd_close <= d_fwd_next:
                    # Ring start arrives first: close the ring.
                    finish_ring(d_end, ring_coords, d_fwd_close)
                    break
                else:
                    # Segment j comes before our ring start.  If j is from a
                    # different source ring and the gap to it is large, only
                    # connect if the ring can't otherwise close without a very
                    # long boundary arc (> 1/5 perimeter).  This prevents an
                    # artifact segment from a source ring outside the domain
                    # from stealing an adjacent exterior segment.
                    if j_mls != current_mls:
                        next_start_xy = line_strings[j].coords[0]
                        gap = math.hypot(
                            ring_coords[-1][0] - next_start_xy[0],
                            ring_coords[-1][1] - next_start_xy[1],
                        )
                        if gap > threshold and d_fwd_close < perimeter / 5:
                            break  # close cannot be long; discard, leave j
                    # Connect to segment j, then continue building the ring.
                    ring_coords += arc_corners(d_end, d_next)
                    ring_coords += list(line_strings[j].coords)
                    remaining_segs.pop(pos)
                    remaining_segs_d.pop(pos)
                    d_end = d_next_end
                    current_mls = j_mls  # track mls of most recently added seg

        # Build LinearRings, skipping degenerate coordinate lists.
        # The minimum useful input is 3 distinct points. A 3-coord list
        # where coords[0] == coords[-1] represents only 2 unique points and
        # would produce a degenerate ring that isn't valid.
        return [
            sgeom.LinearRing(coords)
            for coords in completed
            if len(coords) > 3 or (len(coords) == 3 and coords[0] != coords[-1])
        ]

    def _rings_to_multi_polygon(self, rings, is_ccw):
        exterior_rings = []
        interior_rings = []
        for ring in rings:
            if _ring_is_ccw(ring) != is_ccw:
                interior_rings.append(ring)
            else:
                exterior_rings.append(ring)

        # When ALL exterior rings are smaller than the largest interior ring,
        # the source polygon wraps around the projection boundary (e.g. a
        # global polygon in SouthPolarStereo). Reclassify those small exteriors
        # as interiors so the box-difference path below inverts them correctly.
        if interior_rings and exterior_rings:
            interior_areas = [sgeom.Polygon(r).area for r in interior_rings]
            max_interior_area = max(interior_areas)
            reclassified = [r for r in exterior_rings
                            if sgeom.Polygon(r).area < max_interior_area]
            if len(reclassified) == len(exterior_rings):
                for r in reclassified:
                    exterior_rings.remove(r)
                    interior_rings.append(r)

        polygon_bits = []

        # Turn all the exterior rings into polygon definitions,
        # "slurping up" any interior rings they contain.
        for exterior_ring in exterior_rings:
            polygon = sgeom.Polygon(exterior_ring)
            prep_polygon = prep(polygon)
            holes = []
            for interior_ring in interior_rings[:]:
                if prep_polygon.contains(interior_ring):
                    holes.append(interior_ring)
                    interior_rings.remove(interior_ring)
                elif polygon.crosses(interior_ring):
                    # Likely that we have an invalid geometry such as
                    # that from #509 or #537.
                    holes.append(interior_ring)
                    interior_rings.remove(interior_ring)
            polygon_bits.append((exterior_ring.coords,
                                 [ring.coords for ring in holes]))

        # Any left over "interior" rings need "inverting" with respect
        # to the boundary.
        if interior_rings:
            boundary_poly = self.domain
            x3, y3, x4, y4 = boundary_poly.bounds
            bx = (x4 - x3) * 0.1
            by = (y4 - y3) * 0.1
            x3 -= bx
            y3 -= by
            x4 += bx
            y4 += by

            interior_polys = []

            for ring in interior_rings:
                polygon = shapely.make_valid(sgeom.Polygon(ring))
                if not polygon.is_empty:
                    if isinstance(polygon, sgeom.Polygon):
                        interior_polys.append(polygon)
                    elif isinstance(polygon, sgeom.MultiPolygon):
                        interior_polys.extend(polygon.geoms)
                    elif isinstance(polygon, sgeom.GeometryCollection):
                        for geom in polygon.geoms:
                            if isinstance(geom, sgeom.Polygon):
                                interior_polys.append(geom)
                            elif isinstance(geom, sgeom.MultiPolygon):
                                interior_polys.extend(geom.geoms)
                    else:
                        # make_valid may produce some linestrings.  Ignore these
                        continue

                    x1, y1, x2, y2 = polygon.bounds
                    bx = (x2 - x1) * 0.1
                    by = (y2 - y1) * 0.1
                    x1 -= bx
                    y1 -= by
                    x2 += bx
                    y2 += by

                    x3 = min(x1, x3)
                    x4 = max(x2, x4)
                    y3 = min(y1, y3)
                    y4 = max(y2, y4)

            box = sgeom.box(x3, y3, x4, y4, ccw=is_ccw)

            if interior_polys:
                # Invert any valid interior polygons
                multi_poly = shapely.make_valid(sgeom.MultiPolygon(interior_polys))
                polygon = box.difference(multi_poly)

                # Intersect the inverted polygon with the boundary
                polygon = boundary_poly.intersection(polygon)

                if not polygon.is_empty:
                    if isinstance(polygon, sgeom.MultiPolygon):
                        polygon_bits.extend(polygon.geoms)
                    elif isinstance(polygon, sgeom.Polygon):
                        polygon_bits.append(polygon)
                    elif isinstance(polygon, sgeom.GeometryCollection):
                        for geom in polygon.geoms:
                            if isinstance(geom, sgeom.Polygon):
                                polygon_bits.append(geom)
                            elif isinstance(geom, sgeom.MultiPolygon):
                                polygon_bits.extend(geom.geoms)

        return sgeom.MultiPolygon(polygon_bits)

    def quick_vertices_transform(self, vertices, src_crs):
        """
        Where possible, return a vertices array transformed to this CRS from
        the given vertices array of shape ``(n, 2)`` and the source CRS.

        Note
        ----
            This method may return None to indicate that the vertices cannot
            be transformed quickly, and a more complex geometry transformation
            is required (see :meth:`cartopy.crs.Projection.project_geometry`).

        """
        if vertices.size == 0:
            return vertices

        if self == src_crs:
            x = vertices[:, 0]
            y = vertices[:, 1]
            # Extend the limits a tiny amount to allow for precision mistakes
            epsilon = 1.e-10
            x_limits = (self.x_limits[0] - epsilon, self.x_limits[1] + epsilon)
            y_limits = (self.y_limits[0] - epsilon, self.y_limits[1] + epsilon)
            if (x.min() >= x_limits[0] and x.max() <= x_limits[1] and
                    y.min() >= y_limits[0] and y.max() <= y_limits[1]):
                return vertices


class _RectangularProjection(Projection, metaclass=ABCMeta):
    """
    The abstract superclass of projections with a rectangular domain which
    is symmetric about the origin.

    """
    _wrappable = True

    def __init__(self, proj4_params, half_width, half_height, globe=None):
        self._half_width = half_width
        self._half_height = half_height
        super().__init__(proj4_params, globe=globe)

    @property
    def boundary(self):
        w, h = self._half_width, self._half_height
        return sgeom.LinearRing([(-w, -h), (-w, h), (w, h), (w, -h), (-w, -h)])

    @property
    def x_limits(self):
        return (-self._half_width, self._half_width)

    @property
    def y_limits(self):
        return (-self._half_height, self._half_height)


class _CylindricalProjection(_RectangularProjection, metaclass=ABCMeta):
    """
    The abstract class which denotes cylindrical projections where we
    want to allow x values to wrap around.

    """
    _wrappable = True


def _ellipse_boundary(semimajor=2, semiminor=1, easting=0, northing=0, n=201):
    """
    Define a projection boundary using an ellipse.

    This type of boundary is used by several projections.

    """

    t = np.linspace(0, -2 * np.pi, n)  # Clockwise boundary.
    coords = np.vstack([semimajor * np.cos(t), semiminor * np.sin(t)])
    coords += ([easting], [northing])
    return coords


class PlateCarree(_CylindricalProjection):
    def __init__(self, central_longitude=0.0, globe=None):
        globe = globe or Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS)
        proj4_params = [('proj', 'eqc'), ('lon_0', central_longitude),
                        ('to_meter', math.radians(1) * (
                            globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS)),
                        ('vto_meter', 1)]
        x_max = 180
        y_max = 90
        # Set the threshold around 0.5 if the x max is 180.
        self.threshold = x_max / 360
        super().__init__(proj4_params, x_max, y_max, globe=globe)

    def _bbox_and_offset(self, other_plate_carree):
        """
        Return a pair of (xmin, xmax) pairs and an offset which can be used
        for identification of whether data in ``other_plate_carree`` needs
        to be transformed to wrap appropriately.

        >>> import cartopy.crs as ccrs
        >>> src = ccrs.PlateCarree(central_longitude=10)
        >>> bboxes, offset = ccrs.PlateCarree()._bbox_and_offset(src)
        >>> print(bboxes)
        [[-180, -170.0], [-170.0, 180]]
        >>> print(offset)
        10.0

        The returned values are longitudes in ``other_plate_carree``'s
        coordinate system.

        Warning
        -------
            The two CRSs must be identical in every way, other than their
            central longitudes. No checking of this is done.

        """
        self_lon_0 = self.proj4_params['lon_0']
        other_lon_0 = other_plate_carree.proj4_params['lon_0']

        lon_0_offset = other_lon_0 - self_lon_0

        lon_lower_bound_0 = self.x_limits[0]
        lon_lower_bound_1 = (other_plate_carree.x_limits[0] + lon_0_offset)

        if lon_lower_bound_1 < self.x_limits[0]:
            lon_lower_bound_1 += np.diff(self.x_limits)[0]

        lon_lower_bound_0, lon_lower_bound_1 = sorted(
            [lon_lower_bound_0, lon_lower_bound_1])

        bbox = [[lon_lower_bound_0, lon_lower_bound_1],
                [lon_lower_bound_1, lon_lower_bound_0]]

        bbox[1][1] += self.x_limits[1] - self.x_limits[0]

        return bbox, lon_0_offset

    def quick_vertices_transform(self, vertices, src_crs):
        return_value = super().quick_vertices_transform(vertices, src_crs)

        # Optimise the PlateCarree -> PlateCarree case where no
        # wrapping or interpolation needs to take place.
        if return_value is None and isinstance(src_crs, PlateCarree):
            self_params = self.proj4_params.copy()
            src_params = src_crs.proj4_params.copy()
            self_params.pop('lon_0'), src_params.pop('lon_0')

            xs, ys = vertices[:, 0], vertices[:, 1]

            potential = (self_params == src_params and
                         self.y_limits[0] <= ys.min() and
                         self.y_limits[1] >= ys.max())
            if potential:
                mod = np.diff(src_crs.x_limits)[0]
                bboxes, proj_offset = self._bbox_and_offset(src_crs)
                x_lim = xs.min(), xs.max()
                for poly in bboxes:
                    # Arbitrarily choose the number of moduli to look
                    # above and below the -180->180 range. If data is beyond
                    # this range, we're not going to transform it quickly.
                    for i in [-1, 0, 1, 2]:
                        offset = mod * i - proj_offset
                        if ((poly[0] + offset) <= x_lim[0] and
                                (poly[1] + offset) >= x_lim[1]):
                            return_value = vertices + [[-offset, 0]]
                            break
                    if return_value is not None:
                        break

        return return_value


class TransverseMercator(Projection):
    """
    A Transverse Mercator projection.

    """
    _wrappable = True

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 scale_factor=1.0, globe=None, approx=False):
        """
        Parameters
        ----------
        central_longitude: optional
            The true longitude of the central meridian in degrees.
            Defaults to 0.
        central_latitude: optional
            The true latitude of the planar origin in degrees. Defaults to 0.
        false_easting: optional
            X offset from the planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from the planar origin in metres. Defaults to 0.
        scale_factor: optional
            Scale factor at the central meridian. Defaults to 1.

        globe: optional
            An instance of :class:`cartopy.crs.Globe`. If omitted, a default
            globe is created.

        approx: optional
            Whether to use Proj's approximate projection (True), or the new
            Extended Transverse Mercator code (False). Defaults to True, but
            will change to False in the next release.

        """
        proj4_params = [('proj', 'tmerc'), ('lon_0', central_longitude),
                        ('lat_0', central_latitude), ('k', scale_factor),
                        ('x_0', false_easting), ('y_0', false_northing),
                        ('units', 'm')]
        if approx:
            proj4_params += [('approx', None)]
        super().__init__(proj4_params, globe=globe)

        self.threshold = 1e4

    @property
    def boundary(self):
        x0, x1 = self.x_limits
        y0, y1 = self.y_limits
        return sgeom.LinearRing([(x0, y0), (x0, y1),
                                 (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        return (-2e7, 2e7)

    @property
    def y_limits(self):
        return (-1e7, 1e7)


class OSGB(TransverseMercator):
    def __init__(self, approx=False):
        super().__init__(central_longitude=-2, central_latitude=49,
                         scale_factor=0.9996012717,
                         false_easting=400000, false_northing=-100000,
                         globe=Globe(datum='OSGB36', ellipse='airy'),
                         approx=approx)

    @property
    def boundary(self):
        w = self.x_limits[1] - self.x_limits[0]
        h = self.y_limits[1] - self.y_limits[0]
        return sgeom.LinearRing([(0, 0), (0, h), (w, h), (w, 0), (0, 0)])

    @property
    def x_limits(self):
        return (0, 7e5)

    @property
    def y_limits(self):
        return (0, 13e5)


class OSNI(TransverseMercator):
    def __init__(self, approx=False):
        globe = Globe(semimajor_axis=6377340.189,
                      semiminor_axis=6356034.447938534)
        super().__init__(central_longitude=-8, central_latitude=53.5,
                         scale_factor=1.000035,
                         false_easting=200000, false_northing=250000,
                         globe=globe, approx=approx)

    @property
    def boundary(self):
        w = self.x_limits[1] - self.x_limits[0]
        h = self.y_limits[1] - self.y_limits[0]
        return sgeom.LinearRing([(0, 0), (0, h), (w, h), (w, 0), (0, 0)])

    @property
    def x_limits(self):
        return (18814.9667, 386062.3293)

    @property
    def y_limits(self):
        return (11764.8481, 464720.9559)


class UTM(Projection):
    """
    Universal Transverse Mercator projection.

    """

    def __init__(self, zone, southern_hemisphere=False, globe=None):
        """
        Parameters
        ----------
        zone
            The numeric zone of the UTM required.
        southern_hemisphere: optional
            Set to True if the zone is in the southern hemisphere. Defaults to
            False.
        globe: optional
            An instance of :class:`cartopy.crs.Globe`. If omitted, a default
            globe is created.

        """
        proj4_params = [('proj', 'utm'),
                        ('units', 'm'),
                        ('zone', zone)]
        if southern_hemisphere:
            proj4_params.append(('south', None))
        super().__init__(proj4_params, globe=globe)
        self.threshold = 1e2

    @property
    def boundary(self):
        x0, x1 = self.x_limits
        y0, y1 = self.y_limits
        return sgeom.LinearRing([(x0, y0), (x0, y1),
                                 (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        easting = 5e5
        # allow 50% overflow
        return (0 - easting / 2, 2 * easting + easting / 2)

    @property
    def y_limits(self):
        northing = 1e7
        # allow 50% overflow
        return (0 - northing, 2 * northing + northing / 2)


class EuroPP(UTM):
    """
    UTM Zone 32 projection for EuroPP domain.

    Ellipsoid is International 1924, Datum is ED50.

    """

    def __init__(self):
        globe = Globe(ellipse='intl')
        super().__init__(32, globe=globe)

    @property
    def x_limits(self):
        return (-1.4e6, 2e6)

    @property
    def y_limits(self):
        return (4e6, 7.9e6)


class Mercator(Projection):
    """
    A Mercator projection.

    """
    _wrappable = True

    def __init__(self, central_longitude=0.0,
                 min_latitude=-80.0, max_latitude=84.0,
                 globe=None, latitude_true_scale=None,
                 false_easting=0.0, false_northing=0.0, scale_factor=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The central longitude. Defaults to 0.
        min_latitude: optional
            The maximum southerly extent of the projection. Defaults
            to -80 degrees.
        max_latitude: optional
            The maximum northerly extent of the projection. Defaults
            to 84 degrees.
        globe: A :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.
        latitude_true_scale: optional
            The latitude where the scale is 1. Defaults to 0 degrees.
        false_easting: optional
            X offset from the planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from the planar origin in metres. Defaults to 0.
        scale_factor: optional
            Scale factor at natural origin. Defaults to unused.

        Notes
        -----
        Only one of ``latitude_true_scale`` and ``scale_factor`` should
        be included.
        """
        proj4_params = [('proj', 'merc'),
                        ('lon_0', central_longitude),
                        ('x_0', false_easting),
                        ('y_0', false_northing),
                        ('units', 'm')]

        # If it's None, we don't pass it to Proj4, in which case its default
        # of 0.0 will be used.
        if latitude_true_scale is not None:
            proj4_params.append(('lat_ts', latitude_true_scale))

        if scale_factor is not None:
            if latitude_true_scale is not None:
                raise ValueError('It does not make sense to provide both '
                                 '"scale_factor" and "latitude_true_scale". ')
            else:
                proj4_params.append(('k_0', scale_factor))

        super().__init__(proj4_params, globe=globe)

        # Need to have x/y limits defined for the initial hash which
        # gets used within transform_points for caching
        self._x_limits = self._y_limits = None
        # Calculate limits.
        minlon, maxlon = self._determine_longitude_bounds(central_longitude)
        limits = self.transform_points(self.as_geodetic(),
                                       np.array([minlon, maxlon]),
                                       np.array([min_latitude, max_latitude]))
        self._x_limits = tuple(limits[..., 0])
        self._y_limits = tuple(limits[..., 1])
        self.threshold = min(np.diff(self.x_limits)[0] / 720,
                             np.diff(self.y_limits)[0] / 360)

    def __eq__(self, other):
        res = super().__eq__(other)
        if hasattr(other, "_y_limits") and hasattr(other, "_x_limits"):
            res = res and self._y_limits == other._y_limits and \
                self._x_limits == other._x_limits
        return res

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.proj4_init, self._x_limits, self._y_limits))

    @property
    def boundary(self):
        x0, x1 = self.x_limits
        y0, y1 = self.y_limits
        return sgeom.LinearRing([(x0, y0), (x0, y1),
                                 (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


# Define a specific instance of a Mercator projection, the Google mercator.
Mercator.GOOGLE = Mercator(min_latitude=-85.0511287798066,
                           max_latitude=85.0511287798066,
                           globe=Globe(ellipse=None,
                                       semimajor_axis=WGS84_SEMIMAJOR_AXIS,
                                       semiminor_axis=WGS84_SEMIMAJOR_AXIS,
                                       nadgrids='@null'))
# Deprecated form
GOOGLE_MERCATOR = Mercator.GOOGLE


class LambertCylindrical(_RectangularProjection):
    def __init__(self, central_longitude=0.0, globe=None):
        globe = globe or Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS)
        proj4_params = [('proj', 'cea'), ('lon_0', central_longitude),
                        ('to_meter', math.radians(1) * (
                            globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS))]
        super().__init__(proj4_params, 180, math.degrees(1), globe=globe)


class LambertConformal(Projection):
    """
    A Lambert Conformal conic projection.

    """

    def __init__(self, central_longitude=-96.0, central_latitude=39.0,
                 false_easting=0.0, false_northing=0.0,
                 standard_parallels=(33, 45),
                 globe=None, cutoff=-30):
        """
        Parameters
        ----------
        central_longitude: optional
            The central longitude. Defaults to -96.
        central_latitude: optional
            The central latitude. Defaults to 39.
        false_easting: optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from planar origin in metres. Defaults to 0.
        standard_parallels: optional
            Standard parallel latitude(s). Defaults to (33, 45).
        globe: optional
            A :class:`cartopy.crs.Globe`. If omitted, a default globe is
            created.
        cutoff: optional
            Latitude of map cutoff.
            The map extends to infinity opposite the central pole
            so we must cut off the map drawing before then.
            A value of 0 will draw half the globe. Defaults to -30.

        """
        proj4_params = [('proj', 'lcc'),
                        ('lon_0', central_longitude),
                        ('lat_0', central_latitude),
                        ('x_0', false_easting),
                        ('y_0', false_northing)]

        n_parallels = len(standard_parallels)

        if not 1 <= n_parallels <= 2:
            raise ValueError('1 or 2 standard parallels must be specified. '
                             f'Got {n_parallels} ({standard_parallels})')

        proj4_params.append(('lat_1', standard_parallels[0]))
        if n_parallels == 2:
            proj4_params.append(('lat_2', standard_parallels[1]))

        super().__init__(proj4_params, globe=globe)

        # Compute whether this projection is at the "north pole" or the
        # "south pole" (after the central lon/lat have been taken into
        # account).
        if n_parallels == 1:
            plat = 90 if standard_parallels[0] > 0 else -90
        else:
            # Which pole are the parallels closest to? That is the direction
            # that the cone converges.
            if abs(standard_parallels[0]) > abs(standard_parallels[1]):
                poliest_sec = standard_parallels[0]
            else:
                poliest_sec = standard_parallels[1]
            plat = 90 if poliest_sec > 0 else -90

        self.cutoff = cutoff
        n = 91
        lons = np.empty(n + 2)
        lats = np.full(n + 2, float(cutoff))
        lons[0] = lons[-1] = 0
        lats[0] = lats[-1] = plat
        if plat == 90:
            # Ensure clockwise
            lons[1:-1] = np.linspace(central_longitude + 180 - 0.001,
                                     central_longitude - 180 + 0.001, n)
        else:
            lons[1:-1] = np.linspace(central_longitude - 180 + 0.001,
                                     central_longitude + 180 - 0.001, n)

        points = self.transform_points(self.as_geodetic(), lons, lats)

        self._boundary = sgeom.LinearRing(points)
        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]

        self.threshold = 1e5

    def __eq__(self, other):
        res = super().__eq__(other)
        if hasattr(other, "cutoff"):
            res = res and self.cutoff == other.cutoff
        return res

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.proj4_init, self.cutoff))

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class LambertZoneII(Projection):
    """
    Lambert zone II (extended) projection (https://epsg.io/27572), a
    legacy projection that covers hexagonal France and Corsica.

    """
    def __init__(self):
        crs = pyproj.CRS.from_epsg(27572)
        super().__init__(crs.to_wkt())

        # Projected bounds from https://epsg.io/27572
        self.bounds = [-5242.32, 1212512.16, 1589155.51, 2706796.21]


class LambertAzimuthalEqualArea(Projection):
    """
    A Lambert Azimuthal Equal-Area projection.

    """
    _wrappable = True

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 globe=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The central longitude. Defaults to 0.
        central_latitude: optional
            The central latitude. Defaults to 0.
        false_easting: optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: optional
            A :class:`cartopy.crs.Globe`. If omitted, a default globe is
            created.

        """
        proj4_params = [('proj', 'laea'),
                        ('lon_0', central_longitude),
                        ('lat_0', central_latitude),
                        ('x_0', false_easting),
                        ('y_0', false_northing)]

        super().__init__(proj4_params, globe=globe)

        a = float(self.ellipsoid.semi_major_metre or WGS84_SEMIMAJOR_AXIS)

        # Find the antipode, and shift it a small amount in latitude to
        # approximate the extent of the projection:
        lon = central_longitude + 180
        sign = np.sign(central_latitude) or 1
        lat = -central_latitude + sign * 0.01
        x, max_y = self.transform_point(lon, lat, self.as_geodetic())

        coords = _ellipse_boundary(a * 1.9999, max_y - false_northing,
                                   false_easting, false_northing, 61)
        self._boundary = sgeom.polygon.LinearRing(coords.T)
        mins = np.min(coords, axis=1)
        maxs = np.max(coords, axis=1)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]
        self.threshold = np.diff(self._x_limits)[0] * 1e-3

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class Miller(_RectangularProjection):
    _handles_ellipses = False

    def __init__(self, central_longitude=0.0, globe=None):
        if globe is None:
            globe = Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS, ellipse=None)

        a = globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS

        proj4_params = [('proj', 'mill'), ('lon_0', central_longitude)]
        # See Snyder, 1987. Eqs (11-1) and (11-2) substituting maximums of
        # (lambda-lambda0)=180 and phi=90 to get limits.
        super().__init__(proj4_params, a * np.pi, a * 2.303412543376391,
                         globe=globe)


class RotatedPole(_CylindricalProjection):
    """
    A rotated latitude/longitude projected coordinate system
    with cylindrical topology and projected distance.

    Coordinates are measured in projection metres.

    The class uses proj to perform an ob_tran operation, using the
    pole_longitude to set a lon_0 then performing two rotations based on
    pole_latitude and central_rotated_longitude.
    This is equivalent to setting the new pole to a location defined by
    the pole_latitude and pole_longitude values in the GeogCRS defined by
    globe, then rotating this new CRS about it's pole using the
    central_rotated_longitude value.

    """

    def __init__(self, pole_longitude=0.0, pole_latitude=90.0,
                 central_rotated_longitude=0.0, globe=None):
        """
        Parameters
        ----------
        pole_longitude: optional
            Pole longitude position, in unrotated degrees. Defaults to 0.
        pole_latitude: optional
            Pole latitude position, in unrotated degrees. Defaults to 0.
        central_rotated_longitude: optional
            Longitude rotation about the new pole, in degrees. Defaults to 0.
        globe: optional
            An optional :class:`cartopy.crs.Globe`. Defaults to a "WGS84"
            datum.

        """
        globe = globe or Globe(semimajor_axis=WGS84_SEMIMAJOR_AXIS)
        proj4_params = [('proj', 'ob_tran'), ('o_proj', 'latlon'),
                        ('o_lon_p', central_rotated_longitude),
                        ('o_lat_p', pole_latitude),
                        ('lon_0', 180 + pole_longitude),
                        ('to_meter', math.radians(1) * (
                            globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS))]
        super().__init__(proj4_params, 180, 90, globe=globe)


class Gnomonic(Projection):
    _handles_ellipses = False

    def __init__(self, central_latitude=0.0,
                 central_longitude=0.0, globe=None):
        proj4_params = [('proj', 'gnom'), ('lat_0', central_latitude),
                        ('lon_0', central_longitude)]
        super().__init__(proj4_params, globe=globe)
        self._max = 5e7
        self.threshold = 1e5

    @property
    def boundary(self):
        return sgeom.Point(0, 0).buffer(self._max).exterior

    @property
    def x_limits(self):
        return (-self._max, self._max)

    @property
    def y_limits(self):
        return (-self._max, self._max)


class Stereographic(Projection):
    _wrappable = True

    def __init__(self, central_latitude=0.0, central_longitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 true_scale_latitude=None,
                 scale_factor=None, globe=None):
        proj4_params = [('proj', 'stere'), ('lat_0', central_latitude),
                        ('lon_0', central_longitude),
                        ('x_0', false_easting), ('y_0', false_northing)]

        if true_scale_latitude is not None:
            if central_latitude not in (-90., 90.):
                warnings.warn('"true_scale_latitude" parameter is only used '
                              'for polar stereographic projections. Consider '
                              'the use of "scale_factor" instead.',
                              stacklevel=2)
            proj4_params.append(('lat_ts', true_scale_latitude))

        if scale_factor is not None:
            if true_scale_latitude is not None:
                raise ValueError('It does not make sense to provide both '
                                 '"scale_factor" and "true_scale_latitude". '
                                 'Ignoring "scale_factor".')
            else:
                proj4_params.append(('k_0', scale_factor))

        super().__init__(proj4_params, globe=globe)

        # TODO: Let the globe return the semimajor axis always.
        a = float(self.ellipsoid.semi_major_metre or WGS84_SEMIMAJOR_AXIS)
        b = float(self.ellipsoid.semi_minor_metre or WGS84_SEMIMINOR_AXIS)

        # Note: The magic number has been picked to maintain consistent
        # behaviour with a wgs84 globe. There is no guarantee that the scaling
        # should even be linear.
        x_axis_offset = 5e7 / WGS84_SEMIMAJOR_AXIS
        y_axis_offset = 5e7 / WGS84_SEMIMINOR_AXIS
        self._x_limits = (-a * x_axis_offset + false_easting,
                          a * x_axis_offset + false_easting)
        self._y_limits = (-b * y_axis_offset + false_northing,
                          b * y_axis_offset + false_northing)
        coords = _ellipse_boundary(self._x_limits[1], self._y_limits[1],
                                   false_easting, false_northing, 91)
        self._boundary = sgeom.LinearRing(coords.T)
        self.threshold = np.diff(self._x_limits)[0] * 1e-3

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class NorthPolarStereo(Stereographic):
    def __init__(self, central_longitude=0.0, true_scale_latitude=None,
                 globe=None):
        super().__init__(
            central_latitude=90,
            central_longitude=central_longitude,
            true_scale_latitude=true_scale_latitude,  # None is +90
            globe=globe)


class SouthPolarStereo(Stereographic):
    def __init__(self, central_longitude=0.0, true_scale_latitude=None,
                 globe=None):
        super().__init__(
            central_latitude=-90,
            central_longitude=central_longitude,
            true_scale_latitude=true_scale_latitude,  # None is -90
            globe=globe)


class Orthographic(Projection):
    _handles_ellipses = False

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 azimuth=0.0, globe=None):
        proj4_params = [('proj', 'ortho'), ('lon_0', central_longitude),
                        ('lat_0', central_latitude), ('alpha', azimuth)]
        if pyproj.__proj_version__ < '9.5.0' and azimuth != 0.0:
            warnings.warn(
                'Setting azimuth is not supported with PROJ versions < 9.5.0. '
                'Assuming azimuth=0. '
                f'Current PROJ version: {pyproj.__proj_version__}')
        super().__init__(proj4_params, globe=globe)

        # TODO: Let the globe return the semimajor axis always.
        a = float(self.ellipsoid.semi_major_metre or WGS84_SEMIMAJOR_AXIS)

        # To stabilise the projection of geometries, we reduce the boundary by
        # a tiny fraction at the cost of the extreme edges.
        coords = _ellipse_boundary(a * 0.99999, a * 0.99999, n=61)
        self._boundary = sgeom.polygon.LinearRing(coords.T)
        mins = np.min(coords, axis=1)
        maxs = np.max(coords, axis=1)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]
        self.threshold = np.diff(self._x_limits)[0] * 0.02

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class _WarpedRectangularProjection(Projection, metaclass=ABCMeta):
    _wrappable = True

    def __init__(self, proj4_params, central_longitude,
                 false_easting=None, false_northing=None, globe=None):
        if false_easting is not None:
            proj4_params += [('x_0', false_easting)]
        if false_northing is not None:
            proj4_params += [('y_0', false_northing)]
        super().__init__(proj4_params, globe=globe)

        # Obtain boundary points
        minlon, maxlon = self._determine_longitude_bounds(central_longitude)
        n = 91
        lon = np.empty(2 * n + 1)
        lat = np.empty(2 * n + 1)
        lon[:n] = minlon
        lat[:n] = np.linspace(-90, 90, n)
        lon[n:2 * n] = maxlon
        lat[n:2 * n] = np.linspace(90, -90, n)
        lon[-1] = minlon
        lat[-1] = -90
        points = self.transform_points(self.as_geodetic(), lon, lat)

        self._boundary = sgeom.LinearRing(points)

        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class Aitoff(_WarpedRectangularProjection):
    """
    An Aitoff projection.

    This projection is a modified azimuthal equidistant projection, balancing
    shape and scale distortion. There are no standard lines and only the
    central point is free of distortion.

    """

    _handles_ellipses = False

    def __init__(self, central_longitude=0, false_easting=None,
                 false_northing=None, globe=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        false_easting: float, optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: float, optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

            .. note::
                This projection does not handle elliptical globes.

        """
        proj_params = [('proj', 'aitoff'),
                       ('lon_0', central_longitude)]
        super().__init__(proj_params, central_longitude,
                         false_easting=false_easting,
                         false_northing=false_northing,
                         globe=globe)
        self.threshold = 1e5


class _Eckert(_WarpedRectangularProjection, metaclass=ABCMeta):
    """
    An Eckert projection.

    This class implements all the methods common to the Eckert family of
    projections.

    """

    _handles_ellipses = False

    def __init__(self, central_longitude=0, false_easting=None,
                 false_northing=None, globe=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        false_easting: float, optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: float, optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

            .. note::
                This projection does not handle elliptical globes.

        """
        proj4_params = [('proj', self._proj_name),
                        ('lon_0', central_longitude)]
        super().__init__(proj4_params, central_longitude,
                         false_easting=false_easting,
                         false_northing=false_northing,
                         globe=globe)
        self.threshold = 1e5


class EckertI(_Eckert):
    """
    An Eckert I projection.

    This projection is pseudocylindrical, but not equal-area. Both meridians
    and parallels are straight lines. Its equal-area pair is :class:`EckertII`.

    """
    _proj_name = 'eck1'


class EckertII(_Eckert):
    """
    An Eckert II projection.

    This projection is pseudocylindrical, and equal-area. Both meridians and
    parallels are straight lines. Its non-equal-area pair with equally-spaced
    parallels is :class:`EckertI`.

    """
    _proj_name = 'eck2'


class EckertIII(_Eckert):
    """
    An Eckert III projection.

    This projection is pseudocylindrical, but not equal-area. Parallels are
    equally-spaced straight lines, while meridians are elliptical arcs up to
    semicircles on the edges. Its equal-area pair is :class:`EckertIV`.

    """
    _proj_name = 'eck3'


class EckertIV(_Eckert):
    """
    An Eckert IV projection.

    This projection is pseudocylindrical, and equal-area. Parallels are
    unequally-spaced straight lines, while meridians are elliptical arcs up to
    semicircles on the edges. Its non-equal-area pair with equally-spaced
    parallels is :class:`EckertIII`.

    It is commonly used for world maps.

    """
    _proj_name = 'eck4'


class EckertV(_Eckert):
    """
    An Eckert V projection.

    This projection is pseudocylindrical, but not equal-area. Parallels are
    equally-spaced straight lines, while meridians are sinusoidal arcs. Its
    equal-area pair is :class:`EckertVI`.

    """
    _proj_name = 'eck5'


class EckertVI(_Eckert):
    """
    An Eckert VI projection.

    This projection is pseudocylindrical, and equal-area. Parallels are
    unequally-spaced straight lines, while meridians are sinusoidal arcs. Its
    non-equal-area pair with equally-spaced parallels is :class:`EckertV`.

    It is commonly used for world maps.

    """
    _proj_name = 'eck6'


class EqualEarth(_WarpedRectangularProjection):
    """
    An Equal Earth projection.

    This projection is pseudocylindrical, and equal area. Parallels are
    unequally-spaced straight lines, while meridians are equally-spaced arcs.

    It is intended for world maps.

    Note
    ----
    To use this projection, you must be using Proj 5.2.0 or newer.

    References
    ----------
    Bojan Šavrič, Tom Patterson & Bernhard Jenny (2018)
    The Equal Earth map projection,
    International Journal of Geographical Information Science,
    DOI: 10.1080/13658816.2018.1504949

    """

    def __init__(self, central_longitude=0, false_easting=None,
                 false_northing=None, globe=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        false_easting: float, optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: float, optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

        """
        proj_params = [('proj', 'eqearth'), ('lon_0', central_longitude)]
        super().__init__(proj_params, central_longitude,
                         false_easting=false_easting,
                         false_northing=false_northing,
                         globe=globe)
        self.threshold = 1e5


class Hammer(_WarpedRectangularProjection):
    """
    A Hammer projection.

    This projection is a modified `.LambertAzimuthalEqualArea` projection,
    similar to `.Aitoff`, and intended to reduce distortion in the outer
    meridians compared to `.Mollweide`. There are no standard lines and only
    the central point is free of distortion.

    """

    _handles_ellipses = False

    def __init__(self, central_longitude=0, false_easting=None,
                 false_northing=None, globe=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        false_easting: float, optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: float, optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

            .. note::
                This projection does not handle elliptical globes.

        """
        proj4_params = [('proj', 'hammer'),
                        ('lon_0', central_longitude)]
        super().__init__(proj4_params, central_longitude,
                         false_easting=false_easting,
                         false_northing=false_northing,
                         globe=globe)
        self.threshold = 1e5


class Mollweide(_WarpedRectangularProjection):
    """
    A Mollweide projection.

    This projection is pseudocylindrical, and equal area. Parallels are
    unequally-spaced straight lines, while meridians are elliptical arcs up to
    semicircles on the edges. Poles are points.

    It is commonly used for world maps, or interrupted with several central
    meridians.

    """

    _handles_ellipses = False

    def __init__(self, central_longitude=0, globe=None,
                 false_easting=None, false_northing=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        false_easting: float, optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: float, optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

            .. note::
                This projection does not handle elliptical globes.

        """
        proj4_params = [('proj', 'moll'), ('lon_0', central_longitude)]
        super().__init__(proj4_params, central_longitude,
                         false_easting=false_easting,
                         false_northing=false_northing,
                         globe=globe)
        self.threshold = 1e5


class Robinson(_WarpedRectangularProjection):
    """
    A Robinson projection.

    This projection is pseudocylindrical, and a compromise that is neither
    equal-area nor conformal. Parallels are unequally-spaced straight lines,
    and meridians are curved lines of no particular form.

    It is commonly used for "visually-appealing" world maps.

    """

    _handles_ellipses = False

    def __init__(self, central_longitude=0, globe=None,
                 false_easting=None, false_northing=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        false_easting: float, optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: float, optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

            .. note::
                This projection does not handle elliptical globes.

        """
        proj4_params = [('proj', 'robin'), ('lon_0', central_longitude)]
        super().__init__(proj4_params, central_longitude,
                         false_easting=false_easting,
                         false_northing=false_northing,
                         globe=globe)
        self.threshold = 1e4

    def transform_point(self, x, y, src_crs, trap=True):
        """
        Capture and handle any input NaNs, else invoke parent function,
        :meth:`_WarpedRectangularProjection.transform_point`.

        Needed because input NaNs can trigger a fatal error in the underlying
        implementation of the Robinson projection.

        Note
        ----
            Although the original can in fact translate (nan, lat) into
            (nan, y-value), this patched version doesn't support that.

        """
        if np.isnan(x) or np.isnan(y):
            result = (np.nan, np.nan)
        else:
            result = super().transform_point(x, y, src_crs, trap=trap)
        return result

    def transform_points(self, src_crs, x, y, z=None, trap=False):
        """
        Capture and handle NaNs in input points -- else as parent function,
        :meth:`_WarpedRectangularProjection.transform_points`.

        Needed because input NaNs can trigger a fatal error in the underlying
        implementation of the Robinson projection.

        Note
        ----
            Although the original can in fact translate (nan, lat) into
            (nan, y-value), this patched version doesn't support that.
            Instead, we invalidate any of the points that contain a NaN.

        """
        input_point_nans = np.isnan(x) | np.isnan(y)
        if z is not None:
            input_point_nans |= np.isnan(z)
        handle_nans = np.any(input_point_nans)
        if handle_nans:
            # Remove NaN points from input data to avoid the error.
            x[input_point_nans] = 0.0
            y[input_point_nans] = 0.0
            if z is not None:
                z[input_point_nans] = 0.0
        result = super().transform_points(src_crs, x, y, z, trap=trap)
        if handle_nans:
            # Result always has shape (N, 3).
            # Blank out each (whole) point where we had a NaN in the input.
            result[input_point_nans] = np.nan
        return result


class InterruptedGoodeHomolosine(Projection):
    """
    Composite equal-area projection emphasizing either land or
    ocean features.

    Original Reference:
        Goode, J. P., 1925: The Homolosine Projection: A new device for
        portraying the Earth's surface entire. Annals of the
        Association of American Geographers, 15:3, 119-125,
        DOI: 10.1080/00045602509356949

    A central_longitude value of -160 is recommended for the oceanic view.

    """
    _wrappable = True

    def __init__(self, central_longitude=0, globe=None, emphasis='land'):
        """
        Parameters
        ----------
        central_longitude : float, optional
            The central longitude, by default 0
        globe : :class:`cartopy.crs.Globe`, optional
            If omitted, a default Globe object is created, by default None
        emphasis : str, optional
            Options 'land' and 'ocean' are available, by default 'land'
        """

        if emphasis == 'land':
            proj4_params = [('proj', 'igh'), ('lon_0', central_longitude)]
            super().__init__(proj4_params, globe=globe)

        elif emphasis == 'ocean':
            proj4_params = [('proj', 'igh_o'), ('lon_0', central_longitude)]
            super().__init__(proj4_params, globe=globe)

        else:
            msg = '`emphasis` needs to be either \'land\' or \'ocean\''
            raise ValueError(msg)

        minlon, maxlon = self._determine_longitude_bounds(central_longitude)
        epsilon = 1e-10

        # Obtain boundary points
        n = 31
        if emphasis == 'land':
            top_interrupted_lons = (-40.0,)
            bottom_interrupted_lons = (80.0, -20.0, -100.0)
        elif emphasis == 'ocean':
            top_interrupted_lons = (-90.0, 60.0)
            bottom_interrupted_lons = (90.0, -60.0)
        lons = np.empty(
            (2 + 2 * len(top_interrupted_lons + bottom_interrupted_lons)) * n +
            1)
        lats = np.empty(
            (2 + 2 * len(top_interrupted_lons + bottom_interrupted_lons)) * n +
            1)
        end = 0

        # Left boundary
        lons[end:end + n] = minlon
        lats[end:end + n] = np.linspace(-90, 90, n)
        end += n

        # Top boundary
        for lon in top_interrupted_lons:
            lons[end:end + n] = lon - epsilon + central_longitude
            lats[end:end + n] = np.linspace(90, 0, n)
            end += n
            lons[end:end + n] = lon + epsilon + central_longitude
            lats[end:end + n] = np.linspace(0, 90, n)
            end += n

        # Right boundary
        lons[end:end + n] = maxlon
        lats[end:end + n] = np.linspace(90, -90, n)
        end += n

        # Bottom boundary
        for lon in bottom_interrupted_lons:
            lons[end:end + n] = lon + epsilon + central_longitude
            lats[end:end + n] = np.linspace(-90, 0, n)
            end += n
            lons[end:end + n] = lon - epsilon + central_longitude
            lats[end:end + n] = np.linspace(0, -90, n)
            end += n

        # Close loop
        lons[-1] = minlon
        lats[-1] = -90

        points = self.transform_points(self.as_geodetic(), lons, lats)
        self._boundary = sgeom.LinearRing(points)

        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]

        self.threshold = 2e4

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class _Satellite(Projection):
    def __init__(self, projection, satellite_height=35785831,
                 central_longitude=0.0, central_latitude=0.0,
                 false_easting=0, false_northing=0, globe=None,
                 sweep_axis=None):
        proj4_params = [('proj', projection), ('lon_0', central_longitude),
                        ('lat_0', central_latitude), ('h', satellite_height),
                        ('x_0', false_easting), ('y_0', false_northing),
                        ('units', 'm')]
        if sweep_axis:
            proj4_params.append(('sweep', sweep_axis))
        super().__init__(proj4_params, globe=globe)

    def _set_boundary(self, coords):
        self._boundary = sgeom.LinearRing(coords.T)
        mins = np.min(coords, axis=1)
        maxs = np.max(coords, axis=1)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]
        self.threshold = np.diff(self._x_limits)[0] * 0.02

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class Geostationary(_Satellite):
    """
    A view appropriate for satellites in Geostationary Earth orbit.

    Perspective view looking directly down from above a point on the equator.

    In this projection, the projected coordinates are scanning angles measured
    from the satellite looking directly downward, multiplied by the height of
    the satellite.

    """

    def __init__(self, central_longitude=0.0, satellite_height=35785831,
                 false_easting=0, false_northing=0, globe=None,
                 sweep_axis='y'):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        satellite_height: float, optional
            The height of the satellite. Defaults to 35785831 metres
            (true geostationary orbit).
        false_easting:
            X offset from planar origin in metres. Defaults to 0.
        false_northing:
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.
        sweep_axis: 'x' or 'y', optional. Defaults to 'y'.
            Controls which axis is scanned first, and thus which angle is
            applied first. The default is appropriate for Meteosat, while
            'x' should be used for GOES.
        """

        super().__init__(
            projection='geos',
            satellite_height=satellite_height,
            central_longitude=central_longitude,
            central_latitude=0.0,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe,
            sweep_axis=sweep_axis)

        # TODO: Let the globe return the semimajor axis always.
        a = float(self.ellipsoid.semi_major_metre or WGS84_SEMIMAJOR_AXIS)
        b = float(self.ellipsoid.semi_minor_metre or WGS84_SEMIMINOR_AXIS)
        h = float(satellite_height)

        # To find the bound we trace around where the line from the satellite
        # is tangent to the surface. This involves trigonometry on a sphere
        # centered at the satellite. The two scanning angles form two legs of
        # triangle on this sphere--the hypotenuse "c" (angle arc) is controlled
        # by distance from center to the edge of the ellipse being seen.

        # This is one of the angles in the spherical triangle and used to
        # rotate around and "scan" the boundary
        angleA = np.linspace(0, -2 * np.pi, 91)  # Clockwise boundary.

        # Convert the angle around center to the proper value to use in the
        # parametric form of an ellipse
        th = np.arctan(a / b * np.tan(angleA))

        # Given the position on the ellipse, what is the distance from center
        # to the ellipse--and thus the tangent point
        r = np.hypot(a * np.cos(th), b * np.sin(th))
        sat_dist = a + h

        # Using this distance, solve for sin and tan of c in the triangle that
        # includes the satellite, Earth center, and tangent point--we need to
        # figure out the location of this tangent point on the elliptical
        # cross-section through the Earth towards the satellite, where the
        # major axis is a and the minor is r. With the ellipse centered on the
        # Earth and the satellite on the y-axis (at y = a + h = sat_dist), the
        # equation for an ellipse and some calculus gives us the tangent point
        # (x0, y0) as:
        # y0 = a**2 / sat_dist
        # x0 = r * np.sqrt(1 - a**2 / sat_dist**2)
        # which gives:
        # sin_c = x0 / np.hypot(x0, sat_dist - y0)
        # tan_c = x0 / (sat_dist - y0)
        # A bit of algebra combines these to give directly:
        sin_c = r / np.sqrt(sat_dist ** 2 - a ** 2 + r ** 2)
        tan_c = r / np.sqrt(sat_dist ** 2 - a ** 2)

        # Using Napier's rules for right spherical triangles R2 and R6,
        # (See https://en.wikipedia.org/wiki/Spherical_trigonometry), we can
        # solve for arc angles b and a, which are our x and y scanning angles,
        # respectively.
        coords = np.vstack([np.arctan(np.cos(angleA) * tan_c),  # R6
                            np.arcsin(np.sin(angleA) * sin_c)])  # R2

        # Need to multiply scanning angles by satellite height to get to the
        # actual native coordinates for the projection.
        coords *= h
        coords += np.array([[false_easting], [false_northing]])
        self._set_boundary(coords)


class NearsidePerspective(_Satellite):
    """
    Perspective view looking directly down from above a point on the globe.

    In this projection, the projected coordinates are x and y measured from
    the origin of a plane tangent to the Earth directly below the perspective
    point (e.g. a satellite).

    """

    _handles_ellipses = False

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 satellite_height=35785831,
                 false_easting=0, false_northing=0, globe=None):
        """
        Parameters
        ----------
        central_longitude: float, optional
            The central longitude. Defaults to 0.
        central_latitude: float, optional
            The central latitude. Defaults to 0.
        satellite_height: float, optional
            The height of the satellite. Defaults to 35785831 meters
            (true geostationary orbit).
        false_easting:
            X offset from planar origin in metres. Defaults to 0.
        false_northing:
            Y offset from planar origin in metres. Defaults to 0.
        globe: :class:`cartopy.crs.Globe`, optional
            If omitted, a default globe is created.

            .. note::
                This projection does not handle elliptical globes.

        """
        super().__init__(
            projection='nsper',
            satellite_height=satellite_height,
            central_longitude=central_longitude,
            central_latitude=central_latitude,
            false_easting=false_easting,
            false_northing=false_northing,
            globe=globe)

        # TODO: Let the globe return the semimajor axis always.
        a = self.ellipsoid.semi_major_metre or WGS84_SEMIMAJOR_AXIS

        h = float(satellite_height)
        max_x = a * np.sqrt(h / (2 * a + h))
        coords = _ellipse_boundary(max_x, max_x,
                                   false_easting, false_northing, 61)
        self._set_boundary(coords)


class AlbersEqualArea(Projection):
    """
    An Albers Equal Area projection

    This projection is conic and equal-area, and is commonly used for maps of
    the conterminous United States.

    """

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 standard_parallels=(20.0, 50.0), globe=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The central longitude. Defaults to 0.
        central_latitude: optional
            The central latitude. Defaults to 0.
        false_easting: optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from planar origin in metres. Defaults to 0.
        standard_parallels: optional
            The one or two latitudes of correct scale. Defaults to (20, 50).
        globe: optional
            A :class:`cartopy.crs.Globe`. If omitted, a default globe is
            created.

        """
        proj4_params = [('proj', 'aea'),
                        ('lon_0', central_longitude),
                        ('lat_0', central_latitude),
                        ('x_0', false_easting),
                        ('y_0', false_northing)]
        if standard_parallels is not None:
            try:
                proj4_params.append(('lat_1', standard_parallels[0]))
                try:
                    proj4_params.append(('lat_2', standard_parallels[1]))
                except IndexError:
                    pass
            except TypeError:
                proj4_params.append(('lat_1', standard_parallels))

        super().__init__(proj4_params, globe=globe)

        # bounds
        minlon, maxlon = self._determine_longitude_bounds(central_longitude)
        n = 103
        lons = np.empty(2 * n + 1)
        lats = np.empty(2 * n + 1)
        tmp = np.linspace(minlon, maxlon, n)
        lons[:n] = tmp
        lats[:n] = 90
        lons[n:-1] = tmp[::-1]
        lats[n:-1] = -90
        lons[-1] = lons[0]
        lats[-1] = lats[0]

        points = self.transform_points(self.as_geodetic(), lons, lats)

        self._boundary = sgeom.LinearRing(points)
        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]

        self.threshold = 1e5

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class AzimuthalEquidistant(Projection):
    """
    An Azimuthal Equidistant projection

    This projection provides accurate angles about and distances through the
    central position. Other angles, distances, or areas may be distorted.
    """
    _wrappable = True

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 globe=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The true longitude of the central meridian in degrees.
            Defaults to 0.
        central_latitude: optional
            The true latitude of the planar origin in degrees.
            Defaults to 0.
        false_easting: optional
            X offset from the planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from the planar origin in metres. Defaults to 0.
        globe: optional
            An instance of :class:`cartopy.crs.Globe`. If omitted, a default
            globe is created.

        """
        proj4_params = [('proj', 'aeqd'), ('lon_0', central_longitude),
                        ('lat_0', central_latitude),
                        ('x_0', false_easting), ('y_0', false_northing)]
        super().__init__(proj4_params, globe=globe)

        # TODO: Let the globe return the semimajor axis always.
        a = float(self.ellipsoid.semi_major_metre or WGS84_SEMIMAJOR_AXIS)
        b = float(self.ellipsoid.semi_minor_metre or a)

        coords = _ellipse_boundary(a * np.pi, b * np.pi,
                                   false_easting, false_northing, 61)
        self._boundary = sgeom.LinearRing(coords.T)
        mins = np.min(coords, axis=1)
        maxs = np.max(coords, axis=1)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]

        self.threshold = 1e5

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class Sinusoidal(Projection):
    """
    A Sinusoidal projection.

    This projection is equal-area.

    """

    def __init__(self, central_longitude=0.0, false_easting=0.0,
                 false_northing=0.0, globe=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The central longitude. Defaults to 0.
        false_easting: optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from planar origin in metres. Defaults to 0.
        globe: optional
            A :class:`cartopy.crs.Globe`. If omitted, a default globe is
            created.

        """
        proj4_params = [('proj', 'sinu'),
                        ('lon_0', central_longitude),
                        ('x_0', false_easting),
                        ('y_0', false_northing)]
        super().__init__(proj4_params, globe=globe)

        # Obtain boundary points
        minlon, maxlon = self._determine_longitude_bounds(central_longitude)
        points = []
        n = 91
        lon = np.empty(2 * n + 1)
        lat = np.empty(2 * n + 1)
        lon[:n] = minlon
        lat[:n] = np.linspace(-90, 90, n)
        lon[n:2 * n] = maxlon
        lat[n:2 * n] = np.linspace(90, -90, n)
        lon[-1] = minlon
        lat[-1] = -90
        points = self.transform_points(self.as_geodetic(), lon, lat)

        self._boundary = sgeom.LinearRing(points)
        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]
        self.threshold = max(np.abs(self.x_limits + self.y_limits)) * 1e-5

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


# MODIS data products use a Sinusoidal projection of a spherical Earth
# https://modis-land.gsfc.nasa.gov/GCTP.html
Sinusoidal.MODIS = Sinusoidal(globe=Globe(ellipse=None,
                                          semimajor_axis=6371007.181,
                                          semiminor_axis=6371007.181))


class EquidistantConic(Projection):
    """
    An Equidistant Conic projection.

    This projection is conic and equidistant, and the scale is true along all
    meridians and along one or two specified standard parallels.
    """

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 standard_parallels=(20.0, 50.0), globe=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The central longitude. Defaults to 0.
        central_latitude: optional
            The true latitude of the planar origin in degrees. Defaults to 0.
        false_easting: optional
            X offset from planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from planar origin in metres. Defaults to 0.
        standard_parallels: optional
            The one or two latitudes of correct scale. Defaults to (20, 50).
        globe: optional
            A :class:`cartopy.crs.Globe`. If omitted, a default globe is
            created.

        """
        proj4_params = [('proj', 'eqdc'),
                        ('lon_0', central_longitude),
                        ('lat_0', central_latitude),
                        ('x_0', false_easting),
                        ('y_0', false_northing)]
        if standard_parallels is not None:
            try:
                proj4_params.append(('lat_1', standard_parallels[0]))
                try:
                    proj4_params.append(('lat_2', standard_parallels[1]))
                except IndexError:
                    pass
            except TypeError:
                proj4_params.append(('lat_1', standard_parallels))

        super().__init__(proj4_params, globe=globe)

        # bounds
        n = 103
        lons = np.empty(2 * n + 1)
        lats = np.empty(2 * n + 1)
        minlon, maxlon = self._determine_longitude_bounds(central_longitude)
        tmp = np.linspace(minlon, maxlon, n)
        lons[:n] = tmp
        lats[:n] = 90
        lons[n:-1] = tmp[::-1]
        lats[n:-1] = -90
        lons[-1] = lons[0]
        lats[-1] = lats[0]

        points = self.transform_points(self.as_geodetic(), lons, lats)

        self._boundary = sgeom.LinearRing(points)
        mins = np.min(points, axis=0)
        maxs = np.max(points, axis=0)
        self._x_limits = mins[0], maxs[0]
        self._y_limits = mins[1], maxs[1]

        self.threshold = 1e5

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class ObliqueMercator(Projection):
    """
    An Oblique Mercator projection.

    """
    _wrappable = True

    def __init__(self, central_longitude=0.0, central_latitude=0.0,
                 false_easting=0.0, false_northing=0.0,
                 scale_factor=1.0, azimuth=0.0, globe=None):
        """
        Parameters
        ----------
        central_longitude: optional
            The true longitude of the central meridian in degrees.
            Defaults to 0.
        central_latitude: optional
            The true latitude of the planar origin in degrees. Defaults to 0.
        false_easting: optional
            X offset from the planar origin in metres. Defaults to 0.
        false_northing: optional
            Y offset from the planar origin in metres. Defaults to 0.
        scale_factor: optional
            Scale factor at the central meridian. Defaults to 1.
        azimuth: optional
            Azimuth of centerline clockwise from north at the center point of
            the centre line. Defaults to 0.
        globe: optional
            An instance of :class:`cartopy.crs.Globe`. If omitted, a default
            globe is created.

        Notes
        -----
        The 'Rotated Mercator' projection can be achieved using Oblique
        Mercator with `azimuth` ``=90``.

        """

        if np.isclose(azimuth, 90):
            # Exactly 90 causes coastline 'folding'.
            azimuth -= 1e-3

        proj4_params = [('proj', 'omerc'), ('lonc', central_longitude),
                        ('lat_0', central_latitude), ('k', scale_factor),
                        ('x_0', false_easting), ('y_0', false_northing),
                        ('alpha', azimuth), ('units', 'm')]

        super().__init__(proj4_params, globe=globe)

        # Couple limits to those of Mercator - delivers acceptable plots, and
        #  Mercator has been through much more scrutiny.
        mercator = Mercator(
            central_longitude=central_longitude,
            globe=globe,
            false_easting=false_easting,
            false_northing=false_northing,
            scale_factor=scale_factor,
        )
        self._x_limits = mercator.x_limits
        self._y_limits = mercator.y_limits
        self.threshold = mercator.threshold

    @property
    def boundary(self):
        x0, x1 = self.x_limits
        y0, y1 = self.y_limits
        return sgeom.LinearRing([(x0, y0), (x0, y1),
                                 (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits

class Spilhaus(Projection):
    """
    Spilhaus World Ocean Map in a Square.

    This is a projection based on Adams World in a Square II projection with the
    two major antipodal areas on land: South China
    (115°E and 30°N) and Argentina (65°W and 30°S).
    See https://storymaps.arcgis.com/stories/756bcae18d304a1eac140f19f4d5cb3d
    """
    def __init__(self, rotation=45, false_easting=0.0, false_northing=0.0, globe=None):
        """
        Parameters
        ----------
        rotation : optional
            Clockwise rotation of the map in degrees. Defaults to 45.
        false_easting : optional
            X offset from the planar origin in metres. Defaults to 0.0.
        false_northing : optional
            Y offset from the planar origin in metres. Defaults to 0.0.

        """
        proj4_params = [('proj', 'spilhaus'),
                        ('rot',rotation),
                        ('x_0', false_easting),
                        ('y_0', false_northing)]

        super().__init__(proj4_params, globe=globe)
        # The boundary on https://epsg.io/54099 are wrong
        # The following bounds are calculated based on
        #[-65.00000012, -29.99999981]
        # and [115.00000024,  30.00000036]
        self.bounds = [
            -11802684.083372328,
            11802683.949222516,
            -11801129.925928915,
            11801129.925928915
        ]


def _ring_is_ccw(ring):
    """Return True if *ring* is counter-clockwise via the shoelace formula.

    More robust than ``shapely.LinearRing.is_ccw`` for self-intersecting
    rings (e.g. when a projected polygon spans exactly ±180° longitude and
    boundary attachment creates a degenerate "tail").  The signed-area sum
    cancels self-intersecting tail contributions that would fool Shapely's
    winding test.
    """
    coords = shapely.get_coordinates(ring)
    x, y = coords[:, 0], coords[:, 1]
    return bool(np.dot(x[:-1], y[1:]) > np.dot(x[1:], y[:-1]))


def epsg(code):
    """
    Return the projection which corresponds to the given EPSG code.

    The EPSG code must correspond to a "projected coordinate system",
    so EPSG codes such as 4326 (WGS-84) which define a "geodetic coordinate
    system" will not work.

    Note
    ----
        The conversion is performed by pyproj.CRS.

    """
    import cartopy._epsg
    return cartopy._epsg._EPSGProjection(code)
