# (C) British Crown Copyright 2011 - 2012, Met Office
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
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.


"""
The crs module defines Coordinate Reference Systems and the transformations
between them.

"""
from abc import ABCMeta, abstractproperty
import math
import warnings

import numpy as np
import shapely.geometry as sgeom
from shapely.geometry.polygon import LinearRing
from shapely.prepared import prep

from cartopy._crs import CRS, Geocentric, Geodetic, PROJ4_RELEASE
import cartopy.trace


__document_these__ = ['CRS', 'Geocentric', 'Geodetic']


class RotatedGeodetic(CRS):
    """
    Defines a rotated latitude/longitude coordinate system with spherical
    topology and geographical distance.

    Coordinates are measured in degrees.

    """
    def __init__(self, pole_longitude, pole_latitude, ellipse='WGS84',
                 datum='WGS84'):
        """
        Create a RotatedGeodetic CRS.

        Args:

            * pole_longitude - Pole longitude position, in unrotated degrees.
            * pole_latitude - Pole latitude position, in unrotated degrees.

        Kwargs:

            * ellipse      - Ellipsoid definition.
            * datum        - Datum definition.

        """
        proj4_params = {'proj': 'ob_tran', 'o_proj': 'latlon', 'o_lon_p': 0,
                        'o_lat_p': pole_latitude,
                        'lon_0': 180 + pole_longitude,
                        'to_meter': math.radians(1), 'ellps': ellipse,
                        'datum': datum}
        super(RotatedGeodetic, self).__init__(proj4_params)


class Projection(CRS):
    """
    Defines a projected coordinate system with flat topology and Euclidean
    distance.

    """
    __metaclass__ = ABCMeta

    _method_map = {
        'LineString': '_project_line_string',
        'LinearRing': '_project_linear_ring',
        'Polygon': '_project_polygon',
        'MultiLineString': '_project_multiline',
        'MultiPolygon': '_project_multipolygon',
    }

    def __eq__(self, other):
        # XXX handle params that have been set to the default value on one,
        # but not the other?
        return (isinstance(self, type(other)) and
                self.proj4_params == other.proj4_params)

    def __ne__(self, other):
        return not self == other

    @abstractproperty
    def boundary(self):
        pass

    @abstractproperty
    def threshold(self):
        pass

    @abstractproperty
    def x_limits(self):
        pass

    @abstractproperty
    def y_limits(self):
        pass

    @property
    def cw_boundary(self):
        try:
            boundary = self._cw_boundary
        except AttributeError:
            boundary = sgeom.LineString(self.boundary)
            self._cw_boundary = boundary
        return boundary

    @property
    def ccw_boundary(self):
        try:
            boundary = self._ccw_boundary
        except AttributeError:
            boundary = sgeom.LineString(list(self.boundary.coords)[::-1])
            self._ccw_boundary = boundary
        return boundary

    @property
    def domain(self):
        try:
            domain = self._domain
        except AttributeError:
            domain = self._domain = sgeom.Polygon(self.boundary)
        return domain

    def _as_mpl_axes(self):
        import cartopy.mpl.geoaxes as geoaxes
        return geoaxes.GeoAxes, {'map_projection': self}

    def project_geometry(self, geometry, src_crs=None):
        """
        Projects the given geometry into this projection.

        :param geometry: The geometry to (re-)project.
        :param src_crs: The source CRS, or geodetic CRS if None.
        :rtype: Shapely geometry.

        If src_crs is None, the source CRS is assumed to be a geodetic
        version of the target CRS.

        """
        if src_crs is None:
            src_crs = self.as_geodetic()
        elif not isinstance(src_crs, CRS):
            raise TypeError('Source CRS must be an instance of CRS'
                            ' or one of its subclasses, or None.')
        geom_type = geometry.geom_type
        method_name = self._method_map.get(geom_type)
        if not method_name:
            raise ValueError('Unsupported geometry '
                             'type {!r}'.format(geom_type))
        return getattr(self, method_name)(geometry, src_crs)

    def _project_line_string(self, geometry, src_crs):
        return cartopy.trace.project_linear(geometry, src_crs, self)

    def _project_linear_ring(self, linear_ring, src_crs):
        """
        Projects the given LinearRing from the src_crs into this CRS and
        returns the resultant LinearRing or MultiLineString.

        """
        # 1) Resolve the initial lines into projected segments
        # 1abc
        # def23ghi
        # jkl41
        multi_line_string = cartopy.trace.project_linear(linear_ring,
                                                         src_crs, self)

        # 2) Simplify the segments where appropriate.
        if len(multi_line_string) > 1:
            # Stitch together segments which are close to continuous.
            # This is important when:
            # 1) The first source point projects into the map and the
            # ring has been cut by the boundary.
            # Continuing the example from above this gives:
            #   def23ghi
            #   jkl41abc
            # 2) The cut ends of segments are too close to reliably
            # place into an order along the boundary.
            line_strings = list(multi_line_string)
            any_modified = False
            i = 0
            while i < len(line_strings):
                modified = False
                j = 0
                while j < len(line_strings):
                    if i != j and np.allclose(line_strings[i].coords[0],
                                              line_strings[j].coords[-1],
                                              atol=self.threshold):
                        last_coords = list(line_strings[j].coords)
                        first_coords = list(line_strings[i].coords)[1:]
                        combo = sgeom.LineString(last_coords + first_coords)
                        if j < i:
                            i, j = j, i
                        del line_strings[j], line_strings[i]
                        line_strings.append(combo)
                        modified = True
                        any_modified = True
                        break
                    else:
                        j += 1
                if not modified:
                    i += 1
            if any_modified:
                multi_line_string = sgeom.MultiLineString(line_strings)

        # 3) Check for a single resulting ring.
        if (len(multi_line_string) == 1 and
                len(multi_line_string[0].coords) > 3 and
                np.allclose(multi_line_string[0].coords[0],
                            multi_line_string[0].coords[-1])):
            result_geometry = LinearRing(multi_line_string[0].coords[:-1])
        else:
            result_geometry = multi_line_string

        return result_geometry

    def _project_multiline(self, geometry, src_crs):
        geoms = []
        for geom in geometry.geoms:
            r = self._project_line_string(geom, src_crs)
            if r:
                geoms.extend(r.geoms)
        if geoms:
            return sgeom.MultiLineString(geoms)
        else:
            return []

    def _project_multipolygon(self, geometry, src_crs):
        geoms = []
        for geom in geometry.geoms:
            r = self._project_polygon(geom, src_crs)
            if r:
                geoms.extend(r.geoms)
        if geoms:
            result = sgeom.MultiPolygon(geoms)
        else:
            result = sgeom.MultiPolygon()
        return result

    def _project_polygon(self, polygon, src_crs):
        """
        Returns the projected polygon(s) derived from the given polygon.

        """
        # Determine orientation of polygon.
        # TODO: Consider checking the internal rings have the opposite
        # orientation to the external rings?
        if src_crs.is_geodetic():
            is_ccw = True
        else:
            is_ccw = polygon.exterior.is_ccw

        # Project the polygon exterior/interior rings.
        # Each source ring will result in either a ring, or one or more
        # lines.
        rings = []
        multi_lines = []
        for src_ring in [polygon.exterior] + list(polygon.interiors):
            geometry = self._project_linear_ring(src_ring, src_crs)
            if geometry.geom_type == 'LinearRing':
                rings.append(geometry)
            else:
                multi_lines.append(geometry)

        # Convert any lines to rings by attaching them to the boundary.
        if multi_lines:
            rings.extend(self._attach_lines_to_boundary(multi_lines, is_ccw))

        # Resolve all the inside vs. outside rings, and convert to the
        # final MultiPolygon.
        return self._rings_to_multi_polygon(rings, is_ccw)

    def _attach_lines_to_boundary(self, multi_line_strings, is_ccw):
        """
        Returns a list of LinearRings by attaching the ends of the given lines
        to the boundary, paying attention to the traversal directions of the
        lines and boundary.

        """
        # Accumulate all the boundary and segment end points, along with
        # their distance along the boundary.
        edge_things = []

        # Get the boundary as a LineString of the correct orientation
        # so we can compute distances along it.
        if is_ccw:
            boundary = self.ccw_boundary
        else:
            boundary = self.cw_boundary

        def boundary_distance(xy):
            return boundary.project(sgeom.Point(*xy))

        # Squash all the LineStrings into a single list.
        line_strings = []
        for multi_line_string in multi_line_strings:
            line_strings.extend(multi_line_string)

        # Record the positions of all the segment ends
        for i, line_string in enumerate(line_strings):
            first_dist = boundary_distance(line_string.coords[0])
            thing = _Thing(first_dist, False,
                           (i, 'first', line_string.coords[0]))
            edge_things.append(thing)
            last_dist = boundary_distance(line_string.coords[-1])
            thing = _Thing(last_dist, False,
                           (i, 'last', line_string.coords[-1]))
            edge_things.append(thing)

        # Record the positions of all the boundary vertices
        for xy in list(boundary.coords)[:-1]:
            point = sgeom.Point(*xy)
            dist = boundary.project(point)
            thing = _Thing(dist, True, point)
            edge_things.append(thing)

        # Order everything as if walking around the boundary.
        # NB. We make line end-points take precedence over boundary points
        # to ensure that end-points are still found and followed when they
        # coincide.
        edge_things.sort(key=lambda thing: (thing.distance, thing.kind))
        debug = 0
        if debug:
            print
            print 'Edge things'
            for thing in edge_things:
                print '   ', thing

        to_do = {i: line_string for i, line_string in enumerate(line_strings)}
        done = []
        while to_do:
            i, line_string = to_do.popitem()
            if debug:
                import sys
                sys.stdout.write('+')
                sys.stdout.flush()
                print
                print 'Processing: %s, %s' % (i, line_string)
            filter_fn = lambda t: (t.kind or
                                   t.data[0] != i or
                                   t.data[1] != 'last')
            edge_things = filter(filter_fn, edge_things)

            while True:
                # Find the distance of the last point
                d_last = boundary_distance(line_string.coords[-1])
                if debug:
                    print '   d_last:', d_last
                next_thing = _find_gt(edge_things, d_last)
                if debug:
                    print '   next_thing:', next_thing
                if next_thing.kind:
                    if debug:
                        print '   adding boundary point'
                    boundary_point = next_thing.data
                    combined_coords = (list(line_string.coords) +
                                       [(boundary_point.x, boundary_point.y)])
                    line_string = sgeom.LineString(combined_coords)
                    # XXX
                    #edge_things.remove(next_thing)
                elif next_thing.data[0] == i:
                    if debug:
                        print '   close loop'
                    done.append(line_string)
                    break
                else:
                    if debug:
                        print '   adding line'
                    j = next_thing.data[0]
                    line_to_append = line_strings[j]
                    # XXX pelson: I think this if statement can be removed
                    if j in to_do:
                        del to_do[j]
                    coords_to_append = list(line_to_append.coords)
                    if next_thing.data[1] == 'last':
                        coords_to_append = coords_to_append[::-1]
                    line_string = sgeom.LineString((list(line_string.coords) +
                                                    coords_to_append))

        # filter out any non-valid linear rings
        done = filter(lambda linear_ring: len(linear_ring.coords) > 2, done)

        # XXX Is the last point in each ring actually the same as the first?
        linear_rings = [LinearRing(line) for line in done]

        if debug:
            print '   DONE'

        return linear_rings

    def _rings_to_multi_polygon(self, rings, is_ccw):
        exterior_rings = []
        interior_rings = []
        for ring in rings:
            if ring.is_ccw != is_ccw:
                interior_rings.append(ring)
            else:
                exterior_rings.append(ring)

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
            for ring in interior_rings:
                polygon = sgeom.Polygon(ring)
                if polygon.is_valid:
                    x1, y1, x2, y2 = polygon.bounds
                    bx = (x2 - x1) * 0.1
                    by = (y2 - y1) * 0.1
                    x1 -= bx
                    y1 -= by
                    x2 += bx
                    y2 += by
                    box = sgeom.box(min(x1, x3), min(y1, y3),
                                    max(x2, x4), max(y2, y4))

                    # Invert the polygon
                    polygon = box.difference(polygon)

                    # Intersect the inverted polygon with the boundary
                    polygon = boundary_poly.intersection(polygon)

                    if not polygon.is_empty:
                        polygon_bits.append(polygon)

        if polygon_bits:
            multi_poly = sgeom.MultiPolygon(polygon_bits)
        else:
            multi_poly = sgeom.MultiPolygon()
        return multi_poly


class _RectangularProjection(Projection):
    """
    The abstract superclass of projections with a rectangular domain which
    is symmetric about the origin.

    """
    def __init__(self, proj4_params, half_width, half_height):
        self._half_width = half_width
        self._half_height = half_height
        super(_RectangularProjection, self).__init__(proj4_params)

    @property
    def boundary(self):
        # XXX Should this be a LinearRing?
        w, h = self._half_width, self._half_height
        return sgeom.LineString([(-w, -h), (-w, h), (w, h), (w, -h), (-w, -h)])

    @property
    def x_limits(self):
        return (-self._half_width, self._half_width)

    @property
    def y_limits(self):
        return (-self._half_height, self._half_height)


class _CylindricalProjection(_RectangularProjection):
    """
    The abstract class which denotes cylindrical projections where we
    want to allow x values to wrap around.

    """


class PlateCarree(_CylindricalProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'eqc', 'lon_0': central_longitude,
                        'a': math.degrees(1)}
        super(PlateCarree, self).__init__(proj4_params, 180, 90)

    @property
    def threshold(self):
        return 0.5


class TransverseMercator(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'tmerc', 'lon_0':
                        central_longitude, 'a': math.degrees(1)}
        super(TransverseMercator, self).__init__(proj4_params, 180, 90)

    @property
    def threshold(self):
        return 0.5


# XXX Could become a subclass of TransverseMercator if it exposed enough
# parameters?
class OSGB(Projection):
    def __init__(self):
        proj4_params = {'proj': 'tmerc', 'lat_0': 49, 'lon_0': -2,
                        'k': 0.9996012717, 'x_0': 400000, 'y_0': -100000,
                        'ellps': 'airy', 'datum': 'OSGB36',
                        'units': 'm', 'no_defs': ''}
        super(OSGB, self).__init__(proj4_params)

    @property
    def threshold(self):
        return 1e4

    @property
    def boundary(self):
        # XXX Should this be a LinearRing?
        w, h = 7e5, 13e5
        return sgeom.LineString([(0, 0), (0, h), (w, h), (w, 0), (0, 0)])

    @property
    def x_limits(self):
        return (0, 7e5)

    @property
    def y_limits(self):
        return (0, 13e5)


class OSNI(Projection):
    def __init__(self):
        proj4_params = {'proj': 'tmerc', 'lat_0': 53.5, 'lon_0': -8,
                        'k': 1.000035, 'x_0': 200000, 'y_0': 250000,
                        'a': 6377340.189, 'b': 6356034.447938534,
                        'units': 'm', 'no_defs': ''}
        super(OSNI, self).__init__(proj4_params)

    @property
    def threshold(self):
        return 1e4

    @property
    def boundary(self):
        x0, x1 = self.x_limits
        w = x1 - x0
        y0, y1 = self.y_limits
        h = y1 - y0
        # XXX Should this be a LinearRing?
        return sgeom.LineString([(0, 0), (0, h), (w, h), (w, 0), (0, 0)])

    @property
    def x_limits(self):
        return (18814.9667, 386062.3293)

    @property
    def y_limits(self):
        return (11764.8481, 464720.9559)


class EuroPP(Projection):
    """
    UTM Zone 32 projection for EuroPP domain.

    Ellipsoid is International 1924, Datum is ED50.

    """
    def __init__(self):
        proj4_params = {'proj': 'tmerc',
                        'lat_0': 50, 'lon_0': 9,
                        'k': 0.9996,
                        'x_0': 1750000, 'y_0': 1500000,
                        'zone': 32,
                        'ellps': 'intl',
                        'units': 'm',
                        'towgs84': '-87,-98,-121',
                        'no_defs': ''}
        super(EuroPP, self).__init__(proj4_params)

    @property
    def boundary(self):
        w, h = 3.19e6, 3.8e6
        return shapely.geometry.LineString([(0, 0), (0, h), (w, h),
                                            (w, 0), (0, 0)])

    @property
    def x_limits(self):
        return (0, 3.19e6)

    @property
    def y_limits(self):
        return (0, 3.8e6)

    @property
    def threshold(self):
        return 1e4


class Mercator(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'merc', 'lon_0': central_longitude,
                        'a': math.degrees(1)}
        super(Mercator, self).__init__(proj4_params, 180, 180)

    @property
    def threshold(self):
        return 0.5


class LambertCylindrical(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'cea', 'lon_0': central_longitude,
                        'a': math.degrees(1)}
        super(LambertCylindrical, self).__init__(proj4_params, 180,
                                                 math.degrees(1))

    @property
    def threshold(self):
        return 0.5


class Miller(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'mill', 'lon_0': central_longitude,
                        'a': math.degrees(1)}
        # XXX How can we derive the vertical limit of 131.98?
        super(Miller, self).__init__(proj4_params, 180, 131.98)

    @property
    def threshold(self):
        return 0.5


class RotatedPole(_CylindricalProjection):
    def __init__(self, pole_longitude=0.0, pole_latitude=90.0):
        proj4_params = {'proj': 'ob_tran', 'o_proj': 'latlon', 'o_lon_p': 0,
                        'o_lat_p': pole_latitude,
                        'lon_0': 180 + pole_longitude,
                        'to_meter': math.radians(1)
                        }
        super(RotatedPole, self).__init__(proj4_params, 180, 90)

    @property
    def threshold(self):
        return 0.5


class Gnomonic(Projection):
    def __init__(self, central_latitude=0.0):
        proj4_params = {'proj': 'gnom', 'lat_0': central_latitude}
        super(Gnomonic, self).__init__(proj4_params)
        self._max = 5e7

    @property
    def boundary(self):
        return sgeom.Point(0, 0).buffer(self._max).exterior

    @property
    def threshold(self):
        return 1e5

    @property
    def x_limits(self):
        return (-self._max, self._max)

    @property
    def y_limits(self):
        return (-self._max, self._max)


class Stereographic(Projection):
    def __init__(self, central_latitude=0.0, central_longitude=0.0):
        proj4_params = {'proj': 'stere', 'lat_0': central_latitude,
                        'lon_0': central_longitude}
        super(Stereographic, self).__init__(proj4_params)
        self._max = 5e7

    @property
    def boundary(self):
        return sgeom.Point(0, 0).buffer(self._max).exterior

    @property
    def threshold(self):
        return 1e5

    @property
    def x_limits(self):
        return (-self._max, self._max)

    @property
    def y_limits(self):
        return (-self._max, self._max)


class NorthPolarStereo(Stereographic):
    def __init__(self, central_longitude=0.0):
        super(NorthPolarStereo, self).__init__(
            central_latitude=90,
            central_longitude=central_longitude)


class SouthPolarStereo(Stereographic):
    def __init__(self, central_longitude=0.0):
        super(SouthPolarStereo, self).__init__(
            central_latitude=-90,
            central_longitude=central_longitude)


class Orthographic(Projection):
    def __init__(self, central_longitude=0.0, central_latitude=0.0):
        proj4_params = {'proj': 'ortho', 'lon_0':
                        central_longitude, 'lat_0': central_latitude}
        super(Orthographic, self).__init__(proj4_params)
        self._max = 6.4e6

    @property
    def boundary(self):
        return sgeom.Point(0, 0).buffer(self._max).exterior

    @property
    def threshold(self):
        return 1e5

    @property
    def x_limits(self):
        return (-self._max, self._max)

    @property
    def y_limits(self):
        return (-self._max, self._max)


class _WarpedRectangularProjection(Projection):
    def __init__(self, proj4_params, central_longitude):
        super(_WarpedRectangularProjection, self).__init__(proj4_params)

        # Obtain boundary points
        points = []
        n = 91
        geodetic_crs = self.as_geodetic()
        for lat in np.linspace(-90, 90, n):
            points.append(
                self.transform_point(180 + central_longitude,
                                     lat, geodetic_crs)
            )
        for lat in np.linspace(90, -90, n):
            points.append(
                self.transform_point(-180 + central_longitude,
                                     lat, geodetic_crs)
            )
        points.append(
            self.transform_point(180 + central_longitude, -90, geodetic_crs))

        self._boundary = sgeom.LineString(points[::-1])

        x = [p[0] for p in points]
        y = [p[1] for p in points]
        self._x_limits = min(x), max(x)
        self._y_limits = min(y), max(y)

    @property
    def boundary(self):
        return self._boundary

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class Mollweide(_WarpedRectangularProjection):
    def __init__(self, central_longitude=0):
        proj4_params = {'proj': 'moll', 'lon_0': central_longitude}
        super(Mollweide, self).__init__(proj4_params, central_longitude)

    @property
    def threshold(self):
        return 1e5


class Robinson(_WarpedRectangularProjection):
    def __init__(self, central_longitude=0):
        # Warn when using Robinson with proj4 4.8 due to discontinuity at
        # 40 deg N introduced by incomplete fix to issue #113 (see
        # https://trac.osgeo.org/proj/ticket/113).
        import re
        match = re.search(r"\d\.\d", PROJ4_RELEASE)
        if match is not None:
            proj4_version = float(match.group())
            if proj4_version >= 4.8:
                warnings.warn('The Robinson projection from Proj.4 versions '
                              '4.8.0 and later contains a discontinuity at '
                              '40 deg latitude. Use this projection with '
                              'caution.')
        else:
            warnings.warn('Cannot determine Proj.4 version. The Robinson '
                          'projection may be unreliable and should be used '
                          'with caution.')

        proj4_params = {'proj': 'robin', 'lon_0': central_longitude}
        super(Robinson, self).__init__(proj4_params, central_longitude)

    @property
    def threshold(self):
        return 1e4


class InterruptedGoodeHomolosine(Projection):
    def __init__(self, central_longitude=0):
        proj4_params = {'proj': 'igh', 'lon_0': central_longitude}
        super(InterruptedGoodeHomolosine, self).__init__(proj4_params)

        # Obtain boundary points
        points = []
        n = 31
        geodetic_crs = self.as_geodetic()

        # Right boundary
        for lat in np.linspace(-90, 90, n):
            points.append(self.transform_point(180 + central_longitude,
                                               lat, geodetic_crs))

        # Top boundary
        interrupted_lons = (-40.0,)
        delta = 0.001
        for lon in interrupted_lons:
            for lat in np.linspace(90, 0, n):
                points.append(self.transform_point(lon + delta +
                                                   central_longitude,
                                                   lat, geodetic_crs))
            for lat in np.linspace(0, 90, n):
                points.append(self.transform_point(lon - delta +
                                                   central_longitude,
                                                   lat, geodetic_crs))

        # Left boundary
        for lat in np.linspace(90, -90, n):
            points.append(self.transform_point(-180 + central_longitude,
                                               lat, geodetic_crs))

        # Bottom boundary
        interrupted_lons = (-100.0, -20.0, 80.0)
        delta = 0.001
        for lon in interrupted_lons:
            for lat in np.linspace(-90, 0, n):
                points.append(self.transform_point(lon - delta +
                                                   central_longitude,
                                                   lat, geodetic_crs))
            for lat in np.linspace(0, -90, n):
                points.append(self.transform_point(lon + delta +
                                                   central_longitude,
                                                   lat, geodetic_crs))

        # Close loop
        points.append(self.transform_point(180 + central_longitude, -90,
                                           geodetic_crs))

        self._boundary = sgeom.LineString(points[::-1])

        x = [p[0] for p in points]
        y = [p[1] for p in points]
        self._x_limits = min(x), max(x)
        self._y_limits = min(y), max(y)

    @property
    def boundary(self):
        return self._boundary

    @property
    def threshold(self):
        return 2e4

    @property
    def x_limits(self):
        return self._x_limits

    @property
    def y_limits(self):
        return self._y_limits


class _Thing(object):
    def __init__(self, distance, kind, data):
        self.distance = distance
        self.kind = kind
        self.data = data

    def __repr__(self):
        return '_Thing(%r, %r, %s)' % (self.distance, self.kind, self.data)


def _find_gt(a, x):
    for v in a:
        # TODO: Fix the problem of co-incident boundary & line points
        #if v.distance >= x:
        if v.distance > x:
            return v
    return a[0]
