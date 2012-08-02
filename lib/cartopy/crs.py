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

import numpy
import shapely.geometry as sgeom

from cartopy._crs import CRS, Geocentric, Geodetic
import cartopy.trace


class Projection(CRS):
    """
    Defines a projected coordinate system with flat topology and Euclidean
    distance.

    """
    __metaclass__ = ABCMeta
    
    _method_map = {
        'LineString': '_project_linear',
        'LinearRing': '_project_linear',
        'Polygon': '_project_polygon',
        'MultiLineString': '_project_multiline',
        'MultiPolygon': '_project_multipolygon',
    }
    
    def __eq__(self, other):
        # XXX handle params that have been set to the default value on one, but not the other? 
        return isinstance(self, type(other)) and self.proj4_params == other.proj4_params
    
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
    def domain(self):
        try:
            domain = self._domain
        except AttributeError:
            domain = self._domain = sgeom.Polygon(self.boundary)
        return domain

    def _as_mpl_axes(self):
        import cartopy.mpl_integration.geoaxes as geoaxes
        return geoaxes.GenericProjectionAxes, {'map_projection': self}
    
    def _as_mpl_transform(self, axes=None):
        import cartopy.mpl_integration.geoaxes as geoaxes
        if not isinstance(axes, geoaxes.GenericProjectionAxes):
            raise ValueError('Axes should be an instance of GenericProjectionAxes, got %s' % type(axes))
        if axes.projection == self:
            return axes.transData
        else:
            return geoaxes.InterProjectionTransform(self, axes.projection) + axes.transData
        
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
            raise ValueError('Unsupported geometry type {!r}'.format(geom_type))
        return getattr(self, method_name)(geometry, src_crs)

    def _project_linear(self, geometry, src_crs):
        return cartopy.trace.project_linear(geometry, src_crs, self)

    def _project_multiline(self, geometry, src_crs):
        geoms = []
        for geom in geometry.geoms:
            r = self._project_linear(geom, src_crs)
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
        return sgeom.MultiPolygon(geoms)

    def _project_polygon(self, polygon, src_crs):
        """
        Returns the projected geometry derived from the given geometry.

        TODO Handle polygons containing holes.

        """
        # Convert initial lines into split, projected, contiguous segments.
        # 12341 -> (def23ghi, jkl41abc)
        #
        # 1---2             d-e-f-2       1-a-b-c
        # |   |     =>            |       |
        # 4---3             i-h-g-3       4-l-k-j
        #

        # XXX This shouldn't really be here?
        # The simple 2-dimensional determination of orientation provided
        # by Shapely/GEOS doesn't apply when considering geodetic coordinates.
        # In fact, there is no generic concept of orientation applicable
        # to geodetic coordinates. One could assume the polygon should
        # be oriented to select the smaller section of the surface, but
        # even this breaks down for polygons which divide the surface
        # into two pieces of equal area.
#        if not src_crs.is_geodetic():
        polygon = sgeom.polygon.orient(polygon, -1)
        split, multi_line_string = self._polygon_to_multi_line_string(polygon, src_crs)

        # Patch up the edges if the polygon was split.
        if split:
            multi_line_string = self._attach_to_boundary(multi_line_string)

        polygons = self._convert_to_polygons(multi_line_string)

        return polygons

    def _polygon_to_multi_line_string(self, polygon, src_crs):
        # 1) Resolve the initial lines into projected segments
        # 1abc
        # def23ghi
        # jkl41
        multi_line_string = cartopy.trace.project_linear(polygon.exterior,
                                                         src_crs, self)

        # 2) Stitch together the two ends (if appropriate)
        # i.e. Does the first point of the first line match the last point of
        # the last line?
        # def23ghi
        # jkl41abc
        split = len(multi_line_string) > 1
        if split:
            # XXX Clumsy! (NB. multi_line_string[-1] causes a MemoryFault from shapely)
            line_strings = list(multi_line_string)
            if numpy.allclose(line_strings[0].coords[0], line_strings[-1].coords[-1]):
                last_coords = list(line_strings[-1].coords)
                first_coords = list(line_strings[0].coords)[1:]
                line_strings[-1] = sgeom.LineString(last_coords + first_coords)
                multi_line_string = sgeom.MultiLineString(line_strings[1:])
        elif len(multi_line_string) == 1:
            # If the end-points line up exactly with the map boundary we might not
            # have a split - but we still need to stitch the ends back together.
            split = not numpy.allclose(multi_line_string[0].coords[0], multi_line_string[0].coords[-1])

        return split, multi_line_string

    def _attach_to_boundary(self, multi_line_string):
        """
        Returns a new MultiLineString by attaching the ends of the given lines
        to the boundary, paying attention to the traversal directions of the
        lines and boundary.

        """
        # Accumulate all the boundary and segment end points, along with
        # their distance along the boundary.
        edge_things = []

        # Convert the boundary to a LineString so we can compute distances
        # along it. (It doesn't work with a LinearRing)
        boundary = sgeom.LineString(self.boundary)

        def boundary_distance(xy):
            return boundary.project(sgeom.Point(*xy))

        # Record the positions of all the segment ends
        for i, line_string in enumerate(multi_line_string):
            first_dist = boundary_distance(line_string.coords[0])
            thing = _Thing(first_dist, False, (i, 'first', line_string.coords[0]))
            edge_things.append(thing)
            last_dist = boundary_distance(line_string.coords[-1])
            thing = _Thing(last_dist, False, (i, 'last', line_string.coords[-1]))
            edge_things.append(thing)

        # Record the positions of all the boundary vertices
        for xy in list(boundary.coords)[:-1]:
            point = sgeom.Point(*xy)
            dist = boundary.project(point)
            thing = _Thing(dist, True, point)
            edge_things.append(thing)

        # Order everything as if walking around the boundary.
        # NB. We make line end-points take precedence over boundary points
        # to ensure that end-points are still found and followed when they coincide.
        edge_things.sort(key=lambda thing: (thing.distance, thing.kind))
        debug = 0
        if debug:
            print
            print 'Edge things'
            for thing in edge_things:
                print '   ', thing

        to_do = {i: line_string for i, line_string in enumerate(multi_line_string)}
        done = []
        while to_do:
            i, line_string = to_do.popitem()
            if debug:
                import sys
                sys.stdout.write('+')
                sys.stdout.flush()
                print
                print 'Processing: %s, %s' % (i, line_string)
            edge_things = filter(lambda t: t.kind or t.data[0] != i or t.data[1] != 'last', edge_things)

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
                    combined_coords = list(line_string.coords) + [(boundary_point.x, boundary_point.y)]
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
                    line_to_append = multi_line_string[j]
                    # XXX pelson: this was failing, before I added the if statement, but I don't understand how it could be failing... 
                    if j in to_do:
                        del to_do[j]
                    coords_to_append = list(line_to_append.coords)
                    if next_thing.data[1] == 'last':
                        coords_to_append = coords_to_append[::-1]
                    line_string = sgeom.LineString(list(line_string.coords) + coords_to_append)
        multi_line_string = done

        if debug:
            print '   DONE'

        return multi_line_string

    def _convert_to_polygons(self, multi_line_string):
        # XXX Optimise
        polygons = []
        for line_string in multi_line_string:
            polygon = sgeom.Polygon(line_string)
            # XXX pelson: handle null geometries properly
            if polygon._geom is None:
                continue
            if polygon.exterior.is_ccw:
                # Compute the overall domain of the boundary and the inside-out polygon.
                x1, y1, x2, y2 = polygon.bounds
                bx = (x2 - x1) * 0.1
                by = (y2 - y1) * 0.1
                x1 -= bx
                y1 -= by
                x2 += bx
                y2 += by
                boundary_poly = sgeom.Polygon(self.boundary)
                x3, y3, x4, y4 = boundary_poly.bounds
                bx = (x4 - x3) * 0.1
                by = (y4 - y3) * 0.1
                x3 -= bx
                y3 -= by
                x4 += bx
                y4 += by
                box = sgeom.box(min(x1, x3), min(y1, y3), max(x2, x4), max(y2, y4))

                # Invert the polygon
                polygon = box.difference(polygon)

                # Intersect the inverted polygon with the boundary
                polygon = boundary_poly.intersection(polygon)

                if not polygon.is_empty:
                    polygons.append(polygon)
            else:
                polygons.append(polygon)

        if polygons:
            polygons = sgeom.MultiPolygon(polygons)

        return polygons


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


class PlateCarree(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'eqc', 'lon_0': central_longitude, 'a': math.degrees(1)}
        super(PlateCarree, self).__init__(proj4_params, 180, 90)

    @property
    def threshold(self):
        return 0.5


class TransverseMercator(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'tmerc', 'lon_0': central_longitude, 'a': math.degrees(1)}
        super(TransverseMercator, self).__init__(proj4_params, 180, 90)

    @property
    def threshold(self):
        return 0.5


# XXX Could become a subclass of TransverseMercator if it exposed enough parameters?
class OSGB(Projection):
    def __init__(self):
        proj4_params = {'proj': 'tmerc', 'lat_0': 49, 'lon_0': -2, 'k': 0.9996012717,
            'x_0': 400000, 'y_0': -100000, 'ellps': 'airy', 'datum': 'OSGB36',
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


class Mercator(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'merc', 'lon_0': central_longitude, 'a': math.degrees(1)}
        super(Mercator, self).__init__(proj4_params, 180, 180)

    @property
    def threshold(self):
        return 0.5


class LambertCylindrical(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'cea', 'lon_0': central_longitude, 'a': math.degrees(1)}
        super(LambertCylindrical, self).__init__(proj4_params, 180, math.degrees(1))

    @property
    def threshold(self):
        return 0.5


class Miller(_RectangularProjection):
    def __init__(self, central_longitude=0.0):
        proj4_params = {'proj': 'mill', 'lon_0': central_longitude, 'a': math.degrees(1)}
        # XXX How can we derive the vertical limit of 131.98?
        super(Miller, self).__init__(proj4_params, 180, 131.98)

    @property
    def threshold(self):
        return 0.5


class RotatedPole(_RectangularProjection):
    def __init__(self, pole_longitude=0.0, pole_latitude=90.0):
        proj4_params = {'proj': 'ob_tran', 'o_proj': 'latlon', 'o_lon_p': 0,
            'o_lat_p': pole_latitude, 'lon_0': 180 + pole_longitude,
            #'to_meter': math.degrees(1)}
            'to_meter': math.radians(1)}
            #}
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
    def __init__(self, central_latitude=0.0):
        proj4_params = {'proj': 'stere', 'lat_0': central_latitude}
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
    def __init__(self):
        super(NorthPolarStereo, self).__init__(90)


class SouthPolarStereo(Stereographic):
    def __init__(self):
        super(SouthPolarStereo, self).__init__(-90)


class Orthographic(Projection):
    def __init__(self, central_longitude=0.0, central_latitude=0.0):
        proj4_params = {'proj': 'ortho', 'lon_0': central_longitude, 'lat_0': central_latitude}
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
        for lat in numpy.linspace(-90, 90, n):
            points.append(self.transform_point(180 + central_longitude, lat, geodetic_crs))
        for lat in numpy.linspace(90, -90, n):
            points.append(self.transform_point(-180 + central_longitude, lat, geodetic_crs))
        points.append(self.transform_point(180 + central_longitude, -90, geodetic_crs))

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
        proj4_params = {'proj': 'robin', 'lon_0': central_longitude}
        super(Robinson, self).__init__(proj4_params, central_longitude)

    @property
    def threshold(self):
        return 1e5
    

class InterruptedGoodeHomolosine(Projection):
    def __init__(self, central_longitude=0):
        proj4_params = {'proj': 'igh', 'lon_0': central_longitude}
        super(InterruptedGoodeHomolosine, self).__init__(proj4_params)

    @property
    def threshold(self):
        return 1e5


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
