"""
Combines the shapefile access of the Python Shapefile Library with the
geometry representation of shapely.

>>> import os.path
>>> import cartopy.io.shapereader as shapereader
>>> filename = os.path.join(os.path.dirname(shapereader.__file__), 'data', 'Devon')
>>> reader = shapereader.Reader(filename)
>>> len(reader)
1
>>> list(reader.records()) #doctest: +ELLIPSIS
[<Record: <shapely.geometry.multipolygon.MultiPolygon object at ...>, {'PMS_REGION': 14, 'SHAPE_AREA': 6597719517.55, 'OBJECTID': 15, 'COUNTRY': 'ENGLAND', 'SNAC_GOR': 'South West', 'COUNTY_STR': 'Devon', 'SHAPE_LEN': 570341.652865}, <fields>>]
>>> list(reader.geometries()) #doctest: +ELLIPSIS
[<shapely.geometry.multipolygon.MultiPolygon object at ...>]

"""
from shapely.geometry import MultiLineString, MultiPolygon, Point, Polygon
import shapefile


__version__ = '0.1'

__all__ = ['Reader', 'Record']


def _create_point(shape):
    return Point(shape.points[0])


def _create_polyline(shape):
    parts = list(shape.parts) + [None]
    bounds = zip(parts[:-1], parts[1:])
    lines = [shape.points[slice(lower, upper)] for lower, upper in bounds]
    return MultiLineString(lines)


def _create_polygon(shape):
    # Partition the shapefile rings into outer rings/polygons (clockwise) and
    # inner rings/holes (anti-clockwise).
    parts = list(shape.parts) + [None]
    bounds = zip(parts[:-1], parts[1:])
    outer_polygons_and_holes = []
    inner_polygons = []
    for lower, upper in bounds:
        polygon = Polygon(shape.points[slice(lower, upper)])
        if polygon.exterior.is_ccw:
            inner_polygons.append(polygon)
        else:
            outer_polygons_and_holes.append((polygon, []))

    # Find the appropriate outer ring for each inner ring.
    # aka. Group the holes with their containing polygons.
    for inner_polygon in inner_polygons:
        for outer_polygon, holes in outer_polygons_and_holes:
            if outer_polygon.contains(inner_polygon):
                holes.append(inner_polygon.exterior.coords)
                break

    polygon_defns = [(outer_polygon.exterior.coords, holes) for outer_polygon, holes in outer_polygons_and_holes]
    return MultiPolygon(polygon_defns)


def _make_geometry(geometry_factory, shape):
    geometry = None
    if shape.shapeType != shapefile.NULL:
        geometry = geometry_factory(shape)
    return geometry


# The mapping from shapefile shapeType values to geometry creation functions.
GEOMETRY_FACTORIES = {
    shapefile.POINT: _create_point,
    shapefile.POLYLINE: _create_polyline,
    shapefile.POLYGON: _create_polygon,
}


class Record(object):
    """
    A single logical entry from a shapefile, combining the attributes with
    their associated geometry.

        attributes - A dictionary mapping attribute names to attribute values.
        bounds     - A tuple of (minx, miny, maxx, maxy).
        fields     - A list of field definitions, as per the Python Shapefile Library.
        geometry   - A shapely.geometry instance or None if it was a null shape.

    """
    def __init__(self, shape, geometry_factory, attributes, fields):
        self._shape = shape
        self._geometry_factory = geometry_factory

        if hasattr(shape, 'bbox'):
            self.bounds = tuple(shape.bbox)

        self.attributes = attributes
        self.fields = fields

    def __repr__(self):
        return '<Record: %r, %r, <fields>>' % (self.geometry, self._attributes)

    def __str__(self):
        return 'Record(%s, %s, <fields>)' % (self.geometry, self._attributes)

    def __getattr__(self, name):
        if name == 'bounds':
            value = self.bounds = self.geometry().bounds
        elif name == 'geometry':
            value = self.geometry = _make_geometry(self._geometry_factory, self._shape)
        else:
            value = object.__getattribute__(name)
        return value


class Reader(object):
    """
    Provides iterator based access to the contents of a shapefile.
    
    The shapefile geometry is expressed as ``shapely.geometry`` instances.

    """
    def __init__(self, filename):
        # Validate the filename/shapefile
        self._reader = reader = shapefile.Reader(filename)
        if reader.shp is None or reader.shx is None or reader.dbf is None:
            raise ValueError("Incomplete shapefile definition in '%s'." % filename)

        # Figure out how to make appropriate shapely geometry instances
        shapeType = reader.shapeType
        self._geometry_factory = GEOMETRY_FACTORIES.get(shapeType)
        if self._geometry_factory is None:
            raise ValueError('Unsupported shape type: %s' % shapeType)

        self.fields = self._reader.fields

    def __len__(self):
        return self._reader.numRecords

    def geometries(self):
        """Returns an iterator of shapely geometries."""
        geometry_factory = self._geometry_factory
        for i in xrange(self._reader.numRecords):
            shape = self._reader.shape(i)
            yield _make_geometry(geometry_factory, shape)

    def records(self):
        """Returns an iterator of Record instances."""
        geometry_factory = self._geometry_factory
        # Ignore the "DeletionFlag" field which always comes first
        fields = self._reader.fields[1:]
        field_names = [field[0] for field in fields]
        for i in xrange(self._reader.numRecords):
            shape_record = self._reader.shapeRecord(i)
            attributes = dict(zip(field_names, shape_record.record))
            yield Record(shape_record.shape, geometry_factory, attributes, fields)
