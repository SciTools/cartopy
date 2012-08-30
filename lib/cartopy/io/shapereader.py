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
import os


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
        return '<Record: %r, %r, <fields>>' % (self.geometry, self.attributes)

    def __str__(self):
        return 'Record(%s, %s, <fields>)' % (self.geometry, self.attributes)

    def __getattr__(self, name):
        if name == 'bounds':
            value = self.bounds = self.geometry().bounds
        elif name == 'geometry':
            value = self.geometry = _make_geometry(self._geometry_factory, self._shape)
        else:
            value = object.__getattribute__(self, name)
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


def natural_earth(resolution='110m', category='physical', name='coastline', data_dir=None):
    """
    Returns the path to the requested natural earth shapefile, downloading and unziping if necessary.
    
    """
    import glob
    
    if data_dir is None:
        dname = os.path.dirname
        # be more clever in the data directory so that users can define a setting.
        data_dir = os.path.join(dname(dname(__file__)), 'data', 'shapefiles', 'natural_earth')
        
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    full_name = '%s-%s' % (resolution, name)
    shape_dir = os.path.join(data_dir, full_name)
    if not os.path.exists(shape_dir):
        os.makedirs(shape_dir)
        
    # find the only shapefile in the directory. This is because NE have inconsistent zip file naming conventions.
    glob_pattern = os.path.join(data_dir, full_name, '*.shp')
    shapefiles = glob.glob(glob_pattern)
    
    if not shapefiles:
        # download the zip file
        import urllib2
        import cStringIO as StringIO
        from zipfile import ZipFile
        # note the repeated http. That is intentional
        file_url = ('http://www.naturalearthdata.com/http//www.naturalearthdata.com/'
                    'download/%s/%s/%s.zip' % (resolution, category, full_name))
        
        shapefile_online = urllib2.urlopen(file_url)
        zfh = ZipFile(StringIO.StringIO(shapefile_online.read()), 'r')
        zfh.extractall(shape_dir)
        
        shapefiles = glob.glob(glob_pattern)
    
    if len(shapefiles) != 1:
        raise ValueError('%s shapefiles were found, expecting just one to match %s' % (len(shapefiles), glob_pattern))
    
    return shapefiles[0]


def mpl_axes_plot(axes, geometries):
    """Plot lines on the given axes, given the geometries."""
    # TODO: This interface should be exposed nicely on the geoaxes itself.
    import matplotlib.collections as mcollections
    import cartopy.mpl_integration.patch as patch
    
    paths = []
    for geom in geometries:            
        paths.extend(patch.geos_to_path(axes.projection.project_geometry(geom)))            
    axes.add_collection(mcollections.PathCollection(paths, facecolor='none'), autolim=False)
     

if __name__ == '__main__':
    coastlines = natural_earth(resolution='110m', category='physical', name='coastline')
    for record in Reader(coastlines).records():
        print record.attributes
        
    
    # XXX TODO: Turn into a tutorial
    coastlines = natural_earth(resolution='110m', category='cultural', name='admin-0-countries')
    cntry_size = [(record.attributes['NAME'], int(record.attributes['POP_EST'])) for record in Reader(coastlines).records()]
    
    # return the countries, grouped alphabetically, sorted by size.
    import itertools
    cntry_size.sort(key=lambda (name, population): (name[0], population))
    for k, g in itertools.groupby(cntry_size, key=lambda item: item[0][0]):
        print k, list(g)