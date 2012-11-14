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
Combines the shapefile access of pyshp with the
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

from cartopy.io import DownloadableItem
from cartopy import config


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


def natural_earth(resolution='110m', category='physical', name='coastline'):
    """
    Returns the path to the requested natural earth shapefile, 
    downloading and unziping if necessary.
    
    """
    # get hold of the DownloadableItem (typically a NEShpDownloader instance)
    # which we can then simply call its path method to get the appropriate
    # shapefile (it will download if necessary)
    ne_downloader = DownloadableItem.from_config(('shapefiles', 
                                                  'natural_earth',
                                                  resolution, category, name))
    format_dict = {'config': config, # XXX consider not passing this through
                   'category': category, 'name': name, 'resolution': resolution}
    return ne_downloader.path(format_dict)
    
    
# XXX NEW STUFF.... 





class NEShpDownloader(DownloadableItem):
    """
    Specialises :class:`cartopy.io.DownloadableItem` to download the zipped
    Natural Earth shapefiles and extract them to the defined location 
    (typically user configurable). 
    
    The keys which should be passed through when using the ``format_dict``
    are typically ``category``, ``resolution`` and ``name``.    
     
    """
    FORMAT_KEYS = ('config', 'resolution', 'category', 'resolution', 'name')    
    
    # define the NaturalEarth url template (note the repeat of http is 
    # intentional and part of the NE [Natural Earth] scheme).
    _NE_URL_TEMPLATE = ('http://www.naturalearthdata.com/'
                        'http//www.naturalearthdata.com/download/'
                        '{resolution}/{category}/ne_{resolution}_{name}.zip')
    
    def __init__(self, 
                 url_template=_NE_URL_TEMPLATE, 
                 target_path_template=None, 
                 pre_downloaded_path_template='',
                 ):
        # adds some NE defaults to the __init__ of a DownloadableItem
        DownloadableItem.__init__(self, url_template, 
                                             target_path_template, 
                                             pre_downloaded_path_template)
    
    def zip_file_contents(self, format_dict):
        """
        Returns a generator of the filenames to be found in the downloaded 
        natural earth zip file.
        
        """
        for ext in ['.shp', '.dbf', '.shx']:
            yield ('ne_{resolution}_{name}'
                   '{extension}'.format(extension=ext, **format_dict))
        
    def acquire_resource(self, target_path, format_dict):
        """
        Downloads the zip file and extracts the files listed in
        :meth:`zip_file_contents` to the target path.
        
        """
        import cStringIO as StringIO
        from zipfile import ZipFile
        
        target_dir = os.path.dirname(target_path)
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        
        url = self.url(format_dict)
        
        shapefile_online = self._urlopen(url)
        zfh = ZipFile(StringIO.StringIO(shapefile_online.read()), 'r')
        
        for member_path in self.zip_file_contents(format_dict):
            ext = os.path.splitext(member_path)[1]
            target = os.path.splitext(target_path)[0] + ext
            member = zfh.getinfo(member_path)
            with open(target, 'wb') as fh:
                fh.write(zfh.open(member).read())

        shapefile_online.close()
        zfh.close()
        
        return target_path

    @staticmethod
    def default_downloader():
        """
        Returns a generic, standard, NEShpDownloader instance. 
        
        Typically, a user will not need to call this staticmethod.
        
        To find the path template of the NEShpDownloader:
        
            >>> ne_dnldr = NEShpDownloader.default_downloader()
            >>> print ne_dnldr.target_path_template
            hello world
        
        """ 
        ne_path_template = os.path.join('{config[data_dir]}', 'shapefiles',
                                         'natural_earth', '{category}', 
                                         '{resolution}_{name}.shp')
        return NEShpDownloader(target_path_template=ne_path_template)


# add a generic Natural Earth shapefile downloader to the config dictionary's 
# 'downloads' section.
config['downloads'].setdefault(('shapefiles', 'natural_earth'),
                               NEShpDownloader.default_downloader())


# XXX cartopy's shapefiles are out of date and the new ones cause problems. Temporarily
# use the download mechanism to point to the old files::
config['downloads'][('shapefiles', 'natural_earth')
                    ].target_path_template = os.path.join('{config[data_dir]}', 
                                             'shapefiles',
                                             'natural_earth', 
                                             '{resolution}-{name}', 
                                             '{resolution}_{name}.shp')