# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import os.path

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest
import shapely.geometry as sgeom

import cartopy.io.shapereader as shp


class TestLakes:
    def setup_class(self):
        LAKES_PATH = os.path.join(os.path.dirname(__file__),
                                  'lakes_shapefile', 'ne_110m_lakes.shp')
        self.reader = shp.Reader(LAKES_PATH)
        names = [record.attributes['name'] for record in self.reader.records()]
        # Choose a nice small lake
        self.lake_name = 'Lago de\rNicaragua'
        self.lake_index = names.index(self.lake_name)
        self.test_lake_geometry = \
            list(self.reader.geometries())[self.lake_index]
        self.test_lake_record = list(self.reader.records())[self.lake_index]

    def test_geometry(self):
        lake_geometry = self.test_lake_geometry
        assert lake_geometry.type == 'Polygon'

        # force an orientation due to potential reader differences
        # with pyshp 2.2.0 forcing a specific orientation.
        polygon = sgeom.polygon.orient(lake_geometry, -1)

        expected = np.array([(-84.85548682324658, 11.147898667846633),
                             (-85.29013729525353, 11.176165676310276),
                             (-85.79132117383625, 11.509737046754324),
                             (-85.8851655748783, 11.900100816287136),
                             (-85.5653401354239, 11.940330918826362),
                             (-85.03684526237491, 11.5216484643976),
                             (-84.85548682324658, 11.147898667846633),
                             (-84.85548682324658, 11.147898667846633)])

        assert_array_almost_equal(expected, polygon.exterior.coords)

        assert len(polygon.interiors) == 0

    def test_record(self):
        lake_record = self.test_lake_record
        assert lake_record.attributes.get('name') == self.lake_name
        expected = sorted(['admin', 'featurecla', 'min_label', 'min_zoom',
                           'name', 'name_alt', 'scalerank'])
        actual = sorted(lake_record.attributes.keys())
        assert actual == expected
        assert lake_record.geometry == self.test_lake_geometry

    @pytest.mark.skipif(shp._HAS_FIONA,
                        reason="Fiona reader doesn't support lazy loading.")
    def test_bounds(self):
        # tests that a file which has a record with a bbox can
        # use the bbox without first creating the geometry
        record = next(self.reader.records())
        assert not record._geometry, \
            'The geometry was loaded before it was needed.'
        assert len(record._bounds) == 4
        assert record._bounds == record.bounds
        assert not record._geometry, \
            'The geometry was loaded in order to create the bounds.'


@pytest.mark.filterwarnings("ignore:Downloading")
@pytest.mark.natural_earth
class TestRivers:
    def setup_class(self):
        RIVERS_PATH = shp.natural_earth(resolution='110m',
                                        category='physical',
                                        name='rivers_lake_centerlines')
        self.reader = shp.Reader(RIVERS_PATH)
        names = [record.attributes['name'] for record in self.reader.records()]
        # Choose a nice small river
        self.river_name = 'Peace'
        self.river_index = names.index(self.river_name)
        self.test_river_geometry = \
            list(self.reader.geometries())[self.river_index]
        self.test_river_record = list(self.reader.records())[self.river_index]

    def test_geometry(self):
        geometry = self.test_river_geometry
        assert geometry.type == 'LineString'

        linestring = geometry
        coords = linestring.coords
        assert round(abs(coords[0][0] - -124.83563045947423), 7) == 0
        assert round(abs(coords[0][1] - 56.75692352968272), 7) == 0
        assert round(abs(coords[1][0] - -124.20045039940291), 7) == 0
        assert round(abs(coords[1][1] - 56.243492336646824), 7) == 0

    def test_record(self):
        records = list(self.reader.records())
        assert len(records) == len(self.reader)

        # Choose a nice small lake
        river_record = records[self.river_index]
        expected_attributes = {'featurecla': 'River',
                               'min_label': 3.1,
                               'min_zoom': 2.1,
                               'name': self.river_name,
                               'name_en': self.river_name,
                               'scalerank': 2}
        for key, value in river_record.attributes.items():
            if key in expected_attributes:
                assert value == expected_attributes[key]
        assert river_record.geometry == self.test_river_geometry
