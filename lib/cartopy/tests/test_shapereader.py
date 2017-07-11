# (C) British Crown Copyright 2011 - 2017, Met Office
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

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_array_almost_equal

import cartopy.io.shapereader as shp


LAKES_PATH = shp.natural_earth(resolution='110m',
                               category='physical',
                               name='lakes')
RIVERS_PATH = shp.natural_earth(resolution='110m',
                                category='physical',
                                name='rivers_lake_centerlines')


class TestLakes(object):
    def setup_class(self):
        self.reader = shp.Reader(LAKES_PATH)

    def _assert_geometry(self, geometry):
        assert geometry.type == 'MultiPolygon'
        assert len(geometry) == 1

        polygon = geometry[0]

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

    def test_geometry(self):
        geometries = list(self.reader.geometries())
        assert len(geometries) == len(self.reader)

        # Choose a nice small lake
        lake = geometries[14]
        self._assert_geometry(lake)

    def test_record(self):
        records = list(self.reader.records())
        assert len(records) == len(self.reader)

        # Choose a nice small lake
        lake_record = records[14]
        assert lake_record.attributes.get('name') == 'Lago de\rNicaragua'
        assert sorted(lake_record.attributes.keys()) == \
            sorted(['admin', 'featurecla', 'scalerank', 'name_alt', 'name'])
        lake = lake_record.geometry
        self._assert_geometry(lake)

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


class TestRivers(object):
    def setup_class(self):
        self.reader = shp.Reader(RIVERS_PATH)

    def _assert_geometry(self, geometry):
        assert geometry.type == 'MultiLineString'
        assert len(geometry) == 1

        linestring = geometry[0]
        coords = linestring.coords
        assert round(abs(coords[0][0] - -113.823382738076), 7) == 0
        assert round(abs(coords[0][1] - 58.7102151556671), 7) == 0
        assert round(abs(coords[1][0] - -113.71351864302348), 7) == 0
        assert round(abs(coords[1][1] - 58.669261583075794), 7) == 0

    def test_geometry(self):
        geometries = list(self.reader.geometries())
        assert len(geometries) == len(self.reader)

        # Choose a nice small river
        river = geometries[6]
        self._assert_geometry(river)

    def test_record(self):
        records = list(self.reader.records())
        assert len(records) == len(self.reader)

        # Choose a nice small lake
        river_record = records[6]
        expected_attributes = {'featurecla': 'River',
                               'scalerank': 2,
                               'name': 'Peace'}
        for key, value in river_record.attributes.items():
            if key == 'name_alt':
                # This value changed between pyshp 1.2.10 and 1.2.11, test it
                # as a special case, it should be an empty string once the
                # leading/trailing space is removed:
                assert not len(value.strip())
            else:
                assert value == expected_attributes[key]
        river = river_record.geometry
        self._assert_geometry(river)
