import os.path
import unittest

from shapely.geometry import MultiPolygon, Polygon

import cartopy.io.shapereader as shp

LAKES_PATH = shp.natural_earth(resolution='110m', category='physical', name='lakes')
RIVERS_PATH = shp.natural_earth(resolution='110m', category='physical', name='rivers-lake-centerlines')


class TestLakes(unittest.TestCase):
    def setUp(self):
        self.reader = shp.Reader(LAKES_PATH)

    def _assert_geometry(self, geometry):
        self.assertEqual(geometry.type, 'MultiPolygon')
        self.assertEqual(len(geometry), 1)

        polygon = geometry[0]
        coords = polygon.exterior.coords
        self.assertAlmostEqual(coords[0][0], -80.706437751)
        self.assertAlmostEqual(coords[0][1], 26.7889594589)
        self.assertAlmostEqual(coords[1][0], -80.9324446276)
        self.assertAlmostEqual(coords[1][1], 26.82327261)
        self.assertAlmostEqual(coords[2][0], -80.919706387)
        self.assertAlmostEqual(coords[2][1], 27.0689165309)
        self.assertAlmostEqual(coords[3][0], -80.6936995104)
        self.assertAlmostEqual(coords[3][1], 27.034629218)
        self.assertAlmostEqual(coords[4][0], -80.706437751)
        self.assertAlmostEqual(coords[4][1], 26.7889594589)
        self.assertAlmostEqual(coords[5][0], -80.706437751)
        self.assertAlmostEqual(coords[5][1], 26.7889594589)

        self.assertEqual(len(polygon.interiors), 0)

    def test_geometry(self):
        geometries = list(self.reader.geometries())
        self.assertEqual(len(geometries), len(self.reader))

        # Choose a nice small lake
        lake = geometries[14]
        self._assert_geometry(lake)

    def test_record(self):
        records = list(self.reader.records())
        self.assertEqual(len(records), len(self.reader))

        # Choose a nice small lake
        lake_record = records[14]
        self.assertEqual(
            lake_record.attributes,
            {'FeatureCla': 'Lake', 'Name2': ' ' * 254, 'ScaleRank': 1, 'Name1': 'Lake Okeechobee'})
        lake = lake_record.geometry
        self._assert_geometry(lake)


class TestRivers(unittest.TestCase):
    def setUp(self):
        self.reader = shp.Reader(RIVERS_PATH)

    def _assert_geometry(self, geometry):
        self.assertEqual(geometry.type, 'MultiLineString')
        self.assertEqual(len(geometry), 1)

        linestring = geometry[0]
        coords = linestring.coords
        self.assertAlmostEqual(coords[0][0], -113.823382738076)
        self.assertAlmostEqual(coords[0][1], 58.7102151556671)
        self.assertAlmostEqual(coords[1][0], -113.71351864302348)
        self.assertAlmostEqual(coords[1][1], 58.669261583075794)

    def test_geometry(self):
        geometries = list(self.reader.geometries())
        self.assertEqual(len(geometries), len(self.reader))

        # Choose a nice small river
        river = geometries[6]
        self._assert_geometry(river)

    def test_record(self):
        records = list(self.reader.records())
        self.assertEqual(len(records), len(self.reader))

        # Choose a nice small lake
        river_record = records[6]
        self.assertEqual(
            river_record.attributes,
            {'FeatureCla': 'River', 'Name2': ' ' * 254, 'ScaleRank': 2, 'Name1': 'Peace'})
        river = river_record.geometry
        self._assert_geometry(river)


if __name__ == '__main__':
    unittest.main()
