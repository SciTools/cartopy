# (C) British Crown Copyright 2011 - 2016, Met Office
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

import unittest

import numpy as np
import shapely.geometry as sgeom
import shapely.wkt


import cartopy.crs as ccrs


class TestBoundary(unittest.TestCase):
    def test_no_polygon_boundary_reversal(self):
        # Check that polygons preserve their clockwise or counter-clockwise
        # ordering when they are attached to the boundary.
        # Failure to do so will result in invalid polygons (their boundaries
        # cross-over).
        polygon = sgeom.Polygon([(-10, 30), (10, 60), (10, 50)])
        projection = ccrs.Robinson(170.5)
        multi_polygon = projection.project_geometry(polygon)
        for polygon in multi_polygon:
            self.assertTrue(polygon.is_valid)

    def test_polygon_boundary_attachment(self):
        # Check the polygon is attached to the boundary even when no
        # intermediate point for one of the crossing segments would normally
        # exist.
        polygon = sgeom.Polygon([(-10, 30), (10, 60), (10, 50)])
        projection = ccrs.Robinson(170.6)
        # This will raise an exception if the polygon/boundary intersection
        # fails.
        multi_polygon = projection.project_geometry(polygon)

    def test_out_of_bounds(self):
        # Check that a polygon that is completely out of the map boundary
        # doesn't produce an empty result.
        projection = ccrs.TransverseMercator(central_longitude=0)

        polys = [
            # All valid
            ([(86, -1), (86, 1), (88, 1), (88, -1)], 1),
            # One out of backwards projection range
            ([(86, -1), (86, 1), (130, 1), (88, -1)], 1),
            # An out of backwards projection range segment
            ([(86, -1), (86, 1), (130, 1), (130, -1)], 1),
            # All out of backwards projection range
            ([(120, -1), (120, 1), (130, 1), (130, -1)], 0),
        ]

        # Try all four combinations of valid/NaN vs valid/NaN.
        for coords, expected_polys in polys:
            polygon = sgeom.Polygon(coords)
            multi_polygon = projection.project_geometry(polygon)
            self.assertEqual(len(multi_polygon), expected_polys)


class TestMisc(unittest.TestCase):
    def test_misc(self):
        projection = ccrs.TransverseMercator(central_longitude=-90)
        polygon = sgeom.Polygon([(-10, 30), (10, 60), (10, 50)])
        multi_polygon = projection.project_geometry(polygon)

    def test_small(self):
        projection = ccrs.Mercator()
        polygon = sgeom.Polygon([
            (-179.7933201090486079, -16.0208822567412312),
            (-180.0000000000000000, -16.0671326636424396),
            (-179.9173693847652942, -16.5017831356493616),
        ])
        multi_polygon = projection.project_geometry(polygon)
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].exterior.coords), 4)

    def test_former_infloop_case(self):
        # test a polygon which used to get stuck in an infinite loop
        # see https://github.com/SciTools/cartopy/issues/60
        coords = [(260.625, 68.90383337092122), (360.0, 79.8556091996901),
                  (360.0, 77.76848175458498), (0.0, 88.79068047337279),
                  (210.0, 90.0), (135.0, 88.79068047337279),
                  (260.625, 68.90383337092122)]
        geom = sgeom.Polygon(coords)

        target_projection = ccrs.PlateCarree()
        source_crs = ccrs.Geodetic()

        multi_polygon = target_projection.project_geometry(geom, source_crs)
        # check the result is non-empty
        self.assertFalse(multi_polygon.is_empty)

    def test_project_previous_infinite_loop(self):
        mstring1 = shapely.wkt.loads(
            'MULTILINESTRING ('
            '(-179.9999990464349651 -80.2000000000000171, '
            '-179.5000000001111005 -80.2000000000000171, '
            '-179.5000000001111005 -79.9000000000000199, '
            '-179.9999995232739138 -79.9499999523163041, '
            '-179.8000000001110550 -80.0000000000000000, '
            '-179.8000000001110550 -80.0999999999999943, '
            '-179.9999999047436177 -80.0999999999999943), '
            '(179.9999995231628702 -79.9499999523163041, '
            '179.5000000000000000 -79.9000000000000199, '
            '179.5000000000000000 -80.0000000000000000, '
            '179.9999995231628702 -80.0499999523162842, '
            '179.5000000000000000 -80.0999999999999943, '
            '179.5000000000000000 -80.2000000000000171, '
            '179.9999990463256836 -80.2000000000000171))')
        mstring2 = shapely.wkt.loads(
            'MULTILINESTRING ('
            '(179.9999996185302678 -79.9999999904632659, '
            '179.5999999999999943 -79.9899999999999949, '
            '179.5999999999999943 -79.9399999999999977, '
            '179.9999996185302678 -79.9599999809265114), '
            '(-179.9999999047436177 -79.9600000000000080, '
            '-179.9000000001110777 -79.9600000000000080, '
            '-179.9000000001110777 -80.0000000000000000, '
            '-179.9999999047436177 -80.0000000000000000))')
        multi_line_strings = [mstring1, mstring2]

        src = ccrs.PlateCarree()
        src._attach_lines_to_boundary(multi_line_strings, True)

    def test_3pt_poly(self):
        projection = ccrs.OSGB()
        polygon = sgeom.Polygon([(-1000, -1000),
                                 (-1000, 200000),
                                 (200000, -1000)])
        multi_polygon = projection.project_geometry(polygon, ccrs.OSGB())
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].exterior.coords), 4)

    def test_self_intersecting_1(self):
        # Geometry comes from a matplotlib contourf (see #537)
        wkt = ('POLYGON ((366.22000122 -9.71489298, '
               '366.73212393 -9.679999349999999, '
               '366.77412634 -8.767753000000001, '
               '366.17762962 -9.679999349999999, '
               '366.22000122 -9.71489298), '
               '(366.22000122 -9.692636309999999, '
               '366.32998657 -9.603356099999999, '
               '366.74765799 -9.019999500000001, '
               '366.5094086 -9.63175386, '
               '366.22000122 -9.692636309999999))')
        geom = shapely.wkt.loads(wkt)
        source, target = ccrs.RotatedPole(198.0, 39.25), ccrs.EuroPP()
        projected = target.project_geometry(geom, source)
        # Before handling self intersecting interiors, the area would be
        # approximately 13262233761329.
        area = projected.area
        self.assertTrue(2.2e9 < area < 2.3e9,
                        msg='Got area {}, expecting ~2.2e9'.format(area))

    def test_self_intersecting_2(self):
        # Geometry comes from a matplotlib contourf (see #509)
        wkt = ('POLYGON ((343 20, 345 23, 342 25, 343 22, '
               '340 25, 341 25, 340 25, 343 20), (343 21, '
               '343 22, 344 23, 343 21))')
        geom = shapely.wkt.loads(wkt)
        source = target = ccrs.RotatedPole(193.0, 41.0)
        projected = target.project_geometry(geom, source)
        # Before handling self intersecting interiors, the area would be
        # approximately 64808.
        self.assertTrue(7.9 < projected.area < 8.1)

    def test_tiny_point_between_boundary_points(self):
        # Geometry comes from #259.
        target = ccrs.Orthographic(0, -75)
        source = ccrs.PlateCarree()
        wkt = 'POLYGON ((132 -40, 133 -6, 125.3 1, 115 -6, 132 -40))'
        geom = shapely.wkt.loads(wkt)

        target = ccrs.Orthographic(central_latitude=90., central_longitude=0)
        source = ccrs.PlateCarree()
        projected = target.project_geometry(geom, source)
        area = projected.area
        # Before fixing, this geometry used to fill the whole disk. Approx
        # 1.2e14.
        self.assertTrue(81330 < area < 81340,
                        msg='Got area {}, expecting ~81336'.format(area))


class TestQuality(unittest.TestCase):
    def setUp(self):
        projection = ccrs.RotatedPole(pole_longitude=177.5,
                                      pole_latitude=37.5)
        polygon = sgeom.Polygon([
            (177.5, -57.38460319),
            (180.0, -57.445077),
            (175.0, -57.19913331),
        ])
        self.multi_polygon = projection.project_geometry(polygon)
        # from cartopy.tests.mpl import show
        # show(projection, self.multi_polygon)

    def test_split(self):
        # Start simple ... there should be two projected polygons.
        self.assertEqual(len(self.multi_polygon), 2)

    def test_repeats(self):
        # Make sure we don't have repeated points at the boundary, because
        # they mess up the linear extrapolation to the boundary.

        # Make sure there aren't any repeated points.
        xy = np.array(self.multi_polygon[0].exterior.coords)
        same = (xy[1:] == xy[:-1]).all(axis=1)
        self.assertFalse(any(same), 'Repeated points in projected geometry.')

    def test_symmetry(self):
        # Make sure the number of points added on the way towards the
        # boundary is similar to the number of points added on the way away
        # from the boundary.

        # Identify all the contiguous sets of non-boundary points.
        xy = np.array(self.multi_polygon[0].exterior.coords)
        boundary = np.logical_or(xy[:, 1] == 90, xy[:, 1] == -90)
        regions = (boundary[1:] != boundary[:-1]).cumsum()
        regions = np.insert(regions, 0, 0)

        # For each region, check if the number of increasing steps is roughly
        # equal to the number of decreasing steps.
        for i in range(boundary[0], regions.max(), 2):
            indices = np.where(regions == i)
            x = xy[indices, 0]
            delta = np.diff(x)
            num_incr = np.count_nonzero(delta > 0)
            num_decr = np.count_nonzero(delta < 0)
            self.assertLess(abs(num_incr - num_decr), 3,
                            'Too much asymmetry.')


class PolygonTests(unittest.TestCase):
    def _assert_bounds(self, bounds, x1, y1, x2, y2, delta=1):
        self.assertAlmostEqual(bounds[0], x1, delta=delta)
        self.assertAlmostEqual(bounds[1], y1, delta=delta)
        self.assertAlmostEqual(bounds[2], x2, delta=delta)
        self.assertAlmostEqual(bounds[3], y2, delta=delta)


class TestWrap(PolygonTests):
    # Test that Plate Carree projection "does the right thing"(tm) with
    # source data tha extends outside the [-180, 180] range.
    def test_plate_carree_no_wrap(self):
        proj = ccrs.PlateCarree()
        poly = sgeom.box(0, 0, 10, 10)
        multi_polygon = proj.project_geometry(poly, proj)
        # Check the structure
        self.assertEqual(len(multi_polygon), 1)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, 0, 0, 10, 10)

    def test_plate_carree_partial_wrap(self):
        proj = ccrs.PlateCarree()
        poly = sgeom.box(170, 0, 190, 10)
        multi_polygon = proj.project_geometry(poly, proj)
        # Check the structure
        self.assertEqual(len(multi_polygon), 2)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, 170, 0, 180, 10)
        polygon = multi_polygon[1]
        self._assert_bounds(polygon.bounds, -180, 0, -170, 10)

    def test_plate_carree_wrap(self):
        proj = ccrs.PlateCarree()
        poly = sgeom.box(200, 0, 220, 10)
        multi_polygon = proj.project_geometry(poly, proj)
        # Check the structure
        self.assertEqual(len(multi_polygon), 1)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, -160, 0, -140, 10)


def ring(minx, miny, maxx, maxy, ccw):
    box = sgeom.box(minx, miny, maxx, maxy, ccw)
    return np.array(box.exterior.coords)


class TestHoles(PolygonTests):
    def test_simple(self):
        proj = ccrs.PlateCarree()
        poly = sgeom.Polygon(ring(-40, -40, 40, 40, True),
                             [ring(-20, -20, 20, 20, False)])
        multi_polygon = proj.project_geometry(poly)
        # Check the structure
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].interiors), 1)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, -40, -47, 40, 47)
        self._assert_bounds(polygon.interiors[0].bounds, -20, -21, 20, 21)

    def test_wrapped_poly_simple_hole(self):
        proj = ccrs.PlateCarree(-150)
        poly = sgeom.Polygon(ring(-40, -40, 40, 40, True),
                             [ring(-20, -20, 20, 20, False)])
        multi_polygon = proj.project_geometry(poly)
        # Check the structure
        self.assertEqual(len(multi_polygon), 2)
        self.assertEqual(len(multi_polygon[0].interiors), 1)
        self.assertEqual(len(multi_polygon[1].interiors), 0)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, 110, -47, 180, 47)
        self._assert_bounds(polygon.interiors[0].bounds, 130, -21, 170, 21)
        polygon = multi_polygon[1]
        self._assert_bounds(polygon.bounds, -180, -43, -170, 43)

    def test_wrapped_poly_wrapped_hole(self):
        proj = ccrs.PlateCarree(-180)
        poly = sgeom.Polygon(ring(-40, -40, 40, 40, True),
                             [ring(-20, -20, 20, 20, False)])
        multi_polygon = proj.project_geometry(poly)
        # Check the structure
        self.assertEqual(len(multi_polygon), 2)
        self.assertEqual(len(multi_polygon[0].interiors), 0)
        self.assertEqual(len(multi_polygon[1].interiors), 0)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, 140, -47, 180, 47)
        polygon = multi_polygon[1]
        self._assert_bounds(polygon.bounds, -180, -47, -140, 47)

    def test_inverted_poly_simple_hole(self):
        proj = ccrs.NorthPolarStereo()
        poly = sgeom.Polygon([(0, 0), (-90, 0), (-180, 0), (-270, 0)],
                             [[(0, -30), (90, -30), (180, -30), (270, -30)]])
        multi_polygon = proj.project_geometry(poly)
        # Check the structure
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].interiors), 1)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, -2.4e7, -2.4e7, 2.4e7, 2.4e7, 1e6)
        self._assert_bounds(polygon.interiors[0].bounds,
                            - 1.2e7, -1.2e7, 1.2e7, 1.2e7, 1e6)

    def test_inverted_poly_clipped_hole(self):
        proj = ccrs.NorthPolarStereo()
        poly = sgeom.Polygon([(0, 0), (-90, 0), (-180, 0), (-270, 0)],
                             [[(-135, -60), (-45, -60),
                               (45, -60), (135, -60)]])
        multi_polygon = proj.project_geometry(poly)
        # Check the structure
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].interiors), 1)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, -5.0e7, -5.0e7, 5.0e7, 5.0e7, 1e6)
        self._assert_bounds(polygon.interiors[0].bounds,
                            - 1.2e7, -1.2e7, 1.2e7, 1.2e7, 1e6)
        self.assertAlmostEqual(polygon.area, 7.30e15, delta=1e13)

    def test_inverted_poly_removed_hole(self):
        proj = ccrs.NorthPolarStereo(globe=ccrs.Globe(ellipse='WGS84'))
        poly = sgeom.Polygon([(0, 0), (-90, 0), (-180, 0), (-270, 0)],
                             [[(-135, -75), (-45, -75),
                               (45, -75), (135, -75)]])
        multi_polygon = proj.project_geometry(poly)
        # Check the structure
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].interiors), 1)
        # Check the rough shape
        polygon = multi_polygon[0]
        self._assert_bounds(polygon.bounds, -5.0e7, -5.0e7, 5.0e7, 5.0e7, 1e6)
        self._assert_bounds(polygon.interiors[0].bounds,
                            - 1.2e7, -1.2e7, 1.2e7, 1.2e7, 1e6)
        self.assertAlmostEqual(polygon.area, 7.34e15, delta=1e13)

    def test_multiple_interiors(self):
        exterior = ring(0, 0, 12, 12, True)
        interiors = [ring(1, 1, 2, 2, False), ring(1, 8, 2, 9, False)]

        poly = sgeom.Polygon(exterior, interiors)

        target = ccrs.PlateCarree()
        source = ccrs.Geodetic()

        assert len(list(target.project_geometry(poly, source))) == 1


if __name__ == '__main__':
    unittest.main()
