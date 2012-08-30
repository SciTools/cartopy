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


import unittest

from shapely import geometry

import cartopy.crs as ccrs


class TestBoundary(unittest.TestCase):
    def test_cuts(self):
        # Check that fragments do not start or end with one of the original ... ?
        linear_ring = geometry.polygon.LinearRing([(-10, 30), (10, 60), (10, 50)])
        projection = ccrs.Robinson(170.5)
        multi_line_string = projection.project_geometry(linear_ring)
        from cartopy.tests import show
        #show(projection, multi_line_string)

        # The original ring should have been split into multiple pieces.
        self.assertGreater(len(multi_line_string), 1)

        def assert_intersection_with_boundary(segment_coords):
            # Double the length of the segment.
            start = segment_coords[0]
            end = segment_coords[1]
            end = [end[i] + 2 * (end[i] - start[i]) for i in (0, 1)]
            extended_segment = geometry.LineString([start, end])
            # And see if it crosses the boundary.
            intersection = extended_segment.intersection(projection.boundary)
            self.assertFalse(intersection.is_empty, 'Bad topology near boundary')

        # Each line resulting from the split should start and end with a
        # segment that crosses the boundary when extended to double length.
        # (This is important when considering polygon rings which need to be
        # attached to the boundary.)
        for line_string in multi_line_string:
            coords = list(line_string.coords)
            self.assertGreaterEqual(len(coords), 2)
            assert_intersection_with_boundary(coords[1::-1])
            assert_intersection_with_boundary(coords[-2:])

    def test_out_of_bounds(self):
        # Check that a ring that is completely out of the map boundary
        # produces an empty result.
        # XXX Check efficiency?
        projection = ccrs.TransverseMercator(central_longitude=0)

        rings = [
            # All valid
            [(86, -1), (86, 1), (88, 1), (88, -1)],
            # One NaN
            [(86, -1), (86, 1), (130, 1), (88, -1)],
            # A NaN segment
            [(86, -1), (86, 1), (130, 1), (130, -1)],
            # All NaN
            [(120, -1), (120, 1), (130, 1), (130, -1)],
        ]

        # Try all four combinations of valid/NaN vs valid/NaN.
        for coords in rings:
            linear_ring = geometry.polygon.LinearRing(coords)
            multi_line_string = projection.project_geometry(linear_ring)
            self.assertEqual(len(multi_line_string), 0)


class TestMisc(unittest.TestCase):
    def test_misc(self):
        projection = ccrs.TransverseMercator(central_longitude=-90)
        linear_ring = geometry.polygon.LinearRing([(-10, 30), (10, 60), (10, 50)])
        multi_line_string = projection.project_geometry(linear_ring)
        from cartopy.tests import show
        #show(projection, multi_line_string)

    def test_small(self):
        # What happens when a small (i.e. < threshold) feature crosses the
        # boundary?
        projection = ccrs.Mercator()
        linear_ring = geometry.polygon.LinearRing([
            (-179.9173693847652942, -16.5017831356493616),
            (-180.0000000000000000, -16.0671326636424396),
            (-179.7933201090486079, -16.0208822567412312),
            ])
        multi_line_string = projection.project_geometry(linear_ring)
        from cartopy.tests import show
        #show(projection, multi_line_string)
        self.assertEqual(len(multi_line_string), 1)
        #self.assertGreater(len(multi_line_string[0].coords), 2)


if __name__ == '__main__':
    unittest.main()
