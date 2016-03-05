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

import shapely.geometry as sgeom
import numpy as np

import cartopy.crs as ccrs


class TestBoundary(unittest.TestCase):
    def test_cuts(self):
        # Check that fragments do not start or end with one of the
        # original ... ?
        linear_ring = sgeom.LinearRing([(-10, 30), (10, 60), (10, 50)])
        projection = ccrs.Robinson(170.5)
        rings, multi_line_string = projection.project_geometry(linear_ring)

        # The original ring should have been split into multiple pieces.
        self.assertGreater(len(multi_line_string), 1)
        self.assertFalse(rings)

        def assert_intersection_with_boundary(segment_coords):
            # Double the length of the segment.
            start = segment_coords[0]
            end = segment_coords[1]
            end = [end[i] + 2 * (end[i] - start[i]) for i in (0, 1)]
            extended_segment = sgeom.LineString([start, end])
            # And see if it crosses the boundary.
            intersection = extended_segment.intersection(projection.boundary)
            self.assertFalse(intersection.is_empty,
                             'Bad topology near boundary')

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
            ([(86, 1), (86, -1), (88, -1), (88, 1)], -1),
            # One NaN
            ([(86, 1), (86, -1), (130, -1), (88, 1)], 1),
            # A NaN segment
            ([(86, 1), (86, -1), (130, -1), (130, 1)], 1),
            # All NaN
            ([(120, 1), (120, -1), (130, -1), (130, 1)], 0),
        ]

        # Try all four combinations of valid/NaN vs valid/NaN.
        for coords, expected_n_lines in rings:
            linear_ring = sgeom.LinearRing(coords)
            rings, mlinestr = projection.project_geometry(linear_ring)
            if expected_n_lines == -1:
                self.assertTrue(rings)
                self.assertFalse(mlinestr)
            else:
                self.assertEqual(len(mlinestr), expected_n_lines)
                if expected_n_lines == 0:
                    self.assertTrue(mlinestr.is_empty)


class TestMisc(unittest.TestCase):
    def test_small(self):
        # What happens when a small (i.e. < threshold) feature crosses the
        # boundary?
        projection = ccrs.Mercator()
        linear_ring = sgeom.LinearRing([
            (-179.9173693847652942, -16.5017831356493616),
            (-180.0000000000000000, -16.0671326636424396),
            (-179.7933201090486079, -16.0208822567412312),
        ])
        rings, multi_line_string = projection.project_geometry(linear_ring)
        # There should be one, and only one, returned ring.
        self.assertIsInstance(multi_line_string, sgeom.MultiLineString)
        self.assertEqual(len(multi_line_string), 0)
        self.assertEqual(len(rings), 1)

        # from cartopy.tests.mpl import show
        # show(projection, multi_line_string)

    def test_three_points(self):
        # The following LinearRing when projected from PlateCarree() to
        # PlateCarree(180.0) results in three points all in close proximity.
        # If an attempt is made to form a LinearRing from the three points
        # by combining the first and last an exception will be raised.
        # Check that this object can be projected without error.
        coords = [(0.0, -45.0),
                  (0.0, -44.99974961593933),
                  (0.000727869825138, -45.0),
                  (0.0, -45.000105851567454),
                  (0.0, -45.0)]
        linear_ring = sgeom.LinearRing(coords)
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.PlateCarree(180.0)
        try:
            _ = target_proj.project_geometry(linear_ring, src_proj)
        except ValueError:
            self.fail("Failed to project LinearRing.")

    def test_stitch(self):
        # The following LinearRing wanders in/out of the map domain
        # but importantly the "vertical" lines at 0'E and 360'E are both
        # chopped by the map boundary. This results in their ends being
        # *very* close to each other and confusion over which occurs
        # first when navigating around the boundary.
        # Check that these ends are stitched together to avoid the
        # boundary ordering ambiguity.
        # NB. This kind of polygon often occurs with MPL's contouring.
        coords = [(0.0, -70.70499926182919),
                  (0.0, -71.25),
                  (0.0, -72.5),
                  (0.0, -73.49076371657017),
                  (360.0, -73.49076371657017),
                  (360.0, -72.5),
                  (360.0, -71.25),
                  (360.0, -70.70499926182919),
                  (350, -73),
                  (10, -73)]
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(80)

        linear_ring = sgeom.LinearRing(coords)
        rings, mlinestr = target_proj.project_geometry(linear_ring, src_proj)
        self.assertEqual(len(mlinestr), 1)
        self.assertEqual(len(rings), 0)

        # Check the stitch works in either direction.
        linear_ring = sgeom.LinearRing(coords[::-1])
        rings, mlinestr = target_proj.project_geometry(linear_ring, src_proj)
        self.assertEqual(len(mlinestr), 1)
        self.assertEqual(len(rings), 0)

    def test_at_boundary(self):
        # Check that a polygon is split and recombined correctly
        # as a result of being on the boundary, determined by tolerance.

        exterior = np.array(
            [[177.5, -79.912],
             [178.333, -79.946],
             [181.666, -83.494],
             [180.833, -83.570],
             [180., -83.620],
             [178.438, -83.333],
             [178.333, -83.312],
             [177.956, -83.888],
             [180., -84.086],
             [180.833, -84.318],
             [183., -86.],
             [183., -78.],
             [177.5, -79.912]])
        tring = sgeom.LinearRing(exterior)

        tcrs = ccrs.PlateCarree()
        scrs = ccrs.PlateCarree()

        rings, mlinestr = tcrs._project_linear_ring(tring, scrs)

        # Number of linearstrings
        self.assertEqual(len(mlinestr), 4)
        self.assertFalse(rings)

        # Test area of smallest Polygon that contains all the points in the
        # geometry.
        self.assertAlmostEqual(mlinestr.convex_hull.area, 2347.75623076)


if __name__ == '__main__':
    unittest.main()
