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

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy
from shapely import geometry

import cartopy.crs as ccrs


class TestBoundary(unittest.TestCase):
    def test_no_polygon_boundary_reversal(self):
        # Check that polygons preserve their clockwise or counter-clockwise
        # ordering when they are attached to the boundary.
        # Failure to do so will result in invalid polygons (their boundaries
        # cross-over).
        polygon = geometry.Polygon([(-10, 30), (10, 60), (10, 50)])
        projection = ccrs.Robinson(170.5)
        multi_polygon = projection.project_geometry(polygon)
        for polygon in multi_polygon:
            self.assertTrue(polygon.is_valid)

    def test_polygon_boundary_attachment(self):
        # Check the polygon is attached to the boundary even when no
        # intermediate point for one of the crossing segments would normally
        # exist.
        polygon = geometry.Polygon([(-10, 30), (10, 60), (10, 50)])
        projection = ccrs.Robinson(170.6)
        # This will raise an exception if the polygon/boundary intersection
        # fails.
        multi_polygon = projection.project_geometry(polygon)

    def test_out_of_bounds(self):
        # Check that a polygon that is completely out of the map boundary
        # produces an empty result.
        # XXX Check efficiency?
        projection = ccrs.TransverseMercator(central_longitude=0)

        polys = [
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
        for coords in polys:
            polygon = geometry.Polygon(coords)
            multi_polygon = projection.project_geometry(polygon)
            self.assertEqual(len(multi_polygon), 0)


class TestMisc(unittest.TestCase):
    def test_misc(self):
        projection = ccrs.TransverseMercator(central_longitude=-90)
        polygon = geometry.Polygon([(-10, 30), (10, 60), (10, 50)])
        multi_polygon = projection.project_geometry(polygon)

    def test_small(self):
        projection = ccrs.Mercator()
        polygon = geometry.Polygon([
            (-179.9173693847652942, -16.5017831356493616),
            (-180.0000000000000000, -16.0671326636424396),
            (-179.7933201090486079, -16.0208822567412312),
            ])
        multi_polygon = projection.project_geometry(polygon)
        self.assertEqual(len(multi_polygon), 1)
        self.assertEqual(len(multi_polygon[0].exterior.coords), 5)


class TestQuality(unittest.TestCase):
    def setUp(self):
        projection = ccrs.RotatedPole(pole_longitude=177.5,
                                                     pole_latitude=37.5)
        polygon = geometry.Polygon([
            (175.0, -57.19913331),
            (177.5, -57.38460319),
            (180.0, -57.445077),
            ])
        self.multi_polygon = projection.project_geometry(polygon)
        #from cartopy.tests import show
        #show(projection, self.multi_polygon)

    def test_no_split(self):
        # Start simple ... there should only be one projected polygon.
        self.assertEqual(len(self.multi_polygon), 1)

    def test_repeats(self):
        # Make sure we don't have repeated points at the boundary, because
        # they mess up the linear extrapolation to the boundary.

        # Make sure there aren't any repeated points.
        xy = numpy.array(self.multi_polygon[0].exterior.coords)
        same = (xy[1:] == xy[:-1]).all(axis=1)
        self.assertFalse(any(same), 'Repeated points in projected geometry.')

    def test_symmetry(self):
        # Make sure the number of points added on the way towards the
        # boundary is similar to the number of points added on the way away
        # from the boundary.

        # Identify all the contiguous sets of non-boundary points.
        xy = numpy.array(self.multi_polygon[0].exterior.coords)
        boundary = numpy.logical_or(xy[:, 1] == 90, xy[:, 1] == -90)
        regions = (boundary[1:] != boundary[:-1]).cumsum()
        regions = numpy.insert(regions, 0, 0)

        # For each region, check if the number of increasing steps is roughly
        # equal to the number of decreasing steps.
        for i in range(boundary[0], regions.max(), 2):
            indices = numpy.where(regions == i)
            x = xy[indices, 0]
            delta = numpy.diff(x)
            num_incr = numpy.count_nonzero(delta > 0)
            num_decr = numpy.count_nonzero(delta < 0)
            self.assertLess(abs(num_incr - num_decr), 3, 'Too much assymmetry.')


if __name__ == '__main__':
    unittest.main()
