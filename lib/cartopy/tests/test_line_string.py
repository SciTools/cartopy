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


import itertools
import time
import unittest

import numpy
from shapely import geometry

import cartopy.crs as ccrs


class TestLineString(unittest.TestCase):
    def test_out_of_bounds(self):
        # Check that a line that is completely out of the map boundary produces
        # a valid LineString
        projection = ccrs.TransverseMercator(central_longitude=0)

        # For both start & end, define a point that results in well-defined
        # projection coordinates and one that results in NaN.
        start_points = [(86, 0), (130, 0)]
        end_points = [(88, 0), (120, 0)]

        # Try all four combinations of valid/NaN vs valid/NaN.
        for start, end in itertools.product(start_points, end_points):
            line_string = geometry.LineString([start, end])
            multi_line_string = projection.project_geometry(line_string)
            if start[0] == 130 and end[0] == 120:
                expected = 0
            else:
                expected = 1
            self.assertEqual(len(multi_line_string), expected,
                             'Unexpected line when working from {} '
                             'to {}'.format(start, end))

    def test_simple_fragment_count(self):
        projection = ccrs.PlateCarree()

        tests = [
            ([(150, 0), (-150, 0)], 2),
            ([(10, 0), (90, 0), (180, 0), (-90, 0), (-10, 0)], 2),
            ([(-10, 0), (10, 0)], 1),
            ([(-45, 0), (45, 30)], 1),
        ]

        for coords, pieces in tests:
            line_string = geometry.LineString(coords)
            multi_line_string = projection.project_geometry(line_string)
            #from cartopy.tests import show
            #show(projection, multi_line_string)
            self.assertEqual(len(multi_line_string), pieces)

    def test_split(self):
        projection = ccrs.Robinson(170.5)
        line_string = geometry.LineString([(-10, 30), (10, 60)])
        multi_line_string = projection.project_geometry(line_string)
        from cartopy.tests import show
        #show(projection, multi_line_string)
        self.assertEqual(len(multi_line_string), 2)

    def test_out_of_domain_efficiency(self):
        # Check we're efficiently dealing with lines that project
        # outside the map domain.
        # Because the south pole projects to an *enormous* circle
        # (radius ~ 1e23) this will take a *long* time to project if the
        # within-domain exactness criteria are used.
        line_string = geometry.LineString([(0, -90), (2, -90)])
        tgt_proj = ccrs.NorthPolarStereo()
        src_proj = ccrs.PlateCarree()
        cutoff_time = time.time() + 1
        tgt_proj.project_geometry(line_string, src_proj)
        self.assertLess(time.time(), cutoff_time, 'Projection took too long')


class FakeProjection(ccrs.PlateCarree):
    def __init__(self, left_offset=0, right_offset=0):
        self.left_offset = left_offset
        self.right_offset = right_offset

        self._half_width = 180
        self._half_height = 90
        ccrs.PlateCarree.__init__(self)

    @property
    def boundary(self):
        # XXX Should this be a LinearRing?
        w, h = self._half_width, self._half_height
        from shapely import geometry
        return geometry.LineString([(-w + self.left_offset, -h),
                                    (-w + self.left_offset, h),
                                    (w - self.right_offset, h),
                                    (w - self.right_offset, -h),
                                    (-w + self.left_offset, -h)])


class TestBisect(unittest.TestCase):
    # A bunch of tests to check the bisection algorithm is robust for a
    # variety of simple and/or pathological cases.

    def test_repeated_point(self):
        projection = FakeProjection()
        line_string = geometry.LineString([(10, 0), (10, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_interior_repeated_point(self):
        projection = FakeProjection()
        line_string = geometry.LineString([(0, 0), (10, 0), (10, 0), (20, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 4)

    def test_circular_repeated_point(self):
        projection = FakeProjection()
        line_string = geometry.LineString([(0, 0), (360, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_short(self):
        projection = FakeProjection()
        line_string = geometry.LineString([(0, 0), (1e-12, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_empty(self):
        projection = FakeProjection(right_offset=10)
        line_string = geometry.LineString([(175, 0), (175, 10)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 0)

    def test_simple_run_in(self):
        projection = FakeProjection(right_offset=10)
        line_string = geometry.LineString([(160, 0), (175, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_simple_wrap(self):
        projection = FakeProjection()
        line_string = geometry.LineString([(160, 0), (-160, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 2)
        self.assertEqual(len(multi_line_string[0].coords), 2)
        self.assertEqual(len(multi_line_string[1].coords), 2)

    def test_simple_run_out(self):
        projection = FakeProjection(left_offset=10)
        line_string = geometry.LineString([(-175, 0), (-160, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_point_on_boundary(self):
        projection = FakeProjection()
        line_string = geometry.LineString([(180, 0), (-160, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

        # Add a small offset to the left-hand boundary to make things
        # even more pathological.
        projection = FakeProjection(left_offset=5)
        line_string = geometry.LineString([(180, 0), (-160, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_nan_start(self):
        projection = ccrs.TransverseMercator(central_longitude=-90)
        line_string = geometry.LineString([(10, 50), (-10, 30)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        for line_string in multi_line_string:
            for coord in line_string.coords:
                self.assertFalse(any(numpy.isnan(coord)),
                                 'Unexpected NaN in projected coords.')

    def test_nan_end(self):
        projection = ccrs.TransverseMercator(central_longitude=-90)
        line_string = geometry.LineString([(-10, 30), (10, 50)])
        multi_line_string = projection.project_geometry(line_string)
        from cartopy.tests import show
        #show(projection, multi_line_string)
        self.assertEqual(len(multi_line_string), 1)
        for line_string in multi_line_string:
            for coord in line_string.coords:
                self.assertFalse(any(numpy.isnan(coord)), 
                                 'Unexpected NaN in projected coords.')


class TestMisc(unittest.TestCase):
    def test_misc(self):
        projection = ccrs.TransverseMercator(central_longitude=-90)
        line_string = geometry.LineString([(10, 50), (-10, 30)])
        multi_line_string = projection.project_geometry(line_string)
        from cartopy.tests import show
        #show(projection, multi_line_string)
        for line_string in multi_line_string:
            for coord in line_string.coords:
                self.assertFalse(any(numpy.isnan(coord)),
                                 'Unexpected NaN in projected coords.')

    def test_something(self):
        projection = ccrs.RotatedPole(pole_longitude=177.5,
                                      pole_latitude=37.5)
        line_string = geometry.LineString([(0, 0), (1e-14, 0)])
        multi_line_string = projection.project_geometry(line_string)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string[0].coords), 2)

    def test_global_boundary(self):
        linear_ring = geometry.LineString([(-180, -180), (-180, 180),
                                           (180, 180), (180, -180)])
        pc = ccrs.PlateCarree()
        merc = ccrs.Mercator()
        multi_line_string = pc.project_geometry(linear_ring, merc)
        assert len(multi_line_string) > 0

        # check the identity transform
        multi_line_string = merc.project_geometry(linear_ring, merc)
        assert len(multi_line_string) > 0


class TestSymmetry(unittest.TestCase):
    @unittest.expectedFailure
    def test_curve(self):
        # Obtain a simple, curved path.
        projection = ccrs.PlateCarree()
        coords = [(-0.08, 51.53), (132.00, 43.17)]  # London to Vladivostock
        line_string = geometry.LineString(coords)
        multi_line_string = projection.project_geometry(line_string)

        # Compute the reverse path.
        line_string = geometry.LineString(coords[::-1])
        multi_line_string2 = projection.project_geometry(line_string)

        # Make sure that they generated the same points.
        # (Although obviously they will be in the opposite order!)
        self.assertEqual(len(multi_line_string), 1)
        self.assertEqual(len(multi_line_string2), 1)
        coords = multi_line_string[0].coords
        coords2 = multi_line_string2[0].coords
        numpy.testing.assert_allclose(coords, coords2[::-1],
                                      err_msg='Asymmetric curve generation')


if __name__ == '__main__':
    unittest.main()
