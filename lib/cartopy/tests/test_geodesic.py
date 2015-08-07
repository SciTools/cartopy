# (C) British Crown Copyright 2015, Met Office
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

from __future__ import (absolute_import, division, print_function)

import unittest

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal
from nose.tools import assert_equal

from cartopy import geodesic


class TestGeodesic(unittest.TestCase):
    
    def test_dir(self):
        geod = geodesic.Geodesic()
        geod_dir = geod.direct([40, 50], 38, 2000)
        assert_almost_equal(geod_dir, np.array([[40.01717933, 50.01416786, 
                            38.01316149]]), decimal=3)
                            
    def test_inverse(self):
        geod = geodesic.Geodesic()
        geod_inv = geod.inverse([40, 50], [45, 59])
        assert_almost_equal(geod_inv, np.array([[1.05216248e+06, 1.59036831e+01,
                            1.99877316e+01]]), decimal=3)

    def test_circle(self):
        geod = geodesic.Geodesic()
        geod_circle = geod.circle(40,50,500000,n_samples=3)
        assert_almost_equal(geod_circle, np.array([[40., 54.49349757], 
                            [34.23766162, 47.60355349], 
                            [45.76233838, 47.60355349]]), decimal=3)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
