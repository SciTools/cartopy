# (C) British Crown Copyright 2015 - 2018, Met Office
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

import matplotlib
from matplotlib.path import Path
import pytest
import shapely.geometry as sgeom

import cartopy.mpl.patch as cpatch


class Test_path_to_geos(object):
    def test_empty_polyon(self):
        p = Path([[0, 0], [0, 0], [0, 0], [0, 0],
                  [1, 2], [1, 2], [1, 2], [1, 2]],
                 codes=[1, 2, 2, 79,
                        1, 2, 2, 79])
        geoms = cpatch.path_to_geos(p)
        assert [type(geom) for geom in geoms] == [sgeom.Point, sgeom.Point]
        assert len(geoms) == 2

    @pytest.mark.skipif(matplotlib.__version__ < '2.2.0',
                        reason='Paths may not be closed with old Matplotlib.')
    def test_non_polygon_loop(self):
        p = Path([[0, 10], [170, 20], [-170, 30], [0, 10]],
                 codes=[1, 2, 2, 2])
        geoms = cpatch.path_to_geos(p)
        assert [type(geom) for geom in geoms] == [sgeom.MultiLineString]
        assert len(geoms) == 1

    def test_polygon_with_interior_and_singularity(self):
        # A geometry with two interiors, one a single point.
        p = Path([[0, -90], [200, -40], [200, 40], [0, 40], [0, -90],
                  [126, 26], [126, 26], [126, 26], [126, 26], [126, 26],
                  [114, 5], [103, 8], [126, 12], [126, 0], [114, 5]],
                 codes=[1, 2, 2, 2, 79, 1, 2, 2, 2, 79, 1, 2, 2, 2, 79])
        geoms = cpatch.path_to_geos(p)
        assert [type(geom) for geom in geoms] == [sgeom.Polygon, sgeom.Point]
        assert len(geoms[0].interiors) == 1
