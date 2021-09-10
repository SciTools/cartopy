# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from matplotlib.path import Path
import shapely.geometry as sgeom

import cartopy.mpl.patch as cpatch


class Test_path_to_geos:
    def test_empty_polygon(self):
        p = Path(
            [
                [0, 0], [0, 0], [0, 0], [0, 0],
                [1, 2], [1, 2], [1, 2], [1, 2],
                # The vertex for CLOSEPOLY should be ignored.
                [2, 3], [2, 3], [2, 3], [42, 42],
                # Very close points should be treated the same.
                [193.75, -14.166664123535156], [193.75, -14.166664123535158],
                [193.75, -14.166664123535156], [193.75, -14.166664123535156],
            ],
            codes=[1, 2, 2, 79] * 4)
        geoms = cpatch.path_to_geos(p)
        assert [type(geom) for geom in geoms] == [sgeom.Point] * 4
        assert len(geoms) == 4

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
