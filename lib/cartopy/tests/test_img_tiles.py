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

from nose.tools import assert_equal, assert_raises
from numpy.testing import assert_array_almost_equal as assert_arr_almost
import shapely.geometry

import cartopy.io.img_tiles as cimgt


def test_google_wts():
    gt = cimgt.GoogleTiles()

    extent = [-15, 0.1, 50, 60]
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]],
                                              [extent[2], extent[1]],
                                              [extent[2], extent[3]],
                                              [extent[0], extent[3]],
                                              [extent[0], extent[1]]])

    with assert_raises(AssertionError):
        list(gt.find_images(target_domain, -1))
    assert_equal(tuple(gt.find_images(target_domain, 0)),
                 ((0, 0, 0),))
    assert_equal(tuple(gt.find_images(target_domain, 2)),
                 ((1, 1, 2), (2, 1, 2)))

    assert_equal(list(gt.subtiles((0, 0, 0))),
                 [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)])
    assert_equal(list(gt.subtiles((1, 0, 1))),
                 [(2, 0, 2), (2, 1, 2), (3, 0, 2), (3, 1, 2)])

    with assert_raises(AssertionError):
        gt.tileextent((0, 1, 0))

    assert_arr_almost(gt.tileextent((0, 0, 0)),
                      (-20037508.34278306, -20037508.34279542,
                       19929239.11334565, -19929239.11334564))
    assert_arr_almost(gt.tileextent((2, 0, 2)),
                      (0.0, 10018754.171391532,
                       19929239.113345649, 9964619.5566728264))
    assert_arr_almost(gt.tileextent((0, 2, 2)),
                      (-20037508.342783064, -10018754.171391532,
                       3.7252902984619141e-09, -9964619.556672819))
    assert_arr_almost(gt.tileextent((2, 2, 2)),
                      (0.0, 10018754.171391532,
                       3.7252902984619141e-09, -9964619.556672819))
    assert_arr_almost(gt.tileextent((8, 9, 4)),
                      (0.0, 2504688.542847883, -2491154.8891682029,
                       -4982309.7783364095))  # zoom 4, contains cape town.


def test_quadtree_wts():
    qt = cimgt.QuadtreeTiles()

    extent = [-15, 0.1, 50, 60]
    target_domain = shapely.geometry.Polygon([[extent[0], extent[1]],
                                              [extent[2], extent[1]],
                                              [extent[2], extent[3]],
                                              [extent[0], extent[3]],
                                              [extent[0], extent[1]]])

    with assert_raises(ValueError):
        list(qt.find_images(target_domain, 0))

    assert_equal(qt.tms_to_quadkey((1, 1, 1)), '1')
    assert_equal(qt.quadkey_to_tms('1'), (1, 1, 1))

    assert_equal(qt.tms_to_quadkey((8, 9, 4)), '1220')
    assert_equal(qt.quadkey_to_tms('1220'), (8, 9, 4))

    assert_equal(tuple(qt.find_images(target_domain, 1)), ('0', '1'))
    assert_equal(tuple(qt.find_images(target_domain, 2)), ('03', '12'))

    assert_equal(list(qt.subtiles('0')), ['00', '01', '02', '03'])
    assert_equal(list(qt.subtiles('11')), ['110', '111', '112', '113'])

    with assert_raises(ValueError):
        qt.tileextent('4')

    assert_arr_almost(qt.tileextent(''),
                      (-20037508.34278306, -20037508.34279542,
                       19929239.11334565, -19929239.11334564))
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((2, 0, 2), google=True)),
                      (0.0, 10018754.171391532,
                       19929239.113345649, 9964619.5566728264))
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((0, 2, 2), google=True)),
                      (-20037508.342783064, -10018754.171391532,
                       3.7252902984619141e-09, -9964619.556672819))
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((0, 1, 2), google=True)),
                      (-20037508.342783064, -10018754.171391532,
                       9964619.5566728264, 3.7252902984619141e-09))
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((2, 2, 2), google=True)),
                      (0.0, 10018754.171391532,
                       3.7252902984619141e-09, -9964619.556672819))
    assert_arr_almost(qt.tileextent(qt.tms_to_quadkey((8, 9, 4), google=True)),
                      (0.0, 2504688.542847883,
                       -2491154.8891682029, -4982309.7783364095))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
