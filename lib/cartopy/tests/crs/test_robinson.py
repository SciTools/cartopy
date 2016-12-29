# (C) British Crown Copyright 2013 - 2016, Met Office
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
'''
Tests for Robinson projection.

For now, mostly tests the workaround for a specific problem.
Problem report in : https://github.com/SciTools/cartopy/issues/23
Fix covered in : https://github.com/SciTools/cartopy/pull/277

'''

from __future__ import (absolute_import, division, print_function)

import sys
import unittest

from nose.tools import assert_true, assert_false, assert_equal
import numpy as np
from numpy.testing import assert_array_almost_equal as assert_arr_almost_eq

import cartopy.crs as ccrs


_NAN = float('nan')
_CRS_PC = ccrs.PlateCarree()
_CRS_ROB = ccrs.Robinson()

# Increase tolerance if using older proj.4 releases
_TOL = -1 if ccrs.PROJ4_VERSION < (4, 9) else 7


def test_transform_point():
    # this way has always worked
    result = _CRS_ROB.transform_point(35.0, 70.0, _CRS_PC)
    assert_arr_almost_eq(result, (2376187.27182751, 7275317.81573085), _TOL)

    # this always did something, but result has altered
    result = _CRS_ROB.transform_point(_NAN, 70.0, _CRS_PC)
    assert_true(np.all(np.isnan(result)))

    # this used to crash + is now fixed
    result = _CRS_ROB.transform_point(35.0, _NAN, _CRS_PC)
    assert_true(np.all(np.isnan(result)))


def test_transform_points():
    # these always worked
    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([35.0]),
                                       np.array([70.0]))
    assert_arr_almost_eq(result,
                         [[2376187.27182751, 7275317.81573085, 0]], _TOL)

    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([35.0]),
                                       np.array([70.0]),
                                       np.array([0.0]))
    assert_arr_almost_eq(result,
                         [[2376187.27182751, 7275317.81573085, 0]], _TOL)

    # this always did something, but result has altered
    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([_NAN]),
                                       np.array([70.0]))
    assert_true(np.all(np.isnan(result)))

    # this used to crash + is now fixed
    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([35.0]),
                                       np.array([_NAN]))
    assert_true(np.all(np.isnan(result)))

    # multipoint case
    x = np.array([10.0, 21.0, 0.0, 77.7, _NAN, 0.0])
    y = np.array([10.0, _NAN, 10.0, 77.7, 55.5, 0.0])
    z = np.array([10.0, 0.0, 0.0, _NAN, 55.5, 0.0])
    expect_result = np.array(
        [[9.40422591e+05, 1.06952091e+06, 1.00000000e+01],
         [11.1, 11.2, 11.3],
         [0.0, 1069520.91213902, 0.0],
         [22.1, 22.2, 22.3],
         [33.1, 33.2, 33.3],
         [0.0, 0.0, 0.0]])
    result = _CRS_ROB.transform_points(_CRS_PC, x, y, z)
    assert_equal(result.shape, (6, 3))
    assert_true(np.all(np.isnan(result[[1, 3, 4], :])))
    result[[1, 3, 4], :] = expect_result[[1, 3, 4], :]
    assert_false(np.any(np.isnan(result)))
    assert_true(np.allclose(result, expect_result))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
