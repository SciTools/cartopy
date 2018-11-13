# (C) British Crown Copyright 2013 - 2018, Met Office
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
"""
Tests for Robinson projection.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal
import pytest

import cartopy.crs as ccrs


_CRS_PC = ccrs.PlateCarree()
_CRS_ROB = ccrs.Robinson()

# Increase tolerance if using older proj releases
_TOL = -1 if ccrs.PROJ4_VERSION < (4, 9) else 7
_LIMIT_TOL = -1  # if ccrs.PROJ4_VERSION < (5, 2, 0) else 7


def check_proj_params(crs, other_args):
    expected = other_args | {'proj=robin', 'no_defs'}
    proj_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == proj_params


def test_default():
    robin = ccrs.Robinson()
    other_args = {'a=6378137.0', 'lon_0=0'}
    check_proj_params(robin, other_args)

    assert_almost_equal(robin.x_limits,
                        [-17005833.3305252, 17005833.3305252])
    assert_almost_equal(robin.y_limits,
                        [-8625154.6651000, 8625154.6651000], _LIMIT_TOL)


def test_offset():
    crs = ccrs.Robinson()
    crs_offset = ccrs.Robinson(false_easting=1234, false_northing=-4321)
    other_args = {'a=6378137.0', 'lon_0=0', 'x_0=1234', 'y_0=-4321'}
    check_proj_params(crs_offset, other_args)
    assert tuple(np.array(crs.x_limits) + 1234) == crs_offset.x_limits
    assert tuple(np.array(crs.y_limits) - 4321) == crs_offset.y_limits


@pytest.mark.parametrize('lon', [-10.0, 10.0])
def test_central_longitude(lon):
    robin = ccrs.Robinson(central_longitude=lon)
    other_args = {'a=6378137.0', 'lon_0={}'.format(lon)}
    check_proj_params(robin, other_args)

    assert_almost_equal(robin.x_limits,
                        [-17005833.3305252, 17005833.3305252],
                        decimal=5)
    assert_almost_equal(robin.y_limits,
                        [-8625154.6651000, 8625154.6651000], _LIMIT_TOL)


def test_transform_point():
    """
    Mostly tests the workaround for a specific problem.
    Problem report in: https://github.com/SciTools/cartopy/issues/23
    Fix covered in: https://github.com/SciTools/cartopy/pull/277
    """

    # this way has always worked
    result = _CRS_ROB.transform_point(35.0, 70.0, _CRS_PC)
    assert_array_almost_equal(result, (2376187.27182751, 7275317.81573085),
                              _TOL)

    # this always did something, but result has altered
    result = _CRS_ROB.transform_point(np.nan, 70.0, _CRS_PC)
    assert np.all(np.isnan(result))

    # this used to crash + is now fixed
    result = _CRS_ROB.transform_point(35.0, np.nan, _CRS_PC)
    assert np.all(np.isnan(result))


def test_transform_points():
    """
    Mostly tests the workaround for a specific problem.
    Problem report in: https://github.com/SciTools/cartopy/issues/23
    Fix covered in: https://github.com/SciTools/cartopy/pull/277
    """

    # these always worked
    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([35.0]),
                                       np.array([70.0]))
    assert_array_almost_equal(result,
                              [[2376187.27182751, 7275317.81573085, 0]], _TOL)

    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([35.0]),
                                       np.array([70.0]),
                                       np.array([0.0]))
    assert_array_almost_equal(result,
                              [[2376187.27182751, 7275317.81573085, 0]], _TOL)

    # this always did something, but result has altered
    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([np.nan]),
                                       np.array([70.0]))
    assert np.all(np.isnan(result))

    # this used to crash + is now fixed
    result = _CRS_ROB.transform_points(_CRS_PC,
                                       np.array([35.0]),
                                       np.array([np.nan]))
    assert np.all(np.isnan(result))

    # multipoint case
    x = np.array([10.0, 21.0, 0.0, 77.7, np.nan, 0.0])
    y = np.array([10.0, np.nan, 10.0, 77.7, 55.5, 0.0])
    z = np.array([10.0, 0.0, 0.0, np.nan, 55.5, 0.0])
    expect_result = np.array(
        [[9.40422591e+05, 1.06952091e+06, 1.00000000e+01],
         [11.1, 11.2, 11.3],
         [0.0, 1069520.91213902, 0.0],
         [22.1, 22.2, 22.3],
         [33.1, 33.2, 33.3],
         [0.0, 0.0, 0.0]])
    result = _CRS_ROB.transform_points(_CRS_PC, x, y, z)
    assert result.shape == (6, 3)
    assert np.all(np.isnan(result[[1, 3, 4], :]))
    result[[1, 3, 4], :] = expect_result[[1, 3, 4], :]
    assert not np.any(np.isnan(result))
    assert np.allclose(result, expect_result)
