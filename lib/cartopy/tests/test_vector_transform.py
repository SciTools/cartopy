# (C) British Crown Copyright 2013, Met Office
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
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

import cartopy.vector_transform as vec_trans
import cartopy.crs as ccrs


def test_interpolate_data_extent():
    # Interpolation to a grid with extents of the input data.
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    s = np.array([2, 4, 2, 1.2, 3, 1.2])

    expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.]])
    expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                [7.5, 7.5, 7.5, 7.5, 7.5],
                                [10., 10., 10., 10., 10]])
    expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                [np.nan, 2.5, 3.5, 2.5, np.nan],
                                [2., 3., 4., 3., 2.]])

    x_grid, y_grid, s_grid = vec_trans._interpolate_to_grid(x, y, s, 5, 3)

    assert_array_equal(x_grid, expected_x_grid)
    assert_array_equal(y_grid, expected_y_grid)
    assert_array_almost_equal(s_grid, expected_s_grid)


def test_interpolate_explicit_extent():
    # Interpolation to a grid with explicit extents.
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    s = np.array([2, 4, 2, 1.2, 3, 1.2])

    expected_x_grid = np.array([[-5., 0., 5., 10.],
                                [-5., 0., 5., 10.]])
    expected_y_grid = np.array([[7.5, 7.5, 7.5, 7.5],
                                [10., 10., 10., 10]])
    expected_s_grid = np.array([[2.5, 3.5, 2.5, np.nan],
                                [3., 4., 3., 2.]])

    extent = (-5, 10, 7.5, 10)
    x_grid, y_grid, s_grid = vec_trans._interpolate_to_grid(
        x, y, s, 4, 2, target_extent=extent)

    assert_array_equal(x_grid, expected_x_grid)
    assert_array_equal(y_grid, expected_y_grid)
    assert_array_almost_equal(s_grid, expected_s_grid)


def test_scalar_to_grid_no_transform():
    # Transform and regrid scalar (with no projection transform).
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    s = np.array([2, 4, 2, 1.2, 3, 1.2])

    expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.]])
    expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                [7.5, 7.5, 7.5, 7.5, 7.5],
                                [10., 10., 10., 10., 10]])
    expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                [np.nan, 2.5, 3.5, 2.5, np.nan],
                                [2., 3., 4., 3., 2.]])

    src_crs = target_crs = ccrs.PlateCarree()
    x_grid, y_grid, s_grid = vec_trans.scalar_to_grid(
        src_crs, target_crs, x, y, s, (5, 3), target_extent=None)

    assert_array_equal(x_grid, expected_x_grid)
    assert_array_equal(y_grid, expected_y_grid)
    assert_array_almost_equal(s_grid, expected_s_grid)


def test_scalar_to_grid_with_transform():
    # Transform and regrid scalar.
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    s = np.array([2, 4, 2, 1.2, 3, 1.2])

    target_crs = ccrs.PlateCarree()
    src_crs = ccrs.NorthPolarStereo()

    input_coords = [src_crs.transform_point(xp, yp, target_crs)
                    for xp, yp in zip(x, y)]
    x_nps = np.array([ic[0] for ic in input_coords])
    y_nps = np.array([ic[1] for ic in input_coords])

    expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.]])
    expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                [7.5, 7.5, 7.5, 7.5, 7.5],
                                [10., 10., 10., 10., 10]])
    expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                [np.nan, 2.5, 3.5, 2.5, np.nan],
                                [2., 3., 4., 3., 2.]])

    x_grid, y_grid, s_grid = vec_trans.scalar_to_grid(
        src_crs, target_crs, x_nps, y_nps, s, (5, 3))

    assert_array_almost_equal(x_grid, expected_x_grid)
    assert_array_almost_equal(y_grid, expected_y_grid)
    assert_array_almost_equal(s_grid, expected_s_grid)


def test_vector_to_grid_no_transform():
    # Transform and regrid vector (with no projection transform).
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    u = np.array([2, 4, 2, 1.2, 3, 1.2])
    v = np.array([5.5, 4, 5.5, 1.2, .3, 1.2])

    expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.]])
    expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                [7.5, 7.5, 7.5, 7.5, 7.5],
                                [10., 10., 10., 10., 10]])
    expected_u_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                [np.nan, 2.5, 3.5, 2.5, np.nan],
                                [2., 3., 4., 3., 2.]])
    expected_v_grid = np.array([[np.nan, .8, .3, .8, np.nan],
                                [np.nan, 2.675, 2.15, 2.675, np.nan],
                                [5.5, 4.75, 4., 4.75, 5.5]])

    src_crs = target_crs = ccrs.PlateCarree()
    x_grid, y_grid, u_grid, v_grid = vec_trans.vector_to_grid(
        src_crs, target_crs, x, y, u, v, (5, 3))

    assert_array_equal(x_grid, expected_x_grid)
    assert_array_equal(y_grid, expected_y_grid)
    assert_array_almost_equal(u_grid, expected_u_grid)
    assert_array_almost_equal(v_grid, expected_v_grid)


def test_vector_to_grid_with_transform():
    # Transform and regrid vector.
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    u = np.array([2, 4, 2, 1.2, 3, 1.2])
    v = np.array([5.5, 4, 5.5, 1.2, .3, 1.2])

    target_crs = ccrs.PlateCarree()
    src_crs = ccrs.NorthPolarStereo()

    input_coords = [src_crs.transform_point(xp, yp, target_crs)
                    for xp, yp in zip(x, y)]
    x_nps = np.array([ic[0] for ic in input_coords])
    y_nps = np.array([ic[1] for ic in input_coords])
    u_nps, v_nps = src_crs.transform_vectors(target_crs, x, y, u, v)

    expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.],
                                [-10., -5., 0., 5., 10.]])
    expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                [7.5, 7.5, 7.5, 7.5, 7.5],
                                [10., 10., 10., 10., 10]])
    expected_u_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                [np.nan, 2.5, 3.5, 2.5, np.nan],
                                [2., 3., 4., 3., 2.]])
    expected_v_grid = np.array([[np.nan, .8, .3, .8, np.nan],
                                [np.nan, 2.675, 2.15, 2.675, np.nan],
                                [5.5, 4.75, 4., 4.75, 5.5]])

    x_grid, y_grid, u_grid, v_grid = vec_trans.vector_to_grid(
        src_crs, target_crs, x_nps, y_nps, u_nps, v_nps, (5, 3))

    assert_array_almost_equal(x_grid, expected_x_grid)
    assert_array_almost_equal(y_grid, expected_y_grid)
    # Vector transforms are somewhat approximate, so we are more lenient
    # with the returned values since we have transformed twice.
    assert_array_almost_equal(u_grid, expected_u_grid, decimal=4)
    assert_array_almost_equal(v_grid, expected_v_grid, decimal=4)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
