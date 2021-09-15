# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import numpy as np
from numpy.testing import assert_array_equal
import pytest

import cartopy.img_transform as img_trans
import cartopy.crs as ccrs


@pytest.mark.parametrize('xmin, xmax', [
    (-90, 0), (-90, 90), (-90, None),
    (0, 90), (0, None),
    (None, 0), (None, 90), (None, None)])
@pytest.mark.parametrize('ymin, ymax', [
    (-45, 0), (-45, 45), (-45, None),
    (0, 45), (0, None),
    (None, 0), (None, 45), (None, None)])
def test_mesh_projection_extent(xmin, xmax, ymin, ymax):
    proj = ccrs.PlateCarree()
    nx = 4
    ny = 2

    target_x, target_y, extent = img_trans.mesh_projection(
        proj, nx, ny,
        x_extents=(xmin, xmax),
        y_extents=(ymin, ymax))

    if xmin is None:
        xmin = proj.x_limits[0]
    if xmax is None:
        xmax = proj.x_limits[1]
    if ymin is None:
        ymin = proj.y_limits[0]
    if ymax is None:
        ymax = proj.y_limits[1]
    assert_array_equal(extent, [xmin, xmax, ymin, ymax])
    assert_array_equal(np.diff(target_x, axis=1), (xmax - xmin) / nx)
    assert_array_equal(np.diff(target_y, axis=0), (ymax - ymin) / ny)


def test_gridding_data_std_range():
    # Data which exists inside the standard projection bounds i.e.
    # [-180, 180].
    target_prj = ccrs.PlateCarree()
    # create 3 data points
    lats = np.array([65, 10, -45])
    lons = np.array([-90, 0, 90])
    data = np.array([1, 2, 3])
    data_trans = ccrs.Geodetic()

    target_x, target_y, extent = img_trans.mesh_projection(target_prj, 8, 4)

    image = img_trans.regrid(data, lons, lats, data_trans, target_prj,
                             target_x, target_y,
                             mask_extrapolated=True)

    # The expected image. n.b. on a map the data is reversed in the y axis.
    expected = np.array([[3, 3, 3, 3, 3, 3, 3, 3],
                         [3, 1, 2, 2, 2, 3, 3, 3],
                         [1, 1, 1, 2, 2, 2, 3, 1],
                         [1, 1, 1, 1, 1, 1, 1, 1]], dtype=np.float64)

    expected_mask = np.array(
        [[True, True, True, True, True, True, True, True],
         [True, False, False, False, False, False, False, True],
         [True, False, False, False, False, False, False, True],
         [True, True, True, True, True, True, True, True]])

    assert_array_equal([-180, 180, -90, 90], extent)
    assert_array_equal(expected, image)
    assert_array_equal(expected_mask, image.mask)


def test_gridding_data_outside_projection():
    # Data which exists outside the standard projection e.g. [0, 360] rather
    # than [-180, 180].
    target_prj = ccrs.PlateCarree()
    # create 3 data points
    lats = np.array([65, 10, -45])
    lons = np.array([120, 180, 240])
    data = np.array([1, 2, 3])
    data_trans = ccrs.Geodetic()

    target_x, target_y, extent = img_trans.mesh_projection(target_prj, 8, 4)

    image = img_trans.regrid(data, lons, lats, data_trans, target_prj,
                             target_x, target_y,
                             mask_extrapolated=True)

    # The expected image. n.b. on a map the data is reversed in the y axis.
    expected = np.array(
        [[3, 3, 3, 3, 3, 2, 2, 2],
         [3, 3, 3, 3, 1, 1, 2, 2],
         [3, 3, 3, 3, 1, 1, 1, 2],
         [3, 3, 3, 1, 1, 1, 1, 1]], dtype=np.float64)

    expected_mask = np.array(
        [[True, True, True, True, True, True, True, True],
         [False, False, True, True, True, True, False, False],
         [False, False, True, True, True, True, False, False],
         [True, True, True, True, True, True, True, True]])

    assert_array_equal([-180, 180, -90, 90], extent)
    assert_array_equal(expected, image)
    assert_array_equal(expected_mask, image.mask)
