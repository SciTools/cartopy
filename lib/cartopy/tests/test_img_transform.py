# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from cartopy.tests.conftest import (
    _HAS_PYKDTREE_OR_SCIPY,
    requires_pykdtree,
    requires_scipy,
)


if not _HAS_PYKDTREE_OR_SCIPY:
    pytest.skip("pykdtree or scipy are required", allow_module_level=True)
import cartopy.crs as ccrs
import cartopy.img_transform as img_trans


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
    expected = np.array([[1, 1, 2, 2, 3, 3, 3, 3],
                         [1, 1, 2, 2, 2, 3, 3, 3],
                         [1, 1, 1, 2, 2, 2, 3, 3],
                         [1, 1, 1, 2, 2, 2, 3, 3]], dtype=np.float64)

    expected_mask = np.array(
        [[True, False, False, False, False, False, False, True],
         [True, False, False, False, False, False, False, True],
         [True, False, False, False, False, False, False, True],
         [True, False, False, False, False, False, False, True]])

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
        [[False, False, True, True, True, True, False, False],
         [False, False, True, True, True, True, False, False],
         [False, False, True, True, True, True, False, False],
         [False, False, True, True, True, True, False, False]])

    assert_array_equal([-180, 180, -90, 90], extent)
    assert_array_equal(expected, image)
    assert_array_equal(expected_mask, image.mask)


@pytest.mark.parametrize("target_prj",
                         (ccrs.Mollweide(), ccrs.Orthographic()))
@pytest.mark.parametrize("use_scipy", (pytest.param(True, marks=requires_scipy),
                                       pytest.param(False, marks=requires_pykdtree)))
def test_regridding_with_invalid_extent(target_prj, use_scipy, monkeypatch):
    # tests that when a valid extent results in invalid points in the
    # transformed coordinates, the regridding does not error.

    # create 3 data points
    lats = np.array([65, 10, -45])
    lons = np.array([-170, 10, 170])
    data = np.array([1, 2, 3])
    data_trans = ccrs.Geodetic()

    target_x, target_y, extent = img_trans.mesh_projection(target_prj, 8, 4)

    if use_scipy:
        monkeypatch.setattr(img_trans, "_is_pykdtree", False)
        import scipy.spatial
        monkeypatch.setattr(img_trans, "_kdtreeClass", scipy.spatial.cKDTree)
    _ = img_trans.regrid(data, lons, lats, data_trans, target_prj,
                         target_x, target_y)


@pytest.fixture
def warp_cache():
    img_trans._warp_mapping.cache_clear()
    yield img_trans._warp_mapping
    img_trans._warp_mapping.cache_clear()


def _warp_array_uncached(array, target_proj, source_proj, target_res,
                         source_extent, target_extent, mask_extrapolated):
    source_x, source_y, _ = img_trans.mesh_projection(
        source_proj, array.shape[1], array.shape[0],
        x_extents=source_extent[:2], y_extents=source_extent[2:])
    target_x, target_y, extent = img_trans.mesh_projection(
        target_proj, *target_res,
        x_extents=target_extent[:2], y_extents=target_extent[2:])
    new_array = img_trans.regrid(array, source_x, source_y, source_proj,
                                 target_proj, target_x, target_y,
                                 mask_extrapolated)
    return new_array, extent


def _assert_same_result(actual, expected):
    assert type(actual) is type(expected)
    assert actual.dtype == expected.dtype
    assert_array_equal(np.ma.getdata(actual), np.ma.getdata(expected))
    assert_array_equal(np.ma.getmaskarray(actual),
                       np.ma.getmaskarray(expected))


@pytest.mark.parametrize('mask_extrapolated', [False, True])
@pytest.mark.parametrize('target_prj', [ccrs.Mercator(), ccrs.Robinson(),
                                        ccrs.Orthographic(30, 40)])
def test_warp_array_reuses_mapping(target_prj, mask_extrapolated, warp_cache):
    source_prj = ccrs.PlateCarree()
    source_extent = (-180, 180, -60, 80)
    rng = np.random.default_rng(0)
    nan_data = rng.random((18, 36))
    nan_data[::4, ::5] = np.nan
    frames = [
        rng.random((18, 36)),
        nan_data,
        np.ma.masked_greater(rng.random((18, 36)), 0.8),
        rng.integers(0, 256, (18, 36, 4), dtype=np.uint8),
    ]

    for frame in frames:
        expected, expected_extent = _warp_array_uncached(
            frame, target_prj, source_prj, (30, 20), source_extent,
            [None] * 4, mask_extrapolated)
        result, extent = img_trans.warp_array(
            frame, target_prj, source_prj, target_res=(30, 20),
            source_extent=source_extent, mask_extrapolated=mask_extrapolated)
        _assert_same_result(result, expected)
        assert extent == expected_extent

    info = warp_cache.cache_info()
    assert (info.misses, info.hits) == (1, len(frames) - 1)


@pytest.mark.parametrize('change', [
    dict(target_proj=ccrs.Mollweide()),
    dict(source_proj=ccrs.PlateCarree(central_longitude=180)),
    dict(target_res=(20, 30)),
    dict(source_extent=(-170, 190, -90, 90)),
    dict(source_extent=np.array([-180, 180, -90, 90], dtype=np.float32)),
    dict(target_extent=(-1e7, 1e7, -5e6, 5e6)),
    dict(mask_extrapolated=True),
    dict(array=np.arange(17 * 36.).reshape(17, 36)),
])
def test_warp_array_mapping_cache_key(change, warp_cache):
    kwargs = dict(array=np.arange(18 * 36.).reshape(18, 36),
                  target_proj=ccrs.Robinson(),
                  source_proj=ccrs.PlateCarree(),
                  target_res=(30, 20),
                  source_extent=(-180, 180, -90, 90),
                  target_extent=(-1.5e7, 1.5e7, -8e6, 8e6),
                  mask_extrapolated=False)
    img_trans.warp_array(**kwargs)
    img_trans.warp_array(**kwargs)
    assert warp_cache.cache_info().misses == 1

    kwargs.update(change)
    result, extent = img_trans.warp_array(**kwargs)
    assert warp_cache.cache_info().misses == 2
    expected, expected_extent = _warp_array_uncached(**kwargs)
    _assert_same_result(result, expected)
    assert extent == expected_extent


def test_warp_array_mapping_cache_is_not_modified(warp_cache):
    data = np.arange(18 * 36.).reshape(18, 36)
    kwargs = dict(target_proj=ccrs.Robinson(), source_proj=ccrs.PlateCarree(),
                  target_res=(30, 20), source_extent=(-180, 180, -90, 90))
    expected, expected_extent = img_trans.warp_array(data, **kwargs)

    # The corners of the target grid are outside of Robinson, so masked.
    assert np.ma.is_masked(expected)

    result, extent = img_trans.warp_array(data, **kwargs)
    result[:] = -1
    result.mask = True
    extent[0] = 0
    result, extent = img_trans.warp_array(data, **kwargs)
    _assert_same_result(result, expected)
    assert extent == expected_extent


def test_warp_array_mapping_cache_size(warp_cache):
    data = np.zeros((9, 18))
    maxsize = warp_cache.cache_info().maxsize
    for i in range(maxsize + 2):
        img_trans.warp_array(data, ccrs.Robinson(), target_res=(10 + i, 10))
    assert warp_cache.cache_info().currsize == maxsize
