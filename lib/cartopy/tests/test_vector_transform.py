# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import pytest


try:
    import scipy  # noqa: F401
except ImportError:
    pytest.skip("scipy is required for vector transforms", allow_module_level=True)


import cartopy.crs as ccrs
import cartopy.vector_transform as vec_trans


def _sample_plate_carree_coordinates():
    x = np.array([-10, 0, 10, -9, 0, 9])
    y = np.array([10, 10, 10, 5, 5, 5])
    return x, y


def _sample_plate_carree_scalar_field():
    return np.array([2, 4, 2, 1.2, 3, 1.2])


def _sample_plate_carree_vector_field():
    u = np.array([2, 4, 2, 1.2, 3, 1.2])
    v = np.array([5.5, 4, 5.5, 1.2, .3, 1.2])
    return u, v


class Test_interpolate_to_grid:

    @classmethod
    def setup_class(cls):
        cls.x, cls.y = _sample_plate_carree_coordinates()
        cls.s = _sample_plate_carree_scalar_field()

    def test_data_extent(self):
        # Interpolation to a grid with extents of the input data.
        expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                    [-10., -5., 0., 5., 10.],
                                    [-10., -5., 0., 5., 10.]])
        expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                    [7.5, 7.5, 7.5, 7.5, 7.5],
                                    [10., 10., 10., 10., 10]])
        expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                    [np.nan, 2.5, 3.5, 2.5, np.nan],
                                    [2., 3., 4., 3., 2.]])

        x_grid, y_grid, s_grid = vec_trans._interpolate_to_grid(
            5, 3, self.x, self.y, self.s)

        assert_array_equal(x_grid, expected_x_grid)
        assert_array_equal(y_grid, expected_y_grid)
        assert_array_almost_equal(s_grid, expected_s_grid)

    def test_explicit_extent(self):
        # Interpolation to a grid with explicit extents.
        expected_x_grid = np.array([[-5., 0., 5., 10.],
                                    [-5., 0., 5., 10.]])
        expected_y_grid = np.array([[7.5, 7.5, 7.5, 7.5],
                                    [10., 10., 10., 10]])
        expected_s_grid = np.array([[2.5, 3.5, 2.5, np.nan],
                                    [3., 4., 3., 2.]])

        extent = (-5, 10, 7.5, 10)
        x_grid, y_grid, s_grid = vec_trans._interpolate_to_grid(
            4, 2, self.x, self.y, self.s, target_extent=extent)

        assert_array_equal(x_grid, expected_x_grid)
        assert_array_equal(y_grid, expected_y_grid)
        assert_array_almost_equal(s_grid, expected_s_grid)

    def test_multiple_fields(self):
        # Interpolation of multiple fields in one go.
        expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                    [-10., -5., 0., 5., 10.],
                                    [-10., -5., 0., 5., 10.]])
        expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                    [7.5, 7.5, 7.5, 7.5, 7.5],
                                    [10., 10., 10., 10., 10]])
        expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                    [np.nan, 2.5, 3.5, 2.5, np.nan],
                                    [2., 3., 4., 3., 2.]])

        x_grid, y_grid, s_grid1, s_grid2, s_grid3 = \
            vec_trans._interpolate_to_grid(5, 3, self.x, self.y,
                                           self.s, self.s, self.s)

        assert_array_equal(x_grid, expected_x_grid)
        assert_array_equal(y_grid, expected_y_grid)
        assert_array_almost_equal(s_grid1, expected_s_grid)
        assert_array_almost_equal(s_grid2, expected_s_grid)
        assert_array_almost_equal(s_grid3, expected_s_grid)


class Test_vector_scalar_to_grid:

    @classmethod
    def setup_class(cls):
        cls.x, cls.y = _sample_plate_carree_coordinates()
        cls.u, cls.v = _sample_plate_carree_vector_field()
        cls.s = _sample_plate_carree_scalar_field()

    def test_no_transform(self):
        # Transform and regrid vector (with no projection transform).
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
        x_grid, y_grid, u_grid, v_grid = vec_trans.vector_scalar_to_grid(
            src_crs, target_crs, (5, 3), self.x, self.y, self.u, self.v)

        assert_array_equal(x_grid, expected_x_grid)
        assert_array_equal(y_grid, expected_y_grid)
        assert_array_almost_equal(u_grid, expected_u_grid)
        assert_array_almost_equal(v_grid, expected_v_grid)

    def test_with_transform(self):
        # Transform and regrid vector.
        target_crs = ccrs.PlateCarree()
        src_crs = ccrs.NorthPolarStereo()

        input_coords = [src_crs.transform_point(xp, yp, target_crs)
                        for xp, yp in zip(self.x, self.y)]
        x_nps = np.array([ic[0] for ic in input_coords])
        y_nps = np.array([ic[1] for ic in input_coords])
        u_nps, v_nps = src_crs.transform_vectors(target_crs, self.x, self.y,
                                                 self.u, self.v)

        expected_x_grid = np.array([[-10., -5., 0., 5., 10.],
                                    [-10., -5., 0., 5., 10.],
                                    [-10., -5., 0., 5., 10.]])
        expected_y_grid = np.array([[5., 5., 5., 5., 5.],
                                    [7.5, 7.5, 7.5, 7.5, 7.5],
                                    [10., 10., 10., 10., 10]])
        expected_u_grid = np.array([[np.nan, np.nan, np.nan, np.nan, np.nan],
                                    [np.nan, 2.3838, 3.5025, 2.6152, np.nan],
                                    [2, 3.0043, 4, 2.9022, 2]])
        expected_v_grid = np.array([[np.nan, np.nan, np.nan, np.nan, np.nan],
                                    [np.nan, 2.6527, 2.1904, 2.4192, np.nan],
                                    [5.5, 4.6483, 4, 4.47, 5.5]])

        x_grid, y_grid, u_grid, v_grid = vec_trans.vector_scalar_to_grid(
            src_crs, target_crs, (5, 3), x_nps, y_nps, u_nps, v_nps)

        assert_array_almost_equal(x_grid, expected_x_grid)
        assert_array_almost_equal(y_grid, expected_y_grid)
        # Vector transforms are somewhat approximate, so we are more lenient
        # with the returned values since we have transformed twice.
        assert_array_almost_equal(u_grid, expected_u_grid, decimal=4)
        assert_array_almost_equal(v_grid, expected_v_grid, decimal=4)

    def test_with_scalar_field(self):
        # Transform and regrid vector (with no projection transform) with an
        # additional scalar field.
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
        expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                    [np.nan, 2.5, 3.5, 2.5, np.nan],
                                    [2., 3., 4., 3., 2.]])

        src_crs = target_crs = ccrs.PlateCarree()
        x_grid, y_grid, u_grid, v_grid, s_grid = \
            vec_trans.vector_scalar_to_grid(src_crs, target_crs, (5, 3),
                                            self.x, self.y,
                                            self.u, self.v, self.s)

        assert_array_equal(x_grid, expected_x_grid)
        assert_array_equal(y_grid, expected_y_grid)
        assert_array_almost_equal(u_grid, expected_u_grid)
        assert_array_almost_equal(v_grid, expected_v_grid)
        assert_array_almost_equal(s_grid, expected_s_grid)

    def test_with_scalar_field_non_ndarray_data(self):
        # Transform and regrid vector (with no projection transform) with an
        # additional scalar field which is not a ndarray.
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
        expected_s_grid = np.array([[np.nan, 2., 3., 2., np.nan],
                                    [np.nan, 2.5, 3.5, 2.5, np.nan],
                                    [2., 3., 4., 3., 2.]])

        src_crs = target_crs = ccrs.PlateCarree()
        x_grid, y_grid, u_grid, v_grid, s_grid = \
            vec_trans.vector_scalar_to_grid(src_crs, target_crs, (5, 3),
                                            list(self.x), list(self.y),
                                            list(self.u), list(self.v),
                                            list(self.s))

        assert_array_equal(x_grid, expected_x_grid)
        assert_array_equal(y_grid, expected_y_grid)
        assert_array_almost_equal(u_grid, expected_u_grid)
        assert_array_almost_equal(v_grid, expected_v_grid)
        assert_array_almost_equal(s_grid, expected_s_grid)


def _griddata_per_field(points, fields, xi):
    # The previous implementation, which triangulated once per field.
    from scipy.interpolate import griddata
    return tuple(griddata(points, s.ravel(), xi, method='linear')
                 for s in fields)


def _assert_identical(actual, expected):
    assert len(actual) == len(expected)
    for a, e in zip(actual, expected):
        assert_array_equal(np.ma.getdata(a), np.ma.getdata(e))
        assert_array_equal(np.ma.getmaskarray(a), np.ma.getmaskarray(e))


class Test_interpolate_fields:

    @classmethod
    def setup_class(cls):
        rng = np.random.default_rng(2481)
        # A regular grid has four co-circular points in every cell, so the
        # triangulation is not unique and depends on the full point set.
        x, y = np.meshgrid(np.linspace(-0.1, 1.1, 37),
                           np.linspace(-0.2, 1.2, 23))
        cls.regular = np.column_stack([x.ravel(), y.ravel()])
        cls.scattered = rng.uniform(-0.1, 1.1, size=(500, 2))
        qx, qy = np.meshgrid(np.linspace(-0.3, 1.3, 31),
                             np.linspace(-0.3, 1.3, 29))
        # Some target points can not be projected back to the source
        # CRS; they come through as masked non-finite coordinates.
        qx = np.ma.masked_invalid(np.where(qx > 1.2, np.nan, qx))
        cls.xi = (qx, qy)

    def _fields(self, n):
        s = np.sin(7 * np.random.default_rng(n).random(n)) * 10
        with_nan = s.copy()
        with_nan[::17] = np.nan
        masked = np.ma.masked_greater(s, 5)
        return (s, s.astype(np.float32), (s * 100).astype(int), with_nan,
                masked)

    @pytest.mark.parametrize('points', ['regular', 'scattered'])
    def test_same_as_griddata(self, points):
        points = getattr(self, points)
        fields = self._fields(len(points))
        actual = vec_trans._interpolate_fields(points, fields, self.xi)
        expected = _griddata_per_field(points, fields, self.xi)
        _assert_identical(actual, expected)

    def test_no_fields(self, monkeypatch):
        def fail(points):
            raise AssertionError('no triangulation is needed')

        monkeypatch.setattr(vec_trans, 'Delaunay', fail)
        assert vec_trans._interpolate_fields(self.scattered, (),
                                             self.xi) == ()

    def test_single_triangulation(self, monkeypatch):
        calls = []
        delaunay = vec_trans.Delaunay

        def counting_delaunay(points):
            calls.append(points.shape)
            return delaunay(points)

        monkeypatch.setattr(vec_trans, 'Delaunay', counting_delaunay)
        x, y = _sample_plate_carree_coordinates()
        u, v = _sample_plate_carree_vector_field()
        s = _sample_plate_carree_scalar_field()
        vec_trans.vector_scalar_to_grid(ccrs.PlateCarree(),
                                        ccrs.NorthPolarStereo(), (5, 3),
                                        x, y, u, v, s, s)
        assert calls == [(6, 2)]


@pytest.mark.parametrize('target_crs', [
    ccrs.PlateCarree(),
    ccrs.LambertConformal(central_longitude=-10),
    ccrs.Robinson(central_longitude=180),
], ids=['same_crs', 'lambert', 'robinson_180'])
def test_vector_scalar_to_grid_unchanged(target_crs, monkeypatch):
    # Regridding a regular lat/lon grid with several fields, masked
    # values and NaNs must give exactly the same result as interpolating
    # each field separately.
    lon = np.linspace(-80, 40, 61)
    lat = np.linspace(85, 20, 34)
    x, y = np.meshgrid(lon, lat)
    u = np.ma.masked_less(10 * np.sin(np.deg2rad(4 * x)) *
                          np.cos(np.deg2rad(4 * y)), -8)
    v = 5 * np.sin(np.deg2rad(2 * x + 3 * y))
    v[11:14, 28:31] = np.nan
    speed = np.hypot(u, v)
    src_crs = ccrs.PlateCarree()
    extent = (-50, 15, 50, 70)
    if target_crs != src_crs:
        extent = target_crs.transform_points(
            src_crs, np.array(extent[:2]), np.array(extent[2:]))
        extent = (*extent[:, 0], *extent[:, 1])

    def regrid():
        return vec_trans.vector_scalar_to_grid(
            src_crs, target_crs, (30, 20), lon, lat, u, v, speed, -speed,
            target_extent=extent)

    actual = regrid()
    monkeypatch.setattr(vec_trans, '_interpolate_fields',
                        _griddata_per_field)
    expected = regrid()
    _assert_identical(actual, expected)
    assert np.isfinite(actual[2]).any()
