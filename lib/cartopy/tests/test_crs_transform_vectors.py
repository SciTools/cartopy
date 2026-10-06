# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

import cartopy.crs as ccrs


class TestTransformVectors:

    def test_transform(self):
        # Test some simple vectors to make sure they are transformed
        # correctly.
        rlons = np.array([-90., 0, 90., 180.])
        rlats = np.array([0., 0., 0., 0.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        # transform grid eastward vectors
        ut, vt = target_proj.transform_vectors(src_proj,
                                               rlons,
                                               rlats,
                                               np.ones([4]),
                                               np.zeros([4]))
        assert_array_almost_equal(ut, np.array([0, 1, 0, -1]), decimal=2)
        assert_array_almost_equal(vt, np.array([-1, 0, 1, 0]), decimal=2)
        # transform grid northward vectors
        ut, vt = target_proj.transform_vectors(src_proj,
                                               rlons,
                                               rlats,
                                               np.zeros([4]),
                                               np.ones([4]))
        assert_array_almost_equal(ut, np.array([1, 0, -1, 0]), decimal=2)
        assert_array_almost_equal(vt, np.array([0, 1, 0, -1]), decimal=2)
        # transform grid north-eastward vectors
        ut, vt = target_proj.transform_vectors(src_proj,
                                               rlons,
                                               rlats,
                                               np.ones([4]),
                                               np.ones([4]))
        assert_array_almost_equal(ut, np.array([1, 1, -1, -1]), decimal=2)
        assert_array_almost_equal(vt, np.array([-1, 1, 1, -1]), decimal=2)

    def test_transform_and_inverse(self):
        # Check a full circle transform back to the native projection.
        x = np.arange(-60, 42.5, 2.5)
        y = np.arange(30, 72.5, 2.5)
        x2d, y2d = np.meshgrid(x, y)
        u = np.cos(np.deg2rad(y2d))
        v = np.cos(2. * np.deg2rad(x2d))
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        proj_xyz = target_proj.transform_points(src_proj, x2d, y2d)
        xt, yt = proj_xyz[..., 0], proj_xyz[..., 1]
        ut, vt = target_proj.transform_vectors(src_proj, x2d, y2d, u, v)
        utt, vtt = src_proj.transform_vectors(target_proj, xt, yt, ut, vt)
        assert_array_almost_equal(u, utt, decimal=4)
        assert_array_almost_equal(v, vtt, decimal=4)

    def test_invalid_input_domain(self):
        # If an input coordinate is outside the input projection domain
        # we should be able to handle it correctly.
        rlon = np.array([270.])
        rlat = np.array([0.])
        u = np.array([1.])
        v = np.array([0.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([0]), decimal=2)
        assert_array_almost_equal(vt, np.array([-1]), decimal=2)

    def test_invalid_x_domain(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection x-domain it should be handled correctly as long as
        # it is not a corner point.
        rlon = np.array([180.])
        rlat = np.array([0.])
        u = np.array([1.])
        v = np.array([0.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([-1]), decimal=2)
        assert_array_almost_equal(vt, np.array([0.]), decimal=2)

    def test_invalid_y_domain(self):
        # If the point we need to calculate the vector angle falls outside the
        # source projection y-domain it should be handled correctly as long as
        # it is not a corner point.
        rlon = np.array([0.])
        rlat = np.array([90.])
        u = np.array([0.])
        v = np.array([1.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(src_proj, rlon, rlat, u, v)
        assert_array_almost_equal(ut, np.array([0.]), decimal=2)
        assert_array_almost_equal(vt, np.array([1.]), decimal=2)

    @pytest.mark.parametrize('u, v, expected', [
        pytest.param(1., 1., [-1., -1.], id='xy corner'),
        pytest.param(1., -1., [-1., 1.], id='x corner'),
        pytest.param(-1., 1., [1., -1.], id='y corner')])
    def test_domain_corner(self, u, v, expected):
        # (180, 90) is on both the antimeridian and the pole. At the pole
        # itself, longitude is degenerate: only the latitude direction is
        # well defined. That is still enough to get a correct answer, so
        # no warning is expected. "North" here is just "north" at the
        # prime meridian, flipped 180 degrees, since the target is
        # centred on the pole.
        rlon = np.array([180.])
        rlat = np.array([90.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Stereographic(central_latitude=90,
                                         central_longitude=0)
        ut, vt = target_proj.transform_vectors(
            src_proj, rlon, rlat, np.array([u]), np.array([v]))
        assert_array_almost_equal([ut[0], vt[0]], expected, decimal=2)

    def test_domain_corner_unrecoverable(self):
        # An Orthographic projection only shows one hemisphere, so its own
        # domain edge sits at the antipode of its centre. A point there has
        # no direction in which a probe stays inside that domain, so
        # neither Jacobian column can be estimated and a warning is still
        # expected.
        rlon = np.array([180.])
        rlat = np.array([0.])
        u = np.array([1.])
        v = np.array([1.])
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.Orthographic()
        with pytest.warns(UserWarning, match='source domain corners'):
            target_proj.transform_vectors(src_proj, rlon, rlat, u, v)

    def test_transform_linear(self):
        # Vector transforms should be linear in the input components:
        # transforming (1, 0) + (0, 1) should give the same answer as
        # transforming (1, 1) directly.
        src_proj = ccrs.RotatedPole(pole_longitude=180, pole_latitude=45.0)
        target_proj = ccrs.PlateCarree()
        x = np.array([0.])
        y = np.array([0.])
        zero = np.array([0.])
        one = np.array([1.])

        u10, v10 = target_proj.transform_vectors(src_proj, x, y, one, zero)
        u01, v01 = target_proj.transform_vectors(src_proj, x, y, zero, one)
        u11, v11 = target_proj.transform_vectors(src_proj, x, y, one, one)

        assert_array_almost_equal(u11, u10 + u01, decimal=8)
        assert_array_almost_equal(v11, v10 + v01, decimal=8)

    def test_transform_plate_carree_near_pole(self):
        # Regression test for issue 1179. Close to the pole, a
        # PlateCarree vector transformed into a projection that behaves
        # normally there (NorthPolarStereo) should come back out almost
        # unchanged. https://github.com/SciTools/cartopy/issues/1179
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.NorthPolarStereo()
        lon = np.array([0.])
        lat = np.array([89.])
        u = np.array([-3.])
        v = np.array([0.1])

        ut, vt = target_proj.transform_vectors(src_proj, lon, lat, u, v)

        assert_array_almost_equal(ut, u, decimal=3)
        assert_array_almost_equal(vt, v, decimal=3)

    def test_transform_analytic_north_polar_stereo(self):
        # On a NorthPolarStereo centred at longitude 0, "east" at
        # longitude lon should point exactly lon degrees around from the
        # target's x-axis, at any latitude. This checks that directly,
        # instead of just checking the result is self consistent.
        src_proj = ccrs.PlateCarree()
        target_proj = ccrs.NorthPolarStereo()
        lons = np.array([0., 30., 60., 90., 135., 180., -90., -45.])
        for lat in (10., 45., 85.):
            lats = np.full_like(lons, lat)
            ut, vt = target_proj.transform_vectors(
                src_proj, lons, lats, np.ones_like(lons), np.zeros_like(lons))
            angle = np.degrees(np.arctan2(vt, ut))
            assert_array_almost_equal(
                (angle - lons + 180) % 360 - 180, np.zeros_like(lons),
                decimal=3)

    def test_transform_target_domain_cut(self):
        # Regression test: the source's own domain edge is not the only
        # place a finite-difference probe can cross a coordinate cut.
        # Here the target's cut (its antimeridian) falls in the middle
        # of the source's domain, since the two have different central
        # longitudes. A vector close to it should keep a magnitude of 1
        # and a rotation close to the true local shear (about 12
        # degrees here), not 0 or 90 degrees when this isn't handled correctly.
        src_proj = ccrs.PlateCarree(central_longitude=180)
        target_proj = ccrs.Robinson(central_longitude=0)
        offsets = np.array([-0.5, -1e-2, -1e-4, -1e-6, -1e-8,
                            1e-8, 1e-6, 1e-4, 1e-2, 0.5])
        lons = 180. + offsets
        lats = np.full_like(lons, 30.)
        native = src_proj.transform_points(ccrs.PlateCarree(), lons, lats)
        u = np.ones_like(lons)
        v = np.zeros_like(lons)

        ut, vt = target_proj.transform_vectors(
            src_proj, native[..., 0], native[..., 1], u, v)

        assert_array_almost_equal(np.hypot(ut, vt), np.ones_like(lons),
                                  decimal=6)
        angle = np.abs(np.degrees(np.arctan2(vt, ut)))
        assert np.all((angle > 11.) & (angle < 13.))

    @pytest.mark.parametrize('target, exact_direction', [
        pytest.param(ccrs.NorthPolarStereo(), True, id='conformal'),
        pytest.param(ccrs.Sinusoidal(), False, id='non-conformal')])
    def test_transform_roundtrip_shear(self, target, exact_direction):
        # Magnitude round-trips exactly for any target, since a rotation
        # never changes a vector's length. Direction only round-trips
        # exactly for a conformal target. A non-conformal target has
        # local shear that a pure rotation cannot represent, so some
        # direction is lost. See the note on transform_vectors.
        src = ccrs.PlateCarree()
        lon, lat = np.array([90.]), np.array([65.])
        u, v = np.array([1.]), np.array([2.])

        xyz = target.transform_points(src, lon, lat)
        ut, vt = target.transform_vectors(src, lon, lat, u, v)
        urt, vrt = src.transform_vectors(target, xyz[..., 0], xyz[..., 1],
                                         ut, vt)

        assert_array_almost_equal(np.hypot(urt, vrt), np.hypot(u, v),
                                  decimal=10)
        if exact_direction:
            assert_array_almost_equal([urt[0], vrt[0]], [u[0], v[0]],
                                      decimal=8)
        else:
            assert not np.allclose([urt[0], vrt[0]], [u[0], v[0]], atol=1e-6)
