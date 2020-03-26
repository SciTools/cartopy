# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from __future__ import (absolute_import, division, print_function)

import copy
from io import BytesIO
import pickle

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_equal
from numpy.testing import assert_array_almost_equal as assert_arr_almost_eq
try:
    import pyepsg
except ImportError:
    pyepsg = None
import pytest
import shapely.geometry as sgeom

import cartopy.crs as ccrs


class TestCRS(object):
    def test_hash(self):
        stereo = ccrs.Stereographic(90)
        north = ccrs.NorthPolarStereo()
        assert stereo == north
        assert not stereo != north
        assert hash(stereo) == hash(north)

        assert ccrs.Geodetic() == ccrs.Geodetic()

    @pytest.mark.parametrize('approx', [True, False])
    def test_osni(self, approx):
        osni = ccrs.OSNI(approx=approx)
        ll = ccrs.Geodetic()

        # results obtained by nearby.org.uk.
        lat, lon = np.array([54.5622169298669, -5.54159863617957],
                            dtype=np.double)
        east, north = np.array([359000, 371000], dtype=np.double)

        assert_arr_almost_eq(osni.transform_point(lon, lat, ll),
                             np.array([east, north]),
                             -1)
        assert_arr_almost_eq(ll.transform_point(east, north, osni),
                             np.array([lon, lat]),
                             3)

    def _check_osgb(self, osgb):
        ll = ccrs.Geodetic()

        # results obtained by streetmap.co.uk.
        lat, lon = np.array([50.462023, -3.478831], dtype=np.double)
        east, north = np.array([295131, 63511], dtype=np.double)

        # note the handling of precision here...
        assert_arr_almost_eq(np.array(osgb.transform_point(lon, lat, ll)),
                             np.array([east, north]),
                             1)
        assert_arr_almost_eq(ll.transform_point(east, north, osgb),
                             [lon, lat],
                             2)

        r_lon, r_lat = ll.transform_point(east, north, osgb)
        r_inverted = np.array(osgb.transform_point(r_lon, r_lat, ll))
        assert_arr_almost_eq(r_inverted, [east, north], 3)

        r_east, r_north = osgb.transform_point(lon, lat, ll)
        r_inverted = np.array(ll.transform_point(r_east, r_north, osgb))
        assert_arr_almost_eq(r_inverted, [lon, lat])

    @pytest.mark.parametrize('approx', [True, False])
    def test_osgb(self, approx):
        self._check_osgb(ccrs.OSGB(approx=approx))

    @pytest.mark.network
    @pytest.mark.skipif(pyepsg is None, reason='requires pyepsg')
    def test_epsg(self):
        uk = ccrs.epsg(27700)
        assert uk.epsg_code == 27700
        assert_almost_equal(uk.x_limits, (-118365.7406176, 751581.5647514),
                            decimal=3)
        assert_almost_equal(uk.y_limits, (-5268.1704980, 1272227.7987656),
                            decimal=2)
        assert_almost_equal(uk.threshold, 8699.47, decimal=2)
        self._check_osgb(uk)

    @pytest.mark.network
    @pytest.mark.skipif(pyepsg is None, reason='requires pyepsg')
    def test_epsg_compound_crs(self):
        projection = ccrs.epsg(5973)
        assert projection.epsg_code == 5973

    def test_europp(self):
        europp = ccrs.EuroPP()
        proj4_init = europp.proj4_init
        # Transverse Mercator, UTM zone 32,
        assert '+proj=utm' in proj4_init
        assert '+zone=32' in proj4_init
        # International 1924 ellipsoid.
        assert '+ellps=intl' in proj4_init

    def test_transform_points_nD(self):
        rlons = np.array([[350., 352., 354.], [350., 352., 354.]])
        rlats = np.array([[-5., -0., 1.], [-4., -1., 0.]])

        src_proj = ccrs.RotatedGeodetic(pole_longitude=178.0,
                                        pole_latitude=38.0)
        target_proj = ccrs.Geodetic()
        res = target_proj.transform_points(x=rlons, y=rlats,
                                           src_crs=src_proj)
        unrotated_lon = res[..., 0]
        unrotated_lat = res[..., 1]

        # Solutions derived by proj direct.
        solx = np.array([[-16.42176094, -14.85892262, -11.90627520],
                         [-16.71055023, -14.58434624, -11.68799988]])
        soly = np.array([[46.00724251, 51.29188893, 52.59101488],
                         [46.98728486, 50.30706042, 51.60004528]])
        assert_arr_almost_eq(unrotated_lon, solx)
        assert_arr_almost_eq(unrotated_lat, soly)

    def test_transform_points_1D(self):
        rlons = np.array([350., 352., 354., 356.])
        rlats = np.array([-5., -0., 5., 10.])

        src_proj = ccrs.RotatedGeodetic(pole_longitude=178.0,
                                        pole_latitude=38.0)
        target_proj = ccrs.Geodetic()
        res = target_proj.transform_points(x=rlons, y=rlats,
                                           src_crs=src_proj)
        unrotated_lon = res[..., 0]
        unrotated_lat = res[..., 1]

        # Solutions derived by proj direct.
        solx = np.array([-16.42176094, -14.85892262,
                         -12.88946157, -10.35078336])
        soly = np.array([46.00724251, 51.29188893,
                         56.55031485, 61.77015703])

        assert_arr_almost_eq(unrotated_lon, solx)
        assert_arr_almost_eq(unrotated_lat, soly)

    def test_transform_points_xyz(self):
        # Test geodetic transforms when using z value
        rx = np.array([2574.32516e3])
        ry = np.array([837.562e3])
        rz = np.array([5761.325e3])

        src_proj = ccrs.Geocentric()
        target_proj = ccrs.Geodetic()

        res = target_proj.transform_points(x=rx, y=ry, z=rz,
                                           src_crs=src_proj)

        glat = res[..., 0]
        glon = res[..., 1]
        galt = res[..., 2]

        # Solution generated by pyproj
        solx = np.array([18.0224043189])
        soly = np.array([64.9796515089])
        solz = np.array([5048.03893734])

        assert_arr_almost_eq(glat, solx)
        assert_arr_almost_eq(glon, soly)
        assert_arr_almost_eq(galt, solz)

    def test_globe(self):
        # Ensure the globe affects output.
        rugby_globe = ccrs.Globe(semimajor_axis=9000000,
                                 semiminor_axis=9000000,
                                 ellipse=None)
        footy_globe = ccrs.Globe(semimajor_axis=1000000,
                                 semiminor_axis=1000000,
                                 ellipse=None)

        rugby_moll = ccrs.Mollweide(globe=rugby_globe)
        footy_moll = ccrs.Mollweide(globe=footy_globe)

        rugby_pt = rugby_moll.transform_point(10, 10, ccrs.Geodetic())
        footy_pt = footy_moll.transform_point(10, 10, ccrs.Geodetic())

        assert_arr_almost_eq(rugby_pt, (1400915, 1741319), decimal=0)
        assert_arr_almost_eq(footy_pt, (155657, 193479), decimal=0)

    def test_project_point(self):
        point = sgeom.Point([0, 45])
        multi_point = sgeom.MultiPoint([point, sgeom.Point([180, 45])])

        pc = ccrs.PlateCarree()
        pc_rotated = ccrs.PlateCarree(central_longitude=180)

        result = pc_rotated.project_geometry(point, pc)
        assert_arr_almost_eq(result.xy, [[-180.], [45.]])

        result = pc_rotated.project_geometry(multi_point, pc)
        assert isinstance(result, sgeom.MultiPoint)
        assert len(result) == 2
        assert_arr_almost_eq(result[0].xy, [[-180.], [45.]])
        assert_arr_almost_eq(result[1].xy, [[0], [45.]])

    def test_utm(self):
        utm30n = ccrs.UTM(30)
        ll = ccrs.Geodetic()
        lat, lon = np.array([51.5, -3.0], dtype=np.double)
        east, north = np.array([500000, 5705429.2], dtype=np.double)
        assert_arr_almost_eq(utm30n.transform_point(lon, lat, ll),
                             [east, north],
                             decimal=1)
        assert_arr_almost_eq(ll.transform_point(east, north, utm30n),
                             [lon, lat],
                             decimal=1)
        utm38s = ccrs.UTM(38, southern_hemisphere=True)
        lat, lon = np.array([-18.92, 47.5], dtype=np.double)
        east, north = np.array([763316.7, 7906160.8], dtype=np.double)
        assert_arr_almost_eq(utm38s.transform_point(lon, lat, ll),
                             [east, north],
                             decimal=1)
        assert_arr_almost_eq(ll.transform_point(east, north, utm38s),
                             [lon, lat],
                             decimal=1)


@pytest.fixture(params=[
    [ccrs.PlateCarree, {}],
    [ccrs.PlateCarree, dict(
        central_longitude=1.23)],
    [ccrs.NorthPolarStereo, dict(
        central_longitude=42.5,
        globe=ccrs.Globe(ellipse="helmert"))],
])
def proj_to_copy(request):
    cls, kwargs = request.param
    return cls(**kwargs)


def test_pickle(proj_to_copy):
    # check that we can pickle a simple CRS
    fh = BytesIO()
    pickle.dump(proj_to_copy, fh)
    fh.seek(0)
    pickled_prj = pickle.load(fh)
    assert proj_to_copy == pickled_prj


def test_deepcopy(proj_to_copy):
    prj_cp = copy.deepcopy(proj_to_copy)
    assert proj_to_copy.proj4_params == prj_cp.proj4_params
    assert proj_to_copy == prj_cp


def test_PlateCarree_shortcut():
    central_lons = [[0, 0], [0, 180], [0, 10], [10, 0], [-180, 180], [
        180, -180]]

    target = [([[-180, -180], [-180, 180]], 0),
              ([[-180, 0], [0, 180]], 180),
              ([[-180, -170], [-170, 180]], 10),
              ([[-180, 170], [170, 180]], -10),
              ([[-180, 180], [180, 180]], 360),
              ([[-180, -180], [-180, 180]], -360),
              ]

    assert len(target) == len(central_lons)

    for expected, (s_lon0, t_lon0) in zip(target, central_lons):
        expected_bboxes, expected_offset = expected

        src = ccrs.PlateCarree(central_longitude=s_lon0)
        target = ccrs.PlateCarree(central_longitude=t_lon0)

        bbox, offset = src._bbox_and_offset(target)

        assert offset == expected_offset
        assert bbox == expected_bboxes


def test_transform_points_empty():
    """Test CRS.transform_points with empty array."""
    crs = ccrs.Stereographic()
    result = crs.transform_points(ccrs.PlateCarree(),
                                  np.array([]), np.array([]))
    assert_array_equal(result, np.array([], dtype=np.float64).reshape(0, 3))
