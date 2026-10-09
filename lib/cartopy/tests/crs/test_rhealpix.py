# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Tests for the RHEALPix projection.

"""
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest
import shapely

import cartopy.crs as ccrs
from .helpers import check_proj_params


def test_defaults():
    crs = ccrs.RHEALPix()
    expected = {'ellps=WGS84', 'lon_0=0', 'north_square=0', 'south_square=0'}
    check_proj_params('rhealpix', crs, expected)


def test_square_positions():
    crs = ccrs.RHEALPix(north_square=1, south_square=2)
    expected = {'ellps=WGS84', 'lon_0=0', 'north_square=1', 'south_square=2'}
    check_proj_params('rhealpix', crs, expected)


def test_invalid_square_positions():
    with pytest.raises(ValueError, match='north_square must be'):
        ccrs.RHEALPix(north_square=4)
    with pytest.raises(ValueError, match='south_square must be'):
        ccrs.RHEALPix(south_square=-1)


def test_central_longitude():
    crs = ccrs.RHEALPix(north_square=2, central_longitude=-124.8)
    expected = {'ellps=WGS84', 'lon_0=-124.8', 'north_square=2',
                'south_square=0'}
    check_proj_params('rhealpix', crs, expected)


def test_sphere_globe():
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=1000, ellipse=None)
    crs = ccrs.RHEALPix(globe=globe)
    expected = {'a=1000', 'b=1000', 'lon_0=0', 'north_square=0',
                'south_square=0'}
    check_proj_params('rhealpix', crs, expected)

    assert_almost_equal(np.array(crs.x_limits), [-3141.5927, 3141.5927],
                        decimal=4)
    assert_almost_equal(np.array(crs.y_limits), [-2356.1945, 2356.1945],
                        decimal=4)


def test_eccentric_globe():
    # The map is sized by the radius of the sphere with the same surface
    # area as this globe, which is 830.7145, and not by the semi-major axis.
    globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500, ellipse=None)
    crs = ccrs.RHEALPix(globe=globe)
    expected = {'a=1000', 'b=500', 'lon_0=0', 'north_square=0',
                'south_square=0'}
    check_proj_params('rhealpix', crs, expected)

    assert_almost_equal(np.array(crs.x_limits), [-2609.7664, 2609.7664],
                        decimal=4)
    assert_almost_equal(np.array(crs.y_limits), [-1957.3248, 1957.3248],
                        decimal=4)


@pytest.mark.parametrize('globe', [
    None,
    ccrs.Globe(ellipse='sphere'),
    ccrs.Globe(semimajor_axis=1000, semiminor_axis=500, ellipse=None),
    ccrs.Globe(semimajor_axis=1, semiminor_axis=0.9, ellipse=None),
])
def test_boundary_fits_globe(globe):
    # The boundary should hug the projected globe whatever its shape, see
    # https://github.com/SciTools/cartopy/pull/2458.
    crs = ccrs.RHEALPix(globe=globe, north_square=2, south_square=1)
    lons, lats = np.meshgrid(np.linspace(-180, 180, 181),
                             np.linspace(-90, 90, 181))
    points = crs.transform_points(crs.as_geodetic(), lons.ravel(), lats.ravel())
    boundary = shapely.Polygon(crs.boundary)

    # Nothing pokes out of the boundary, and the boundary is no bigger than
    # it needs to be.
    width = crs.x_limits[1] - crs.x_limits[0]
    assert boundary.distance(shapely.MultiPoint(points[:, :2])) == 0
    assert_allclose(crs.x_limits, [points[:, 0].min(), points[:, 0].max()],
                    atol=0.01 * width)
    assert_allclose(crs.y_limits, [points[:, 1].min(), points[:, 1].max()],
                    atol=0.01 * width)
