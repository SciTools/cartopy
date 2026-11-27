# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import pytest
import shapely

import cartopy.crs as ccrs


def test_extents():
    # tests that one can set the extents of a map in a variety of coordinate
    # systems, for a variety of projection
    uk = [-12.5, 4, 49, 60]
    uk_crs = ccrs.Geodetic()

    ax = plt.axes(projection=ccrs.PlateCarree(), label='pc')
    ax.set_extent(uk, crs=uk_crs)
    # enable to see what is going on (and to make sure it is a plot of the uk)
    # ax.coastlines()
    assert_array_almost_equal(ax.viewLim.get_points(),
                              np.array([[-12.5, 49.], [4., 60.]]))

    ax = plt.axes(projection=ccrs.NorthPolarStereo(), label='npstere')
    ax.set_extent(uk, crs=uk_crs)
    # enable to see what is going on (and to make sure it is a plot of the uk)
    # ax.coastlines()
    assert_array_almost_equal(ax.viewLim.get_points(),
                              np.array([[-1034046.22566261, -4765889.76601514],
                                        [333263.47741164, -3345219.0594531]])
                              )

    # given that we know the PolarStereo coordinates of the UK, try using
    # those in a PlateCarree plot
    ax = plt.axes(projection=ccrs.PlateCarree(), label='pc')
    ax.set_extent([-1034046, 333263, -4765889, -3345219],
                  crs=ccrs.NorthPolarStereo())
    # enable to see what is going on (and to make sure it is a plot of the uk)
    # ax.coastlines()
    assert_array_almost_equal(ax.viewLim.get_points(),
                              np.array([[-17.17698577, 48.21879707],
                                        [5.68924381, 60.54218893]])
                              )


def test_get_extent():
    # tests that getting the extents of a map produces something reasonable.
    uk = [-12.5, 4, 49, 60]
    uk_crs = ccrs.PlateCarree()

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(uk, crs=uk_crs)
    assert_array_almost_equal(ax.get_extent(uk_crs), uk)

    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent(uk, crs=uk_crs)
    assert_array_almost_equal(ax.get_extent(uk_crs), uk)

    ax = plt.axes(projection=ccrs.Mercator(min_latitude=uk[2],
                                           max_latitude=uk[3]))
    ax.set_extent(uk, crs=uk_crs)
    assert_array_almost_equal(ax.get_extent(uk_crs), uk, decimal=1)


@pytest.mark.parametrize('projection', [
    ccrs.Orthographic(central_longitude=0, central_latitude=0),
    ccrs.Geostationary(central_longitude=0),
])
def test_get_extent_near_boundary(projection):
    # Close to the horizon, a narrow strip of the map covers a large part of
    # the globe. It must not be dropped from the extent.
    ax = plt.axes(projection=projection)
    x0, x1 = projection.x_limits
    y0, y1 = projection.y_limits
    ax.set_xlim(0, x1 / 2)
    ax.set_ylim(y1 * 0.9, y1)
    lon0, lon1, lat0, lat1 = ax.get_extent(ccrs.PlateCarree())

    # All visible points of the map have to be within the extent.
    x, y = np.meshgrid(np.linspace(0, x1 / 2, 100),
                       np.linspace(y1 * 0.9, y1, 100))
    visible = shapely.contains_xy(shapely.Polygon(projection.boundary), x, y)
    lonlat = ccrs.PlateCarree().transform_points(
        projection, x[visible], y[visible])
    assert np.all(lonlat[:, 0] >= lon0 - 0.1)
    assert np.all(lonlat[:, 0] <= lon1 + 0.1)
    assert np.all(lonlat[:, 1] >= lat0 - 0.1)
    assert np.all(lonlat[:, 1] <= lat1 + 0.1)


def test_get_extent_wrapped_global():
    # The left and right edges of the map lie on the same meridian. The
    # boundary is eroded slightly when calculating the extent, as otherwise
    # the extent would collapse to that meridian.
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_global()
    np.testing.assert_allclose(ax.get_extent(ccrs.PlateCarree()),
                               [-180, 180, -90, 90], atol=1)


def test_domain_extents():
    # Setting the extent to global or the domain limits.
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent((-180, 180, -90, 90))
    assert_array_equal(ax.viewLim.get_points(), [[-180, -90], [180, 90]])
    ax.set_extent((-180, 180, -90, 90), ccrs.PlateCarree())
    assert_array_equal(ax.viewLim.get_points(), [[-180, -90], [180, 90]])

    ax = plt.axes(projection=ccrs.PlateCarree(90))
    ax.set_extent((-180, 180, -90, 90))
    assert_array_equal(ax.viewLim.get_points(), [[-180, -90], [180, 90]])
    ax.set_extent((-180, 180, -90, 90), ccrs.PlateCarree(90))
    assert_array_equal(ax.viewLim.get_points(), [[-180, -90], [180, 90]])

    ax = plt.axes(projection=ccrs.OSGB(approx=False))
    ax.set_extent((0, 7e5, 0, 13e5), ccrs.OSGB(approx=False))
    assert_array_equal(ax.viewLim.get_points(), [[0, 0], [7e5, 13e5]])


def test_update_lim():
    # check that the standard data lim setting works
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.update_datalim([(-10, -10), (-5, -5)])
    assert_array_almost_equal(ax.dataLim.get_points(),
                              np.array([[-10., -10.], [-5., -5.]]))


def test_limits_contour():
    xs, ys = np.meshgrid(np.linspace(250, 350, 15), np.linspace(-45, 45, 20))
    data = np.sin((xs * ys) * 1.e7)

    resulting_extent = np.array([[250 - 180, -45.], [-10. + 180, 45.]])

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.contourf(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)

    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.contour(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)


def test_limits_pcolor():
    xs, ys = np.meshgrid(np.linspace(250, 350, 15), np.linspace(-45, 45, 20))
    data = (np.sin((xs * ys) * 1.e7))[:-1, :-1]

    resulting_extent = np.array([[250 - 180, -45.], [-10. + 180, 45.]])

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.pcolor(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)

    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.pcolormesh(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)


def test_view_lim_autoscaling():
    x = np.linspace(0.12910209, 0.42141822)
    y = np.linspace(0.03739792, 0.33029076)
    x, y = np.meshgrid(x, y)
    ax = plt.axes(projection=ccrs.RotatedPole(37.5, 357.5))
    ax.scatter(x, y, x * y, transform=ccrs.PlateCarree())

    expected = np.array([[86.12433701, 52.51570463],
                         [86.69696603, 52.86372057]])

    assert_array_almost_equal(ax.viewLim.frozen().get_points(), expected,
                              decimal=2)
    plt.draw()
    assert_array_almost_equal(ax.viewLim.frozen().get_points(), expected,
                              decimal=2)
    ax.autoscale_view(tight=False)
    expected_non_tight = np.array([[86, 52.45], [86.8, 52.9]])
    assert_array_almost_equal(ax.viewLim.frozen().get_points(),
                              expected_non_tight, decimal=1)


def test_view_lim_default_global():
    ax = plt.axes(projection=ccrs.PlateCarree())
    # viewLim is set to projection bounds in __clear.
    expected = np.array([[-180, -90], [180, 90]])
    assert_array_almost_equal(ax.viewLim.frozen().get_points(),
                              expected)
    plt.gcf().draw_without_rendering()
    assert_array_almost_equal(ax.viewLim.frozen().get_points(),
                              expected)
