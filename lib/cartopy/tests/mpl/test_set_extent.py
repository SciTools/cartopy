# (C) British Crown Copyright 2011 - 2016, Met Office
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

from __future__ import (absolute_import, division, print_function)

import tempfile

from matplotlib.testing.decorators import cleanup
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

import cartopy.crs as ccrs


@cleanup
def test_extents():
    # tests that one can set the extents of a map in a variety of coordinate
    # systems, for a variety of projection
    uk = [-12.5, 4, 49, 60]
    uk_crs = ccrs.Geodetic()

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(uk, crs=uk_crs)
    # enable to see what is going on (and to make sure it is a plot of the uk)
    # ax.coastlines()
    assert_array_almost_equal(ax.viewLim.get_points(),
                              np.array([[-12.5, 49.], [4., 60.]]))

    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(uk, crs=uk_crs)
    # enable to see what is going on (and to make sure it is a plot of the uk)
    # ax.coastlines()
    assert_array_almost_equal(ax.viewLim.get_points(),
                              np.array([[-1034046.22566261, -4765889.76601514],
                                        [333263.47741164, -3345219.0594531]])
                              )

    # given that we know the PolarStereo coordinates of the UK, try using
    # those in a PlateCarree plot
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-1034046, 333263, -4765889, -3345219],
                  crs=ccrs.NorthPolarStereo())
    # enable to see what is going on (and to make sure it is a plot of the uk)
    # ax.coastlines()
    assert_array_almost_equal(ax.viewLim.get_points(),
                              np.array([[-17.17698577, 48.21879707],
                                        [5.68924381, 60.54218893]])
                              )


@cleanup
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

    ax = plt.axes(projection=ccrs.OSGB())
    ax.set_extent((0, 7e5, 0, 13e5), ccrs.OSGB())
    assert_array_equal(ax.viewLim.get_points(), [[0, 0], [7e5, 13e5]])


def test_update_lim():
    # check that the standard data lim setting works
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.update_datalim([(-10, -10), (-5, -5)])
    assert_array_almost_equal(ax.dataLim.get_points(),
                              np.array([[-10., -10.], [-5., -5.]]))
    plt.close()


def test_limits_contour():
    xs, ys = np.meshgrid(np.linspace(250, 350, 15), np.linspace(-45, 45, 20))
    data = np.sin((xs * ys) * 1.e7)

    resulting_extent = np.array([[250 - 180, -45.], [-10. + 180, 45.]])

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.contourf(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)
    plt.close()

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.contour(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)
    plt.close()


def test_limits_pcolor():
    xs, ys = np.meshgrid(np.linspace(250, 350, 15), np.linspace(-45, 45, 20))
    data = (np.sin((xs * ys) * 1.e7))[:-1, :-1]

    resulting_extent = np.array([[250 - 180, -45.], [-10. + 180, 45.]])

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.pcolor(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)
    plt.close()

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    plt.pcolormesh(xs, ys, data, transform=ccrs.PlateCarree(180))
    assert_array_almost_equal(ax.dataLim, resulting_extent)
    plt.close()


def test_view_lim_autoscaling():
    x = np.linspace(0.12910209, 0.42141822)
    y = np.linspace(0.03739792, 0.33029076)
    x, y = np.meshgrid(x, y)
    ax = plt.axes(projection=ccrs.RotatedPole(37.5, 357.5))
    plt.scatter(x, y, x * y, transform=ccrs.PlateCarree())

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
    plt.close()


def test_view_lim_default_global():
    ax = plt.axes(projection=ccrs.PlateCarree())
    # The view lim should be the default unit bbox until it is drawn.
    assert_array_almost_equal(ax.viewLim.frozen().get_points(),
                              [[0, 0], [1, 1]])
    with tempfile.TemporaryFile() as tmp:
        plt.savefig(tmp)
    expected = np.array([[-180, -90], [180, 90]])
    assert_array_almost_equal(ax.viewLim.frozen().get_points(),
                              expected)
    plt.close()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
