# (C) British Crown Copyright 2011 - 2012, Met Office
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

from matplotlib.testing.decorators import cleanup
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_array_almost_equal

import cartopy.crs as ccrs
from cartopy.tests.mpl import ImageTesting

@cleanup
def test_extents():
    # tests that one can set the extents of a map in a variety of coordinate systems, for a variety
    # of projection
    uk = [-12.5, 4, 49, 60]
    uk_crs = ccrs.Geodetic()
    
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(uk, crs=uk_crs)
#    ax.coastlines() # <- enable to see what is going on (and to make sure it is a plot of the uk)
    assert_array_almost_equal(ax.viewLim.get_points(), np.array([[-12.5, 49.], [4. , 60.]]))
    
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent(uk, crs=uk_crs)
#    ax.coastlines() # <- enable to see what is going on (and to make sure it is a plot of the uk)
    assert_array_almost_equal(ax.viewLim.get_points(), 
                              np.array([[-1034046.22566261, -4765889.76601514],
                                        [  333263.47741164, -3345219.0594531 ]])
                              )
    
    # given that we know the PolarStereo coordinates of the UK, try using those in a PlateCarree plot
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-1034046, 333263, -4765889, -3345219], crs=ccrs.NorthPolarStereo())
#    ax.coastlines() # <- enable to see what is going on (and to make sure it is a plot of the uk)
    assert_array_almost_equal(ax.viewLim.get_points(), 
                              np.array([[-17.17698577,  48.21879707],
                                        [  5.68924381, 60.54218893]])
                              )

@ImageTesting(['set_extent_wrapping'])
def test_wrapping():
    # Tests that set_extent() handles longitudes in 0 to 360 format
    # and that an extent that crosses the boundary is handled appropriately.
    fig = plt.figure(figsize=(10, 6))

    # Extent of Australia region in 0 to 360.
    extent_0_to_360 = (85.0, 220.0, -55.0, 20.0)

    # In PlateCarree(180) Australia region is central and
    # should not wrap. 0 to 360 conversion should be handled implicitly.
    ax = fig.add_subplot(2, 2, 1, projection=ccrs.PlateCarree(180))
    ax.coastlines()
    ax.set_extent(extent_0_to_360, ccrs.PlateCarree())

    # In PlateCarree() the region crosses the boundary so the result
    # should be a -180 to 180 strip covering -55 to 20 deg lat.
    ax = fig.add_subplot(2, 2, 2, projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent(extent_0_to_360, ccrs.PlateCarree())

    # The extent when expressed in the projection of the axes
    # should not require the crs to be specified.
    extent = (-95.0, 40.0, -55.0, 20.0)
    ax = fig.add_subplot(2, 2, 3, projection=ccrs.PlateCarree(180))
    ax.coastlines()
    ax.set_extent(extent)

    # If we do set the crs explicitly we should still get the same result as
    # the previous plot.
    ax = fig.add_subplot(2, 2, 4, projection=ccrs.PlateCarree(180))
    ax.coastlines()
    ax.set_extent(extent, ccrs.PlateCarree(180))

def test_update_lim():
    # check that the standard data lim setting works
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.update_datalim([(-10, -10), (-5, -5)])
    assert_array_almost_equal(ax.dataLim.get_points(), np.array([[-10., -10.], [-5., -5.]]))
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


if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
