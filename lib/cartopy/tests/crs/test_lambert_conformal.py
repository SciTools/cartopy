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


import unittest

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
from cartopy.tests.mpl import ImageTesting


class TestLambertConformal(unittest.TestCase):

    @ImageTesting(['lambert_conformal_default'])
    def test_default(self):
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.coastlines(resolution='110m')
        ax.gridlines()

    @ImageTesting(['lambert_conformal_cutoff'])
    def test_lambert_conformal_cutoff(self):
        ax = plt.axes(projection=ccrs.LambertConformal(cutoff=-80))
        ax.coastlines(resolution='110m')
        ax.gridlines()

    @ImageTesting(['lambert_conformal_inspire'])
    def test_lambert_inspire(self):
        # EPSG Projection 3034 - ETRS89 / ETRS-LCC
        # TODO: 1) Find a reference image.
        # TODO: 2) Get +ellps=GRS80 into the proj4 string somehow.
        #          We currently use WGS84, which is not correct.
        epsg3034 = ccrs.LambertConformal(central_longitude=10,
                                         secant_latitudes=(35, 65),
                                         central_latitude=52,
                                         false_easting=4000000,
                                         false_northing=2800000)

        expects = ["+lat_0=52", " +lat_1=35", "+lat_2=65", "+lon_0=10",
                   "+proj=lcc", "+x_0=4000000", "+y_0=2800000"]  # +ellps=GRS80
        for i in expects:
            self.assert_(i in epsg3034.proj4_init, "expected {}".format(i))

        ax = plt.axes(projection=epsg3034)
        ax.coastlines(resolution='110m')
        ax.gridlines()

    @ImageTesting(['lambert_conformal_south'])
    def test_lambert_south(self):
        # Reference image: http://www.icsm.gov.au/mapping/map_projections.html
        ax = plt.axes(
            projection=ccrs.LambertConformal(central_longitude=140,
                                             secant_latitudes=(-30, -60,),
                                             cutoff=65))
        ax.coastlines(resolution='110m')
        ax.gridlines()

    @ImageTesting(['lambert_conformal_oz'])
    def test_lambert_oz(self):
        # Reference image: http://www.icsm.gov.au/mapping/map_projections.html
        ax = plt.axes(
            projection=ccrs.LambertConformal(central_longitude=140,
                                             secant_latitudes=(-18, -36,),
                                             cutoff=65))
        ax.set_extent((90, 190, -45, 0))
        ax.coastlines(resolution='110m')
        ax.gridlines()

    @ImageTesting(['lambert_conformal_tangental'])
    def test_lambert_tangental(self):
        ax = plt.axes(
            projection=ccrs.LambertConformal(central_longitude=0,
                                             secant_latitudes=(45, 45)))
        ax.coastlines(resolution='110m')
        ax.gridlines()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
