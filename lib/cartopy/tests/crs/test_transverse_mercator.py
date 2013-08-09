# (C) British Crown Copyright 2013, Met Office
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


class TestTransverseMercator(unittest.TestCase):
    @ImageTesting(['tmerc_default'])
    def test_default(self):
        ax = plt.axes(projection=ccrs.TransverseMercator())
        ax.coastlines()
        ax.gridlines()

    @ImageTesting(['osgb'])
    def test_osgb_vals(self):
        proj = ccrs.TransverseMercator(central_longitude=-2,
                                       central_latitude=49,
                                       scale_factor=0.9996012717,
                                       false_easting=400000,
                                       false_northing=-100000,
                                       globe=ccrs.Globe(datum='OSGB36',
                                                        ellipse='airy'))
        ax = plt.axes(projection=proj)
        ax.set_xlim(0, 7e5)
        ax.set_ylim(0, 13e5)
        ax.coastlines()
        ax.gridlines()


class TestOSGB(unittest.TestCase):
    @ImageTesting(['osgb'])
    def test_default(self):
        ax = plt.axes(projection=ccrs.OSGB())
        ax.coastlines()
        ax.gridlines()


class TestOSNI(unittest.TestCase):
    @ImageTesting(['osni'])
    def test_default(self):
        ax = plt.axes(projection=ccrs.OSNI())
        ax.coastlines()
        ax.gridlines()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
