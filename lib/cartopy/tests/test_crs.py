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

import numpy

import cartopy.crs as ccrs


class TestCRS(unittest.TestCase):

    def test_osgb(self):
        osgb = ccrs.OSGB()
        ll = ccrs.Geodetic()
        
        # results obtained by streetmap.co.uk.
        lat, lon = numpy.array([50.462023, -3.478831], dtype=numpy.double)
        east, north = numpy.array([295131, 63511], dtype=numpy.double)
        
        # note the handling of precision here...
        numpy.testing.assert_array_almost_equal(numpy.array(osgb.transform_point(lon, lat, ll)), numpy.array([east, north]), 1)
        numpy.testing.assert_array_almost_equal(ll.transform_point(east, north, osgb), [lon, lat], 2)
        
        r_lon, r_lat = ll.transform_point(east, north, osgb)
        r_inverted = numpy.array(osgb.transform_point(r_lon, r_lat, ll))
        numpy.testing.assert_array_almost_equal(r_inverted, [east, north], 3)
        
        r_east, r_north = osgb.transform_point(lon, lat, ll)
        r_inverted = numpy.array(ll.transform_point(r_east, r_north, osgb))
        numpy.testing.assert_array_almost_equal(r_inverted, [lon, lat])


if __name__ == '__main__':
    unittest.main()