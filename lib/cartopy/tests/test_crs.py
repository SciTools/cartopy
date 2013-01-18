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


from io import BytesIO
import pickle
import unittest

import numpy
from numpy.testing import assert_array_almost_equal as assert_arr_almost_eq


import cartopy.crs as ccrs


class TestCRS(unittest.TestCase):
    def test_hash(self):
        stereo = ccrs.Stereographic(90)
        north = ccrs.NorthPolarStereo()
        self.assertEqual(stereo, north)
        self.assertFalse(stereo != north)
        self.assertNotEqual(hash(stereo), hash(north))

        self.assertEqual(ccrs.Geodetic(), ccrs.Geodetic())

    def test_osgb(self):
        osgb = ccrs.OSGB()
        ll = ccrs.Geodetic()

        # results obtained by streetmap.co.uk.
        lat, lon = numpy.array([50.462023, -3.478831], dtype=numpy.double)
        east, north = numpy.array([295131, 63511], dtype=numpy.double)

        # note the handling of precision here...
        assert_arr_almost_eq(numpy.array(osgb.transform_point(lon, lat, ll)),
                             numpy.array([east, north]),
                             1)
        assert_arr_almost_eq(ll.transform_point(east, north, osgb),
                             [lon, lat],
                             2)

        r_lon, r_lat = ll.transform_point(east, north, osgb)
        r_inverted = numpy.array(osgb.transform_point(r_lon, r_lat, ll))
        assert_arr_almost_eq(r_inverted, [east, north], 3)

        r_east, r_north = osgb.transform_point(lon, lat, ll)
        r_inverted = numpy.array(ll.transform_point(r_east, r_north, osgb))
        assert_arr_almost_eq(r_inverted, [lon, lat])

    def test_europp(self):
        europp = ccrs.EuroPP()
        proj4_init = europp.proj4_init
        # Transverse Mercator, UTM zone 32,
        self.assertTrue('+proj=tmerc' in proj4_init)
        self.assertTrue('+zone=32' in proj4_init)
        # International 1924 ellipsoid.
        self.assertTrue('+ellps=intl' in proj4_init)
        # Scale factor on central meridian
        self.assertTrue('+k=0.9996' in proj4_init)
        # True origin Standard latitude and longitude
        self.assertTrue('+lat_0=50' in proj4_init)
        self.assertTrue('+lon_0=9' in proj4_init)
        # Map co-ordinates of true origin (metres)
        self.assertTrue('+x_0=1750000' in proj4_init)
        self.assertTrue('+y_0=1500000' in proj4_init)


def test_pickle():
    # check that we can pickle a simple CRS
    fh = BytesIO()
    pickle.dump(ccrs.PlateCarree(), fh)
    fh.seek(0)
    pc = pickle.load(fh)
    assert pc == ccrs.PlateCarree()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
