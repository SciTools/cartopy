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
import warnings

import matplotlib.pyplot as plt
from numpy.testing import assert_array_almost_equal

import cartopy.crs as ccrs

class TestGetExtent(unittest.TestCase):
    def test_get_extent(self):
        # Set up known extent around the UK.
        uk = (-12.5, 4., 49., 60.)
        uk_crs = ccrs.Geodetic()
        projection = ccrs.Mollweide()
        ax = plt.axes(projection=projection)
        ax.set_extent(uk, crs=uk_crs)

        # Result of get_extent with no specified coordinate system should be
        # in the projection of the axes (Mollweide in this case).
        assert_array_almost_equal(ax.get_extent(), ax.get_extent(projection))

        # Obtain the extent in a range of projections and check against
        # expected values.

        # Mollweide
        expected = (-963125.8421709834, 308200.26949471474,
                    5768352.350108459, 6876758.993328802)
        extent = ax.get_extent(ccrs.Mollweide())
        assert_array_almost_equal(extent, expected)

        # Geodetic (should raise warning, but still give the PlateCarree result).
        expected_plate_carree = (-14.850129, 4.752041, 49., 60.)
        expected = expected_plate_carree
        with warnings.catch_warnings():
            # Check for warning.
            warnings.simplefilter('error')
            with self.assertRaises(UserWarning):
                ax.get_extent(ccrs.Geodetic())
            # Check result.
            warnings.simplefilter('ignore')
            extent = ax.get_extent(ccrs.Geodetic())
            assert_array_almost_equal(extent, expected)

        # PlateCarree
        expected = expected_plate_carree
        extent = ax.get_extent(ccrs.PlateCarree())
        assert_array_almost_equal(extent, expected)

        # PlateCarree with nonzero central_longitude
        central_longitude = 20
        expected = (expected_plate_carree[0] - central_longitude,
                    expected_plate_carree[1] - central_longitude,
                    expected_plate_carree[2],
                    expected_plate_carree[3])
        extent = ax.get_extent(ccrs.PlateCarree(central_longitude))
        assert_array_almost_equal(extent, expected)

        # NorthPolarStereo
        expected = (-1034046.2256626057, 333263.47741164186,
                    -4765889.766015143, -3311994.6422885)
        extent = ax.get_extent(ccrs.NorthPolarStereo())
        assert_array_almost_equal(extent, expected)


if __name__ == '__main__':
    unittest.main()
