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
import numpy as np
from numpy.testing import assert_array_equal

import cartopy.img_transform as img_trans
import cartopy.crs as ccrs


def test_griding_data():
    target_prj = ccrs.PlateCarree()
    import numpy
    # create 3 data points
    lats = np.array([45,  20, -45], dtype=np.float64)
    lons = numpy.array([-90, 90,   0], dtype=numpy.float64)
    data = numpy.array([1,    2,   3], dtype=numpy.float64)
    data_trans = ccrs.Geodetic()

    target_x, target_y, extent = img_trans.mesh_projection(target_prj, 8, 4)

    image = img_trans.regrid(data, lons, lats, data_trans, target_prj, target_x, target_y)

    # The expected image. n.b. on a map the data is reversed in the y axis.
    expected = np.array([[3., 3., 3., 3., 3., 3., 3., 3.],
                         [1., 1., 3., 3., 3., 2., 2., 2.],
                         [1., 1., 1., 1., 2., 2., 2., 2.],
                         [1., 1., 1., 1., 1., 2., 2., 1.]], dtype=np.float64)

    assert_array_equal([-180,  180,  -90,   90], extent)
    assert_array_equal(expected, image)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
