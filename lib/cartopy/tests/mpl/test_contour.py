# (C) British Crown Copyright 2016, Met Office
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

import matplotlib.pyplot as plt
from matplotlib.testing.decorators import cleanup
import numpy as np
from numpy.testing import assert_array_almost_equal

import cartopy.crs as ccrs


@cleanup
def test_contour_plot_bounds():
    x = np.linspace(-2763217.0, 2681906.0, 200)
    y = np.linspace(-263790.62, 3230840.5, 130)
    data = np.hypot(*np.meshgrid(x, y)) / 2e5

    proj_lcc = ccrs.LambertConformal(central_longitude=-95,
                                     central_latitude=25,
                                     standard_parallels=[25])
    ax = plt.axes(projection=proj_lcc)
    ax.contourf(x, y, data, levels=np.arange(0, 40, 1))
    assert_array_almost_equal(ax.get_extent(),
                              np.array([x[0], x[-1], y[0], y[-1]]))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
