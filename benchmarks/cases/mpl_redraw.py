# (C) British Crown Copyright 2019, Met Office
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

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import io


# No need for anything other than the agg backend, and we don't want
# windows popping up as we are running these tests.
plt.switch_backend('agg')


def create_pc_png():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.stock_img()
    plt.savefig(io.BytesIO(), format='png')
    plt.close(fig)


def time_basic_draw_speed():
    create_pc_png()


def time_second_figure():
    # Successive figures with Axes of the same projection
    # could have various caching mechanisms in place.
    # At the time of writing, there is no
    # noticable performance speedup during the second figure :(
    create_pc_png()
    create_pc_png()
