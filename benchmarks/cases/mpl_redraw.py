# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

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
