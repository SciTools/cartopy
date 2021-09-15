# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt

import cartopy.crs as ccrs


class Gridliner:
    params = [
        (False, True),
        (False, True),
    ]
    param_names = ['draw_labels', 'inline']

    def setup(self, draw_labels, inline):
        self.proj = ccrs.PlateCarree()
        fig, ax = plt.subplots(subplot_kw=dict(projection=self.proj))
        ax.gridlines(draw_labels=draw_labels, x_inline=inline, y_inline=inline,
                     auto_inline=False)
        self.figure = fig

    def time_gridlines(self, draw_labels, inline):
        self.figure.canvas.draw()
