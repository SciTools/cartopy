# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.scalebar import fancy_scalebar


def scalebar_invoqued_as_an_independent_object_geoaxes():
    """Test that we do not trigger matplotlib's line subslice optimization."""
    # This behavior caused lines with > 1000 points and
    # sorted data to disappear

    fig, axes = plt.subplots(1, 2,
                             subplot_kw={'projection':
                                         ccrs.Mercator()})

    projections = [ccrs.Mercator(), ccrs.PlateCarree()]

    axes = axes.ravel()

    for proj, ax in zip(projections, axes):
        ax.projection = proj
        fancy_scalebar(ax,
                       location=(0.5, 0.2),
                       length=1000_000,
                       unit_name='km',
                       angle=0,
                       max_stripes=5,
                       fontsize=8,
                       dy=0.05)

        ax.gridlines(draw_labels=True)
        ax.stock_img()
        ax.coastlines()
    


def scalebar_from_within_geoaxes():
    """Test that we do not trigger matplotlib's line subslice optimization."""
    # This behavior caused lines with > 1000 points and
    # sorted data to disappear

    fig, axes = plt.subplots(1, 2,
                             subplot_kw={'projection':
                                         ccrs.Mercator()})

    projections = [ccrs.Mercator(),
                   ccrs.PlateCarree()]

    axes = axes.ravel()
    
    for proj, ax in zip(projections, axes):
        ax.projection = proj

        ax.set_extent([-60, -35, -40, 10])
        ax.gridlines(draw_labels=True)
        ax.add_scalebar(location=(0.5, 0.5),
                        length=250_000,
                        dy=5,
                        max_stripes=3)
        ax.stock_img()
        ax.coastlines()

    fig.show()

