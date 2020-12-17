# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt

import cartopy.crs as ccrs

from cartopy.mpl.scalebar import fancy_scalebar


def test_scalebar():
    """Test that we do not trigger matplotlib's line subslice optimization."""
    # This behavior caused lines with > 1000 points and
    # sorted data to disappear

    fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.Mercator()})

    projections = [ccrs.Mercator(), ccrs.PlateCarree()]

    axes = axes.ravel()
    try:
        for proj, ax in zip(projections, axes):
            ax.projection = proj
            fancy_scalebar(ax,
                           location=(0.5, 0.2),
                           length=10000,
                           metres_per_unit=1000,
                           unit_name='km',
                           tol=0.01,
                           angle=0,
                           max_stripes=5,
                           fontsize=8,
                           dy=0.05)
        condition = True
    except BaseException:
        condition = False

    assert(condition)


def test_scalebar_within_geoaxes():
    """Test that we do not trigger matplotlib's line subslice optimization."""
    # This behavior caused lines with > 1000 points and
    # sorted data to disappear

    fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.Mercator()})

    projections = [ccrs.Mercator(), ccrs.PlateCarree()]

    axes = axes.ravel()
    try:
        for proj, ax in zip(projections, axes):
            ax.projection = proj

            ax.set_extent([-60, -35, -40, 10])

            ax.add_scalebar(location=(0.5, 0.2),
                            length=250_000,
                            dy=5,
                            max_stripes=3)

            fig.show()

        condition = True
    except BaseException:
        condition = False

    assert(condition)
