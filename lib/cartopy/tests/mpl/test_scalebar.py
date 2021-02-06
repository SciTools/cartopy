# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.


import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs

from cartopy.tests.mpl import MPL_VERSION, ImageTesting

from cartopy.mpl.scalebar import add_scalebar


# The tolerance on these tests are particularly high because of the high number
# of text objects. A new testing strategy is needed for this kind of test.
if MPL_VERSION < "3":
    TOL = 15
else:
    TOL = 0.5

grid_label_tol = grid_label_inline_tol = grid_label_inline_usa_tol = TOL

scalebar_plot = 'scalebar_plot'
scalebar_inline_plot = 'scalebar_inline_plot'


@pytest.mark.natural_earth
@ImageTesting([scalebar_plot], tolerance=grid_label_tol)
def test_scalebar():
    """Test"""

    fig, axes = plt.subplots(1, 2,
                             subplot_kw={'projection':
                                         ccrs.Mercator()})

    projections = [ccrs.Mercator(), ccrs.PlateCarree()]

    axes = axes.ravel()

    for proj, ax in zip(projections, axes):

        ax.projection = proj

        ax.set_title(ax.projection.__class__.__name__)

        add_scalebar(ax,
                     bbox_to_anchor=(0.1, 0.2),
                     length=10_000_000,
                     ruler_unit='km',
                     max_stripes=3,
                     fontsize=8,
                     frameon=True,
                     ruler_unit_fontsize=15,
                     ruler_fontweight='bold',
                     tick_fontweight='bold',
                     dy=0.085)

        ax.gridlines(draw_labels=True)
        ax.stock_img()
        ax.coastlines()

    return axes


@pytest.mark.natural_earth
@ImageTesting([scalebar_inline_plot], tolerance=grid_label_tol)
def test_scalebar2():
    """Test"""

    fig, axes = plt.subplots(1, 2,
                             subplot_kw={'projection':
                                         ccrs.Mercator()})

    projections = [ccrs.Mercator(), ccrs.PlateCarree()]

    axes = axes.ravel()

    for proj, ax in zip(projections, axes):

        ax.projection = proj

        ax.set_title(ax.projection.__class__.__name__)

        ax.add_scalebar(ax,
                        bbox_to_anchor=(0.1, 0.2),
                        length=10_000_000,
                        ruler_unit='km',
                        max_stripes=3,
                        fontsize=8,
                        frameon=True,
                        ruler_unit_fontsize=15,
                        ruler_fontweight='bold',
                        tick_fontweight='bold',
                        dy=0.085)

        ax.gridlines(draw_labels=True)
        ax.stock_img()
        ax.coastlines()

    return axes


test_scalebar2()
