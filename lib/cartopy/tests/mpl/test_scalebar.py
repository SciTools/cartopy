# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.scalebar import fancy_scalebar
import pytest
from cartopy.tests.mpl import MPL_VERSION, ImageTesting


if MPL_VERSION < "3":
    TOL = 15
else:
    TOL = 0.5
    
grid_label_tol = grid_label_inline_tol = grid_label_inline_usa_tol = TOL
grid_label_inline_tol += 1.1
grid_label_image = 'scalebar_plot'
grid_label_inline_image = 'scalebar_inline_plot'


@pytest.mark.natural_earth
@ImageTesting([grid_label_image], tolerance=grid_label_tol)
def test_scalebar():
    """Test the scalebar"""

    fig, axes = plt.subplots(1, 2,
                             subplot_kw={'projection':
                                         ccrs.Mercator()})

    projections = [ccrs.Mercator(), ccrs.PlateCarree()]

    axes = axes.ravel()
    try:
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
        plt.close('all')
        condition = True
    except BaseException as err:
        print(err)
        condition = False

    assert(condition)




@pytest.mark.natural_earth
@ImageTesting([grid_label_inline_image], tolerance=grid_label_tol)
def test_scalebar_within_geoaxes():
    """Test scalebat within the geoaxes"""

    fig, axes = plt.subplots(1, 2,
                             subplot_kw={'projection':
                                         ccrs.Mercator()})

    projections = [ccrs.Mercator(),
                   ccrs.PlateCarree()]

    axes = axes.ravel()
    try:
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

        plt.close('all')

        condition = True
    except BaseException as err:
        print(err)
        condition = False

    assert(condition)
