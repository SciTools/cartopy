"""
Feature Creation
----------------

This example manually instantiates a
:class:`cartopy.feature.NaturalEarthFeature` to access administrative
boundaries (states and provinces).

Note that this example is intended to illustrate the ability to construct
Natural Earth features that cartopy does not necessarily know about
*a priori*.
In this instance however, it would be possible to make use of the
pre-defined :data:`cartopy.feature.STATES` constant.

"""
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature


def main():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([80, 170, -45, 30], crs=ccrs.PlateCarree())

    # Put a background image on for nice sea rendering.
    ax.stock_img()

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth.
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    SOURCE = 'Natural Earth'
    LICENSE = 'public domain'

    # Add our states feature.
    ax.add_feature(states_provinces, edgecolor='gray')
    # Add land feature, overriding the default negative zorder so it shows
    # above the background image.
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='k')

    # Add a text annotation for the license information to the
    # the bottom right corner.
    text = AnchoredText(f'\u00A9 {SOURCE}; license: {LICENSE}',
                        loc=4, prop={'size': 12}, frameon=True)
    ax.add_artist(text)

    plt.show()


if __name__ == '__main__':
    main()
