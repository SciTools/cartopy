"""
This example plots the Natual Earth data along with a copyright note.
"""
__tags__ = ['Lines and polygons']
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText

# Define the license information:
# NOTE: Adapt to your data source
data_source = 'Natural Earth'
# NOTE: adapt to your data source
data_license = 'public domain'


def main():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([80, 170, -45, 30])

    # Put a background image on for nice sea rendering.
    ax.stock_img()

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(states_provinces, edgecolor='gray')

    # Include the license information:
    # Add a text annotation for the license information to the
    # the bottom right corner .
    # ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    text = AnchoredText(r'$\mathcircled{c}$ ' + data_source + '; license: ' +
                        data_license, loc=4, prop={'size': 10},
                        frameon=True)
    ax.add_artist(text)

    plt.show()


if __name__ == '__main__':
    main()
