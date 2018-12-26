"""
Cartopy Logo
------------

The actual code to produce cartopy's logo.

"""
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import InterProjectionTransform
import matplotlib.pyplot as plt
import matplotlib.textpath
import matplotlib.patches
from matplotlib.font_manager import FontProperties
import numpy as np


def add_cartopy_logo(ax):
    # Generate a matplotlib path representing the word "cartopy".
    fp = FontProperties(family='Bitstream Vera Sans', weight='bold')
    logo_path = matplotlib.textpath.TextPath((-175, -35), 'cartopy',
                                             size=1, prop=fp)
    # Scale the letters up to sensible longitude and latitude sizes.
    logo_path._vertices *= np.array([80, 160])

    # Add a background image.
    im = ax.stock_img()

    # Clip the image according to the logo_path.
    pc_to_data = InterProjectionTransform(ccrs.PlateCarree(), ax.projection)
    im.set_clip_path(logo_path, transform=pc_to_data + ax.transData)

    # Add the path as a patch, drawing black outlines around the text.
    patch = matplotlib.patches.PathPatch(logo_path,
                                         facecolor='none', edgecolor='black',
                                         transform=ccrs.PlateCarree())
    ax.add_patch(patch)


def main():
    fig = plt.figure(figsize=[12, 6])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    ax.coastlines()
    ax.gridlines()
    add_cartopy_logo(ax)

    plt.show()


if __name__ == '__main__':
    main()
