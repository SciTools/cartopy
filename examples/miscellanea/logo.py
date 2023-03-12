"""
Cartopy Logo
------------

The actual code to produce cartopy's logo.

"""
from matplotlib.font_manager import FontProperties
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.textpath
import matplotlib.transforms

import cartopy.crs as ccrs


def main():
    fig = plt.figure(figsize=[12, 6])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    ax.coastlines()
    ax.gridlines()

    # generate a matplotlib path representing the word "cartopy"
    fp = FontProperties(family='DejaVu Sans', weight='bold')
    logo_path = matplotlib.textpath.TextPath((-175, -35), 'cartopy',
                                             size=80, prop=fp)
    # scale the letters up to sensible longitude and latitude sizes
    transform = matplotlib.transforms.Affine2D().scale(1, 2).translate(0, 35)

    # add a background image
    im = ax.stock_img()
    # Apply the scale transform and then the map coordinate transform
    plate_carree_transform = (transform +
                              ccrs.PlateCarree()._as_mpl_transform(ax))

    # add the path as a patch, drawing black outlines around the text
    patch = matplotlib.patches.PathPatch(logo_path,
                                         facecolor='none', edgecolor='black',
                                         transform=plate_carree_transform)
    im.set_clip_path(patch)
    ax.add_patch(patch)

    plt.show()


if __name__ == '__main__':
    main()
