"""
Cartopy Favicon
---------------

The actual code to generate cartopy's favicon.

"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.textpath
import matplotlib.patches
from matplotlib.font_manager import FontProperties


def main():
    fig = plt.figure(figsize=[8, 8])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())

    ax.coastlines()
    ax.gridlines()
    ax.stock_img()

    # Generate a matplotlib path representing the character "C".
    fp = FontProperties(family='DejaVu Sans', weight='bold')
    logo_path = matplotlib.textpath.TextPath((-4.5e7, -3.7e7),
                                             'C', size=103250000, prop=fp)

    # Add the path as a patch, drawing black outlines around the text.
    patch = matplotlib.patches.PathPatch(logo_path, facecolor='white',
                                         edgecolor='black', linewidth=10,
                                         transform=ccrs.SouthPolarStereo())
    ax.add_patch(patch)
    plt.show()


if __name__ == '__main__':
    main()
