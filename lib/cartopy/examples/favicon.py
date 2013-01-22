import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.patch as pt
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy

import matplotlib.textpath
import matplotlib.patches
from matplotlib.font_manager import FontProperties


def main():

    plt.figure(figsize=[8, 8])
    ax = plt.axes(projection=ccrs.SouthPolarStereo())

    ax.coastlines()
    ax.gridlines()

    # Add a background image. Note: unlike the Robinson projection, the
    # background image spills outside the boundary due to the fact that
    # points outside the boundary have 1:1 mappings in this projection
    im = ax.stock_img()

    # Clip the image to the current background boundary. This will not
    # be sufficient if zooming or saving
    im.set_clip_path(ax.background_patch.get_path(),
                     transform=ax.background_patch.get_transform())

    # Generate a matplotlib path representing the character "C"
    fp = FontProperties(family='Arial', weight='bold')
    xy = (-4.5e7,-3.7e7)
    logo_path = matplotlib.textpath.TextPath(xy, 'C', size=1, prop=fp)

    # Scale the letters up to sensible X and Y sizes
    XY = [123500000, 103250000]
    logo_path._vertices *= numpy.array(XY)

    # Add the path as a patch, drawing black outlines around the text
    patch = matplotlib.patches.PathPatch(logo_path, facecolor='white',
                                         edgecolor='black', linewidth=10,
                                         transform=ccrs.SouthPolarStereo())
    ax.add_patch(patch)

    plt.show()


if __name__ == '__main__':
    main()
