

"""
Cartopy Logo
------------

The actual code to produce cartopy's logo.

"""
import cartopy.crs as ccrs
from cartopy.examples.logo import add_cartopy_logo
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.path as mpath
from matplotlib.collections import PathCollection


def rect_healpix_tabs(ax):
    tab_width = 2e6
    tab_step = 2e6

    horiz_tab = mpath.Path([[0, 0], [2e6, 2e6],
                            [1e7 - 2e6, 2e6], [1e7, 0]])

    paths = [
         # Side.
         mpath.Path([[2e7, 5e6], [2e7 + tab_width, 3e6], [2.2e7, -3e6], [2e7, -5e6]]),

         mpath.Path(horiz_tab.vertices + [-2e7, 5e6]),
         mpath.Path(horiz_tab.vertices + [0, 5e6]),
         mpath.Path(horiz_tab.vertices + [1e7, 5e6]),

         mpath.Path([1, -1] * horiz_tab.vertices + [-2e7, -5e6]),
         mpath.Path([1, -1] * horiz_tab.vertices + [-1e7, -5e6]),
         mpath.Path([1, -1] * horiz_tab.vertices + [1e7, -5e6])]

    healpix_tabs = PathCollection(
        paths, facecolor='none',
        edgecolor='black', transform=ax.transData, clip_on=False, zorder=-1)
    ax.add_collection(healpix_tabs)

def main():
    fig = plt.figure(figsize=[12, 6])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.RectangularHealpix(central_longitude=0, north_square=1, south_square=2))

    ax.coastlines()
    ax.gridlines()
    add_cartopy_logo(ax)

    rect_healpix_tabs(ax)

    plt.show()


if __name__ == '__main__':
    main()
