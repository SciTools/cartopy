__tags__ = ['Lines and polygons']
import matplotlib.pyplot as plt

import cartopy.crs as ccrs


def main():
    ax = plt.axes(projection=ccrs.Robinson())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()

    ax.stock_img()
    ax.coastlines()

    plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
    plt.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.PlateCarree())
    plt.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic())

    plt.show()


if __name__ == '__main__':
    main()
