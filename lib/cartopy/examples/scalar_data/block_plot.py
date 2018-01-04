"""
Block Plot
----------

Shows how to display a colour-mapped blockplot of data over a world map.

.. _examples-blockplot:
"""
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt


def main():
    # load some sample iris data
    fname = iris.sample_data_path('rotated_pole.nc')
    temperature = iris.load_cube(fname)

    # iris comes complete with a method to put bounds on a simple point
    # coordinate. This is very useful...
    temperature.coord('grid_latitude').guess_bounds()
    temperature.coord('grid_longitude').guess_bounds()

    # turn the iris Cube data structure into numpy arrays
    gridlons = temperature.coord('grid_longitude').contiguous_bounds()
    gridlats = temperature.coord('grid_latitude').contiguous_bounds()
    temperature = temperature.data

    # set up a map
    ax = plt.axes(projection=ccrs.PlateCarree())

    # define the coordinate system that the grid lons and grid lats are on
    rotated_pole = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
    plt.pcolormesh(gridlons, gridlats, temperature, transform=rotated_pole)

    # add coastlines
    ax.coastlines()

    plt.show()


if __name__ == '__main__':
    main()
