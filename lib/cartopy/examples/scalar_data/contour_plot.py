"""
Contour Plot
------------

Shows how to display filled contours of data over a world map.

.. _examples-contour-plot:
"""
import os

import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset

from cartopy import config
import cartopy.crs as ccrs


def main():
    # get the path of the file. It can be found in the repo data directory.
    fname = os.path.join(config["repo_data_dir"],
                         'netcdf', 'HadISST1_SST_update.nc'
                         )

    # extract file variables containing the data and location coordinates
    dataset = netcdf_dataset(fname)
    sst = dataset.variables['sst'][0, :, :]
    lats = dataset.variables['lat'][:]
    lons = dataset.variables['lon'][:]

    # make a colour-filled contour plot over a map
    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.contourf(lons, lats, sst, 60,
                 transform=ccrs.PlateCarree())

    # overlay coastlines for context
    ax.coastlines()
    plt.show()


if __name__ == '__main__':
    main()
