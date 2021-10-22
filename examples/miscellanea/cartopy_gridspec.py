"""
Using Cartopy and AxesGrid toolkit
----------------------------------
This example demonstrates how to use cartopy `GeoAxes` with
matplotlib's gridspec. Creates a list of axes for a figure
with gridspec layout of pre-defined number of rows and columns

Parameters:
rows (int): Number of rows in grid
cols (int): Number of columns in grid
proj (cartopy.crs projection): cartopy.crs projection to be used

Returns:
axlist: List of axes for each individual figure in the grid.

>>> make_grid(2,2,ccrs.PlateCarree())
[<AxesSubplot:>, <AxesSubplot:>, <AxesSubplot:>, <AxesSubplot:>, <AxesSubplot:>]

>>> make_grid(2,2,ccrs.PlateCarree())
[<AxesSubplot:>, <AxesSubplot:>, <AxesSubplot:>, <AxesSubplot:>]
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.gridspec as gs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
import numpy as np


def sample_data_3d(shape):
    """Return `lons`, `lats`, `times` and fake `data`"""
    ntimes, nlats, nlons = shape
    lats = np.linspace(-np.pi / 2, np.pi / 2, nlats)
    lons = np.linspace(0, 2 * np.pi, nlons)
    lons, lats = np.meshgrid(lons, lats)
    wave = 0.75 * (np.sin(2 * lats) ** 8) * np.cos(4 * lons)
    mean = 0.5 * np.cos(2 * lats) * ((np.sin(2 * lats)) ** 2 + 2)

    lats = np.rad2deg(lats)
    lons = np.rad2deg(lons)
    data = wave + mean

    times = np.linspace(-1, 1, ntimes)
    new_shape = data.shape + (ntimes,)
    data = np.rollaxis(data.repeat(ntimes).reshape(new_shape), -1)
    data *= times[:, np.newaxis, np.newaxis]

    return lons, lats, times, data


def make_grid(rows: int = 0, cols: int = 0, proj=ccrs.PlateCarree()) -> list:

    axlist = []
    fig2 = plt.figure(figsize=(20, 10))

    spec2 = gs.GridSpec(ncols=cols, nrows=rows, figure=fig2)
    spec2.update(wspace=0.2, hspace=0.5)

    for x in range(rows):
        for y in range(cols):
            axlist.append(fig2.add_subplot(spec2[x, y], projection=proj))

    return axlist


def main():
    projection = ccrs.PlateCarree()

    lons, lats, times, data = sample_data_3d((6, 73, 145))

    axlist = make_grid(rows=3, cols=2, proj=projection)

    for i, ax in enumerate(axlist):

        ax.set_xticks(np.linspace(-180, 180, 5), crs=projection)
        ax.set_yticks(np.linspace(-90, 90, 5), crs=projection)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

        p = ax.contourf(lons, lats, data[i, ...], transform=projection, cmap="RdBu")
        ax.coastlines()

    plt.colorbar(p, ax=axlist[:], orientation="vertical", shrink=1.0)
    plt.show()


if __name__ == "__main__":
    main()
