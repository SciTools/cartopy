"""
"""
import matplotlib.pyplot as plt
import numpy

import cartopy.crs as ccrs
import cartopy.io.srtm as io_srtm


def main():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm(-4, 50)

    elev = numpy.ma.masked_less_equal(elev, 0, copy=False)

    ax.imshow(numpy.ma.log(elev**2),
                extent=extent,
                transform=crs,
                cmap='Greens',
                )

    ax.coastlines()
    plt.show()


def main2():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm_composite(-6, 50, 4, 3)

    elev = numpy.ma.masked_less_equal(elev, 0, copy=False)

    ax.imshow(numpy.ma.log(elev**2),
              extent=extent,
              transform=crs,
              cmap='Greens',
              )

    ax.coastlines()

    plt.show()

if __name__ == '__main__':
    main2()
