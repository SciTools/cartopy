"""
"""
import matplotlib.pyplot as plt
import numpy

import cartopy.crs as ccrs
import cartopy.io.srtm as io_srtm


def main():
    pc = ccrs.PlateCarree()

    ax = plt.axes(projection=pc)

    img, crs, extent = io_srtm.srtm(-4, 50)
    ax.imshow(numpy.log(img),
                extent=extent,
                transform=crs)


    plt.show()


if __name__ == '__main__':
    main()
