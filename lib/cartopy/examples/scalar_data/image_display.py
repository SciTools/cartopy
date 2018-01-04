"""
Image Display
-------------

Shows how to display an image over a world map.

By specifying the geolocation of the corner points, it can be positioned
precisely.

.. _examples-image-display:
"""
import os

import matplotlib.pyplot as plt

from cartopy import config
import cartopy.crs as ccrs


def main():
    plt.figure(figsize=(8, 12))

    # get the path of the file. It can be found in the repo data directory.
    fname = os.path.join(config["repo_data_dir"],
                         'raster', 'sample', 'Miriam.A2012270.2050.2km.jpg')
    img = plt.imread(fname)

    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.title('Hurricane Miriam from the Aqua/MODIS satellite\n'
              '2012 09/26/2012 20:50 UTC')

    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    # encode the known extent of the image : the minimum and maximum X and Y
    # coordinate values it covers, in the coordinate system it is sampled in.
    img_extent = (-120.67660000000001, -106.32104523100001,
                  13.2301484511245, 30.766899999999502)

    # add the image. Because this image was a tif, the "origin" of the image is
    # in the upper left corner
    ax.imshow(img, origin='upper', extent=img_extent,
              transform=ccrs.PlateCarree())
    ax.coastlines(resolution='50m', color='black', linewidth=1)

    # mark a known place to help us geo-locate ourselves
    ax.plot(-117.1625, 32.715, 'bo', markersize=7, transform=ccrs.Geodetic())
    ax.text(-117, 33, 'San Diego', transform=ccrs.Geodetic())

    plt.show()


if __name__ == '__main__':
    main()
