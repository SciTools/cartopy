"""
Plotting the Aurora Forecast from NOAA on Orthographic Polar Projection
-----------------------------------------------------------------------

The National Oceanic and Atmospheric Administration (NOAA) monitors the
solar wind conditions using the ACE spacecraft orbiting close to the L1
Lagrangian point of the Sun-Earth system. This data is fed into the
OVATION-Prime model to forecast the probability of visible aurora at
various locations on Earth. Every five minutes a new forecast is
published for the coming 30 minutes. The data is provided as a
1024 by 512 grid of probabilities in percent of visible aurora. The
data spaced equally in degrees from -180 to 180 and -90 to 90.

"""
__tags__ = ["Scalar data"]
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

from io import StringIO

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def aurora_forecast():
    """
    Gets the latest Aurora Forecast from http://swpc.noaa.gov.

    Returns
    -------
    img : numpy array
        The pixels of the image in a numpy array.
    img_proj : cartopy CRS
        The rectangular coordinate system of the image.
    img_extent : tuple of floats
        The extent of the image ``(x0, y0, x1, y1)`` referenced in
        the ``img_proj`` coordinate system.
    origin : str
        The origin of the image to be passed through to matplotlib's imshow.

    """

    # GitHub gist to download the example data from
    url = ('https://gist.githubusercontent.com/belteshassar/'
           'c7ea9e02a3e3934a9ddc/raw/aurora-nowcast-map.txt')
    # To plot the current forecast instead, uncomment the following line
    # url = 'http://services.swpc.noaa.gov/text/aurora-nowcast-map.txt'

    response_text = StringIO(urlopen(url).read().decode('utf-8'))
    img = np.loadtxt(response_text)
    img_proj = ccrs.PlateCarree()
    img_extent = (-180, 180, -90, 90)
    return img, img_proj, img_extent, 'lower'


def main():
    fig = plt.figure(figsize=[10, 5])
    # Orthographic projection looks natural and distortion is
    # small around the poles where the aurora is most likely

    # ax1 for Northern Hemisphere
    ax1 = plt.subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))
    # ax2 for Southern Hemisphere
    ax2 = plt.subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))

    # viridis was added in matplotlib 1.5. Try to use it.
    try:
        plt.set_cmap('viridis')
    except ValueError:
        plt.set_cmap('YlGnBu_r')

    img, crs, extent, origin = aurora_forecast()
    for ax in [ax1, ax2]:
        ax.coastlines()
        ax.gridlines()
        ax.imshow(img, vmin=0, vmax=100, transform=crs,
                  extent=extent, origin=origin)

    plt.show()


if __name__ == '__main__':
    main()
