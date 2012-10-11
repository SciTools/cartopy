import os

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.image import imread

ax = plt.axes(projection=ccrs.Robinson())

ax.set_global()
# get the path to an image (in this case, a stock image which ships with cartopy)
fname = os.path.join(os.path.dirname(ccrs.__file__), 'data',
                     'raster', 'natural_earth', '50-natural-earth-1-downsampled.png')
img = imread(fname)
ax.imshow(img, origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -90, 90])

ax.coastlines()

plt.show()