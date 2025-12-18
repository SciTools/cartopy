"""
Extending longitude beyond 360 degrees
======================================

This example demonstrates the use of the over parameter with
Cartopy.  When using a cylindrical projection, setting over to True
engages proj's "+over" switch and enables longitudes to be extended
so that a single map can show more than 360 degrees of the globe.  This
can be useful for ensuring that large structures can be shown in their
entirety, unbroken by the edge of the map.

The underlying data needs to be explicitly extended, for which the
utility routines extend_lons_np and extend_lons_xr are included in
cartopy.util.  This allows for the possibility of non-repeated data,
such as an object's trajectory; a trivial example of this is shown
by plotting Chicago twice with different coloured points.

The script also applies Nightshade, which requires the "over"
parameter to be set, and tissot and coastlines, which work transparently.

Note that due a limitation in the underlying proj library, the longitudes
are limited to [-572.95, 572.95].
"""

import datetime

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
from cartopy.util import extend_lons_np


date = datetime.datetime.now()
#date = datetime.datetime(2025, 3, 20, 9, 1, 0) # NH Spring equinox 2025
#date = datetime.datetime(2025, 6, 21, 2, 42, 0) # NH Summer solstice 2025

# Coordinates of a few places
toulouse = 43.604500, 1.444000
nyalesund = 78.9, 11.9
chicago = 39.162, -84.45689

### Load some Copernicus ocean altimetry data
# ds = open_dataset("nrt_global_allsat_phy_l4_20241125_20241125.nc")
# u = ds.ugosa[0,2:-1:4,1::3]
# v = ds.vgosa[0,2:-1:4,1::3]
# u, v = ds.ugosa[0, ...] , ds.vgosa[0,...]
# speed = np.sqrt(u**2 + v**2)
# x, y = speed.longitude, speed.latitude

### Or just make up some random, pretend meteorological data
incr = 0.25
nlats, nlons = 720, 1440
lat = np.arange(incr/2, 90, incr)
lat = np.concat([-lat[-1::-1], lat])
lon = np.arange(incr/2, 180, incr)
lon = np.concat([-lon[-1::-1], lon])
size = (nlats, nlons)
rfield = np.random.normal(size=size)
rfield = np.concat((rfield, rfield), axis=1)
feature = gaussian_filter(rfield,
            sigma=[nlats/12,nlons/4])[:,int(nlons/2):int(nlons*3/2)]

# var = speed
var = feature

# Extrema:
factor = 1
vmin = np.nanmin(var)/factor
vmax = np.nanmax(var)/factor
# # Make symmetrical
# vmin = -vmax

extend_cbar = "both"

def a_transform(arr):
    """Some random transformation to differentiate data."""
    amax = arr.max()
    return amax/2 - arr

# Design a few test configurations
mapconf1 = dict(title="Standard Mercator", lonmin=-180, lonmax=180, over=False,
                trans=ccrs.PlateCarree, proj=ccrs.Mercator,
                central_longitude=0)
mapconf2 = dict(title="Extended Mercator", lonmin=-390, lonmax=525, over=True,
                trans=ccrs.PlateCarree, proj=ccrs.Mercator,
                central_longitude=0)
mapconf3 = dict(title="Extended Plate Carrée", lonmin=-390, lonmax=525, over=True,
                trans=ccrs.PlateCarree, proj=ccrs.PlateCarree,
                central_longitude=0)

mapconfs = [mapconf1, mapconf2, mapconf3]

for ind, mapconf in enumerate(mapconfs[1:2]):
    print(ind, mapconf)
    lonmin, lonmax = mapconf["lonmin"], mapconf["lonmax"]
    over = mapconf["over"]
    central_longitude = mapconf["central_longitude"]
    proj = mapconf["proj"]
    trans = mapconf["trans"]
    projection = proj(over=over, central_longitude=central_longitude)
    transform = trans(over=over, central_longitude=central_longitude)
    title = mapconf["title"] + " [" + str(lonmin) + "," + str(lonmax) + "]"
    lon_ext, lat_ext, var_ext = extend_lons_np(lon, lat, var, lonmin, lonmax)
    ## Uncomment the following line to highlight the extra data (for Xarrays)
    # var_ext = var_ext.where(
    #     np.abs(var_ext.longitude)<180, a_transform(var_ext))

    fig = plt.figure(figsize=(10,4), facecolor=(0,0,0,0))
    ax = plt.axes(projection=projection)
    ax.pcolormesh(lon_ext, lat_ext, var_ext,
                  vmin=vmin, vmax = vmax, transform=transform)
    ax.coastlines()
    ax.gridlines(draw_labels=True, dms=True,
                 x_inline=False, y_inline=False,
                 crs=ccrs.PlateCarree(over=over))
    ax.add_feature(Nightshade(date, alpha=0.2, delta=0.1, over=over))
    ax.plot(toulouse[1], toulouse[0], "ro", transform=ccrs.Geodetic())
    ax.plot(nyalesund[1], nyalesund[0], "ro", transform=ccrs.Geodetic())
    ax.plot(chicago[1], chicago[0], "ro", transform=ccrs.Geodetic())
    ax.plot(chicago[1]+360, chicago[0], "orange", marker="o", transform=ccrs.Geodetic())
    ax.stock_img()
    ax.tissot(
        lons=np.arange(-590, 580,  200),
        lats=[-75, -60, -45, -30, -10, 20, 50, 65, 80],
        n_samples=40
    )
    ax.set_title(title)

plt.show()
#plt.savefig("eke_cartopy_extended_mercator.png")
