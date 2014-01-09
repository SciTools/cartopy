import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs


def gen_data():
    # Generate some data
    x = np.arange(-60, 45, 5)
    y = np.arange(30, 80, 5)

    x2d, y2d = np.meshgrid(x, y)
    u = 40 * np.cos(np.deg2rad(y2d))
    v = 40 * np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    return dict(x=x, y=y, u=u, v=v, magnitude=mag)


# Get some data with x, y, u, v, components
wind = gen_data()

# Setup our figure
plot_extent = [-60, 40, 30, 70]
plt.figure(figsize=(12, 12))

# plot on the native projection
ax = plt.subplot(321, projection=ccrs.PlateCarree())
ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
ax.coastlines()
ax.barbs(wind['x'], wind['y'], wind['u'], wind['v'], length=4,
         linewidth=.25)
ax.stock_img()
plt.title('Wind barbs, native PlateCarree projection')

# plot on a different projection
ax = plt.subplot(322, projection=ccrs.NorthPolarStereo())
ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
ax.coastlines()
ax.barbs(wind['x'], wind['y'], wind['u'], wind['v'], length=4,
         linewidth=.25, transform=ccrs.PlateCarree())
ax.stock_img()
plt.title('Wind barbs, North Polar Stereographic projection')

# Wind arrow plots
ax = plt.subplot(323, projection=ccrs.PlateCarree())
ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
ax.coastlines()
ax.stock_img()
plt.title('Wind arrows, native PlateCarree projection')
ax.quiver(wind['x'], wind['y'], wind['u'], wind['v'], wind['magnitude'])

ax = plt.subplot(324, projection=ccrs.NorthPolarStereo())
ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
ax.coastlines()
ax.stock_img()
plt.title('Wind arrows, North Polar Stereographic projection')
ax.quiver(wind['x'], wind['y'], wind['u'], wind['v'], wind['magnitude'],
          transform=ccrs.PlateCarree())

# Stream plots
ax = plt.subplot(325, projection=ccrs.PlateCarree())
ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
ax.coastlines()
ax.stock_img()
plt.title('Stream plot, native PlateCarree projection')
ax.streamplot(wind['x'], wind['y'], wind['u'], wind['v'], density=1,
              arrowsize=1,
              linewidth=3 * (wind['magnitude']/wind['magnitude'].max()),
              color=wind['magnitude'])

ax = plt.subplot(326, projection=ccrs.NorthPolarStereo())
ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
ax.coastlines()
ax.stock_img()
plt.title('Stream plot, North Polar Stereographic projection')
ax.streamplot(wind['x'], wind['y'], wind['u'], wind['v'], density=1,
              arrowsize=1, transform=ccrs.PlateCarree(),
              linewidth=3 * (wind['magnitude']/wind['magnitude'].max()),
              color=wind['magnitude'])

plt.show()
