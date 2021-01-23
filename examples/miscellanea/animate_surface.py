"""
Animating a gridded surface
---------------------------

This example demonstrates how to animate
gridded data using `pcolormesh()`.
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import cartopy.crs as ccrs

fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()

x = np.linspace(-80, 80)
xs, ys = np.meshgrid(2 * x + 180, x)
zs = xs + ys
vmin, vmax = np.min(zs), np.max(zs)
mesh = ax.pcolormesh(xs, ys, np.zeros_like(zs), transform=ccrs.PlateCarree(),
                     shading='auto', vmin=vmin, vmax=vmax)

n = 10


def update_mesh(t):
    mesh.set_array(zs.ravel() * t)


ts = [i / n for i in range(n)]
# Go back to the start to make it a smooth repeat
ts += ts[::-1]
ani = FuncAnimation(fig, update_mesh, frames=ts,
                    interval=100)

plt.show()
