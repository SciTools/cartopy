import matplotlib.pyplot as plt
import cartopy.crs as ccrs

plt.figure(figsize=(6, 3))
delta = 0.125
ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5))
#ax.set_global()
ax.coastlines()
ax.gridlines()