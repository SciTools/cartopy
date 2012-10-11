import matplotlib.pyplot as plt
import cartopy.crs as ccrs

plt.figure(figsize=(4.09152901955, 3))
delta = 0.125
ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Miller())
#ax.set_global()
ax.coastlines()
ax.gridlines()