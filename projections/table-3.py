import matplotlib.pyplot as plt
import cartopy.crs as ccrs

plt.figure(figsize=(9.42477796077, 3))
delta = 0.125
ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.LambertCylindrical())
#ax.set_global()
ax.coastlines()
ax.gridlines()