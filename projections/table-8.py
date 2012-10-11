import matplotlib.pyplot as plt
import cartopy.crs as ccrs

plt.figure(figsize=(6, 3))
delta = 0.125
ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.PlateCarree(central_longitude=180))
#ax.set_global()
ax.coastlines()
ax.gridlines()