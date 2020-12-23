"""
Nightshade feature
------------------

Draws a polygon where there is no sunlight for the given datetime.

"""
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade


fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

date = datetime.datetime(1999, 12, 31, 12)

ax.set_title(f'Night time shading for {date}')
ax.stock_img()
ax.add_feature(Nightshade(date, alpha=0.2))
plt.show()
