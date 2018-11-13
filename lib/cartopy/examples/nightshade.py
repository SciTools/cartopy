"""
Nightshade feature
------------------

Shade areas of night terminator

"""
__tags__ = ['Lines and polygons']

import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade


ax = plt.axes(projection=ccrs.PlateCarree())

date = datetime.datetime(1999, 12, 31, 12, tzinfo=datetime.timezone.utc)
plt.title('Night time shading for {}'.format(date))
ax.stock_img()
ax.add_feature(Nightshade(date))
plt.show()
