import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# The code below fails without the bugfix proposed
fig = plt.figure(figsize=(18, 12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()

water_color = cfeature.COLORS['water']
rivers = cfeature.NaturalEarthFeature(scale='10m', category='physical',
        name='rivers_lake_centerlines', edgecolor=water_color, facecolor='none')

ax.add_feature(rivers)
plt.savefig('rivers.png')

