import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as feature

colours = {'land': (240/256., 240/256., 220/256.),
           'water': (152/256., 183/256., 226/256.)}

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-30., 30., 20., 80.])

# set background colour (water)
ax.background_patch.set_facecolor(colours['water'])

# create a new feature
countries = feature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor=colours['land'])

# add features
ax.add_feature(countries)
ax.add_feature(feature.RIVERS)
ax.add_feature(feature.LAKES)

# plot stock image
ax.stock_img()

plt.show()
