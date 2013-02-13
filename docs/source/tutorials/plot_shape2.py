import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as feature


ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-30., 30., 20., 80.])

# set background colour (water)
ax.background_patch.set_facecolor(feature._COLOURS['water'])

# create a new feature
countries = feature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='50m',
    facecolor=feature._COLOURS['land'])

# add features
ax.add_feature(countries)
ax.add_feature(feature.RIVERS)
ax.add_feature(feature.LAKES)

# plot stock image
ax.stock_img()

plt.show()
