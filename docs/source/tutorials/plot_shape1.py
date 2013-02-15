import matplotlib.pyplot as plt
from matplotlib import patches

import cartopy.feature as feature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shp


# Set up axes.
ax = plt.axes(projection=ccrs.PlateCarree())
lon_min, lon_max = (-13., 5.)
lat_min, lat_max = (48., 60.)
ax.set_extent([lon_min, lon_max, lat_min, lat_max])

# Get paths to shapefiles.
ctry_fname = shp.natural_earth(resolution='50m',
                               category='cultural',
                               name='admin_0_countries')

pop_fname = shp.natural_earth(resolution='50m',
                              category='cultural',
                              name='populated_places_simple')

rivers_fname = shp.natural_earth(resolution='10m',
                                 category='physical',
                                 name='rivers_lake_centerlines')

# Define source reference coordinate system.
src_crs = ccrs.PlateCarree()

# Read files.
countries = shp.Reader(ctry_fname).geometries()
population = shp.Reader(pop_fname).records()
rivers = shp.Reader(rivers_fname).geometries()

# Return city name records.
clip_population = [pop for pop in population if
                   lat_min < pop.attributes['latitude'] < lat_max and
                   lon_min < pop.attributes['longitude'] < lon_max]

# Set background color (water).
ax.background_patch.set_facecolor(feature.COLORS['water'])

# Add features - countries and rivers.
ax.add_geometries(countries, crs=src_crs, edgecolor='black',
                  facecolor=feature.COLORS['land'])
ax.add_geometries(rivers, crs=src_crs, edgecolor='blue', facecolor='none',
                  linewidth=1.0, zorder=10)

# Add geometries - city names.
for city in clip_population:
    plt.text(city.geometry.x, city.geometry.y, city.attributes['name'],
             horizontalalignment='right', transform=src_crs)

    city_patch = patches.Circle([city.geometry.x, city.geometry.y],
                                radius=0.1, edgecolor='black',
                                facecolor='red', zorder=2,
                                linewidth=0.8)
    ax.add_patch(city_patch)

plt.show()
