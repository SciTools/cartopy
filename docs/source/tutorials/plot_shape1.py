import matplotlib.pyplot as plt
from matplotlib import patches
from shapely.geometry import MultiPolygon, Polygon

import cartopy.feature as feature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shp


COLOURS = feature._COLOURS
COLOURS.update({'land_ig':(220/256., 220/256., 220/256.)})

# Set up axes.
ax = plt.axes(projection=ccrs.PlateCarree())
lon_min = -13.
lon_max = 5.
lat_min = 48.
lat_max = 60.
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
read_ctry = shp.Reader(ctry_fname)
countries = read_ctry.geometries()
countries = feature.ShapelyFeature(
    countries, crs=src_crs, edgecolor='black',
    facecolor=COLOURS['land_ig'])

read_pop = shp.Reader(pop_fname)
population = read_pop.records()

read_riv = shp.Reader(rivers_fname)
rivers = read_riv.geometries()
rivers = feature.ShapelyFeature(
    rivers, crs=src_crs, edgecolor='blue', facecolor='none',
    linewidth=1.0, zorder=10)

# Return city name records.
pop_uk = [pop for pop in population if
          lat_min < pop.attributes['latitude'] < lat_max and
          lon_min < pop.attributes['longitude'] < lon_max]

# Set background colour (water).
ax.background_patch.set_facecolor(COLOURS['water'])


# Add features - countries and rivers.
ax.add_feature(countries)
ax.add_feature(rivers)

# Add geometries - city names.
for city in pop_uk:
    plt.text(city.geometry.x, city.geometry.y, city.attributes['name'],
              horizontalalignment='right', transform=src_crs)

    city_patch = patches.Circle([city.geometry.x, city.geometry.y],
                                radius=0.1, edgecolor='black',
                                facecolor='red', zorder=2,
                                linewidth=0.8)
    ax.add_patch(city_patch)
    
plt.show()
