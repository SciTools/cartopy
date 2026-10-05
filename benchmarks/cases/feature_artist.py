# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import io

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.feature_artist import FeatureArtist


# No need for anything other than the agg backend, and we don't want
# windows popping up as we are running these tests.
plt.switch_backend('agg')


class DrawNaturalEarthFeature:
    params = [
        ('regional', 'global'),
        ('110m', '50m', '10m'),
    ]
    param_names = ['extent', 'resolution']
    # Every call has to start from an empty path cache.
    number = 1

    def setup(self, extent, resolution):
        feature = cfeature.NaturalEarthFeature(
            'physical', 'land', resolution, edgecolor='black', facecolor='0.8')
        # Read the shapefile outside of the timing.
        list(feature.geometries())

        fig = plt.figure()
        ax = fig.add_subplot(projection=ccrs.LambertConformal(
            central_longitude=-92, central_latitude=29))
        if extent == 'regional':
            # Gulf of Mexico, a small part of the globe spanning land records.
            ax.set_extent([-98, -84, 25, 31], crs=ccrs.PlateCarree())
        else:
            ax.set_global()
        ax.add_feature(feature)
        self.figure = fig

        FeatureArtist._geom_key_to_geometry_cache.clear()
        FeatureArtist._geom_key_to_path_cache.clear()

    def teardown(self, extent, resolution):
        plt.close(self.figure)

    def time_draw(self, extent, resolution):
        self.figure.canvas.draw()


class DrawNaturalEarthFeaturesLazily:
    """
    The example from issue #2102: 10m ocean, land and borders around Italy.

    Like in a script, the geometries are only read when the figure is drawn,
    which also replaces the CRS of each feature with the one from its
    shapefile.

    """
    params = ['PlateCarree', 'Mercator']
    param_names = ['projection']
    number = 1
    repeat = (1, 3, 120.0)
    timeout = 600

    def setup(self, projection):
        features = [
            cfeature.OCEAN.with_scale('10m'),
            cfeature.LAND.with_scale('10m'),
            cfeature.BORDERS.with_scale('10m'),
        ]
        # Download the shapefiles if necessary, but read them in the timing.
        for feature in features:
            shpreader.natural_earth(resolution='10m', category=feature.category,
                                    name=feature.name)
        cfeature._NATURAL_EARTH_GEOM_CACHE.clear()
        FeatureArtist._geom_key_to_geometry_cache.clear()
        FeatureArtist._geom_key_to_path_cache.clear()

        fig = plt.figure(figsize=(15, 10))
        if projection == 'PlateCarree':
            proj = ccrs.PlateCarree()
        else:
            proj = ccrs.Mercator(central_longitude=5, min_latitude=30,
                                 max_latitude=75)
        ax = fig.add_subplot(projection=proj)
        ax.set_extent([6, 19, 36, 48], ccrs.PlateCarree())
        ax.add_feature(features[0], facecolor='#2081C3', zorder=2)
        ax.add_feature(features[1], facecolor='lightgray', zorder=0)
        ax.add_feature(features[2], edgecolor='black', zorder=2)
        self.figure = fig

    def teardown(self, projection):
        plt.close(self.figure)

    def time_savefig(self, projection):
        self.figure.savefig(io.BytesIO(), format='png', bbox_inches='tight',
                            dpi=100)
