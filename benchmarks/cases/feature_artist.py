# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
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
