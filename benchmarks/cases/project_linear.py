# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs


class Suite:
    params = [
        ('PlateCarree', 'NorthPolarStereo', 'Robinson',
         'InterruptedGoodeHomolosine'),
        ('110m', '50m'),
    ]
    param_names = ['projection', 'resolution']

    def setup(self, projection, resolution):
        shpfilename = shpreader.natural_earth(
            resolution=resolution, category='physical', name='ocean')
        reader = shpreader.Reader(shpfilename)
        oceans = list(reader.geometries())
        self.geoms = oceans[0]

        self.projection = getattr(ccrs, projection)()

    def time_project_linear(self, projection, resolution):
        self.projection.project_geometry(self.geoms)
