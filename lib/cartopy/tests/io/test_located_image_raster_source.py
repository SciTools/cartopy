# (C) British Crown Copyright 2018, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

import os

from matplotlib.pyplot import imread
import numpy as np

from cartopy import config
import cartopy.crs as ccrs
from cartopy.io import LocatedImage, LocatedImageRasterSource


RESOLUTION = (30, 30)


class TestLocatedImageRasterSource(object):
    fname = os.path.join(config["repo_data_dir"],
                         'raster', 'natural_earth',
                         '50-natural-earth-1-downsampled.png')
    extent = [-180, 180, -90, 90]
    projection = ccrs.PlateCarree()
    image = LocatedImage(imread(fname), extent)
    source = LocatedImageRasterSource(image, projection)

    def test_valid_projection(self):
        # Validation of projection should always return True as we can regrid
        # the raster source to any projection.
        tgt_proj = ccrs.Mercator()
        result = self.source.validate_projection(tgt_proj)
        assert result is True

    def test_fetch_raster(self):
        located_image, = self.source.fetch_raster(self.projection,
                                                  self.extent,
                                                  RESOLUTION)
        img = np.array(located_image.image)
        extent = located_image.extent
        assert isinstance(located_image, LocatedImage)
        assert img.shape == RESOLUTION + (3,)
        assert extent == self.extent

    def test_fetch_raster_limited_extent(self):
        limited_extent = [-15, 25, 45, 85]
        (_, extent), = self.source.fetch_raster(self.projection,
                                                limited_extent,
                                                RESOLUTION)
        assert extent == limited_extent

    def test_fetch_raster_different_projection(self):
        tgt_proj = ccrs.Orthographic()
        tgt_extent = [-570000, 5100000, 870000, 3500000]
        located_image, = self.source.fetch_raster(tgt_proj,
                                                  tgt_extent,
                                                  RESOLUTION)
        img = np.array(located_image.image)
        extent = located_image.extent
        assert img.shape == RESOLUTION + (3,)
        assert extent == tgt_extent

    def test_float_resolution(self):
        # The resolution (in pixels) should be rounded to ints.
        tgt_proj = ccrs.Orthographic()
        tgt_extent = [-570000, 5100000, 870000, 3500000]
        located_image, = self.source.fetch_raster(tgt_proj,
                                                  tgt_extent,
                                                  [19.5, 39.1])
        img = np.array(located_image.image)
        assert img.shape == (40, 20, 3)
