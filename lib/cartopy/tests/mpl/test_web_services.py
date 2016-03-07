# (C) British Crown Copyright 2014 - 2016, Met Office
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

import unittest

import matplotlib.pyplot as plt

from cartopy.tests.mpl import ImageTesting
import cartopy.crs as ccrs
from cartopy.io.ogc_clients import _OWSLIB_AVAILABLE


@unittest.skipIf(not _OWSLIB_AVAILABLE, 'OWSLib is unavailable.')
@ImageTesting(['wmts'])
def test_wmts():
    ax = plt.axes(projection=ccrs.PlateCarree())
    url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
    # Use a layer which doesn't change over time.
    ax.add_wmts(url, 'MODIS_Water_Mask')


@unittest.skipIf(not _OWSLIB_AVAILABLE, 'OWSLib is unavailable.')
@ImageTesting(['wms'])
def test_wms():
    ax = plt.axes(projection=ccrs.Orthographic())
    url = 'http://vmap0.tiles.osgeo.org/wms/vmap0'
    layer = 'basic'
    ax.add_wms(url, layer)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
