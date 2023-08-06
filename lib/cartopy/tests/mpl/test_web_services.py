# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
from cartopy.io.ogc_clients import _OWSLIB_AVAILABLE


@pytest.mark.filterwarnings("ignore:TileMatrixLimits")
@pytest.mark.network
@pytest.mark.skipif(not _OWSLIB_AVAILABLE, reason='OWSLib is unavailable.')
@pytest.mark.mpl_image_compare(filename='wmts.png', tolerance=0.03)
def test_wmts():
    ax = plt.axes(projection=ccrs.PlateCarree())
    url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
    # Use a layer which doesn't change over time.
    ax.add_wmts(url, 'MODIS_Water_Mask')
    return ax.figure


@pytest.mark.network
@pytest.mark.skipif(not _OWSLIB_AVAILABLE, reason='OWSLib is unavailable.')
def test_wms_tight_layout():
    ax = plt.axes(projection=ccrs.PlateCarree())
    url = 'http://vmap0.tiles.osgeo.org/wms/vmap0'
    layer = 'basic'
    ax.add_wms(url, layer)
    ax.figure.tight_layout()


@pytest.mark.network
@pytest.mark.skipif(not _OWSLIB_AVAILABLE, reason='OWSLib is unavailable.')
@pytest.mark.mpl_image_compare(filename='wms.png', tolerance=0.02)
def test_wms():
    ax = plt.axes(projection=ccrs.Orthographic())
    url = 'http://vmap0.tiles.osgeo.org/wms/vmap0'
    layer = 'basic'
    ax.add_wms(url, layer)
    return ax.figure
