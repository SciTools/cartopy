# (C) British Crown Copyright 2011 - 2012, Met Office
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
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import absolute_import

import cartopy.io.ogc_clients as ogc
from owslib.wms import WebMapService
from owslib.wmts import ContentMetadata, WebMapTileService
import unittest
import cartopy.crs as ccrs
import numpy as np


class test_WMSRasterSource(unittest.TestCase):
    URI = 'http://vmap0.tiles.osgeo.org/wms/vmap0'
    layer = 'basic'
    layers = ['basic', 'ocean']
    projection = ccrs.PlateCarree()

    def test_string_service(self):
        source = ogc.WMSRasterSource(self.URI, self.layer)
        self.assertIsInstance(source.service, WebMapService)
        self.assertIsInstance(source.layers, list)
        self.assertEqual(source.layers, [self.layer])

    def test_wms_service_instance(self):
        service = WebMapService(self.URI)
        source = ogc.WMSRasterSource(service, self.layer)
        self.assertIs(source.service, service)

    def test_multiple_layers(self):
        source = ogc.WMSRasterSource(self.URI, self.layers)
        self.assertEqual(source.layers, self.layers)

    def test_no_layers(self):
        msg = 'One or more layers must be defined.'
        with self.assertRaisesRegexp(ValueError, msg):
            ogc.WMSRasterSource(self.URI, [])

    def test_extra_kwargs_empty(self):
        source = ogc.WMSRasterSource(self.URI, self.layer,
                                     getmap_extra_kwargs={})
        self.assertEqual(source.getmap_extra_kwargs, {})

    def test_extra_kwargs_None(self):
        source = ogc.WMSRasterSource(self.URI, self.layer,
                                     getmap_extra_kwargs=None)
        self.assertEqual(source.getmap_extra_kwargs, {'transparent': True})

    def test_extra_kwargs_non_empty(self):
        kwargs = {'another': 'kwarg'}
        source = ogc.WMSRasterSource(self.URI, self.layer,
                                     getmap_extra_kwargs=kwargs)
        self.assertEqual(source.getmap_extra_kwargs, kwargs)

    def test_supported_projection(self):
        source = ogc.WMSRasterSource(self.URI, self.layer)
        source.validate_projection(self.projection)

    def test_unsupported_projection(self):
        source = ogc.WMSRasterSource(self.URI, self.layer)
        msg = 'was not convertible to a suitable WMS SRS.'
        with self.assertRaisesRegexp(ValueError, msg):
            source.validate_projection(ccrs.Miller())

    def test_fetch_img(self):
        source = ogc.WMSRasterSource(self.URI, self.layer)
        extent = [-10, 10, 40, 60]
        img, extent_out = source.fetch_raster(self.projection, extent,
                                              (30, 30))
        img = np.array(img)
        self.assertEqual(img.shape, (30, 30, 4))
        # No transparency in this image.
        self.assertEqual(img[:, :, 3].min(), 255)
        self.assertEqual(extent, extent_out)


class test_WMTSRasterSource(unittest.TestCase):
    URI = 'http://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
    layer_name = 'VIIRS_CityLights_2012'
    projection = ccrs.PlateCarree()

    def test_string_service(self):
        source = ogc.WMTSRasterSource(self.URI, self.layer_name)
        self.assertIsInstance(source.wmts, WebMapTileService)
        self.assertIsInstance(source.layer, ContentMetadata)
        self.assertEqual(source.layer.name, self.layer_name)

    def test_wmts_service_instance(self):
        service = WebMapTileService(self.URI)
        source = ogc.WMTSRasterSource(service, self.layer_name)
        self.assertIs(source.wmts, service)

    def test_supported_projection(self):
        source = ogc.WMTSRasterSource(self.URI, self.layer_name)
        source.validate_projection(self.projection)

    def test_unsupported_projection(self):
        source = ogc.WMTSRasterSource(self.URI, self.layer_name)
        msg = 'Unable to find tile matrix for projection.'
        with self.assertRaisesRegexp(ValueError, msg):
            source.validate_projection(ccrs.Miller())

    def test_fetch_img(self):
        source = ogc.WMTSRasterSource(self.URI, self.layer_name)
        extent = [-10, 10, 40, 60]
        img, extent_out = source.fetch_raster(self.projection, extent,
                                              (30, 30))
        img = np.array(img)
        self.assertEqual(img.shape, (512, 512, 4))
        # No transparency in this image.
        self.assertEqual(img[:, :, 3].min(), 255)
        self.assertEqual((-180.0, 107.99999999999994,
                          -197.99999999999994, 90.0), extent_out)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
