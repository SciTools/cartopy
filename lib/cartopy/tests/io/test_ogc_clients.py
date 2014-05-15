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
import unittest
import cartopy.crs as ccrs
import numpy as np


class test_WMSRasterSource(unittest.TestCase):
    URI = 'http://vmap0.tiles.osgeo.org/wms/vmap0'
    layer = 'basic'
    layers = ['basic', 'ocean']
    projection = ccrs.PlateCarree()

    def test_string_service(self):
        fetcher = ogc.WMSRasterSource(self.URI, self.layer, self.projection)
        self.assertIsInstance(fetcher.service, WebMapService)
        self.assertIsInstance(fetcher.layers, list)
        self.assertEqual(fetcher.layers, [self.layer])

    def test_wms_service_instance(self):
        service = WebMapService(self.URI)
        fetcher = ogc.WMSRasterSource(service, self.layer)
        self.assertIs(fetcher.service, service)

    def test_multiple_layers(self):
        fetcher = ogc.WMSRasterSource(self.URI, self.layers)
        self.assertEqual(fetcher.layers, self.layers)

    def test_no_layers(self):
        msg = 'One or more layers must be defined.'
        with self.assertRaisesRegexp(ValueError, msg):
            ogc.WMSRasterSource(self.URI, [])

    def test_extra_kwargs_empty(self):
        fetcher = ogc.WMSRasterSource(self.URI, self.layer,
                                      getmap_extra_kwargs={})
        self.assertEqual(fetcher.getmap_extra_kwargs, {})

    def test_extra_kwargs_None(self):
        fetcher = ogc.WMSRasterSource(self.URI, self.layer,
                                      getmap_extra_kwargs=None)
        self.assertEqual(fetcher.getmap_extra_kwargs, {'transparent': True})

    def test_extra_kwargs_non_empty(self):
        kwargs = {'another': 'kwarg'}
        fetcher = ogc.WMSRasterSource(self.URI, self.layer,
                                      getmap_extra_kwargs=kwargs)
        self.assertEqual(fetcher.getmap_extra_kwargs, kwargs)

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


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
