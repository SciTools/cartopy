# (C) British Crown Copyright 2014, Met Office
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
"""
Implements RasterFetcher classes which can retrieve imagery from web services
such as WMS and WMTS.

"""
from __future__ import absolute_import

from io import BytesIO

from PIL import Image
from owslib.wms import WebMapService

from cartopy.io import RasterFetcher
import cartopy.crs as ccrs


# Hardcode some known EPSG codes for now.
_CRS_TO_OGC_SRS = {ccrs.PlateCarree(): 'EPSG:4326'
                   }


class WMSFetcher(RasterFetcher):

    """
    A WMS imagery retriever which can be added to a map.

    .. note::

        No caching of retrieved maps is done with this WMSFetcher.

        To reduce load on the WMS server it is encouraged to tile
        map requests and subsequently stitch them together to recreate
        a single raster, thus allowing for a more aggressive caching scheme,
        but this WMSFetcher does not currently implement WMS tile fetching.

        Whilst not the same service, there is also a WMTSFetcher which makes
        use of tiles and comes with built-in caching for fast repeated
        map retrievals.

    """
    def __init__(self, service, layers, projection, getmap_extra_kwargs=None):
        """
        Parameters
        ----------
        service : string or WebMapService instance
            The WebMapService instance, or URL of a WMS service, from whence
            to retrieve the image.
        layers : string or list of strings
            The name(s) of layers to use from the WMS service.
        projection : :class:`cartopy.crs.Projection` instance
            The projection of the resulting images.
        getmap_extra_kwargs : dict or None
            Extra keywords to pass through to the service's getmap method.
            If None, a dictionary with ``{'transparent': True}`` will
            be defined.

        """
        if isinstance(service, basestring):
            service = WebMapService(service)

        if isinstance(layers, basestring):
            layers = [layers]

        if getmap_extra_kwargs is None:
            getmap_extra_kwargs = {'transparent': True}

        if len(layers) == 0:
            raise ValueError('One or more layers must be defined.')
        for layer in layers:
            if layer not in service.contents:
                raise ValueError('The {} layer does not exist in '
                                 'this service.'.format(layer))

        #: The OWSLib WebMapService instance.
        self.service = service

        #: The name of the layers to fetch.
        self.layers = layers

        #: Extra kwargs passed through to the service's getmap request.
        self.getmap_extra_kwargs = getmap_extra_kwargs

        srs = _CRS_TO_OGC_SRS.get(projection)
        if srs is None:
            raise ValueError('The projection {!r} was not convertible to a '
                             'suitable WMS SRS.'.format(projection))

        for layer in self.layers:
            if srs not in self.service.contents[layer].crsOptions:
                raise ValueError('The SRS {} is not a valid SRS for the '
                                 '{!r} WMS layer.'.format(srs, layer))

        self._srs = srs
        self.projection = projection

    def fetch_raster(self, extent, target_resolution):
        service = self.service
        min_x, max_x, min_y, max_y = extent

        wms_image = service.getmap(layers=self.layers, srs=self._srs,
                                   bbox=(min_x, min_y, max_x, max_y),
                                   size=target_resolution, format='image/png',
                                   **self.getmap_extra_kwargs)

        wms_image = Image.open(BytesIO(wms_image.read()))
        return wms_image, extent
