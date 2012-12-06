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
"""
Provides the :class:`WMS` class for use with :func:`cartopy.mpl_integration.geoaxes.GeoAxes.add_image`.

"""

from PIL import Image, ImageOps
import urllib

import StringIO
from copy import deepcopy
import math
import warnings


class WMS(object):
    """
    Image retrieval factory for use with :func:`cartopy.mpl_integration.geoaxes.GeoAxes.add_image`.
    
    Initialised with a server and retrieval parameters,
    the image is retrieved by the "image_for_domain" method, when the data limits are known.
    
    """
    
    def __init__(self, server, layers, origin='lower', **kwargs):
        """
        Create the factory with the given retrieval parameters.
        
        Args:
        
            * server    -    Server name.
            * layers    -    Comma-separated names of required layers (server-specific).
            
        Kwargs:

            * origin    -    The WMS image origin: 'upper' or 'lower' (default).
        
        Other keyword arguments are passed through to :func:`WmsFactory.retrieve`. 

        """
        if "srs" in kwargs:
            warnings.warn("The srs keyword is forbidden as it is set by GeoAxes.")
        
        self.server = server
        self.layers = layers
        self.origin = origin
        self.kwargs = kwargs
        
    def image_for_domain(self, target_domain, dummy_target_z, srs=None, image_size=None):
        """
        Called at draw-time by :func:`cartopy.mpl_intergration.geoaxes.GeoAxes.draw`.
        
        Args:
        
            * target_domain     -    Data extent, described by a LineString.
            * dummy_target_z    -    Unused.

        Kwargs:

            * srs               -    WMS SRS string.
            * image_size        -    Required image resolution.
        
        """
        # Format the bounds.
        bbox = [math.floor(i) for i in target_domain.bounds]
        
        # Copy the kwargs and add a 'size' keyword if necessary.
        kwargs = deepcopy(self.kwargs)
        if "size" not in kwargs.keys():
            if image_size is not None:
                kwargs["size"] = image_size
            else:
                raise ValueError("Expected image_size")

        # Add a srs keyword.
        if srs is None:
            raise ValueError("Expected a srs")
        kwargs["srs"] = srs

        # Get the image from the server.
        img = self.retrieve(bbox, **kwargs)
        return img, (bbox[0], bbox[2], bbox[1], bbox[3]), 'lower'
    
    def retrieve(self, bbox, styles=None, size=None, transparent=None, srs=None):
        """
        Retrieve a WMS jpeg image from our server.
    
        Args:
        
            * bbox    -    Bounding box for the image, in projection coordinates.
    
        Kwargs:
        
            * styles       - List of styles, as available for the requested layers.
            * size         - Image resolution.
            * transparent  - Useful for "overlay" images, such as political borders.
                             Default is False.
            * srs          - WMS SRS string.
    
        Returns:
    
            * A PIL :class:`~PIL.Image.Image`
    
        """
        if isinstance(styles, basestring):
            styles = [styles]
            
        # Construct request string.
        request = self.server
        if not request.endswith("?"):
            request += "?"

        request += "&request=getmap"
        request += "&version=1.1.1"
        request += "&layers={}".format(self.layers)
        request += "&bbox={}".format(",".join([str(i) for i in bbox]))

        if styles is not None:                
            request += "&styles={}".format(",".join([str(s) for s in styles]))
        else:
            request += "&styles="

        if size is None:
            size = (256, 256)                
        request += "&width={}&height={}".format(size[0], size[1])

        if transparent is not None:                
            request += "&transparent={}".format(int(transparent))

        if srs is not None:                
            request += "&srs={}".format(srs.as_wms_srs())
            
        # Get jpeg from server.
        request += "&format=image/jpeg"
        jpeg_bytes = urllib.urlopen(request).read()
        pil_img = Image.open(StringIO.StringIO(jpeg_bytes))
        if self.origin.lower() == 'upper':
            pil_img = ImageOps.flip(pil_img)
        
        return pil_img
