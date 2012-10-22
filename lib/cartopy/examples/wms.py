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
This example retrieves map images from the World OSM WMS server. 

First, we interrogate the server by typing the following into a browser: 

    http://129.206.228.72/cached/osm?request=getCapabilities

This will return a text file containing the layer names available for use
in the code below (search for the string "layer").

"""


import matplotlib.pyplot as plt

import cartopy.crs as ccrs
from cartopy.io.wms import WMS


def main():
    plt.figure(figsize=(12,10))

    # Create a WMS image retriever, pointing to the World OSM WMS server. 
    server = "http://129.206.228.72/cached/osm?"
    layers = "osm_auto:all"
    origin = 'upper'
    wms = WMS(server, layers, origin=origin)
    
    # Global
    ax = plt.subplot(2,2,1, projection=ccrs.PlateCarree())
    ax.add_image(wms)
    
    # Europe
    ax = plt.subplot(2,2,2, projection=ccrs.PlateCarree(), xlim=[-15, 25], ylim=[30,70])
    ax.add_image(wms)
    
    # Britain
    ax = plt.subplot(2,2,3, projection=ccrs.PlateCarree(), xlim=[-13, 3], ylim=[50,60])
    ax.add_image(wms)
    
    # South-west Britain
    ax = plt.subplot(2,2,4, projection=ccrs.PlateCarree(), xlim=[-7, -3], ylim=[49,53])
    ax.add_image(wms)
    
    plt.show()
    
    
if __name__ == '__main__':
    main()
