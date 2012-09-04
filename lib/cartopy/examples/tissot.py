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
This example demonstrates the way a box is warped when it is defined in a rotated space. 


XXX Keep the meaning, but re-word the following line: 
Notice that the box contains the north pole, adding extra complexity to the underlying transformation. 


"""
# XXX Needs Geodetic functionality from Basemap to work. Consider using http://geographiclib.sourceforge.net

import matplotlib.pyplot as plt

import cartopy
import numpy

def main():
    pc = cartopy.prj.PlateCarree()
    
    ax = plt.axes(projection=cartopy.prj.Mercator())
    ax.coastlines()
    ax.set_global()
#    g = ax.projection.as_geodetic()
    radius = 15 # degrees
    
    for lat in numpy.arange(-80, 90, 20):
        for lon in numpy.arange(-160, 180, 40):
            ax.geod_circle_meters(lon, lat, 500e3, facecolor='blue', alpha=0.7)
    plt.show()

if __name__ == '__main__':
    main()
