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
Notice that the box contains the north pole, adding extra complexity to the underlying transformation. 

"""

import matplotlib.pyplot as plt

import cartopy.crs as ccrs

def main():
    rp = ccrs.RotatedPole(pole_latitude=45, pole_longitude=180)
    pc = ccrs.PlateCarree()
    
    ax = plt.subplot(211, projection=rp)
    ax.bluemarble()
    ax.coastlines()
    x, y = [-44, -44, 45, 45, -44], [-45, 45, 45, -45, -45]
    ax.plot(x, y, marker='o', transform=rp)
    ax.fill(x, y, color='coral', transform=rp, alpha=0.4)
    ax.gridlines(15)   
    
    
    ax = plt.subplot(212, projection=pc)
    ax.bluemarble()
    ax.coastlines()
    ax.plot(x, y, transform=rp)
    ax.fill(x, y, transform=rp, color='coral', alpha=0.4)
    ax.gridlines(15)
    plt.show()

if __name__ == '__main__':
    main()
