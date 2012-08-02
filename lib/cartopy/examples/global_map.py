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
    pc = ccrs.PlateCarree()
    rob = ccrs.Robinson()
    sph = ccrs.Geodetic()
    
    ax = plt.axes(projection=rob)
#    ax = plt.axes(projection=pc)
    
#    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    
    ax.set_global()
    
    from cartopy.examples.waves import sample_data
    
    x, y, data = sample_data((20, 40))
#    ax.contourf(x, y, data, transform=pc, alpha=0.3)
#    ax.contour(x, y, data, transform=pc, alpha=0.3)
#    ax.contourf(x, y, data, 3, transform=pc, alpha=0.3)
#    ax.contourf(x, y, data, 5, transform=pc, alpha=0.3)
    #print 'getting domain'
    #print ax.native_extents()
    #print ax.map_domain(ccrs.PlateCarree())
    #print ax.ll_boundary_poly()
#    ax.stock_img('bluemarble')
#    ax.coastlines()
    ax.gshhs_line()
#    ax.coastlines_land()
    
    plt.plot(-0.08, 51.53, 'o', transform=pc)
    plt.plot([-0.08, 132], [51.53, 43.17], transform=pc)
    plt.plot([-0.08, 132], [51.53, 43.17], transform=sph) 
    
    
    
#    ax.gshhs_line(resolution='coarse', domain=ax.boundary_poly())
    
    plt.show()

if __name__ == '__main__':
    main()
