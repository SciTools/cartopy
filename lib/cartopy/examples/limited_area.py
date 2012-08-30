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
    
#    ax = plt.subplot(211, projection=rp)
#    ax.bluemarble()
#    ax.coastlines()
    # XXX The change in projection is a result of plt.plot making a poly or a line depending on some things...
    x, y = [-44, -44, 45, 45], [-45, 45, 45, -45]
#    x, y = [-44, -44, 45, 45, -44], [-45, 45, 45, -45, -45]
#    # XXX why is the horizontal line on the native projection not straight? BUG!
#    ax.plot(x, y, marker='o', transform=rp)
#    ax.gridlines()   
    
    
#    ax = plt.subplot(212, projection=pc)
    ax = plt.axes(projection=pc)
    ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
#    ax = plt.axes(projection=rp)
    
    x, y = [-20, -20, 20, 20], [-20, 20, 20, -20]
#    x, y = [-40, -40, 20, 20], [-20, 20, 20, -20]
    ax.plot(x, y, transform=rp)
    
    
    #ax.ll_boundary_poly_draw()
    
#    ax.stock_img('bluemarble')
#    ax.coastlines()
#    ax.ll_extents_draw()
    #ax.gshhs_line(outline_color='gray')
    
    plt.show()
    
    return
    
    
#    ax.bluemarble()
#    ax.coastlines()
#    ax.gshhs_line()

    x, y = [-20, -20, 20, 20], [-20, 20, 20, -20]
    ax.plot(x, y, transform=rp)
    x, y = [20, -20], [-20, -20]
    ax.plot(x, y, transform=rp)
#    ax.gridlines()

    print ax.ll_extents()
#    ax.set_xlim(-50, 50)
#    ax.set_ylim(30, 70)
    
    
    x1, y1, x2, y2 =  ax.ll_boundary_poly().bounds
    from matplotlib.patches import Rectangle    
    rect = Rectangle([x1, y1], x2-x1, y2-y1, edgecolor='gray', facecolor='none')
    ax.add_patch(rect)
    
#    ax.gshhs_line(outline_color='gray')
    #ax.gshhs_line()
    
    ax.stock_img('ne_shaded')
    
    plt.draw()
    plt.show()

if __name__ == '__main__':
    main()
