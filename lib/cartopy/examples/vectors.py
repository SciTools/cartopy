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


import matplotlib.pyplot as plt
import numpy as np

import cartopy


def sample_data(shape=(73, 145)):
    """Returns lons, lats, u data and v data of some interesting data on a regular grid."""
    nlats, nlons = shape 
    lats = np.linspace(-np.pi/2, np.pi/2, nlats)
    lons = np.linspace(0, 2*np.pi, nlons)
    lons, lats = np.meshgrid(lons, lats)
        
    vortex_spacing = 0.5
    extra_factor = 2.
    
    a = np.array([1, 0]) * vortex_spacing
    b = np.array([np.cos(np.pi/3), np.sin(np.pi/3)]) * vortex_spacing
    rnv = int(2 * extra_factor / vortex_spacing)
    vortices = [n * a + m * b for n in xrange(-rnv, rnv) for m in xrange(-rnv, rnv)]
    in_range = lambda x_y: -extra_factor < x_y[0] < extra_factor and -extra_factor < x_y[1] < extra_factor
    vortices = filter(in_range, vortices)
   
    
    xs = np.linspace(-1,1,nlons).astype(np.float32)[None, :]
    ys = np.linspace(-1,1,nlats).astype(np.float32)[:, None]
    udata = np.zeros((nlats, nlons), dtype=np.float32)
    vdata = np.zeros((nlats, nlons), dtype=np.float32)
    for (x,y) in vortices:
        rsq = (xs-x)**2 + (ys-y)**2
        udata +=  (ys-y) / rsq
        vdata += -(xs-x) / rsq
    
    lats = np.rad2deg(lats)
    lons = np.rad2deg(lons)    
    
    return lons, lats, udata, vdata


def main():
    # XXX do a vector plot
    lons, lats, udata, vdata = sample_data()
    
#    plt.streamplot(lons, lats, udata, vdata)
    ax = plt.axes(projection=cartopy.prj.PlateCarree())
    ax.streamplot(*sample_data(), transform=cartopy.prj.PlateCarree())
    ax.coastlines()

    plt.show()
    
if __name__ == '__main__':
    main()