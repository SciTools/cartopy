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
    """Returns lons, lats and data of some interesting data on a regular grid."""
    nlats, nlons = shape 
    lats = np.linspace(-np.pi/2, np.pi/2, nlats)
    lons = np.linspace(0, 2*np.pi, nlons)
    lons, lats = np.meshgrid(lons, lats)
    wave = 0.75*(np.sin(2*lats)**8) * np.cos(4*lons)
    mean = 0.5*np.cos(2*lats) * ((np.sin(2*lats))**2 + 2)

    lats = np.rad2deg(lats)
    lons = np.rad2deg(lons)
    data = wave + mean
    
    return lons, lats, data


def main():
#    ax = plt.axes(projection=cartopy.prj.NorthPolarStereo())
    ax = plt.axes(projection=cartopy.prj.PlateCarree())
    plt.contourf(*sample_data()), #transform=cartopy.prj.PlateCarree())
    plt.show()
    

def main():
    ax = plt.axes(projection=cartopy.prj.NorthPolarStereo())
    ax.contourf(*sample_data(), transform=cartopy.prj.PlateCarree())
    plt.show()
    

if __name__ == '__main__':
    main()
