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

import cartopy
from cartopy.examples.waves import sample_data


def main():
    pc = cartopy.prj.PlateCarree()
    rob = cartopy.prj.Robinson()
    igh = cartopy.prj.InterruptedGoodeHomolosine()
    nps = cartopy.prj.NorthPolarStereo()
    
    x, y, z = sample_data()
    slices = slice(20, 55), slice(85, 115)
    x = x.__getitem__(slices)
    y = y.__getitem__(slices)
    z = z.__getitem__(slices)
    
    plt.subplot(2, 2, 1, projection=pc)
    plt.scatter(x, y, c=z, transform=pc)
    
    plt.subplot(2, 2, 2, projection=rob)
    plt.scatter(x, y, c=z, transform=pc)
    
    # XXX Ask the mpl mailing list to find out how you might create a subplot and subsequently modify it's projection.
#    plt.subplot(2, 2, 3, )#projection=igh)
#    plt.scatter(x, y, c=z, transform=pc)
    
    plt.subplot(2, 2, 4, projection=nps)
    plt.scatter(x, y, c=z, transform=pc)
    
    plt.show()

if __name__ == '__main__':
    main()
