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


import cartopy
import cartopy.mpl_integration.patch as pt
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy


import matplotlib.textpath
import matplotlib.patches
from matplotlib.font_manager import FontProperties


def main():
    pc = cartopy.prj.PlateCarree()
    rp = cartopy.prj.Robinson()
        
    plt.figure(figsize=[12, 6])
    ax = plt.axes([0, 0, 1, 1], projection=rp)
    
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_global()
    # create a transform from PlateCarree to the axes' projection
    pc_t = pc._as_mpl_transform(ax)
    
    ax.coastlines()
    ax.gridlines()
    im = ax.bluemarble()
    
    name, pos = 'pycart', (-145, -37.5)
    name, pos = 'cartopy', (-180, -30)    
    logo_path = matplotlib.textpath.TextPath(pos, name, size=1, prop=FontProperties(family='Arial', weight='bold'))
    # put the letters in the right place
    logo_path._vertices *= numpy.array([95, 160])
    
    im.set_clip_path(logo_path, transform=pc_t)

#    # add the letters as patches...
#    # make a dictionary of letters to paths (slightly round about by converting to geos and back)
#    letter_paths = {letter: path for letter, path in zip(name, pt.geos_to_path(pt.path_to_geos(logo_path)))}
#    o_path = letter_paths.pop('o')
#    for path in letter_paths.values():
#        ax.add_patch(matplotlib.patches.PathPatch(path, color='purple', transform=pc_t, zorder=2))
    
    plt.show()
        
if __name__ == '__main__':
    main()