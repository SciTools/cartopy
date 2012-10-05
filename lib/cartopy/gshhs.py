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



# PLEASE NOTE: DUE TO SOME MPL RELATED ISSUES, THE GSHHS SUPPORT HAS BEEN DISABLED. 
#              IT IS ANTICIPATED THAT BY 0.5 THERE WILL BE A CLEAN AND TIDY INTERFACE
#              TO USE THIS USEFUL DATASET. - pelson


 
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.collections import PatchCollection
import matplotlib.cm
import numpy
import os

from shapely.geometry import Polygon

# XXX Make the data dir configurable
project_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
data_dir = os.path.join(project_dir, 'data')

gshhs_data_dir = os.path.join(data_dir, 'gshhs')

fnames = {
          'coarse': os.path.join(gshhs_data_dir, 'gshhs_c.b'),
          'low': os.path.join(gshhs_data_dir, 'gshhs_l.b'),
          'intermediate': os.path.join(gshhs_data_dir, 'gshhs_i.b'),
          'high': os.path.join(gshhs_data_dir, 'gshhs_h.b'),
          'full': os.path.join(gshhs_data_dir, 'gshhs_f.b'),
          }


def read_gshhc(filename, poly=True, domain=None, filter_predicate=None):
    """
    Reads:
    
    Global Self-consistent Hierarchical High-resolution Shorelines
        version 2.0 July 15, 2009
        
    .. seealso:: http://www.soest.hawaii.edu/pwessel/gshhs/README.TXT
    
    XXX: Return internal polygons when appropriate
    
    """
    DEBUG = False
    
    fh = open(filename, 'rb')
    
    #(0, 360, -90, 90)
    if domain is None:
        domain = Polygon([[0, -90], [360, -90], [360, 90], [0, 90], [0, -90]])
    
    extent_w, extent_s, extent_e, extent_n = [v * 1e6 for v in domain.bounds]
    
#    corners = [extent_w, extent_n], [extent_w, extent_s], [extent_e, extent_s], [extent_e, extent_n]
#    
#    poly_extent = Polygon(numpy.array(corners) / 1e6)
    poly_extent = domain
    
    i=-1
    # XXX
    while True:
        i += 1
#    for i in xrange(10000):
#        if i % 10000 == 1: print i
        
        header = numpy.fromfile(fh, dtype='>i4', count=11)
        
        # If no header was received, we are at the end of the file
        if len(header) == 0:
            break
        
        if DEBUG:
            if i not in ([
#                      0, # Europe & Asia
#                      1, # Africa
#                      2, # USA
#                      3, # S.America
#                      4, # Antarctic
                      14, # UK
#                      25, # Ireland
                      ]):
                continue

        flag = header[2]
        crosses_greenwich = (flag >> 16) & 1
               
        
        flag = header[2]
        level = flag & 255

        # ###########################
        # Filter the shapes by extent
        # ###########################
        
        # get the maximum extent in microdegrees
        w, e, south, north = header[3:7]
        in_x_range = False
        in_y_range = False
        
        # handle the case where the uk has an extent of -6230861 1765806 and Ireland has an extent of 349515833 354569167
        # XXX I'm sure this could be done more cleanly
        for off in range(2):
            west = w - 360 * 1e6 * off
            east = e - 360 * 1e6 * off 
            in_x_range = in_x_range or (extent_w <= west <= extent_e or extent_w <= east <= extent_e or (east >= extent_e and west <= extent_w))
            in_y_range = in_y_range or (extent_s <= south <= extent_n or extent_s <= north <= extent_n or (north >= extent_n and south <= extent_s))
            
        if not (in_x_range and in_y_range):
            if DEBUG: print in_x_range, in_y_range, w, e, south, north, extent_w, extent_e
            fh.seek(header[1]*2 * 4, 1)
            continue
        else:
            if DEBUG: print in_x_range, in_y_range, w, e, south, north, extent_w, extent_e
        
        
        points = numpy.fromfile(fh, dtype='>i4', count=header[1]*2) * 1.0e-6
        points = points.reshape(-1, 2)
        
        intersects = False
        for off in range(2):
##            west = points - numpy.array([[360 * off, 0]])
#            east = points - numpy.array([[360 * off, 0]])
            poly_shape = Polygon(points - numpy.array([[360 * off, 0]]))
#            print (points - numpy.array([[360 * off, 0]]))[:10, ...]
#            print corners
#            print 'intersect? ', i, off*360, poly_shape.intersects(poly_extent)

            intersects = intersects or poly_shape.intersects(poly_extent)
        
        if not intersects:
            continue
        
        lons, lats = points[:, 0:1], points[:, 1:2]

        if poly:
            if ( level == 1 and 
                 points.shape[0] > 4
                ):
                yield points
        else:
            yield points
            
#        break
#            yield header, lons, lats
#        if points.shape[0] > 4:
#            yield header, lons, lats
#            yield points
        
#        if crosses_greenwich:
#            # If the greenwich has been crossed, then 360 is added to any number below 0 in this format.
#            # To fix this, identify any points which are more than 180 degrees apart, using this information we can identify
#            # polygon groups and shift them appropriately.  
#            delta = numpy.diff(lons)
#            step = numpy.where(numpy.abs(delta) > 180)[0]
#            step = [0] + list(step+1) + [None]
#            for s1, s2 in zip(step[:-1] , step[1:]):
#                if delta[s1-1] > 180:
#                    lons[s1:s2] -= 360
#            
#        if i == 4:
#            # antarctic
#            lons = numpy.array(list(lons) + [lons[-1], lons[0], lons[0]])
#            lats = numpy.array(list(lats) + [-90, -90, lats[0]])

#        yield header, lons, lats

