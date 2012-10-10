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

import numpy as np
from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
import shapely.geometry

import cartopy.crs as ccrs
import cartopy.mpl_integration.patch as cpatch
import cartopy.io.shapereader
import cartopy.mpl_integration.geoaxes as cgeoaxes
from cartopy.examples.waves import sample_data
import cartopy.mpl_integration.patch

    
from cartopy.tests.mpl import image_comparison



class CallCounter(object):
    """
    Exposes a context manager which can count the number of calls to a specific
    function. (useful for cache checking!)
    
    Internally, the target function is replaced with a new one created by this context
    manager which then increments ``self.count`` every time it is called.
    
    Example usage::
    
        show_counter = CallCounter(plt, 'show')   
        with show_counter:
            plt.show()
            plt.show()
            plt.show()
        
        print show_counter.count    # <--- outputs 3
    
    
    """
    def __init__(self, parent, function_name):
        self.count = 0
        self.parent = parent
        self.function_name = function_name
        self.orig_fn = getattr(parent, function_name)
        
    def __enter__(self):
        def replacement_fn(*args, **kwargs):
            self.count += 1
            return self.orig_fn(*args, **kwargs)
        
        setattr(self.parent, self.function_name, replacement_fn)
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        setattr(self.parent, self.function_name, self.orig_fn)


def test_coastline_loading_cache():
    # a5caae040ee11e72a62a53100fe5edc355304419 added coastline caching.
    # this test ensures it is working...

    
    # count the number of times shapereader.Reader is created.
    shapereader_counter = CallCounter(cartopy.io.shapereader.Reader, '__init__')
       
    with shapereader_counter:
        ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
        ax1.coastlines()
        ax2 = plt.subplot(2, 1, 1, projection=ccrs.Robinson())
        ax2.coastlines()
    
    assert shapereader_counter.count == 1, ('The shapereader Reader class was created '
                                            'more than (actually %s times) - '
                                            ' the caching is not working.' % 
                                            shapereader_counter.count)

    plt.close()


def test_shapefile_transform_cache():
    # a5caae040ee11e72a62a53100fe5edc355304419 added shapefile mpl geometry caching
    # based on geometry object id. This test ensures it is working...
    coastline_path = cartopy.io.shapereader.natural_earth(resolution="50m",
                                               category='physical',
                                               name='coastline')
    geoms = tuple(cartopy.io.shapereader.Reader(coastline_path).geometries())[:10]
    n_geom = len(geoms)
    
    ax = plt.axes(projection=ccrs.Robinson())
    
    project_geometry_counter = CallCounter(ax.projection, 'project_geometry')
    
    with project_geometry_counter:
        c = cartopy.io.shapereader.mpl_axes_plot(ax, geoms)
        c = cartopy.io.shapereader.mpl_axes_plot(ax, geoms)
        c = cartopy.io.shapereader.mpl_axes_plot(ax, geoms[:]) # <- caching will not work for this case
    
    # before the performance enhancement, the count would have been n_calls * n_geom,
    # but should now be just n_geom * n_unique_id_geoms_calls (i.e. in this case 2 * n_calls))
    assert project_geometry_counter.count == n_geom, ('The given geometry was transformed '
                                                      'too many times (expected: %s; got %s) - '
                                                      ' the caching is not working.' % 
                                                      (n_geom * 2, project_geometry_counter.count))
    
    plt.close()


def test_contourf_transform_path_counting():
    ax = plt.axes(projection=ccrs.Robinson())
    plt.draw()

    path_to_geos_counter = CallCounter(cartopy.mpl_integration.patch, 
                                         'path_to_geos')

    with path_to_geos_counter:
        x, y, z = sample_data((30, 60))
        cs = plt.contourf(x, y, z, 5, transform=ccrs.PlateCarree())
        plt.draw()
    
    n_geom = sum([len(c.get_paths()) for c in cs.collections])
    
    # before the performance enhancement, the count would have been 2 * n_geom,
    # but should now be just n_geom
    assert path_to_geos_counter.count == n_geom, ('The given geometry was transfomed '
                                                  'too many times (expected: %s; got %s) - '
                                                  ' the caching is not working.' % 
                                                  (n_geom, path_to_geos_counter.count))
    
    plt.close()
    

if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
