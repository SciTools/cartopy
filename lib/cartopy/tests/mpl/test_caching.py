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
import gc

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
import shapely.geometry

import cartopy.crs as ccrs
import cartopy.io.shapereader
import cartopy.mpl.geoaxes as cgeoaxes
import cartopy.mpl.patch
from cartopy.examples.waves import sample_data


class CallCounter(object):
    """
    Exposes a context manager which can count the number of calls to a specific
    function. (useful for cache checking!)

    Internally, the target function is replaced with a new one created
    by this context manager which then increments ``self.count`` every
    time it is called.

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
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        setattr(self.parent, self.function_name, self.orig_fn)


def test_coastline_loading_cache():
    # a5caae040ee11e72a62a53100fe5edc355304419 added coastline caching.
    # this test ensures it is working...

    # count the number of times shapereader.Reader is created.
    shapereader_counter = CallCounter(cartopy.io.shapereader.Reader,
                                      '__init__')

    with shapereader_counter:
        ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
        ax1.coastlines()
        ax2 = plt.subplot(2, 1, 1, projection=ccrs.Robinson())
        ax2.coastlines()

    msg = ('The shapereader Reader class was created more than (actually %s '
           'times) - the caching is not working.' % shapereader_counter.count)
    assert shapereader_counter.count == 1, msg

    plt.close()


def test_shapefile_transform_cache():
    # a5caae040ee11e72a62a53100fe5edc355304419 added shapefile mpl
    # geometry caching based on geometry object id. This test ensures
    # it is working...
    coastline_path = cartopy.io.shapereader.natural_earth(resolution="50m",
                                                          category='physical',
                                                          name='coastline')
    geoms = cartopy.io.shapereader.Reader(coastline_path).geometries()
    # filter just the first 10 of them
    geoms = tuple(geoms)[:10]
    n_geom = len(geoms)

    ax = plt.axes(projection=ccrs.Robinson())

    project_geometry_counter = CallCounter(ax.projection, 'project_geometry')

    # Capture the size of the cache before our test
    gc.collect()
    initial_cache_size = len(cgeoaxes._GEOMETRY_TO_PATH_CACHE)

    with project_geometry_counter:
        c = ax.add_geometries(geoms, ccrs.Geodetic())
        c = ax.add_geometries(geoms, ccrs.Geodetic())
        c = ax.add_geometries(geoms[:], ccrs.Geodetic())

    # Before the performance enhancement, the count would have been
    # n_calls * n_geom, but should now be just n_geom.
    msg = ('The given geometry was transformed too many times (expected: '
           '%s; got %s) - the caching is not working.'
           '' % (n_geom, project_geometry_counter.count))
    assert project_geometry_counter.count == n_geom, msg

    # Check the cache has an entry for each geometry.
    assert len(cgeoaxes._GEOMETRY_TO_PATH_CACHE) == initial_cache_size + n_geom

    # Check that the cache is empty again once we've dropped all references
    # to the source paths.
    plt.clf()
    del geoms
    gc.collect()
    assert len(cgeoaxes._GEOMETRY_TO_PATH_CACHE) == initial_cache_size

    plt.close()


def test_contourf_transform_path_counting():
    ax = plt.axes(projection=ccrs.Robinson())
    plt.draw()

    path_to_geos_counter = CallCounter(cartopy.mpl.patch, 'path_to_geos')

    with path_to_geos_counter:
        x, y, z = sample_data((30, 60))
        cs = plt.contourf(x, y, z, 5, transform=ccrs.PlateCarree())
        n_geom = sum([len(c.get_paths()) for c in cs.collections])
        del cs, c
        plt.draw()

    # before the performance enhancement, the count would have been 2 * n_geom,
    # but should now be just n_geom
    msg = ('The given geometry was transfomed too many times (expected: %s; '
           'got %s) - the caching is not working.'
           '' % (n_geom, path_to_geos_counter.count))
    assert path_to_geos_counter.count == n_geom, msg

    # Check the cache has an entry for each geometry.
    assert len(cgeoaxes._PATH_TRANSFORM_CACHE) == n_geom

    # Check that the cache is empty again once we've dropped all references
    # to the source paths.
    plt.clf()
    gc.collect()
    assert len(cgeoaxes._PATH_TRANSFORM_CACHE) == 0

    plt.close()


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
