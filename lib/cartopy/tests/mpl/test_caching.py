# (C) British Crown Copyright 2011 - 2016, Met Office
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
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

import gc

import six
import unittest

try:
    from owslib.wmts import WebMapTileService
except ImportError as e:
    WebMapTileService = None
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.path import Path

import cartopy.crs as ccrs
from cartopy.mpl.feature_artist import FeatureArtist
from cartopy.io.ogc_clients import WMTSRasterSource, _OWSLIB_AVAILABLE
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
    # This test ensures it is working.

    # Create coastlines to ensure they are cached.
    ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
    ax1.coastlines()
    plt.draw()
    # Create another instance of the coastlines and count
    # the number of times shapereader.Reader is created.
    counter = CallCounter(cartopy.io.shapereader.Reader, '__init__')
    with counter:
        ax2 = plt.subplot(2, 1, 1, projection=ccrs.Robinson())
        ax2.coastlines()
        plt.draw()

    assert counter.count == 0, ('The shapereader Reader class was created {} '
                                'times, indicating that the caching is not '
                                'working.'.format(counter.count))

    plt.close()


def test_shapefile_transform_cache():
    # a5caae040ee11e72a62a53100fe5edc355304419 added shapefile mpl
    # geometry caching based on geometry object id. This test ensures
    # it is working.
    coastline_path = cartopy.io.shapereader.natural_earth(resolution="50m",
                                                          category='physical',
                                                          name='coastline')
    geoms = cartopy.io.shapereader.Reader(coastline_path).geometries()
    # Use the first 10 of them.
    geoms = tuple(geoms)[:10]
    n_geom = len(geoms)

    ax = plt.axes(projection=ccrs.Robinson())

    # Empty the cache.
    FeatureArtist._geom_key_to_geometry_cache.clear()
    FeatureArtist._geom_key_to_path_cache.clear()
    assert len(FeatureArtist._geom_key_to_geometry_cache) == 0
    assert len(FeatureArtist._geom_key_to_path_cache) == 0

    counter = CallCounter(ax.projection, 'project_geometry')
    with counter:
        ax.add_geometries(geoms, ccrs.PlateCarree())
        ax.add_geometries(geoms, ccrs.PlateCarree())
        ax.add_geometries(geoms[:], ccrs.PlateCarree())
        ax.figure.canvas.draw()

    # Without caching the count would have been
    # n_calls * n_geom, but should now be just n_geom.
    assert counter.count == n_geom, ('The given geometry was transformed too '
                                     'many times (expected: {}; got {}) - the'
                                     ' caching is not working.'
                                     ''.format(n_geom, counter.count))

    # Check the cache has an entry for each geometry.
    assert len(FeatureArtist._geom_key_to_geometry_cache) == n_geom
    assert len(FeatureArtist._geom_key_to_path_cache) == n_geom

    # Check that the cache is empty again once we've dropped all references
    # to the source paths.
    plt.clf()
    del geoms
    gc.collect()
    assert len(FeatureArtist._geom_key_to_geometry_cache) == 0
    assert len(FeatureArtist._geom_key_to_path_cache) == 0

    plt.close()


def test_contourf_transform_path_counting():
    ax = plt.axes(projection=ccrs.Robinson())
    ax.figure.canvas.draw()

    # Capture the size of the cache before our test.
    gc.collect()
    initial_cache_size = len(cgeoaxes._PATH_TRANSFORM_CACHE)

    path_to_geos_counter = CallCounter(cartopy.mpl.patch, 'path_to_geos')
    with path_to_geos_counter:
        x, y, z = sample_data((30, 60))
        cs = plt.contourf(x, y, z, 5, transform=ccrs.PlateCarree())
        n_geom = sum([len(c.get_paths()) for c in cs.collections])
        del cs
        if not six.PY3:
            del c
        ax.figure.canvas.draw()

    # Before the performance enhancement, the count would have been 2 * n_geom,
    # but should now be just n_geom.
    msg = ('The given geometry was transformed too many times (expected: {}; '
           'got {}) - the caching is not working.'
           '').format(n_geom, path_to_geos_counter.count)
    assert path_to_geos_counter.count == n_geom, msg

    # Check the cache has an entry for each geometry.
    assert len(cgeoaxes._PATH_TRANSFORM_CACHE) == initial_cache_size + n_geom

    # Check that the cache is empty again once we've dropped all references
    # to the source paths.
    plt.clf()
    gc.collect()
    assert len(cgeoaxes._PATH_TRANSFORM_CACHE) == initial_cache_size

    plt.close()


@unittest.skipIf(not _OWSLIB_AVAILABLE, 'OWSLib is unavailable.')
def test_wmts_tile_caching():
    image_cache = WMTSRasterSource._shared_image_cache
    image_cache.clear()
    assert len(image_cache) == 0

    url = 'https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi'
    wmts = WebMapTileService(url)
    layer_name = 'MODIS_Terra_CorrectedReflectance_TrueColor'

    source = WMTSRasterSource(wmts, layer_name)

    gettile_counter = CallCounter(wmts, 'gettile')
    crs = ccrs.PlateCarree()
    extent = (-180, 180, -90, 90)
    resolution = (20, 10)
    with gettile_counter:
        source.fetch_raster(crs, extent, resolution)
    n_tiles = 2
    assert gettile_counter.count == n_tiles, ('Too many tile requests - '
                                              'expected {}, got {}.'.format(
                                                  n_tiles,
                                                  gettile_counter.count)
                                              )
    gc.collect()
    assert len(image_cache) == 1
    assert len(image_cache[wmts]) == 1
    tiles_key = (layer_name, '0')
    assert len(image_cache[wmts][tiles_key]) == n_tiles

    # Second time around we shouldn't request any more tiles so the
    # call count will stay the same.
    with gettile_counter:
        source.fetch_raster(crs, extent, resolution)
    assert gettile_counter.count == n_tiles, ('Too many tile requests - '
                                              'expected {}, got {}.'.format(
                                                  n_tiles,
                                                  gettile_counter.count)
                                              )
    gc.collect()
    assert len(image_cache) == 1
    assert len(image_cache[wmts]) == 1
    tiles_key = (layer_name, '0')
    assert len(image_cache[wmts][tiles_key]) == n_tiles

    # Once there are no live references the weak-ref cache should clear.
    del source, wmts, gettile_counter
    gc.collect()
    assert len(image_cache) == 0


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
