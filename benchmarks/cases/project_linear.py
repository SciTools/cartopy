# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import shapely.geometry as sgeom


class Oceans:
    def prepare(self):
        shpfilename = shpreader.natural_earth(
            resolution='50m', category='physical', name='ocean')
        reader = shpreader.Reader(shpfilename)
        oceans = reader.geometries()
        oceans = sgeom.MultiPolygon(oceans)
        self.geoms = oceans


OCEAN = Oceans()


def use_setup(setup_fn):
    # A decorator to create a decorator...
    def decorator(test_func):
        # This decorator attaches the setup function to the test.
        test_func.setup = setup_fn
        return test_func
    return decorator


@use_setup(OCEAN.prepare)
def time_ocean_pc():
    ccrs.PlateCarree().project_geometry(OCEAN.geoms)


@use_setup(OCEAN.prepare)
def time_ocean_np():
    ccrs.NorthPolarStereo().project_geometry(OCEAN.geoms)


@use_setup(OCEAN.prepare)
def time_ocean_rob():
    ccrs.Robinson().project_geometry(OCEAN.geoms)


@use_setup(OCEAN.prepare)
def time_ocean_igh():
    ccrs.InterruptedGoodeHomolosine().project_geometry(OCEAN.geoms)
