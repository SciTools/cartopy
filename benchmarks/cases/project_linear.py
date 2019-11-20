# (C) British Crown Copyright 2019, Met Office
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
