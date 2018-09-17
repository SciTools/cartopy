# (C) British Crown Copyright 2017 - 2018, Met Office
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

import os
import pytest
import cartopy.feature as cfeature
from cartopy import config

small_extent = (-6, -8, 56, 59)
medium_extent = (-20, 20, 20, 60)
large_extent = (-40, 40, 0, 80)

auto_scaler = cfeature.AdaptiveScaler('110m', (('50m', 50), ('10m', 15)))

auto_land = cfeature.NaturalEarthFeature('physical', 'land', auto_scaler)


class TestFeatures(object):
    def test_change_scale(self):
        # Check that features can easily be retrieved with a different
        # scale.
        new_lakes = cfeature.LAKES.with_scale('10m')
        assert new_lakes.scale == '10m'
        assert new_lakes.kwargs == cfeature.LAKES.kwargs
        assert new_lakes.category == cfeature.LAKES.category
        assert new_lakes.name == cfeature.LAKES.name

    def test_scale_from_extent(self):
        # Check that scaler.scale_from_extent returns the appropriate
        # scales.
        small_scale = auto_land.scaler.scale_from_extent(small_extent)
        medium_scale = auto_land.scaler.scale_from_extent(medium_extent)
        large_scale = auto_land.scaler.scale_from_extent(large_extent)
        assert small_scale == '10m'
        assert medium_scale == '50m'
        assert large_scale == '110m'

    def test_intersecting_geometries_small(self):
        # Check that intersecting_geometries will set the scale to
        # '10m' when the extent is small and autoscale is True.
        auto_land.intersecting_geometries(small_extent)
        assert auto_land.scale == '10m'

    def test_intersecting_geometries_medium(self):
        # Check that intersecting_geometries will set the scale to
        # '50m' when the extent is medium and autoscale is True.
        auto_land.intersecting_geometries(medium_extent)
        assert auto_land.scale == '50m'

    def test_intersecting_geometries_large(self):
        # Check that intersecting_geometries will set the scale to
        # '110m' when the extent is large and autoscale is True.
        auto_land.intersecting_geometries(large_extent)
        assert auto_land.scale == '110m'

    @pytest.mark.parametrize("disk_caching", [False, True])
    def test_gshhs(self, disk_caching):

        scale = 'l'
        extent = [-10, 0, 45, 50]
        config['allow_disk_caching'] = disk_caching

        cfeature.GSHHSFeature._geometries_cache = {}
        gshhs = cfeature.GSHHSFeature(scale, level=1)

        if disk_caching:
            cache_file = cfeature.FeatureDiskCaching(gshhs, extent,
                                                     scale=scale,
                                                     level=1).cache_file
            if os.path.exists(cache_file):
                os.remove(cache_file)

        # First read
        geoms0 = list(gshhs.intersecting_geometries(extent))
        assert geoms0

        # From mem cache
        geoms1 = list(gshhs.intersecting_geometries(extent))
        assert len(geoms0) == len(geoms1)
        assert geoms0[0] is geoms1[0]

        # From disk cache
        if disk_caching:
            assert os.path.exists(cache_file)
            cfeature.GSHHSFeature._geometries_cache = {}
            geoms2 = list(gshhs.intersecting_geometries(extent))
            assert len(geoms2) == len(geoms0)
