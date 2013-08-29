#!/usr/bin/env python
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
"""
This module provides a command-line tool for triggering the download of
the data used by various Feature instances.

For detail on how to use this tool, execute it with the `-h` option:

    python download.py -h

"""
from __future__ import print_function

import argparse

from cartopy.feature import Feature, GSHHSFeature, NaturalEarthFeature
from cartopy.crs import PlateCarree
import matplotlib.pyplot as plt


ALL_SCALES = ('110m', '50m', '10m')


FEATURE_DEFN_GROUPS = {
    # Only need one GSHHS resolution because they *all* get downloaded
    # from one file.
    'gshhs': GSHHSFeature(scale='c'),
    'physical': (
        ('physical', 'coastline', ALL_SCALES),
        ('physical', 'land', ALL_SCALES),
        ('physical', 'ocean', ALL_SCALES),
        ('physical', 'rivers_lake_centerlines', ALL_SCALES),
        ('physical', 'lakes', ALL_SCALES),
        ('physical', 'geography_regions_polys', ALL_SCALES),
        ('physical', 'geography_regions_points', ALL_SCALES),
        ('physical', 'geography_marine_polys', ALL_SCALES),
        ('physical', 'glaciated_areas', ALL_SCALES)
    ),
    'cultural': (
        ('cultural', 'admin_0_countries', ALL_SCALES),
        ('cultural', 'admin_0_countries_lakes', ALL_SCALES),
        ('cultural', 'admin_0_sovereignty', ALL_SCALES),
        ('cultural', 'admin_0_boundary_lines_land', ALL_SCALES),

        ('cultural', 'urban_areas', ('50m', '10m')),

        #('cultural', 'roads', '10m'), # ERROR in NE dataset?
        ('cultural', 'roads_north_america', '10m'),
        ('cultural', 'railroads', '10m'),
        ('cultural', 'railroads_north_america', '10m'),
    ),
    'cultural-extra': (
        ('cultural', 'admin_0_map_units', '110m'),
        ('cultural', 'admin_0_scale_rank', '110m'),
        ('cultural', 'admin_0_tiny_countries', '110m'),
        ('cultural', 'admin_0_pacific_groupings', '110m'),
        ('cultural', 'admin_1_states_provinces_shp', '110m'),
        ('cultural', 'admin_1_states_provinces_lines', '110m'),
    ),
}


def download_features(group_names, hold):
    plt.ion()
    ax = plt.axes(projection=PlateCarree())
    ax.set_global()
    for group_name in group_names:
        feature_defns = FEATURE_DEFN_GROUPS[group_name]
        if isinstance(feature_defns, Feature):
            features = [feature_defns]
        else:
            features = []
            for category, name, scales in feature_defns:
                if not isinstance(scales, tuple):
                    scales = (scales,)
                for scale in scales:
                    features.append(NaturalEarthFeature(category, name, scale))
        for feature in features:
            ax.add_feature(feature)
            plt.draw()

    plt.ioff()
    if hold:
        plt.show()


if __name__ == '__main__':
    def group_name(string):
        if string not in FEATURE_DEFN_GROUPS:
            msg = '{!r} is not a valid feature group (choose from {!s})'
            msg = msg.format(string, list(FEATURE_DEFN_GROUPS.keys()))
            raise argparse.ArgumentTypeError(msg)
        return string

    parser = argparse.ArgumentParser(description='Download feature datasets.')
    parser.add_argument('group_names', nargs='*',
                        type=group_name,
                        metavar='GROUP_NAME',
                        help='Feature group name')
    parser.add_argument('--hold', action='store_true',
                        help='keep the matplotlib window open')
    parser.add_argument('--show', action='store_true',
                        help='show the list of valid feature group names')
    args = parser.parse_args()
    if args.show:
        print('Feature group names:')
        for name in sorted(FEATURE_DEFN_GROUPS.keys()):
            print('   ', name)
    elif not args.group_names:
       parser.error('Please supply one or more feature group names.')
    download_features(args.group_names, args.hold)
