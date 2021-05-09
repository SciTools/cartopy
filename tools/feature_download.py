#!/usr/bin/env python
# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
"""
This module provides a command-line tool for triggering the download of
the data used by various Feature instances.

For detail on how to use this tool, execute it with the `-h` option:

    python download.py -h

"""

import argparse

from cartopy import config
from cartopy.feature import Feature, GSHHSFeature, NaturalEarthFeature
from cartopy.io import Downloader


ALL_SCALES = ('110m', '50m', '10m')


FEATURE_DEFN_GROUPS = {
    # Only need one GSHHS resolution because they *all* get downloaded
    # from one file.
    'gshhs': GSHHSFeature(scale='f'),
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

        ('cultural', 'roads', '10m'),
        ('cultural', 'roads_north_america', '10m'),
        ('cultural', 'railroads', '10m'),
        ('cultural', 'railroads_north_america', '10m'),
    ),
    'cultural-extra': (
        ('cultural', 'admin_0_map_units', '110m'),
        ('cultural', 'admin_0_scale_rank', '110m'),
        ('cultural', 'admin_0_tiny_countries', '110m'),
        ('cultural', 'admin_0_pacific_groupings', '110m'),
        ('cultural', 'admin_1_states_provinces', '110m'),
        ('cultural', 'admin_1_states_provinces_lines', '110m'),
    ),
}


def download_features(group_names, dry_run=True):
    for group_name in group_names:
        feature_defns = FEATURE_DEFN_GROUPS[group_name]
        if isinstance(feature_defns, Feature):
            feature = feature_defns
            level = list(feature._levels)[0]
            downloader = Downloader.from_config(('shapefiles', 'gshhs',
                                                 feature._scale, level))
            format_dict = {'config': config, 'scale': feature._scale,
                           'level': level}
            if dry_run:
                print(f'URL: {downloader.url(format_dict)}')
            else:
                downloader.path(format_dict)
                geoms = list(feature.geometries())
                print(f'Feature {feature} length: {len(geoms)}')
        else:
            for category, name, scales in feature_defns:
                if not isinstance(scales, tuple):
                    scales = (scales,)
                for scale in scales:
                    downloader = Downloader.from_config(('shapefiles',
                                                         'natural_earth',
                                                         scale, category,
                                                         name))
                    feature = NaturalEarthFeature(category, name, scale)
                    format_dict = {'config': config, 'category': category,
                                   'name': name, 'resolution': scale}
                    if dry_run:
                        print(f'URL: {downloader.url(format_dict)}')
                    else:
                        downloader.path(format_dict)
                        geoms = list(feature.geometries())
                        print('Feature {}, {}, {} length: {}'
                              ''.format(category, name, scale, len(geoms)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download feature datasets.')
    parser.add_argument('group_names', nargs='+',
                        choices=FEATURE_DEFN_GROUPS,
                        metavar='GROUP_NAME',
                        help='Feature group name: %(choices)s')
    parser.add_argument('--output', '-o',
                        help='save datasets in the specified directory '
                             '(default: user cache directory)')
    parser.add_argument('--dry-run',
                        help='just print the URLs to download',
                        action='store_true')
    parser.add_argument('--ignore-repo-data', action='store_true',
                        help='ignore existing repo data when downloading')
    args = parser.parse_args()

    if args.output:
        config['pre_existing_data_dir'] = args.output
        config['data_dir'] = args.output
    if args.ignore_repo_data:
        config['repo_data_dir'] = config['data_dir']
    download_features(args.group_names, dry_run=args.dry_run)
